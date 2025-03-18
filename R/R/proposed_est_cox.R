#' @description
#' Computes the proposed estimator for mean Time in Range (TIR) using a Cox model.
#'
#' This function estimates TIR using a Cox proportional hazards model, weighting observations
#' based on their predicted survival. If bootstrapping is enabled, confidence intervals are computed.
#'
#' @param data A data frame containing time-to-event data.
#' @param min_time The minimum time threshold for filtering (default = 0).
#' @param max_time The maximum time threshold for filtering (default = 10075, equivalent to 1440 * 7 - 5).
#' @param id_col Column name identifying subjects (default = "patient_id").
#' @param event_col Column name representing event occurrence (default = "event").
#' @param start_col Column name representing the start time (default = "time").
#' @param stop_col Column name representing the stop time (default = "time2").
#' @param formula Right-hand side of the Cox model formula as a character string (default = "var1").
#' @param boot Number of bootstrap iterations (default = NULL, no bootstrapping).
#' @param value_in_range Column name of the binary indicator (0/1) for time in range (default = "value_in_range").
#' @return A list containing TIR estimation results, including the estimate, standard error, and confidence intervals.
#' @export
proposed_est_cox <- function(data, min_time = 0, max_time = (1440 * 7 - 5),
                             id_col = "patient_id", event_col = "event",
                             start_col = "time", stop_col = "time2", formula = "var1",
                             boot = NULL, value_in_range = "value_in_range") {
  # data <- data.table::as.data.table(data)
  # Fit Cox model
  cox_fit <- survival::coxph(as.formula(paste0("survival::Surv(", start_col, ",", stop_col, ",", "event==1", ") ~ ", formula)), data = data, method="breslow")
  # Baseline cumulative hazard
  baseline_hazard <- survival::basehaz(cox_fit, centered = FALSE)
  # Add time and hazard difference
  baseline_hazard <- baseline_hazard |> dplyr::arrange(time) |> dplyr::mutate(hazard_diff = c(0, diff(hazard)))
  # Predict partial hazard
  data[, "predict_partial_hazard"] <- predict(cox_fit, newdata = data, type = "risk")
  # Merge cumulative hazard with dataset
  data <- merge(data, baseline_hazard, by.x = start_col, by.y = "time", all.x = TRUE)
  data <- data |>
    dplyr::mutate(across(everything(), ~ tidyr::replace_na(.x, 0)))
  # Compute lambda_exp_diff
  data <- data |> dplyr::mutate(lambda_exp_diff = predict_partial_hazard * hazard_diff)
  data <- data |> dplyr::group_by(.data[[id_col]]) |>
    dplyr::mutate(cum_lambda_exp_diff = cumsum(lambda_exp_diff), weight = 1 / exp(-cum_lambda_exp_diff)) |> dplyr::ungroup()

  # Calculate TIR
  data <- data[data[[start_col]] <= max_time & data[[start_col]] >= min_time, ]
  TIR <- data |>
    dplyr::group_by(.data[[start_col]]) |>
    dplyr::summarise(weighted_avg = weighted.mean(.data[[value_in_range]], weight), .groups = "drop") |>
    dplyr::summarise(TIR = mean(weighted_avg), .groups = "drop") |>
    dplyr::pull(TIR)

  if (is.null(boot)) {
    return(list(est = TIR))
  } else {
    unique_ids <- unique(data[[id_col]])
    n_ids <- length(unique_ids)
    boot_TIR <- replicate(boot, {
      boot_id <- data.frame(ID = sample((unique_ids), n_ids, replace = T))
      colnames(boot_id) <- id_col
      boot_sample_temp <- boot_id |> dplyr::left_join(data, by = id_col, relationship = "many-to-many")
      count_boot <- NULL
      boot_sample <- boot_sample_temp |>
        dplyr::group_by(.data[[id_col]], .data[[start_col]]) |>
        dplyr::mutate(count_boot = dplyr::row_number()) |>
        dplyr::ungroup() |>
        dplyr::mutate(ID_boot = ifelse(count_boot > 1, paste0(.data[[id_col]], "BOOT", count_boot), .data[[id_col]])) |>
        dplyr::select(-count_boot, -.data[[id_col]], -predict_partial_hazard, -hazard, -hazard_diff, -lambda_exp_diff) |>
        dplyr::rename(patient_id = ID_boot)
      proposed_est_cox(boot_sample, min_time, max_time, "patient_id", event_col, start_col, stop_col, formula, boot = NULL, value_in_range)$est
    })
    return(list(
      est = TIR,
      `std err` = sd(boot_TIR),
      `CI 025` = quantile(boot_TIR, 0.025),
      `CI 975` = quantile(boot_TIR, 0.975)
    ))
  }
}
