#' @description
#' Check if a number is within the specified range.
#' @param x A numeric value to check.
#' @param range A numeric vector of length 2 specifying the lower and upper bounds.
#' @return An integer: 1 if x is within the range, 0 otherwise.
#' @noRd
value_in_range <- function(x, range = c(-.Machine$integer.max, .Machine$integer.max)) {
  if (range[1] <= x & x <= range[2]) {
    return(1)
  } else {
    return(0)
  }
}

#' @description
#' Computes the naive estimator for mean Time in Range (TIR).
#'
#' This function estimates the mean TIR for a given dataset, filtering based on a
#' specified time window. If bootstrapping is enabled, it computes confidence
#' intervals using resampling.
#'
#' @param data A data frame containing long format time-series data.
#' @param min_time The minimum time threshold for filtering (default = 0).
#' @param max_time The maximum time threshold for filtering (default = 10075, equivalent to 1440 * 7 - 5).
#' @param boot Number of bootstrap iterations (default = NULL, no bootstrapping).
#' @param id_col Column name identifying subjects (default = "patient_id").
#' @param time Column name representing time (default = "time").
#' @param value_in_range Column name of the binary indicator (0/1) for time in range (default = "value_in_range").
#' @return A list containing:
#'   \item{est}{Estimated mean TIR.}
#'   \item{std err}{Standard error of the estimate (if bootstrapping is enabled).}
#'   \item{CI 025}{Lower bound of the 95% confidence interval (if bootstrapping is enabled).}
#'   \item{CI 975}{Upper bound of the 95% confidence interval (if bootstrapping is enabled).}
#' @noRd
naive_est <- function(data,
                      min_time = 0, max_time = (1440 * 7 - 5),
                      boot = NULL,
                      id_col = "patient_id", time = "time", value_in_range = "value_in_range") {
  # Filter data based on time range
  data <- data[data[[time]] <= max_time & data[[time]] >= min_time, ]
  # Compute mean TIR per subject, then overall mean TIR
  TIR <- data|>
    dplyr::group_by(.data[[id_col]]) |>
    dplyr::summarise(TIR_i = mean(.data[[value_in_range]], na.rm = TRUE), .groups = "drop") |>
    dplyr::summarise(TIR = mean(TIR_i), .groups = "drop") |>
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
        dplyr::group_by(.data[[id_col]], .data[[time]]) |>
        dplyr::mutate(count_boot = dplyr::row_number()) |>
        dplyr::ungroup() |>
        dplyr::mutate(ID_boot = ifelse(count_boot > 1, paste0(.data[[id_col]], "BOOT", count_boot), .data[[id_col]])) |>
        dplyr::select(-count_boot)
      naive_est(boot_sample, min_time, max_time, boot = NULL, "ID_boot", time, value_in_range)$est
    })
    return(list(
      est = TIR,
      `std err` = sd(boot_TIR),
      `CI 025` = quantile(boot_TIR, 0.025),
      `CI 975` = quantile(boot_TIR, 0.975)
    ))
  }
}

#' @description
#' Computes the proposed estimator for mean Time in Range (TIR) under the
#' assumption of noninformative follow-up duration.
#'
#' This function estimates the mean TIR using a time-stratified averaging
#' approach. If bootstrapping is enabled, it computes confidence intervals
#' using resampling.
#'
#' @param data A data frame containing time-series data.
#' @param min_time The minimum time threshold for filtering (default = 0).
#' @param max_time The maximum time threshold for filtering (default = 10075, equivalent to 1440 * 7 - 5).
#' @param boot Number of bootstrap iterations (default = NULL, no bootstrapping).
#' @param id_col Column name identifying subjects (default = "patient_id").
#' @param time Column name representing time (default = "time").
#' @param value_in_range Column name of the binary indicator (0/1) for time in range (default = "value_in_range").
#' @return A list containing:
#'   \item{est}{Estimated mean TIR.}
#'   \item{std err}{Standard error of the estimate (if bootstrapping is enabled).}
#'   \item{CI 025}{Lower bound of the 95% confidence interval (if bootstrapping is enabled).}
#'   \item{CI 975}{Upper bound of the 95% confidence interval (if bootstrapping is enabled).}
#' @noRd
proposed_est_noninfo <- function(data,
                                 min_time = 0, max_time = (1440 * 7 - 5),
                                 boot = NULL,
                                 id_col = "patient_id", time = "time", value_in_range = "value_in_range") {
  data <- data[data[[time]] <= max_time & data[[time]] >= min_time, ]
  TIR <- data |>
    dplyr::group_by(.data[[time]]) |>
    dplyr::summarise(avg_val = mean(value_in_range, na.rm = TRUE), .groups = "drop") |>
    dplyr::summarise(TIR = mean(avg_val), .groups = "drop") |>
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
        dplyr::group_by(.data[[id_col]], .data[[time]]) |>
        dplyr::mutate(count_boot = dplyr::row_number()) |>
        dplyr::ungroup() |>
        dplyr::mutate(ID_boot = ifelse(count_boot > 1, paste0(.data[[id_col]], "BOOT", count_boot), .data[[id_col]])) |>
        dplyr::select(-count_boot)
      proposed_est_noninfo(boot_sample, min_time, max_time, boot = NULL, "ID_boot", time, value_in_range)$est
    })
    return(list(
      est = TIR,
      `std err` = sd(boot_TIR),
      `CI 025` = quantile(boot_TIR, 0.025),
      `CI 975` = quantile(boot_TIR, 0.975)
    ))
  }
}

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
#' @return A list containing:
#'   \item{est}{Estimated mean TIR.}
#'   \item{std err}{Standard error of the estimate (if bootstrapping is enabled).}
#'   \item{CI 025}{Lower bound of the 95% confidence interval (if bootstrapping is enabled).}
#'   \item{CI 975}{Upper bound of the 95% confidence interval (if bootstrapping is enabled).}
#' @noRd
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


#' @description
#' Rounding number
#' @noRd
#'
round_nested <- function(data, decimals = 3) {
  if (is.list(data)) {
    return(lapply(data, round_nested, decimals))
  } else if (is.numeric(data)) {
    return(round(data, decimals))
  }
  return(data)
}

#' @description
#' Print results
#' @noRd
#' @export
#'
printTIR <- function(est, decimals = 3) {
  DF <- data.frame(round_nested(est, decimals))
  rownames(DF) <- 'meanTIR'
  print(DF)
}
