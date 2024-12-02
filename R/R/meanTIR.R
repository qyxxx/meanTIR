#' @description
#' Indicating whether a number is within the target range
#' @noRd
#'
value_in_range <- function(x, range = c(-.Machine$integer.max, .Machine$integer.max)) {
  if (range[1] <= x & x <= range[2]) {
    return(1)
  } else {
    return(0)
  }
}

#' @description
#' The naive estimator for mean TIR
#' @noRd
#'
naive_est <- function(data, min_time = 0, max_time = (1440 * 7 - 5), boot = NULL, id_col = "patient_id", time = "time", value_in_range = "value_in_range") {
  data <- data[data[[time]] <= max_time & data[[time]] >= min_time, ]
  TIR <- data|>
    dplyr::group_by_at(id_col) |>
    dplyr::summarise(TIR_i = mean(value_in_range, na.rm = TRUE)) |>
    dplyr::summarise(TIR = mean(TIR_i)) |>
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
#' The proposed estimator with noninformative follow-up duration assumption
#' @noRd
#'
proposed_est_noninfo <- function(data, min_time = 0, max_time = (1440 * 7 - 5), boot = NULL, id_col = "patient_id", time = "time", value_in_range = "value_in_range") {
  data <- data[data[[time]] <= max_time & data[[time]] >= min_time, ]
  TIR <- data |>
    dplyr::group_by_at(time) |>
    dplyr::summarise(avg_val = mean(value_in_range, na.rm = TRUE)) |>
    dplyr::summarise(TIR = mean(avg_val)) |>
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
      `CI 025` = quantile(boot_TIR, 0.975)
    ))
  }
}

#' @description
#' The proposed estimator with Cox model
#' @noRd
#'
proposed_est_cox <- function(data, min_time = 0, max_time = (1440 * 7 - 5), id_col = "patient_id", event_col = "event", start_col = "time", stop_col = "time2", formula = "var1", boot = NULL, value_in_range = "value_in_range") {
  data <- as.data.table(data)

  # Fit Cox model
  cox_fit <- coxph(as.formula(paste0("Surv(", start_col, ", ", stop_col, ", ", event_col, ") ~ ", formula)), data = data)
  baseline_hazard <- basehaz(cox_fit, centered = FALSE)
  data[, "predict_partial_hazard"] <- predict(cox_fit, newdata = data, type = "risk")

  # Calculate weights
  data <- merge(data, baseline_hazard, by.x = start_col, by.y = "time", all.x = TRUE)
  data[, "lambda_exp_diff" := c(0, diff(hazard)) * predict_partial_hazard]
  data[, "cum_lambda_exp_diff" := cumsum(lambda_exp_diff), by = id_col]
  data[, "weight" := 1 / exp(-cum_lambda_exp_diff)]

  # Calculate TIR
  data <- data[data[[start_col]] <= max_time & data[[start_col]] >= min_time, ]
  TIR <- data |>
    dplyr::group_by_at(start_col) |>
    dplyr::summarise(weighted_avg = weighted.mean(get(value_in_range), weight)) |>
    dplyr::summarise(TIR = mean(weighted_avg)) |>
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
      proposed_est_cox(boot_sample, min_time, max_time, "ID_boot", event_col, start_col, stop_col, formula, boot = NULL, value_in_range)$est
    })
    return(list(
      est = TIR,
      `std err` = sd(boot_TIR),
      `CI 025` = quantile(boot_TIR, 0.025),
      `CI 025` = quantile(boot_TIR, 0.975)
    ))
  }
}

#' @description
#' Estimation of mean TIR
#'
estTIR <- function(data, method = "proposed", model = "NULL", time = c(0, 1440 * 7 - 5), range = c(70, 180), boot = NULL, id = "patient_id", glucose = "glucose", time_col = "time", period = 5, formula = "var1") {
  data$value_in_range <- sapply(data[[glucose]], value_in_range, range = range)
  data$event <- FALSE
  data[which.max(data[[time_col]]), "event"] <- TRUE
  data$time2 <- data[[time_col]] + period

  if (method == "naive" && model == "NULL") {
    return(naive_est(data, min_time = time[1], max_time = time[2], boot = boot, id_col = id, time = time_col, value_in_range = "value_in_range"))
  } else if (method == "proposed" && model == "NULL") {
    return(proposed_est_noninfo(data, min_time = time[1], max_time = time[2], boot = boot, id_col = id, time = time_col, value_in_range = "value_in_range"))
  } else if (method == "proposed" && model == "cox") {
    return(proposed_est_cox(data, min_time = time[1], max_time = time[2], id_col = id, event_col = "event", start_col = time_col, stop_col = "time2", formula = formula, boot = boot, value_in_range = "value_in_range"))
  } else {
    stop("Error: model not recognized")
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
#'
printTIR <- function(est, decimals = 3) {
  DF <- data.frame(round_nested(est, decimals))
  rownames(DF) <- 'meanTIR'
  print(DF)
}
