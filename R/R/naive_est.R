#' @title Computes the naive estimator for mean Time in Range (TIR).
#' @description
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
#' @return A list containing TIR estimation results, including the estimate, standard error, and confidence intervals.
#' @export
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
