#' @title Estimate Time in Range (TIR)
#' @description
#' Estimates the mean Time in Range (TIR) using different methods.
#'
#' This function estimates TIR using either the naive estimator, the proposed
#' noninformative follow-up estimator, or the Cox model-based estimator.
#'
#' @param data A data frame containing glucose monitoring data.
#' @param method A string specifying the estimation method: "naive" or "proposed".
#' @param model A string specifying the model type (default: "NULL" for non-model-based approaches, or "cox" for Cox model).
#' @param time A numeric vector of length 2 specifying the start and end time window (default: c(0, 10075)).
#' @param range A numeric vector specifying the glucose range for in-range classification (default: c(70, 180)).
#' @param boot Number of bootstrap iterations (default: NULL, no bootstrapping).
#' @param id Column name for subject IDs (default: "patient_id").
#' @param glucose Column name for glucose values (default: "glucose").
#' @param time_col Column name representing time (default: "time").
#' @param period Additional time offset for stop time in Cox model (default: 5).
#' @param formula A formula string for the Cox model (default: "var1").
#' @return A list containing TIR estimation results, including the estimate, standard error, and confidence intervals.
#' @export
estTIR <- function(data, method = "proposed", model = "NULL",
                   time = c(0, 1440 * 7 - 5), range = c(70, 180),
                   boot = NULL, id = "patient_id", glucose = "glucose",
                   time_col = "time", period = 5, formula = "var1") {
  # Compute value_in_range based on the specified range
  data <- data |>
    dplyr::mutate(value_in_range = sapply(.data[[glucose]], value_in_range, range = range))
  # Generate event column for Cox model (1 if last observation for a subject, otherwise 0)
  data <- data |> dplyr::group_by(.data[[id]]) |> dplyr::mutate(event = dplyr::if_else(.data[[time_col]] == max(.data[[time_col]]),1,0)) |> dplyr::ungroup()
  # Compute stop time for Cox model
  data$time2 <- data[[time_col]] + period

  # Initialize `est` variable before assignment
  est <- NULL
  # Use switch() for cleaner method selection
  est <- switch(method,
         "naive" = {
           if (model == "NULL") {
             naive_est(data, min_time = time[1], max_time = time[2],
                       boot = boot, id_col = id, time = time_col,
                       value_in_range = "value_in_range")
           } else stop("Error: Model not recognized for 'naive' method")
         },
         "proposed" = {
           switch(model,
                  "NULL" = proposed_est_noninfo(data, min_time = time[1], max_time = time[2],
                                                boot = boot, id_col = id, time = time_col,
                                                value_in_range = "value_in_range"),
                  "cox" = proposed_est_cox(data, min_time = time[1], max_time = time[2],
                                           id_col = id, event_col = "event",
                                           start_col = time_col, stop_col = "time2",
                                           formula = formula, boot = boot,
                                           value_in_range = "value_in_range"),
                  stop("Error: Model not recognized for 'proposed' method")
           )
         },
         stop("Error: Method not recognized")
  )

  # Ensure `est` is not NULL before assigning the class
  if (is.null(est)) stop("Error: Estimation method failed to return a result")

  # Assign "TIR" class to result
  class(est) <- "TIR"
  return(est)
}
