#' @title Recursively rounds numeric values in nested lists, data frames, and tibbles.
#' @description
#' This function applies rounding to numeric values within a list, including nested structures.
#' If the input is a numeric vector, it rounds it directly. Non-numeric values remain unchanged.
#'
#' @param data A list, data frame, tibble, or numeric vector to be rounded.
#' @param decimals Number of decimal places to round to (default = 3).
#' @return The input data with all numeric values rounded.
#' @export
round_nested <- function(data, decimals = 3) {
  if (is.list(data)) {
    return(lapply(data, round_nested, decimals))
  } else if (is.numeric(data)) {
    return(round(data, decimals))
  }
  return(data)
}
