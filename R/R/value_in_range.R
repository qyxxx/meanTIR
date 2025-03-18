#' @title Check if a number is within the specified range.
#' @description
#' Check if a number is within the specified range.
#' @param x A numeric value to check.
#' @param range A numeric vector of length 2 specifying the lower and upper bounds.
#' @return An integer: 1 if x is within the range, 0 otherwise.
#' @export
value_in_range <- function(x, range = c(-.Machine$integer.max, .Machine$integer.max)) {
  if (range[1] <= x & x <= range[2]) {
    return(1)
  } else {
    return(0)
  }
}
