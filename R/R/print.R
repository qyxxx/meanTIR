#' @title Print TIR Estimation Results
#' @description
#' Prints Time in Range (TIR) estimation results, ensuring numeric values are rounded.
#'
#' @param x An object of class "TIR" (a list containing estimation results).
#' @param decimals Number of decimal places to round to (default = 3).
#' @param ... Additional arguments (for compatibility with generic print method).
#' @export
#' @method print TIR
print.TIR <- function(x, decimals = 3, ...) {
  # Ensure x is of class "TIR"
  if (!inherits(x, "TIR")) {
    stop("Error: Object is not of class 'TIR'")
  }

  # Round numeric values within the list
  rounded_x <- round_nested(x, decimals)

  # Convert to data frame for better printing
  DF <- as.data.frame(t(unlist(rounded_x)))  # Transpose for better display

  # Set row name
  rownames(DF) <- "meanTIR"

  # Print formatted output
  print(DF, row.names = TRUE)
}
