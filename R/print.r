#' @export
#' @title Print
#' @method print SCGLR
#' @S3method print SCGLR
#' @param x SCGLR result to print
#' @param ... unused
print.SCGLR <- function(x, ...) {
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), sep = "","\n")
  cat("\nInertia:\n")
  print.default(x$inertia,print.gap=2)
  cat("\nDeviance:\n")
  print.default(x$deviance,print.gap=2)
}