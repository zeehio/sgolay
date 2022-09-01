#' @useDynLib sgolay cconvolve_circular
convolve_circular <- function(x, y, conj=TRUE) {
  .Call(cconvolve_circular, x, y, conj)
}
