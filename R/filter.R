#' @useDynLib sgolay cfilter
filter <- function(x, y) {
  .Call(cfilter, x, y)
}
