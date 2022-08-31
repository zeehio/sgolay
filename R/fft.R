#' @importFrom stats fft
convolve_prepare <- function(x, conj = FALSE, plan = NULL, impl = "auto") {
  if (conj) {
    do_conj <- Conj
  } else {
    do_conj <- identity
  }
  do_conj(fft(x))
}

#' @importFrom stats fft
convolve_do <- function(fft_x, conj_fft_y) {
  Re(fft(fft_x * conj_fft_y, inverse = TRUE))/length(fft_x)
}
