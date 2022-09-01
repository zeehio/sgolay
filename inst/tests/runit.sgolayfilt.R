

test01_sgolayfilt_vectors <- function() {
  filt <- sgolay::sgolay(p = 2, n = 5)
  x <- runif(500)
  yref <- signal::sgolayfilt(x, filt)
  y_filter <- sgolayfilt(x, filt, engine = "filter")
  y_fft <- sgolayfilt(x, filt, engine = "fft")
  checkEquals(yref, y_filter, msg = "Check sgolayfilt on a vector (filter engine)")
  checkEquals(yref, y_fft, msg = "Check sgolayfilt on a vector (fft engine)")
}

test02_sgolayfilt_matrix_rowwisefalse <- function() {
  filt <- sgolay::sgolay(p = 2, n = 5)
  x <- matrix(runif(1000), nrow = 500, ncol = 2)
  yref <- apply(x, 2, function(xin) signal::sgolayfilt(xin, filt))
  y_filter <- sgolayfilt(x, filt, engine = "filter")
  y_fft <- sgolayfilt(x, filt, engine = "fft")
  checkEquals(yref, y_filter, msg = "Check sgolayfilt on a matrix (rowwise=FALSE, filter engine)")
  checkEquals(yref, y_fft, msg = "Check sgolayfilt on a matrix (rowwise=FALSE, fft engine)")
}


test02_sgolayfilt_matrix_rowwisetrue <- function() {
  filt <- sgolay::sgolay(p = 2, n = 5)
  x <- matrix(runif(1000), nrow = 2, ncol = 500)
  yref <- t(apply(x, 1, function(xin) signal::sgolayfilt(xin, filt)))
  y_filter <- sgolayfilt(x, filt, engine = "filter", rowwise = TRUE)
  y_fft <- sgolayfilt(x, filt, engine = "fft", rowwise = TRUE)
  checkEquals(yref, y_filter, msg = "Check sgolayfilt on a matrix (rowwise=TRUE, filter engine)")
  checkEquals(yref, y_fft, msg = "Check sgolayfilt on a matrix (rowwise=TRUE, fft engine)")
}
