

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

test03_sgolayfilt_plain_matrix_p <- function() {
  # p can be a plain matrix without the "sgolayFilter" class (e.g. after
  # unclass()), not only an object returned directly by sgolay()
  filt_classed <- sgolay::sgolay(p = 2, n = 5)
  filt <- unclass(filt_classed)
  x <- runif(500)
  yref <- signal::sgolayfilt(x, filt_classed)
  y_filter <- sgolayfilt(x, filt, engine = "filter")
  y_fft <- sgolayfilt(x, filt, engine = "fft")
  checkEquals(yref, y_filter, msg = "Check sgolayfilt with a plain matrix p (filter engine)")
  checkEquals(yref, y_fft, msg = "Check sgolayfilt with a plain matrix p (fft engine)")
}

test04_sgolayfilt_even_rows_p_errors <- function() {
  # An even number of rows has no well-defined center row and must be rejected
  even_filt <- matrix(1, nrow = 4, ncol = 4)
  x <- runif(50)
  checkException(sgolayfilt(x, even_filt), msg = "Check sgolayfilt rejects an even-row filter matrix p")
}

test05_sgolayfilt_short_x_errors <- function() {
  # x shorter than the filter length must be rejected with a clear error
  # instead of an obscure out-of-bounds indexing failure
  filt <- sgolay::sgolay(p = 2, n = 5)
  checkException(sgolayfilt(runif(3), filt), msg = "Check sgolayfilt rejects a vector shorter than the filter")
  x <- matrix(runif(2 * 3), nrow = 2, ncol = 3)
  checkException(sgolayfilt(x, filt, rowwise = TRUE), msg = "Check sgolayfilt rejects a rowwise matrix shorter than the filter")
}

test06_sgolayfilt_x_length_equals_filter_length <- function() {
  # x exactly as long as the filter is the smallest valid input and must
  # still work (boundary just above the rejected case in test05)
  filt <- sgolay::sgolay(p = 2, n = 5)
  x <- runif(5)
  yref <- signal::sgolayfilt(x, filt)
  y_filter <- sgolayfilt(x, filt, engine = "filter")
  y_fft <- sgolayfilt(x, filt, engine = "fft")
  checkEquals(yref, y_filter, msg = "Check sgolayfilt when length(x) == n (filter engine)")
  checkEquals(yref, y_fft, msg = "Check sgolayfilt when length(x) == n (fft engine)")
}
