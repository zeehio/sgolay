


test_that("sgolayfilt results are all equal for a vector", {
  filt <- sgolay::sgolay(p = 2, n = 5)
  x <- runif(500)
  yref <- signal::sgolayfilt(x, filt)
  y_filter <- sgolayfilt(x, filt, engine = "filter")
  y_fft <- sgolayfilt(x, filt, engine = "fft")
  expect_equal(yref, y_filter)
  expect_equal(yref, y_fft)
})


test_that("sgolayfilt results are all equal with rowwise=FALSE", {
  filt <- sgolay::sgolay(p = 2, n = 5)
  x <- matrix(runif(1000), nrow = 500, ncol = 2)
  yref <- apply(x, 2, function(xin) signal::sgolayfilt(xin, filt))
  y_filter <- sgolayfilt(x, filt, engine = "filter")
  y_fft <- sgolayfilt(x, filt, engine = "fft")
  expect_equal(yref, y_filter)
  expect_equal(yref, y_fft)
})


test_that("sgolayfilt results are all equal with rowwise=TRUE", {
  filt <- sgolay::sgolay(p = 2, n = 5)
  x <- matrix(runif(1000), nrow = 2, ncol = 500)
  yref <- t(apply(x, 1, function(xin) signal::sgolayfilt(xin, filt)))
  y_filter <- sgolayfilt(x, filt, engine = "filter", rowwise = TRUE)
  y_fft <- sgolayfilt(x, filt, engine = "fft", rowwise = TRUE)
  expect_equal(yref, y_filter)
  expect_equal(yref, y_fft)
})
