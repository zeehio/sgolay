


test_that("sgolayfilt results are all equal", {
  filt <- sgolay::sgolay(p=2, n=5)
  x <- runif(500)
  yref <- signal::sgolayfilt(x, filt)
  y_filter <- sgolayfilt(x, filt, engine = "filter")
  y_fft <- sgolayfilt(x, filt, engine = "fft")
  expect_equal(yref, y_filter)
  expect_equal(yref, y_fft)
})


test_that("sgolayfilt results are all equal", {
  filt <- sgolay::sgolay(p=2, n=5)
  x <- matrix(runif(1000), nrow=500, ncol = 2)
  yref <- apply(x, 2, function(xin) signal::sgolayfilt(xin, filt))
  y_filter <- sgolayfilt(x, filt, engine = "filter")
  y_fft <- sgolayfilt(x, filt, engine = "fft")
  expect_equal(yref, y_filter)
  expect_equal(yref, y_fft)
})
