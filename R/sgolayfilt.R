# filter_sweep <- function(x, filt) {
#   # FIXME: It does not do what I want it to do. Probably I still don't know what I want
#   len <- length(x)
#   n <- length(filt)
#   k <- floor(n/2)
#   #lengtout <- len - 2*k
#   out <- numeric(length(x))
#   for (j in seq_along(filt)) {
#     first_idx <- max(1L, j - k)
#     last_idx <- min(len, j + k)
#     out[first_idx:last_idx] <- out[first_idx:last_idx] + filt[j]*x[first_idx:last_idx]
#   }
#   out[(k + 1L):(len - k)]
# }

convolve_prepare <- function(x, conj = FALSE, plan = NULL, impl = "auto") {
  if (conj) {
    do_conj <- Conj
  } else {
    do_conj <- identity
  }
  do_conj(stats::fft(x))
}

convolve_do <- function(fft_x, conj_fft_y) {
  Re(stats::fft(fft_x * conj_fft_y, inverse = TRUE))/length(fft_x)
}


#' @importFrom signal sgolay
#' @export
signal::sgolay


choose_engine <- function(x, filter_length, orig_engine) {
  engine <- orig_engine
  if (engine == "filter") {
    return(engine)
  }
  if (engine == "auto") {
    if (filter_length > 29) {
      engine <- "fft"
    } else {
      engine <- "filter"
    }
  }
  if (engine == "fft" && anyNA(x)) {
    if (orig_engine == "fft") {
      warning('Switching sgolayfilt engine from "fft" to "filter". The fft engine does not handle missing values')
    }
    engine <- "filter"
  }
  engine
}

#' Savitzky-Golay filtering
#'
#' @param x A numeric matrix or vector
#' @inheritParams signal::sgolayfilt
#' @param rowwise If `TRUE`, Apply the filter by rows instead of by columns
#' @param engine "fft" Uses the Fast Fourier Transform to apply the filter. "filter" uses a convolution. Both give
#' the same results, fft is usually faster on larger filter lengths, and larger matrices.
#'
#' @return A matrix or vector of the same dimensions or length as `x`, with the result of the filter
#' @export
#'
#' @examples
#' x <- runif(300)
#' y <- sgolayfilt(x, p=2, n = 21)
sgolayfilt <- function(x, p = 3, n = p + 3 - p %% 2, m = 0, ts = 1, rowwise = FALSE,
                       engine = c("auto", "fft", "filter")) {
  engine <- match.arg(engine)
  if (inherits(p, "sgolayFilter") || (!is.null(dim(p)) && dim(p) > 1)) {
    filt <- p
  } else {
    filt <- sgolay(p, n, m, ts)
  }
  mode(x) <- "double"
  return_matrix <- TRUE
  if (!is.matrix(x)) {
    x <- matrix(x, ncol = 1)
    return_matrix <- FALSE
    rowwise <- FALSE
  }
  engine <- choose_engine(x = x, filter_length = nrow(filt), orig_engine = engine)
  if (rowwise) {
    num_ser <- nrow(x)
    len <- ncol(x)
  } else {
    num_ser <- ncol(x)
    len <- nrow(x)
  }

  n <- nrow(filt)
  k <- floor(n/2)
  out <- matrix(0, nrow = nrow(x), ncol = ncol(x))

  coefs_for_first_points <- filt[1:k, , drop = FALSE]
  coefs_for_last_points <- filt[(k + 2L):n, , drop = FALSE]

  if (rowwise) {
    out[, 1L:k] <- x[, 1L:n, drop = FALSE] %*% t(coefs_for_first_points)
    out[, (len - k + 1L):len] <- x[, (len - n + 1L):len, drop = FALSE] %*% t(coefs_for_last_points)
  } else {
    out[1L:k, ] <- coefs_for_first_points %*% x[1L:n, , drop = FALSE]
    out[(len - k + 1L):len,] <- coefs_for_last_points %*% x[(len - n + 1L):len, , drop = FALSE]
  }

  if (engine == "fft") {
    conv_coefs <- filt[k + 1L, n:1L]
    center_points_idx <- (k + 1L):(len - k)
    if (rowwise) {
      out[,center_points_idx] <- t(convolve_circular(t(x), conv_coefs)[n:len,])
    } else {
      out[center_points_idx,] <- convolve_circular(x, conv_coefs)[n:len,]
    }
  } else if (engine == "filter") {
    xvec <- numeric(n + len - 1L)
    filt_cent <- filt[k + 1L, n:1L]
    offset <- 2L*n - 1L
    center_points_idx <- (k + 1L):(len - k)
    for (i in seq_len(num_ser)) {
      if (rowwise) {
        xvec[n:length(xvec)] <- x[i,]
        out[i, center_points_idx] <- filter(xvec, filt_cent)[offset:(offset + len - 2*k - 1L)]
      } else {
        xvec[n:length(xvec)] <- x[,i]
        out[center_points_idx, i] <- filter(xvec, filt_cent)[offset:(offset + len - 2*k - 1L)]
      }
    }
  } else {
    stop("Wrong engine. Use fft or filter")
  }
  if (!return_matrix) {
    attr(out, "dim") <- NULL
  }
  out
}
