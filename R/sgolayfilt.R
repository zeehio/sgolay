filter_sweep <- function(x, filt) {
  len <- length(x)
  n <- length(filt)
  k <- floor(n/2)
  how_many <- len - 2*k
  out <- numeric(length(x))
  for (j in seq_along(filt)) {
    first_idx <- max(1L, j - k)
    last_idx <- min(len, j + k)
    out[first_idx:last_idx] <- out[first_idx:last_idx] + filt[j]*x[first_idx:last_idx]
  }
  out[(k + 1L):(len - k)]
}

sgolayfilt_impl <- function(x, filt, rowwise, return_matrix, engine = c("fft", "filter"), filter_embed = "cfilter") {
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

  coefs_for_first_points <- filt[1:k, ] # k x n
  coefs_for_last_points <- filt[(k + 2L):n, ]

  if (rowwise) {
    out[, 1L:k] <- x[, 1L:n, drop = FALSE] %*% t(coefs_for_first_points)
    out[, (len - k + 1L):len] <- x[, (len - n + 1L):len, drop = FALSE] %*% t(coefs_for_last_points)
  } else {
    out[1L:k, ] <- coefs_for_first_points %*% x[1L:n, , drop = FALSE]
    out[(len - k + 1L):len,] <- coefs_for_last_points %*% x[(len - n + 1L):len, , drop = FALSE]
  }

  if (engine == "fft") {
    conv_coefs <- filt[k + 1L, n:1L]
    fft_length <- length(conv_coefs) + len - 1L
    conv_coefs_padded <- c(rev(conv_coefs), rep(0, len - 1L))
    conj_fft_y_prep <- sgolay:::convolve_prepare(conv_coefs_padded, conj = TRUE)
    x_padded <- numeric(fft_length)
    center_points_idx <- (k + 1L):(len - k)
    for (i in seq_len(num_ser)) {
      if (rowwise) {
        x_padded[length(conv_coefs):fft_length] <- x[i,]
      } else {
        x_padded[length(conv_coefs):fft_length] <- x[,i]
      }
      fft_x_prep <- sgolay:::convolve_prepare(x_padded, conj = FALSE)
      center_points <- sgolay:::convolve_do(
        fft_x_prep,
        conj_fft_y_prep
      )[n:len]
      if (rowwise) {
        out[i, center_points_idx] <- center_points
      } else {
        out[center_points_idx, i] <- center_points
      }
    }
  } else if (engine == "filter") {
    if (filter_embed == "cfilter") {
      xvec <- numeric(n + len - 1L)
      filt_cent <- filt[k + 1L, n:1L]
      offset <- 2L*n - 1L
      center_points_idx <- (k + 1L):(len - k)
      for (i in seq_len(num_ser)) {
        if (rowwise) {
          xvec[n:length(xvec)] <- x[i,]
          out[i, center_points_idx] <- .Call(stats:::C_cfilter, xvec, filt_cent, 1L, FALSE)[offset:(offset + len - 2*k - 1L)]
        } else {
          xvec[n:length(xvec)] <- x[,i]
          potato2 <- .Call(stats:::C_cfilter, xvec, filt_cent, 1L, FALSE)[offset:(offset + len - 2*k - 1L)]
          out[center_points_idx, i] <- potato2
        }
      }
    } else if (filter_embed == "embed") {
      center_points_idx <- (k + 1L):(len - k)
      filt_cent <- filt[k + 1L, n:1L]
      embedded_signal <- matrix(0, nrow = len - 2*k, ncol = n)
      for (i in seq_len(num_ser)) {
        if (rowwise) {
          embedded_signal[] <- stats::embed(x[i,], n)
        } else {
          embedded_signal[] <- stats::embed(x[,i], n)
        }
        cosa2 <- embedded_signal %*% filt_cent
        if (rowwise) {
          out[i, center_points_idx] <- cosa2
        } else {
          out[center_points_idx, i] <- cosa2
        }
      }
    } else if (filter_embed == "slider") {
      center_points_idx <- (k + 1L):(len - k)
      filt_cent <- filt[k + 1L, n:1L]
      for (i in seq_len(num_ser)) {
        if (rowwise) {
          out[i, center_points_idx] <- slider::slide_dbl(x[i,], ~ c(. %*% filt_cent), .before = k, .after = k, .complete = TRUE)[(k + 1L):(len - k)]
        } else {
          cosa2 <- slider::slide_dbl(x[,i], ~ sum(. * filt_cent), .before = k, .after = k, .complete = TRUE)
          out[center_points_idx, i] <- cosa2[(k + 1L):(len - k)]
        }
      }
    } else if (filter_embed == "filter_sweep") {
      center_points_idx <- (k + 1L):(len - k)
      filt_cent <- filt[k + 1L, n:1L]
      for (i in seq_len(num_ser)) {
        if (rowwise) {
          out[i, center_points_idx] <- filter_sweep(x[i,], filt_cent)
        } else {
          out[center_points_idx, i] <- filter_sweep(x[,i], filt_cent)
        }
      }
    } else {
      stop("Wrong filter_embed method")
    }
  } else {
    stop("Wrong engine")
  }
  if (!return_matrix) {
    attr(out, "dim") <- NULL
  }
  out
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
    if (is.matrix(x)) {
      if (filter_length > 100) {
        engine <- "fft"
      } else {
        engine <- "filter"
      }
    }
  }
  if (engine == "fft" && anyNA(x)) {
    if (orig_engine == "fft") {
      rlang::warn(
        message = c(
          'Switching sgolayfilt engine from "fft" to "filter"',
          "!" = "The fft engine does not handle missing values.",
          "i" = "Using engine = 'filter' instead."
        )
      )
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
  return_matrix <- TRUE
  if (!is.matrix(x)) {
    x <- matrix(x, ncol = 1)
    return_matrix <- FALSE
    rowwise <- FALSE
  }
  engine <- choose_engine(x = x, filter_length = nrow(filt), orig_engine = engine)
  out <- sgolayfilt_impl(x, filt, rowwise, return_matrix, engine = engine)
  out
}
