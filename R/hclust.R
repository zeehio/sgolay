#  File src/library/stats/R/hclust.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2020 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

## Hierarchical clustering, on raw input data; we will use Euclidean
## distance.  A range of criteria are supported; also there is a
## storage-economic option.
##
## We use the general routine, `hc', which caters for 7 criteria,
## using a half dissimilarity matrix; (BTW, this uses the very efficient
## nearest neighbor chain algorithm, which makes this algorithm of
## O(n^2) computational time, and differentiates it from the less
## efficient -- i.e. O(n^3) -- implementations in all commercial
## statistical packages -- as far as I am aware -- except Clustan.)
##
## Clustering Methods:
##
## 1. Ward's minimum variance or error sum of squares method (using raw d) -> "ward.D"
## 2. single linkage or nearest neighbor method.
## 3. complete linkage or diameter.
## 4. average linkage, group average, or UPGMA method.
## 5. McQuitty's or WPGMA method.
## 6. median, Gower's or WPGMC method.
## 7. centroid or UPGMC method (7).
## 8. Ward's ... "correct" method using d^2 (in Fortran) -> "ward.D2"
##
## Original author: F. Murtagh, May 1992
## R Modifications: Ross Ihaka, Dec 1996
##		    Friedrich Leisch, Apr 1998, Jun 2000
## "ward.D" and "ward.D2" from suggestions by Pierre Legendre,
## by Martin Maechler, mostly in the Fortran part.

#' Hierarchical clustering
#'
#' Experiments on hierarchical clustering
#'
#' @inheritParams stats::hclust
#' @seealso stats::hclust
hclust <- function(d, method="complete", members=NULL)
{
  ## order of METHODS --> i.meth -> Fortran's  iOpt  codes
  METHODS <- c("ward.D", "single", # 1, 2,
               "complete", "average", "mcquitty", # 3, 4, 5,
               "median", "centroid", "ward.D2") # 6, 7, 8
  if(method == "ward") { # do not deprecate earlier than 2015!
    message("The \"ward\" method has been renamed to \"ward.D\"; note new \"ward.D2\"")
    method <- "ward.D"
  }
  i.meth <-  pmatch(method, METHODS)
  if(is.na(i.meth))
    ## TODO: use gettextf() [-> translation string change]
    stop("invalid clustering method", paste("", method))
  if(i.meth == -1)
    stop("ambiguous clustering method", paste("", method))

  n <- as.integer(attr(d, "Size"))
  if (is.null(n))
    stop("invalid dissimilarities")
  if (is.na(n) || n > 65536L) stop("size cannot be NA nor exceed 65536")
  if (n < 2)
    stop("must have n >= 2 objects to cluster")
  len <- as.integer(n*(n-1)/2)
  if (length(d) != len)
    (if (length(d) < len) stop else warning
    )("dissimilarities of improper length")


  storage.mode(d) <- "double"
  hcl <- hclust_impl(n = n, d = d, method = i.meth, members = members)

  ## 2nd step: interpret the information that we now have
  ## as merge, height, and order lists.

  hcass <- hcass_impl(n, hcl$ia, hcl$ib)
  structure(
    list(
      merge = cbind(hcass$iia[1L:(n - 1)], hcass$iib[1L:(n - 1)]),
      height = hcl$crit[1L:(n - 1)],
      order = hcass$order,
      labels = attr(d, "Labels"),
      method = METHODS[i.meth],
      call = match.call(),
      dist.method = attr(d, "method")
    ),
    class = "hclust"
  )
}

hcass_impl <- function(n, ia, ib) {
  out <- .Call(
    c_hcass,
    as.integer(n),
    as.integer(ia),
    as.integer(ib)
  )
  names(out) <- c("order", "iia", "iib")
  out
}


hclust_impl <- function(n, d, method, members) {

  if (is.null(members))
    members <- rep(1., n)
  else if (length(members) != n)
    stop("invalid length of members")

  out <- .Call(
    c_hclust,
    as.integer(n),
    as.double(d),
    as.integer(method),
    as.double(members)
  )
  names(out) <- c("ia", "ib", "crit")
  out
}

