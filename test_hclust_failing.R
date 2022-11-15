x <- matrix(c(1,2,3,4.5,6,7,8), ncol = 1)
rownames(x) <- sprintf("P%d", seq_len(nrow(x)))
peak2peak_dist <- dist(x)

# M1:
# P1, P2, P4
peak2peak_dist[1] <- 20
peak2peak_dist[3] <- 20
peak2peak_dist[8] <- 20

# M2:
# P3, P5, P6
# 12
peak2peak_dist[13:14] <- 20
peak2peak_dist[19] <- 20

# M3, P7
peak2peak_dist


clust <- stats::hclust(peak2peak_dist, method = "complete")
clust2 <- sgolay:::hclust(peak2peak_dist, method = "complete")

n <- as.integer(attr(peak2peak_dist, "Size"))
hcl1 <- .Fortran(
  stats:::C_hclust,
  n = n,
  len = as.integer(n * (n - 1)/2),
  method = as.integer(3L),
  ia = integer(n),
  ib = integer(n),
  crit = double(n),
  members = as.double(rep(1, n)),
  nn = integer(n),
  disnn = double(n),
  diss = peak2peak_dist
)
hcl2 <- sgolay:::hclust_impl(
  n = n,
  d = peak2peak_dist,
  method = 3L,
  members = rep(1, n)
)

clust$height
clust2$height

plot(clust)
clust$merge
as.data.frame(clust$merge)



