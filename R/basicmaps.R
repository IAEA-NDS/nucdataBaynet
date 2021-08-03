get_convolute_rect_matrix <- function(src_x, tar_x, winsize,
                                      src_idx=seq_along(src_x), tar_idx=seq_along(tar_x),
                                      dims = c(length(tar_x), length(src_x))) {
  stopifnot(!is.unsorted(src_x))
  tar_x_min <- tar_x - winsize/2
  tar_x_max <- tar_x + winsize/2
  if (any(tar_x_min < src_x[1]))
    stop("some left limit of a bin is outside the mesh given by src_x")
  if (any(tar_x_max > tail(src_x, n=1)))
    stop("some right limit of a bin is outside the mesh given by src_x")
  loc_tar_idx <- seq_along(tar_x)
  min_idx <- findInterval(tar_x_min, src_x) + 1
  max_idx <- findInterval(tar_x_max, src_x)
  # for all segments of src_x completely inside apply simpson integration rule
  xdiff <- src_x[-1] - head(src_x, n=-1)
  coeff <- xdiff/2
  imi <- sequence(max_idx-min_idx, loc_tar_idx, by=0L)
  jmi <- sequence(max_idx-min_idx, min_idx, by=1)
  x <- (src_x[jmi+1] - src_x[jmi])/2
  # edge correction lower bound
  xmeshdiff <- src_x[min_idx] - src_x[min_idx-1]
  xrealdiff <- src_x[min_idx] - tar_x_min
  ilo <- loc_tar_idx
  jlo <- min_idx
  # coefficients for linear interpolation to get funval at tar_x_min
  c1lo <- (src_x[min_idx] - tar_x_min) / xmeshdiff
  c2lo <- (tar_x_min - src_x[min_idx-1]) / xmeshdiff
  # multiplication to get integral
  c1lo <- c1lo * xrealdiff/2
  c2lo <- (1+c2lo) * xrealdiff/2
  # edge correction upper bound
  xmeshdiff <- src_x[max_idx+1] - src_x[max_idx]
  xrealdiff <- tar_x_max - src_x[max_idx]
  ihi <- loc_tar_idx
  jhi <- max_idx
  # coefficients for linear interpolation to get funval at tar_x_max
  c1hi <- (src_x[max_idx+1] - tar_x_max) / xmeshdiff
  c2hi <- (tar_x_max - src_x[max_idx]) / xmeshdiff
  # multiplication to get integral
  c1hi <- (1+c1hi) * xrealdiff/2
  c2hi <- c2hi * xrealdiff/2
  # put everything together
  convmat <- sparseMatrix(i = tar_idx[c(rep(imi, 2), rep(ilo, 2), rep(ihi, 2))],
               j = src_idx[c(jmi, jmi+1, jlo-1, jlo, jhi, jhi+1)],
               x = c(rep(x, 2), c1lo, c2lo, c1hi, c2hi),
               dims=dims)
  return(convmat)
}


get_convolute_rect_derivative <- function(src_x, src_y, tar_x, winsize) {
  stopifnot(!is.unsorted(src_x))
  tar_x_min <- tar_x - winsize/2
  tar_x_max <- tar_x + winsize/2
  tar_idx <- seq_along(tar_x)
  min_idx <- findInterval(tar_x_min, src_x, rightmost.closed=TRUE)
  max_idx <- findInterval(tar_x_max, src_x, rightmost.closed=TRUE)
  # get the function value
  f1 <- (src_y[min_idx] * (src_x[min_idx+1] - tar_x_min) +
           src_y[min_idx+1] * (tar_x_min - src_x[min_idx])) / (src_x[min_idx+1] - src_x[min_idx])
  f2 <- (src_y[max_idx] * (src_x[max_idx+1] - tar_x_max) +
         src_y[max_idx+1] * (tar_x_max - src_x[max_idx])) / (src_x[max_idx+1] - src_x[max_idx])
  deriv <- (f1 + f2) / 2
  return(deriv)
}
