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
  in_several_segments <- min_idx <= max_idx
  # save the original specifications
  orig_min_idx <- min_idx
  orig_max_idx <- max_idx
  orig_loc_tar_idx <- loc_tar_idx
  orig_tar_x_min <- tar_x_min
  orig_tar_x_max <- tar_x_max
  # remove the target points in one segment and treat them specially
  min_idx <- orig_min_idx[in_several_segments]
  max_idx <- orig_max_idx[in_several_segments]
  loc_tar_idx <- orig_loc_tar_idx[in_several_segments]
  tar_x_min <- tar_x_min[in_several_segments]
  tar_x_max <- tar_x_max[in_several_segments]
  # for all segments of src_x completely inside apply simpson integration rule
  xdiff <- src_x[-1] - head(src_x, n=-1)
  coeff <- xdiff/2
  imi <- sequence(max_idx-min_idx, loc_tar_idx, by=0L)
  jmi <- sequence(max_idx-min_idx, min_idx, by=1)
  x <- (src_x[jmi+1] - src_x[jmi])/2
  # divide to go from area to average value
  x <- x / (tar_x_max - tar_x_min)
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
  c1lo <- c1lo / (tar_x_max - tar_x_min)
  c2lo <- c2lo / (tar_x_max - tar_x_min)
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
  c1hi <- c1hi / (tar_x_max - tar_x_min)
  c2hi <- c2hi / (tar_x_max - tar_x_min)
  # now calculate the coefficients for convolutions being completely in one segment of src_x
  tar_x_min <- orig_tar_x_min[!in_several_segments]
  tar_x_max <- orig_tar_x_max[!in_several_segments]
  ione <- orig_loc_tar_idx[!in_several_segments]
  idx <- orig_max_idx[!in_several_segments]
  dmesh <- src_x[idx+1] - src_x[idx]
  # coeffs due to linear interpolation
  c1one <- ((src_x[idx+1] - tar_x_max) + (src_x[idx+1] - tar_x_min)) / dmesh
  c2one <- ((tar_x_max - src_x[idx]) + (tar_x_min - src_x[idx])) / dmesh
  # integration factors in trapezoidal rule
  c1one <- c1one * (tar_x_max - tar_x_min)/2
  c2one <- c2one * (tar_x_max - tar_x_min)/2
  # apply multiplication to get average value
  c1one <- c1one / (tar_x_max - tar_x_min)
  c2one <- c2one / (tar_x_max - tar_x_min)

  # put everything together
  convmat <- sparseMatrix(
    i = tar_idx[c(rep(imi, 2),
                  rep(ilo, 2),
                  rep(ihi, 2),
                  rep(ione, 2)
                  )],
    j = src_idx[c(jmi, jmi+1,
                  jlo-1, jlo,
                  jhi, jhi+1,
                  idx, idx+1)],
    x = c(rep(x, 2),
          c1lo, c2lo,
          c1hi, c2hi,
          c1one, c2one),
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
  return(list(flo=f1, fhi=f2))
}
