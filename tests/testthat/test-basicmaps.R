test_that("get convolute_rect_matrix yields correct result", {
  set.seed(27)
  src_x <- cumsum(runif(100, min=1, max=5))
  tar_x <- 50
  winsize <- 30
  tar_x_min <- tar_x - winsize / 2
  tar_x_max <- tar_x + winsize / 2
  vals <- seq_along(src_x)
  S <- get_convolute_rect_matrix(src_x, tar_x, winsize)
  aug_src_x <- sort(c(src_x, tar_x-winsize/2, tar_x+winsize/2))
  aug_vals <- approx(src_x, vals, xout=aug_src_x)$y
  xdiff <- aug_src_x[-1] - head(aug_src_x, n=-1)
  c1 <- head(aug_vals, n=-1) * xdiff/2
  c2 <- aug_vals[-1] * xdiff/2
  area <- c1 + c2
  cumarea <- c(0, cumsum(area))
  idcs <- findInterval(tar_x+c(-winsize,winsize)/2, aug_src_x)
  expres <- diff(cumarea[idcs]) / (tar_x_max - tar_x_min)
  res <- as.vector(S %*% vals)
  expect_equal(res, expres)
  # also check with augmented grid
  S <- get_convolute_rect_matrix(aug_src_x, tar_x, winsize)
  res <- as.vector(S %*% aug_vals)
  expect_equal(res, expres)
})


test_that("get_convolute_rect_matrix yields correct result if no segment completely inside", {
  winsize <- 30
  src_x <- c(1, 55, 100)
  tar_x <- 50
  tar_x_min <- tar_x - winsize / 2
  tar_x_max <- tar_x + winsize / 2
  vals <- seq_along(src_x)
  S <- get_convolute_rect_matrix(src_x, tar_x, winsize)
  aug_src_x <- sort(c(src_x, tar_x-winsize/2, tar_x+winsize/2))
  aug_vals <- approx(src_x, vals, xout=aug_src_x)$y
  xdiff <- aug_src_x[-1] - head(aug_src_x, n=-1)
  c1 <- head(aug_vals, n=-1) * xdiff/2
  c2 <- aug_vals[-1] * xdiff/2
  area <- c1 + c2
  cumarea <- c(0, cumsum(area))
  idcs <- findInterval(tar_x+c(-winsize,winsize)/2, aug_src_x)
  expres <- diff(cumarea[idcs]) / (tar_x_max - tar_x_min)
  res <- as.vector(S %*% vals)
  expect_equal(res, expres)
})


test_that("get_convolute_rect_matrix yields correct result if window completely in one segment", {
  winsize <- 50
  src_x <- c(1, 100)
  tar_x <- 50
  tar_x_min <- tar_x - winsize / 2
  tar_x_max <- tar_x + winsize / 2
  vals <- c(1, 100)
  S <- get_convolute_rect_matrix(src_x, tar_x, winsize)
  aug_src_x <- sort(c(src_x, tar_x-winsize/2, tar_x+winsize/2))
  aug_vals <- approx(src_x, vals, xout=aug_src_x)$y
  xdiff <- aug_src_x[-1] - head(aug_src_x, n=-1)
  c1 <- head(aug_vals, n=-1) * xdiff/2
  c2 <- aug_vals[-1] * xdiff/2
  area <- c1 + c2
  cumarea <- c(0, cumsum(area))
  idcs <- findInterval(tar_x+c(-winsize,winsize)/2, aug_src_x)
  expres <- diff(cumarea[idcs]) / (tar_x_max - tar_x_min)
  res <- as.vector(S %*% vals)
  expect_equal(res, expres)
})


test_that("compound convolution map for several points equals sum of individual convolution maps", {
  set.seed(27)
  src_x <- cumsum(runif(30, min=3, max=5))
  tar_x <- c(30, 35, 40)
  winsize <- 20
  Scomp <- get_convolute_rect_matrix(src_x, tar_x, winsize)
  Sadd <- NULL
  for (cur_tar_x in tar_x) {
    Sadd <- rbind(Sadd, get_convolute_rect_matrix(src_x, cur_tar_x, winsize))
  }
  expect_equal(Scomp, Sadd)
})


test_that("convolute_rect_derivative conincides with numerical jacobian", {
  src_x <- c(1, 55, 100)
  tar_x <- 50
  winsize <- 30
  tar_x_min <- tar_x - winsize/2
  tar_x_max <- tar_x + winsize/2
  vals <- seq_along(src_x)
  fun <- function(x) {
    as.vector(get_convolute_rect_matrix(src_x, tar_x, winsize=x) %*% vals)
  }
  Fval <- fun(winsize)
  expres <- as.vector(jacobian(fun, winsize))
  tmpres <- get_convolute_rect_derivative(src_x, vals, tar_x, winsize)
  res <- (tmpres$flo + tmpres$fhi)/2 / (tar_x_max - tar_x_min)
  res <- res - Fval / (tar_x_max - tar_x_min)
  expect_equal(res, expres)
})
