set.seed(98)

convmap_params <- list(
  maptype = "convolution_with_xtrafo_map",
  mapname = "convolution_testmap",
  src_idx = c(4,10,6,9,2),
  tar_idx = c(11, 14, 12),
  src_x = c(1,4,50,27,98),
  tar_x = c(27, 56, 12),
  shiftx_idx = 17,
  scalex_idx = 15,
  winsize_idx = 13
)


test_that("jacobian of convolution_with_xtrafo_map coincides with numerical jacobian of propagate function if no energy derivatives", {
  cur_pars <- convmap_params
  cur_pars[c("shiftx_idx", "scalex_idx", "winsize_idx")] <- NULL
  cur_pars[["winsize"]] <- 10
  cur_pars[["shiftx"]] <- 1
  cur_pars[["scalex"]] <- 0.1
  inp <- runif(20, min=1, max=10)
  inp[convmap_params$shiftx_idx] <- 1.5
  inp[convmap_params$scalex_idx] <- 0.1
  inp[convmap_params$winsize_idx] <- 12
  convmap <- create_convolution_with_xtrafo_map()
  convmap$setup(cur_pars)
  # with.id=TRUE
  expres <- jacobian(convmap$propagate, inp, with.id=TRUE)
  res <- as.matrix(convmap$jacobian(inp, with.id=TRUE))
  dimnames(expres) <- dimnames(res) <- NULL
  expect_equal(res, expres)
  # with.id=FALSE
  expres <- jacobian(convmap$propagate, inp, with.id=FALSE)
  res <- as.matrix(convmap$jacobian(inp, with.id=FALSE))
  dimnames(expres) <- dimnames(res) <- NULL
  expect_equal(res, expres)
})


test_that("jacobian of convolution_with_xtrafo_map coincides with numerical jacobian of propagate function including energy derivatives", {
  cur_pars <- convmap_params
  inp <- runif(20, min=1, max=10)
  inp[convmap_params$shiftx_idx] <- 1.5
  inp[convmap_params$scalex_idx] <- 0.1
  inp[convmap_params$winsize_idx] <- 12
  convmap <- create_convolution_with_xtrafo_map()
  convmap$setup(cur_pars)
  # with.id=FALSE
  expres <- jacobian(convmap$propagate, inp, with.id=FALSE)
  res <- as.matrix(convmap$jacobian(inp, with.id=FALSE))
  dimnames(expres) <- dimnames(res) <- NULL
  expect_equal(res, expres)
  # with.id=TRUE
  expres <- jacobian(convmap$propagate, inp, with.id=TRUE)
  res <- as.matrix(convmap$jacobian(inp, with.id=TRUE))
  dimnames(expres) <- dimnames(res) <- NULL
  expect_equal(res, expres)
})
