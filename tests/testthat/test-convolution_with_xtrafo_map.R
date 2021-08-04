set.seed(98)

convmap_params <- list(
  maptype = "convolution_with_xtrafo_map",
  mapname = "convolution_testmap",
  src_idx = c(4,10,6,9,2),
  tar_idx = c(11, 14, 12),
  src_x = c(1,4,50,27,98),
  tar_x = c(27, 56, 12), # problem: c(28.6, 56, 12), # problem: c(17.8, 56, 12), #c(27, 56, 12),
  shiftx_idx = 17,
  scalex_idx = 15,
  winsize_idx = 13
)


test_that("propagate of convolution_with_xtrafo map with energy transformation coincides with convolution with manual energy trafo", {
  inp <- runif(20, min=1, max=10)
  inp[convmap_params$shiftx_idx] <- 1.5
  inp[convmap_params$scalex_idx] <- 0.1
  inp[convmap_params$winsize_idx] <- 12
  cur_pars <- convmap_params
  cur_pars2 <- convmap_params
  cur_pars2[["shiftx_idx"]] <- NULL
  cur_pars2[["shiftx"]] <- 0
  cur_pars2[["scalex_idx"]] <- NULL
  cur_pars2[["scalex"]] <- 0
  cur_pars2[["tar_x"]] <- inp[cur_pars[["shiftx_idx"]]] +
    (1+inp[cur_pars[["scalex_idx"]]]) * cur_pars[["tar_x"]]
  convmap1 <- create_convolution_with_xtrafo_map()
  convmap1$setup(cur_pars)
  convmap2 <- create_convolution_with_xtrafo_map()
  convmap2$setup(cur_pars2)
  res1 <- convmap1$propagate(inp)
  res2 <- convmap2$propagate(inp)
  expect_equal(res1, res2)
})


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
  set.seed(30)
  cur_pars <- convmap_params
  inp <- runif(20, min=1, max=5)
  inp[convmap_params$shiftx_idx] <- 1.5
  inp[convmap_params$scalex_idx] <- 0.1 # 0.1
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

