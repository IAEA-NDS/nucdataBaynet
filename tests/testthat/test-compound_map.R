normerr_params <- list(
  maptype = "normerr_map",
  src_idx = 1,
  tar_idx = 6:10,
  src_feat = 1L,
  tar_feat = rep(1, 5)
)


linmod_params <- list(
  maptype = "linmod_map",
  src_idx = 2:5,
  tar_idx = 6:10,
  src_x = 11:14,
  tar_x = c(12, 12.5, 13, 13.5, 14)
)


compound_params <- list(
  maptype = "compound_map",
  maps = list(
    normerr_params,
    linmod_params
  )
)


test_that("compound_map propagate without id matrix works correctly", {
  map <- create_compound_map()
  map$setup(compound_params)
  inp <- c(5, 21:24, rep(7,5))
  res <- map$propagate(inp, with.id = FALSE)
  exp_out <- c(rep(0, 5), 5 + c(22, 22.5, 23, 23.5, 24))
  expect_equal(res, exp_out)
})


test_that("compound_map propagate with id matrix works correctly", {
  map <- create_compound_map()
  map$setup(compound_params)
  inp <- c(5, 21:24, rep(7,5))
  res <- map$propagate(inp, with.id = TRUE)
  exp_out <- c(5, 21:24, 5 + c(22, 22.5, 23, 23.5, 24) + 7)
  expect_equal(res, exp_out)
})


test_that("sensitivity matrix without id matrix correctly computed", {
  linmap <- create_linmod_map()
  linmap$setup(linmod_params)
  normerrmap <- create_normerr_map()
  normerrmap$setup(normerr_params)
  compmap <- create_compound_map()
  compmap$setup(compound_params)
  inp <- 1:10
  compS <- compmap$jacobian(inp, with.id=FALSE)
  linS <- linmap$jacobian(inp, with.id=FALSE)
  normerrS <- normerrmap$jacobian(inp, with.id=FALSE)
  # check the linmap submatrix
  S1 <- compS[linmap$get_tar_idx(), linmap$get_src_idx()]
  S2 <- linS[linmap$get_tar_idx(), linmap$get_src_idx()]
  expect_equal(S1, S2)
  # check the normerr submatrix
  S1 <- compS[normerrmap$get_tar_idx(), normerrmap$get_src_idx(), drop=FALSE]
  S2 <- normerrS[normerrmap$get_tar_idx(), normerrmap$get_src_idx(), drop=FALSE]
  expect_equal(S1, S2)
  # check if required blocks are zero
  expect_true(all(compS[compmap$get_src_idx(), compmap$get_src_idx()] == 0))
  expect_true(all(compS[compmap$get_tar_idx(), compmap$get_tar_idx()] == 0))
})


test_that("sensitivity matrix with id matrix correctly computed", {
  linmap <- create_linmod_map()
  linmap$setup(linmod_params)
  normerrmap <- create_normerr_map()
  normerrmap$setup(normerr_params)
  compmap <- create_compound_map()
  compmap$setup(compound_params)
  inp <- 1:10
  compS <- compmap$jacobian(inp, with.id=TRUE)
  linS <- linmap$jacobian(inp, with.id=FALSE)
  normerrS <- normerrmap$jacobian(inp, with.id=FALSE)
  # check the linmap submatrix
  S1 <- compS[linmap$get_tar_idx(), linmap$get_src_idx()]
  S2 <- linS[linmap$get_tar_idx(), linmap$get_src_idx()]
  expect_equal(S1, S2)
  # check the normerr submatrix
  S1 <- compS[normerrmap$get_tar_idx(), normerrmap$get_src_idx(), drop=FALSE]
  S2 <- normerrS[normerrmap$get_tar_idx(), normerrmap$get_src_idx(), drop=FALSE]
  expect_equal(S1, S2)
  # check if required blocks are zero
  diag_idcs <- (0:(ncol(compS)-1))*nrow(compS) + seq_len(ncol(compS))
  expect_true(all(compS[diag_idcs] == 1))
})

# test nested relationships

linmod_params2 <- list(
  maptype = "linmod_map",
  src_idx = 6:10,
  tar_idx = 11:15,
  src_x = c(12, 12.5, 13, 13.5, 14),
  tar_x = c(12, 12.5, 13, 13.5, 14)
)


compound_params2 <- list(
  maptype = "compound_map",
  maps = list(
    linmod_params2,
    linmod_params,
    normerr_params
  )
)


test_that("nested compound map does propagate correctly without id contribution", {
  compmap <- create_compound_map()
  compmap$setup(compound_params2)
  inp <- 1:15
  res <- compmap$propagate(inp, with.id=FALSE)
  linmod_res <- c(4, 4.5, 5, 5.5, 6)
  expres <- c(rep(0, 5), linmod_res, linmod_res)
  expect_equal(res, expres)
})


test_that("nested compound map jacobian without id contribution calculated correctly", {
  compmap <- create_compound_map()
  compmap$setup(compound_params2)
  inp <- 1:15
  S <- compmap$jacobian(inp, with.id=FALSE)

  expect_equal(drop0(S[6:10,1:5]), drop0(S[11:15,1:5]))
  expect_true(all(S[1:15,6:10]==0))
  expect_true(all(S[1:5,1:15])==0)
})


test_that("nested compound map does propagate correctly with id contribution", {
  compmap <- create_compound_map()
  compmap$setup(compound_params2)
  inp <- 1:15
  res <- compmap$propagate(inp, with.id=TRUE)
  linmod_res <- c(4, 4.5, 5, 5.5, 6)
  expres <- c(rep(0, 5), linmod_res, linmod_res)
  expres[1:5] <- inp[1:5]
  expres[6:10] <- expres[6:10] + inp[6:10]
  expres[11:15] <- expres[11:15] + inp[6:10] + inp[11:15]
  expect_equal(res, expres)
})


test_that("nested compound map jacobian with id contribution calculated correctly", {
  compmap <- create_compound_map()
  compmap$setup(compound_params2)
  inp <- 1:15
  S <- compmap$jacobian(inp, with.id=TRUE)
  expect_equal(drop0(S[6:10,1:5]), drop0(S[11:15,1:5]))
  S1 <- as.matrix(S[11:15,11:15])
  S2 <- diag(rep(1,5))
  dimnames(S1) <- NULL
  dimnames(S2) <- NULL
  expect_equal(S1, S2)
  S1 <- as.matrix(S[11:15,6:10])
  S2 <- diag(rep(1,5))
  dimnames(S1) <- NULL
  dimnames(S2) <- NULL
  expect_equal(S1, S2)
  expect_true(all(S[1:5,6:15]==0))
  expect_true(all(diag(S)==1))
})


# test nested relationships including selfmaps

nonlinear_params <- list(
  maptype = "nonlinear_map",
  funname = "exp",
  src_idx = 11:15,
  tar_idx = 11:15
)


compound_params3 <- list(
  maptype = "compound_map",
  maps = list(
    nonlinear_params,
    linmod_params2,
    linmod_params,
    normerr_params
  )
)


test_that("nested map with nonlinear selfmap propagates correctly without id contribution", {
  nestedmap <- create_compound_map()
  nestedmap$setup(compound_params3)
  inp <- 1:15
  linmod_res <- c(4, 4.5, 5, 5.5, 6)
  expres <- c(rep(0, 5), linmod_res, linmod_res)
  expres[11:15] <- exp(expres[11:15])
  res <- nestedmap$propagate(inp, with.id=FALSE)
  expect_equal(res, expres)
})


test_that("nested map with nonlinear selfmap propagates correctly with id contribution", {
  nestedmap <- create_compound_map()
  nestedmap$setup(compound_params3)
  inp <- 1:15
  linmod_res <- c(4, 4.5, 5, 5.5, 6)
  expres <- c(rep(0, 5), linmod_res, linmod_res)
  expres[1:5] <- expres[1:5] + inp[1:5]
  expres[6:10] <- expres[6:10] + inp[6:10]
  expres[11:15] <- exp(expres[11:15] + inp[6:10]) + inp[11:15]
  res <- nestedmap$propagate(inp, with.id=TRUE)
  expect_equal(res, expres)
})


test_that("numeric jacobian and analytic jacobian of nested map with nonlinear selfmap and with id contribution coincide", {
  nestedmap <- create_compound_map()
  nestedmap$setup(compound_params3)
  inp <- 1:15
  expres <- as.matrix(jacobian(nestedmap$propagate, inp, with.id=TRUE))
  res <- as.matrix(nestedmap$jacobian(inp, with.id=TRUE))
  dimnames(expres) <- dimnames(res) <- NULL
  expect_equal(res, expres, tolerance=1e-6)
})


test_that("numeric jacobian and analytic jacobian of nested map with nonlinear selfmap and without id contribution coincide", {
  nestedmap <- create_compound_map()
  nestedmap$setup(compound_params3)
  inp <- 1:15
  expres <- as.matrix(jacobian(nestedmap$propagate, inp, with.id=FALSE))
  res <- as.matrix(nestedmap$jacobian(inp, with.id=FALSE))
  dimnames(expres) <- dimnames(res) <- NULL
  expect_equal(res, expres, tolerance=1e-6)
})
