linint_params <- list(
  maptype = "linearinterpol_map",
  src_idx = 1:4,
  tar_idx = 6:10,
  src_x = c(11, 13, 15, 20),
  tar_x = c(12, 12.5, 13, 13, 13.3)
)

linint_xtrafo_params <- modifyList(linint_params, list(
  maptype = "linearinterpol_with_xtrafo_map",
  shiftx_idx = 11,
  scalex_idx = 12
))


test_that("propagate of linearinterpol_with_xtrafo_map without scale/shift yields identical result as linearinterpol_map", {
  linint_map <- create_linearinterpol_map()
  linint_map$setup(linint_params)
  linint_xtrafo_map <- create_linearinterpol_with_xtrafo_map()
  linint_xtrafo_map$setup(linint_xtrafo_params)
  linint_inp <- 1:10
  linint_xtrafo_inp <- c(linint_inp, c(0, 1))
  # with.id=FALSE
  expres <- linint_map$propagate(linint_inp, with.id=FALSE)
  res <- linint_xtrafo_map$propagate(linint_xtrafo_inp, with.id=FALSE)[1:10]
  expect_equal(res, expres)
  # with.id=TRUE
  expres <- linint_map$propagate(linint_inp, with.id=TRUE)
  res <- linint_xtrafo_map$propagate(linint_xtrafo_inp, with.id=TRUE)[1:10]
  expect_equal(res, expres)
})


test_that("propagate of linearinptol_with_xtrafo_map with shift and scale identical to manually shifted/scaled mesh of linearinterpol_map", {
  shiftx <- 0.1
  scalex <- 1.05
  mod_linint_params <- linint_params
  mod_linint_params$tar_x <- shiftx + linint_params$tar_x * scalex
  mod_linint_map <- create_linearinterpol_map()
  mod_linint_map$setup(mod_linint_params)
  linint_xtrafo_map <- create_linearinterpol_with_xtrafo_map()
  linint_xtrafo_map$setup(linint_xtrafo_params)
  mod_linint_inp <- 1:10
  linint_xtrafo_inp <- c(mod_linint_inp, c(shiftx, scalex))
  # with.id=FALSE
  expres <- mod_linint_map$propagate(mod_linint_inp, with.id=FALSE)
  res <- linint_xtrafo_map$propagate(linint_xtrafo_inp, with.id=FALSE)[1:10]
  expect_equal(res, expres)
  # with.id=TRUE
  expres <- mod_linint_map$propagate(mod_linint_inp, with.id=TRUE)
  res <- linint_xtrafo_map$propagate(linint_xtrafo_inp, with.id=TRUE)[1:10]
  expect_equal(res, expres)
})


test_that("jacobian of linearinterpol_with_xtrafo_map coincides with numerical jacobian of propagate function", {
  shiftx <- 0.1
  scalex <- 1.05
  inp <- c(1:10, shiftx, scalex)
  linint_xtrafo_map <- create_linearinterpol_with_xtrafo_map()
  linint_xtrafo_map$setup(linint_xtrafo_params)
  linint_xtrafo_map$propagate(inp)
  # with.id=TRUE
  expres <- jacobian(linint_xtrafo_map$propagate, inp, with.id=TRUE)
  res <- as.matrix(linint_xtrafo_map$jacobian(inp), with.id=TRUE)
  dimnames(res) <- dimnames(expres) <- NULL
  expect_equal(res, expres)
  # with.id=FALSE
  expres <- jacobian(linint_xtrafo_map$propagate, inp, with.id=FALSE)
  res <- as.matrix(linint_xtrafo_map$jacobian(inp, with.id=FALSE))
  dimnames(res) <- dimnames(expres) <- NULL
  expect_equal(res, expres)
})
