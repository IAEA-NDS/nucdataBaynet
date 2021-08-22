set.seed(77)

lin_params <- list(
  maptype = "linearinterpol_map",
  src_idx = 1:4,
  tar_idx = 6:10,
  src_x = 11:14,
  tar_x = c(12,12.5,13,13.5,14)
)

relmap_params <- list(
  maptype = "relativemap_map",
  basemap = lin_params
)


test_that("propagation works as expected", {
  linmap <- create_map(lin_params)
  relmap <- create_map(relmap_params)
  x <- runif(10)
  # with.id = FALSE
  res <- relmap$propagate(x, with.id=FALSE)
  expres <- linmap$propagate(x, with.id=FALSE)
  expres[linmap$get_tar_idx()] <- expres[linmap$get_tar_idx()] * x[linmap$get_tar_idx()]
  expect_equal(res, expres)
  # with.id = TRUE
  res <- relmap$propagate(x, with.id=TRUE)
  expres <- linmap$propagate(x, with.id=TRUE)
  expres[linmap$get_tar_idx()] <- expres[linmap$get_tar_idx()] * x[linmap$get_tar_idx()]
  expres[linmap$get_tar_idx()] <- expres[linmap$get_tar_idx()] + x[linmap$get_tar_idx()]
  expect_equal(res, expres)
})


test_that("analytic and numeric jacobian of relativemap are equal", {
  relmap <- create_map(relmap_params)
  x <- runif(10)
  # with.id = FALSE
  expres <- jacobian(relmap$propagate, x, with.id=FALSE)
  res <- as.matrix(relmap$jacobian(x, with.id=FALSE))
  dimnames(res) <- dimnames(expres) <- NULL
  expect_equal(res, expres)
  # with.id = TRUE
  expres <- jacobian(relmap$propagate, x, with.id=FALSE)
  diag(expres) <- diag(expres) + 1
  res <- as.matrix(relmap$jacobian(x, with.id=TRUE))
  dimnames(res) <- dimnames(expres) <- NULL
  expect_equal(res, expres)
})
