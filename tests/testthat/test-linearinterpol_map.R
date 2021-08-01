params <- list(
  maptype = "linearinterpol_map",
  src_idx = 1:4,
  tar_idx = 6:10,
  src_x = 11:14,
  tar_x = c(12,12.5,13,13.5,14)
)


test_that("linearinterpol_map propagate without id matrix works correctly", {
  map <- create_linearinterpol_map()
  map$setup(params)
  res <- map$propagate(c(21:24, rep(0,6)), with.id = FALSE)
  expect_equal(res[6:10], c(22, 22.5, 23, 23.5, 24))
  expect_equal(res[1:5], rep(0, 5))
})


test_that("linearinterpol_map propagate with id matrix works correctly", {
  map <- create_linearinterpol_map()
  map$setup(params)
  res <- map$propagate(c(21:24, rep(0,6)), with.id = TRUE)
  expect_equal(res[6:10], c(22, 22.5, 23, 23.5, 24))
  expect_equal(res[1:5], c(21:24, 0))
})


test_that("expansion point of linearinterpol_map does not impact Jacobian matrix", {
  map <- create_linearinterpol_map()
  map$setup(params)
  r1 <- 1:15
  r2 <- 31:45
  S1 <- map$jacobian(r1, with.id=FALSE)
  S2 <- map$jacobian(r2, with.id=FALSE)
  expect_equal(S1, S2)
  S1 <- map$jacobian(r1, with.id=TRUE)
  S2 <- map$jacobian(r2, with.id=TRUE)
  expect_equal(S1, S2)
})
