params <- list(
  maptype = "nonlinear_map",
  src_idx = 1:5,
  tar_idx = 6:10,
  funname = "exp"
)


test_that("propagate function of nonlinear_map works as expected with with.id=FALSE", {
  map <- create_nonlinear_map()
  map$setup(params)
  inp <- 1:10
  res <- map$propagate(inp, with.id=FALSE)
  expres <- c(rep(0,5), exp(1:5))
  expect_equal(res, expres)
})


test_that("propagate function of nonlinear_map works as expected with with.id=TRUE", {
  map <- create_nonlinear_map()
  map$setup(params)
  inp <- 1:10
  res <- map$propagate(inp, with.id=TRUE)
  expres <- c(rep(0,5), exp(1:5))
  expres <- expres + inp
  expect_equal(res, expres)
})
