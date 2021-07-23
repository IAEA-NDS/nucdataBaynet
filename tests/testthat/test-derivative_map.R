params <- list(
  mapname = "derivative_map",
  src_idx = 1:5,
  tar_idx = 6:10,
  src_x = 11:15,
  tar_x = c(12, 12, 13, 13, 14)
)


test_that("map_derivative works as expected with with.id=FALSE", {
  testmap <- create_derivative_map()
  testmap$setup(params)
  inp <- c(1,3,6,10,15, rep(0, 5))
  res <- testmap$propagate(inp, with.id=FALSE)
  expres <- c(0,0,0,0,0,3,3,4,4,5)
  expect_equal(res, expres)
})


test_that("map_derivative works as expected with with.id=TRUE", {
  testmap <- create_derivative_map()
  testmap$setup(params)
  inp <- c(1,3,6,10,15, rep(0, 5))
  res <- testmap$propagate(inp, with.id=TRUE)
  expres <- c(1,3,6,10,15,3,3,4,4,5)
  expect_equal(res, expres)
})


test_that("map_derivative correctly adds noise term if with.id=TRUE", {
  testmap <- create_derivative_map()
  testmap$setup(params)
  inp <- c(1,3,6,10,15, 1:5)
  res <- testmap$propagate(inp, with.id=TRUE)
  expres <- c(1,3,6,10,15,3,3,4,4,5)
  expres[6:10] <- expres[6:10] + 1:5
  expect_equal(res, expres)
})
