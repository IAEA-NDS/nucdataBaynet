params <- list(
  mapname = "derivative2nd_map",
  src_idx = 1:5,
  tar_idx = 6:10,
  src_x = 11:15,
  tar_x = c(12, 12, 13, 13, 14)
)


test_that("map_derivative2nd works as expected with with.id=FALSE", {
  testmap <- create_derivative2nd_map()
  testmap$setup(params)
  # case 1, no curvature
  inp <- 1:10
  res <- testmap$propagate(inp, with.id=FALSE)
  expres <- rep(0, 10)
  expect_equal(res, expres)
  # case 2
  inp <- c(1,0,1,-1,2,0,0,0,0,0)
  res <- testmap$propagate(inp, with.id=FALSE)
  expres <- c(rep(0, 5), 1, 1, -1.5, -1.5, 2.5)
  expect_equal(res, expres)
})


test_that("map_derivative2nd works as expected with with.id=TRUE", {
  testmap <- create_derivative2nd_map()
  testmap$setup(params)
  inp <- c(1,0,1,-1,2,0,0,0,0,0)
  res <- testmap$propagate(inp, with.id=TRUE)
  expres <- c(inp[1:5], 1, 1, -1.5, -1.5, 2.5)
  expect_equal(res, expres)
})


test_that("map_derivative2nd correctly adds noise term if with.id=TRUE", {
  testmap <- create_derivative2nd_map()
  testmap$setup(params)
  inp <- c(1,0,1,-1,2,0,0,0,0,0)
  addnoise <- 1:5
  inp[6:10] <- inp[6:10] + addnoise
  res <- testmap$propagate(inp, with.id=TRUE)
  expres <- c(inp[1:5], 1, 1, -1.5, -1.5, 2.5)
  expres[6:10] <- expres[6:10] + addnoise
  expect_equal(res, expres)
})
