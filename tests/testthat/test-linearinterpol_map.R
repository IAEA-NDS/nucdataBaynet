params <- list(
  maptype = "linearinterpol_map",
  src_idx = 1:4,
  tar_idx = 6:10,
  src_x = 11:14,
  tar_x = c(12,12.5,13,13.5,14)
)


test_that("index shuffling of src_idx and tar_idx produce the expected permutated result", {
  perm <- sample(10)
  invperm <- match(1:10, perm)
  perm_params <- with(params, list(
    maptype = maptype,
    src_x = src_x,
    tar_x = tar_x,
    src_idx = perm[src_idx],
    tar_idx = perm[tar_idx]
  ))
  map <- create_linearinterpol_map()
  perm_map <- create_linearinterpol_map()
  map$setup(params)
  perm_map$setup(perm_params)
  inp <- 1:10
  perm_inp <- inp[invperm]
  # with.id = FALSE
  res1 <- map$propagate(inp, with.id=FALSE)
  res2 <- perm_map$propagate(perm_inp, with.id=FALSE)
  expect_equal(res1, res2[perm])
  # with.id = TRUE
  res1 <- map$propagate(inp, with.id=TRUE)
  res2 <- perm_map$propagate(perm_inp, with.id=TRUE)
  expect_equal(res1, res2[perm])
})



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


test_that("setup and set_tar_x fail if some tar_x outside mesh if zero_outside not TRUE", {
  cur_params <- params
  cur_params[["tar_x"]][2] <- 300
  map <- create_linearinterpol_map()
  expect_error(map$setup(cur_params))
  cur_params <- params
  map <- create_linearinterpol_map()
  map$setup(cur_params)
  tar_x <- map$get_tar_x()
  tar_x[3] <- 300
  expect_error(map$set_tar_x(tar_x))
})


test_that("setup and set_tar_x do not fail if some tar_x outside mesh and zero_outside TRUE", {
  cur_params <- params
  cur_params[["tar_x"]][2] <- 300
  cur_params[["zero_outside"]] <- TRUE
  map <- create_linearinterpol_map()
  map$setup(cur_params)
  inp <- 1:10
  # one element erased
  res <- map$propagate(1:10, with.id=FALSE)
  expect_true(res[map$get_tar_idx()[2]] == 0)
  res <- map$propagate(1:10, with.id=TRUE)
  expect_true((res-inp)[map$get_tar_idx()[2]] == 0)
  # all tar_x outside src_x mesh
  map$set_tar_x(c(-10, 20, -30))
  res <- map$propagate(1:10, with.id=FALSE)
  expect_true(all(res[map$get_tar_idx()] == 0))
  res <- map$propagate(1:10, with.id=TRUE)
  expect_true((all((res-inp) == 0)))
})
