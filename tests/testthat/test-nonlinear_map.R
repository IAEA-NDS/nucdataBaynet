# for reproducible "random" numbers
set.seed(31)

params <- list(
  maptype = "nonlinear_map",
  src_idx = 1:5,
  tar_idx = 6:10,
  funname = "exp"
)


test_that("index shuffling of src_idx and tar_idx produce the expected permutated result", {
  perm <- sample(10)
  invperm <- match(1:10, perm)
  perm_params <- with(params, list(
    src_idx = perm[src_idx],
    tar_idx = perm[tar_idx],
    funname = funname
  ))
  map <- create_nonlinear_map()
  perm_map <- create_nonlinear_map()
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


test_that("propagate function of nonlinear_map works as expected with exp nonlinearity", {
  map <- create_nonlinear_map()
  map$setup(params)
  inp <- 1:10
  # with.id = FALSE
  res <- map$propagate(inp, with.id=FALSE)
  expres <- c(rep(0,5), exp(1:5))
  expect_equal(res, expres)
  # with.id = TRUE
  res <- map$propagate(inp, with.id=TRUE)
  expres <- c(rep(0,5), exp(1:5))
  expres <- expres + inp
  expect_equal(res, expres)
})


test_that("propagate function of nonlinear_map works as expected with relu nonlinearity", {
  cur_params <- modifyList(params, list(funname = "relu"))
  map <- create_nonlinear_map()
  map$setup(cur_params)
  inp <- c(1,-2,3,-4,5, 6:10)
  # with.id = FALSE
  res <- map$propagate(inp, with.id=FALSE)
  expres <- c(rep(0,5), c(1,0,3,0,5))
  expect_equal(res, expres)
  # with.id = TRUE
  res <- map$propagate(inp, with.id=TRUE)
  expres <- c(rep(0,5), c(1,0,3,0,5)) + inp
  expect_equal(res, expres)
  # with.id = TRUE
})


test_that("jacobian function of nonlinear_map equivalent to numeric jacobian with exp nonlinearity", {
  map <- create_nonlinear_map()
  map$setup(params)
  inp <- 1:10
  # with.id = FALSE
  expres <- jacobian(map$propagate, inp, with.id=FALSE)
  res <- as.matrix(map$jacobian(inp, with.id=FALSE))
  dimnames(res) <- NULL
  expect_equal(res, expres)
  # with.id = TRUE
  expres <- jacobian(map$propagate, inp, with.id=TRUE)
  res <- as.matrix(map$jacobian(inp, with.id=TRUE))
  dimnames(res) <- NULL
  expect_equal(res, expres)
})


test_that("jacobian function of nonlinear_map equivalent to numeric jacobian with relu nonlinearity", {
  cur_params <- modifyList(params, list(funname = "relu"))
  map <- create_nonlinear_map()
  map$setup(cur_params)
  inp <- c(1,-2,3,-4,5, 6:10)
  # with.id = FALSE
  expres <- jacobian(map$propagate, inp, with.id=FALSE)
  res <- as.matrix(map$jacobian(inp, with.id=FALSE))
  dimnames(res) <- NULL
  expect_equal(res, expres)
  # with.id = TRUE
  expres <- jacobian(map$propagate, inp, with.id=TRUE)
  res <- as.matrix(map$jacobian(inp, with.id=TRUE))
  dimnames(res) <- NULL
  expect_equal(res, expres)
})
