set.seed(59)

yref <- runif(3)
pref <- runif(3)
S <- matrix(runif(9), 3, 3)

linearmap_def <- list(
  maptype = "linearmap_map",
  mapname = "linearmap",
  src_idx = c(5,1,7),
  tar_idx = c(2,6,3),
  S = S, yref = yref, pref = pref
)


test_that("linearmap_map returns correct info", {
  m <- linearmap_def
  linmap <- create_map(m)
  expect_equal(linmap$getType(), "linearmap_map")
  expect_equal(linmap$get_src_idx(), m$src_idx)
  expect_equal(linmap$get_tar_idx(), m$tar_idx)
})


test_that("linearmap_map propagates correctly", {
  m <- linearmap_def
  linmap <- create_map(m)
  x <- runif(10)
  # with.id = TRUE
  expres <- x
  expres[m$tar_idx] <- expres[m$tar_idx] + m$yref + m$S %*% (x[m$src_idx] - m$pref)
  expres <- as.vector(expres)
  res <- linmap$propagate(x, with.id=TRUE)
  expect_equal(res, expres)
  # with.id = FALSE
  expres <- rep(0, 10)
  expres[m$tar_idx] <- m$yref + m$S %*% (x[m$src_idx] - m$pref)
  res <- linmap$propagate(x, with.id=FALSE)
  expect_equal(res, expres)
})


test_that("linearmap_map returns correct jacobian matrix", {
  m <- linearmap_def
  linmap <- create_map(m)
  x <- runif(10)
  # with.id = TRUE
  res <- linmap$jacobian(x, with.id=TRUE)
  expres <- Diagonal(n=length(x), x=1)
  expres[m$tar_idx, m$src_idx] <- m$S
  expect_equal(res, expres)
  # with.id = FALSE
  diag(expres) <- diag(expres) - 1
  res <- linmap$jacobian(x, with.id=FALSE)
  expect_equal(res, expres)
})
