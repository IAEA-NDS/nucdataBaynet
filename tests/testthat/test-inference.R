# for reproducibility of tests involving random numbers
set.seed(23)

# system definition

sysdt <- data.table(
  idx = 1:10,
  type = c(rep("mod", 5), rep("exp", 5)),
  en = c(seq(0, 10, length=5), seq(1, 9, length=5)),
  data = c(1:10),
  unc = c(rep(10, 5), rep(1, 5)),
  obs = c(rep(NA,5), rep(4,5))
)

sysdt2 <- data.table(
  idx = 1:10,
  type = c(rep("mod", 5), rep("exp", 5)),
  en = c(seq(0, 10, length=5), seq(1, 9, length=5)),
  data = c(1:10),
  unc = c(rep(10, 5), rep(0, 3), rep(1, 2)),
  obs = c(rep(NA,8), rep(4,2))
)


params <- list(
  maptype = "linearinterpol_map",
  src_idx = sysdt[type=="mod",idx],
  tar_idx = sysdt[type=="exp",idx],
  src_x = sysdt[type=="mod",en],
  tar_x = sysdt[type=="exp",en]
)


# tests starting

test_that("Bayesian network inference provides correct posterior estimate", {
  linearinterpol_map <- create_linearinterpol_map()
  linearinterpol_map$setup(params)
  zprior <- sysdt[, data]
  U <- Diagonal(n=nrow(sysdt), x=sysdt$unc)
  # prepare standard Bayesian update
  obs <- sysdt$obs
  obsmask <- which(!is.na(obs))
  x0 <- sysdt$data[-obsmask]
  A <- U[-obsmask, -obsmask]
  y <- sysdt$obs[obsmask] - sysdt$data[obsmask]
  B <- U[obsmask, obsmask]
  S <- linearinterpol_map$jacobian(zprior)[obsmask, -obsmask]
  # perform the update
  A1 <- solve(solve(A) + t(S)%*%solve(B)%*%S)
  x1 <- as.vector(A1 %*% (t(S) %*% solve(B) %*% y + solve(A) %*% x0))
  # prepare Bayesian network update
  res <- glsalgo(linearinterpol_map, zprior, U, obs)
  res <- res[-obsmask]
  expect_equal(x1, res)
})


test_that("Deterministic relationships are handled correctly", {
  cursysdt <- copy(sysdt)
  cursysdt[1, unc := 0]
  linearinterpol_map <- create_linearinterpol_map()
  linearinterpol_map$setup(params)
  zprior <- sysdt[, data]
  U <- Diagonal(n=nrow(cursysdt), x=cursysdt$unc)
  # prepare standard Bayesian update
  obs <- cursysdt$obs
  obsmask <- which(!is.na(obs))
  x0 <- zprior[-obsmask]
  effobs <- obs[obsmask] - cursysdt$data[obsmask]
  A <- U[-obsmask, -obsmask]
  y <- effobs
  B <- U[obsmask, obsmask]
  S <- linearinterpol_map$jacobian(zprior)[obsmask, -obsmask]
  # perform the update
  A1 <- A - A %*% t(S) %*% solve(S%*%A%*%t(S) + B) %*% S %*% A
  x1 <- as.vector(x0 + A %*% t(S) %*% solve(S%*%A%*%t(S) + B) %*% (effobs - S %*% x0))
  # prepare Bayesian network update
  res <- glsalgo(linearinterpol_map, zprior, U, obs)
  res <- res[-obsmask]
  expect_equal(x1, res)
})


test_that("Parts of posterior covariance matrix associated with independent variables correctly computed", {
  cursysdt <- copy(sysdt)
  linearinterpol_map <- create_linearinterpol_map()
  linearinterpol_map$setup(params)
  zprior <- sysdt[, data]
  obs <- cursysdt$obs
  U <- Diagonal(n=nrow(cursysdt), x=cursysdt$unc)
  zpost <- glsalgo(linearinterpol_map, zprior, U, obs)
  # prepare expected result
  is_indep <- is.na(obs)
  S <- linearinterpol_map$jacobian(zpost, with.id=TRUE)
  expres <- solve(t(S[,is_indep]) %*% solve(U) %*% S[,is_indep])
  # compare to obtained result
  res <- get_posterior_cov(linearinterpol_map, zpost, U, obs, 1:3, 1:3)
  expect_equal(res, expres[1:3,1:3])
  res <- get_posterior_cov(linearinterpol_map, zpost, U, obs, c(1,3,5), 2:4)
  expect_equal(res, expres[c(1,3,5),2:4])
})


test_that("Parts of posterior covariance matrix associated with noise nodes of observed variables correctly computed", {
  cursysdt <- copy(sysdt)
  linearinterpol_map <- create_linearinterpol_map()
  linearinterpol_map$setup(params)
  zprior <- sysdt[, data]
  obs <- cursysdt$obs
  U <- Diagonal(n=nrow(cursysdt), x=cursysdt$unc)
  zpost <- glsalgo(linearinterpol_map, zprior, U, obs)
  # prepare expected result
  is_indep <- is.na(obs)
  is_obs <- !is.na(obs)
  S <- linearinterpol_map$jacobian(zpost, with.id=TRUE)
  expres <- U - U %*% t(S[is_obs,]) %*% solve(S[is_obs,] %*% U %*% t(S[is_obs,])) %*% S[is_obs,] %*% U
  # compare to obtained result
  res <- get_posterior_cov(linearinterpol_map, zpost, U, obs, c(1,3,5,7,8), c(2,6,4))
  expect_equal(res, expres[c(1,3,5,7,8), c(2,6,4)])
})


test_that("Parts of posterior covariance matrix associated with observed variables correctly computed for ret.dep=TRUE", {
  cursysdt <- copy(sysdt)
  linearinterpol_map <- create_linearinterpol_map()
  linearinterpol_map$setup(params)
  zprior <- sysdt[, data]
  obs <- cursysdt$obs
  U <- Diagonal(n=nrow(cursysdt), x=cursysdt$unc)
  zpost <- glsalgo(linearinterpol_map, zprior, U, obs)
  res <- get_posterior_cov(linearinterpol_map, zpost, U, obs, c(1,3,5,7,8), c(2,7,4), ret.dep=TRUE)
  # prepare expected result
  is_indep <- is.na(obs)
  is_obs <- !is.na(obs)
  S <- linearinterpol_map$jacobian(zpost, with.id=TRUE)
  expres <- U - U %*% t(S[is_obs,]) %*% solve(S[is_obs,] %*% U %*% t(S[is_obs,])) %*% S[is_obs,] %*% U
  expres <- expres[c(1,3,5,7,8),c(2,7,4)]
  expres[c(4,5),] <- 0
  expres[,2] <- 0
  expect_equal(res, expres)
})


test_that("Parts of posterior covariance matrix correctly computed for deterministic nodes", {
  cursysdt <- copy(sysdt)
  cursysdt[idx==1 & type=="mod", unc:=0]
  linearinterpol_map <- create_linearinterpol_map()
  linearinterpol_map$setup(params)
  zprior <- sysdt[, data]
  obs <- cursysdt$obs
  U <- Diagonal(n=nrow(cursysdt), x=cursysdt$unc)
  zpost <- glsalgo(linearinterpol_map, zprior, U, obs)
  # prepare expected result
  is_indep <- is.na(obs)
  is_obs <- !is.na(obs)
  S <- linearinterpol_map$jacobian(zpost, with.id=TRUE)
  globexpres <- U - U %*% t(S[is_obs,]) %*% solve(S[is_obs,] %*% U %*% t(S[is_obs,])) %*% S[is_obs,] %*% U
  # compare to obtained results
  res <- as.matrix(get_posterior_cov(linearinterpol_map, zpost, U, obs, c(1,3,5,7,8), 2:4))
  expres <- as.matrix(globexpres[c(1,3,5,7,8), 2:4])
  expect_equal(res, expres)
  res <- as.matrix(get_posterior_cov(linearinterpol_map, zpost, U, obs, 1, 2:6))
  expres <- as.matrix(globexpres[1, 2:6, drop=FALSE])
  dimnames(res) <- dimnames(expres) <- NULL
  expect_equal(res, expres)
})


test_that("Parts of posterior covariance matrix for deterministic nodes correctly computed for ret.dep=TRUE", {
  cursysdt <- copy(sysdt2)
  linearinterpol_map <- create_linearinterpol_map()
  linearinterpol_map$setup(params)
  zprior <- sysdt[, data]
  obs <- cursysdt$obs
  U <- Diagonal(n=nrow(cursysdt), x=cursysdt$unc)
  zpost <- glsalgo(linearinterpol_map, zprior, U, obs)
  # prepare expected result
  is_indep <- is.na(obs)
  is_obs <- !is.na(obs)
  is_fixed <- is.na(obs) & diag(U) == 0
  S <- linearinterpol_map$jacobian(zpost, with.id=TRUE)
  globexpres <- U - U %*% t(S[is_obs,]) %*% solve(S[is_obs,] %*% U %*% t(S[is_obs,])) %*% S[is_obs,] %*% U
  # compare to obtained results
  expres <- as.matrix(S %*% globexpres %*% t(S))
  res <- as.matrix(get_posterior_cov(linearinterpol_map, zpost, U, obs, 1:10, 1:10, ret.dep=TRUE))
  expect_equal(res, expres)
})


test_that("Posterior samples are reasonably consistent with posterior distribution", {
  cursysdt <- copy(sysdt)
  linearinterpol_map <- create_linearinterpol_map()
  linearinterpol_map$setup(params)
  zprior <- sysdt[, data]
  obs <- cursysdt$obs
  U <- Diagonal(n=nrow(cursysdt), x=cursysdt$unc)
  zpost <- glsalgo(linearinterpol_map, zprior, U, obs)
  postcov <- get_posterior_cov(linearinterpol_map, zpost, U, obs,
                               seq_along(zpost), seq_along(zpost))
  num <- 1e4
  smpl <- get_posterior_sample(linearinterpol_map, zpost, U, obs, num)
  smplcov <- cov.wt(t(smpl), method="unbiased")$cov
  sigma <- sqrt(2/num) * diag(postcov)
  zscores <- (diag(postcov) - diag(smplcov)) / sigma
  expect_true(all(zscores < 3))
})
