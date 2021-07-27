# system definition

sysdt <- data.table(
  idx = 1:10,
  type = c(rep("mod", 5), rep("exp", 5)),
  en = c(seq(0, 10, length=5), seq(1, 9, length=5)),
  data = c(rep(0, 10)),
  unc = c(rep(10, 5), rep(1, 5)),
  obs = c(rep(NA,5), rep(4,5))
)

params <- list(
  mapname = "linmod_map",
  src_idx = sysdt[type=="mod",idx],
  tar_idx = sysdt[type=="exp",idx],
  src_x = sysdt[type=="mod",en],
  tar_x = sysdt[type=="exp",en]
)

# tests starting

test_that("Bayesian network inference provides correct posterior estimate", {
  linmod_map <- create_linmod_map()
  linmod_map$setup(params)
  zprior <- sysdt[, data]
  U <- Diagonal(n=nrow(sysdt), x=sysdt$unc)
  # prepare standard Bayesian update
  obs <- sysdt$obs
  obsmask <- which(!is.na(obs))
  x0 <- sysdt$data[-obsmask]
  A <- U[-obsmask, -obsmask]
  y <- sysdt$obs[obsmask]
  B <- U[obsmask, obsmask]
  S <- linmod_map$jacobian(zprior)[obsmask, -obsmask]
  # perform the update
  A1 <- solve(solve(A) + t(S)%*%solve(B)%*%S)
  x1 <- as.vector(A1 %*% (t(S) %*% solve(B) %*% y + solve(A) %*% x0))
  # prepare Bayesian network update
  res <- gls(linmod_map, zprior, U, obs)
  res <- res[-obsmask]
  expect_equal(x1, res)
})