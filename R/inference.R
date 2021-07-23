
gls <- function(map, zprior, U, obs, zref=zprior) {

  obsmask <- which(!is.na(obs))

  # prepare the Taylor expansion
  zref[obsmask] <- 0
  yref <- map$propagate(zref, with.id=TRUE)
  S <- map$jacobian(zref, with.id=TRUE)
  S <- drop0(S)

  # prepare the statistical model description
  # and auxiliary quantities
  invU <- solve(U, sparse=TRUE)

  v <- zprior
  v[obsmask] <- obs[obsmask]
  vref <- zref
  vref[obsmask] <- yref[obsmask]

  Sred <- S[obsmask,-obsmask]
  Sred1 <- S[,-obsmask]
  invPostU_red <- forceSymmetric(t(Sred1) %*% invU %*% Sred1)

  s1 <- solve(U, v - vref)
  s2 <- t(Sred1) %*% s1
  s3 <- solve(invPostU_red, s2)

  zpost <- rep(0, length(zprior))
  zpost[-obsmask] <- zref[-obsmask] + s3

  # calculate posterior of experimental values
  zpost[obsmask] <- zprior[obsmask] - Sred %*% zpost[-obsmask]
  zpost
}

