gls <- function(map, zprior, U, obs, zref=zprior) {

  obsmask <- !is.na(obs)

  # prepare the Taylor expansion
  zref[obsmask] <- 0
  yref <- map$propagate(zref, with.id=TRUE)
  S <- map$jacobian(zref, with.id=TRUE)
  S <- drop0(S)

  # prepare the statistical model description
  # and auxiliary quantities
  adjustable <- diag(U) > 0
  isindep <- !obsmask & adjustable
  stopifnot(all(diag(U)[obsmask] > 0))

  Ured <- U[adjustable,adjustable]
  invUred <- solve(Ured, sparse=TRUE)

  v <- zprior
  v[obsmask] <- obs[obsmask] - zprior[obsmask]
  vref <- zref
  vref[obsmask] <- yref[obsmask]
  # select only adjustable or observed
  v <- v[adjustable]
  vref <- vref[adjustable]

  Sred <- S[obsmask, isindep]
  Sred1 <- S[adjustable, isindep]
  invPostU_red <- forceSymmetric(t(Sred1) %*% invUred %*% Sred1)

  s1 <- solve(Ured, v - vref)
  s2 <- t(Sred1) %*% s1
  s3 <- solve(invPostU_red, s2)

  zpost <- zprior
  zpost[isindep] <- zref[isindep] + s3
  # calculate posterior of experimental values
  zpost[obsmask] <- obs[obsmask] - (yref[obsmask] + Sred %*% (zpost[isindep] - zref[isindep]))
  zpost
}
