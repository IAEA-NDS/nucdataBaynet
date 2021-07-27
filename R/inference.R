gls <- function(map, zprior, U, obs, zref=zprior, ret.list=FALSE) {

  obsmask <- !is.na(obs)

  # prepare the Taylor expansion
  zref[obsmask] <- 0
  S <- map$jacobian(zref, with.id=TRUE)
  S <- drop0(S)
  yref <- map$propagate(zref, with.id=TRUE)

  # prepare the statistical model description
  # and auxiliary quantities
  adjustable <- diag(U) > 0
  isindep <- !obsmask & adjustable
  stopifnot(all(diag(U)[obsmask] > 0))

  Ured <- U[adjustable,adjustable]

  v <- zprior
  v[obsmask] <- obs[obsmask] - zprior[obsmask]
  vref <- zref
  vref[obsmask] <- yref[obsmask]
  # select only adjustable or observed
  v <- v[adjustable]
  vref <- vref[adjustable]

  Sred <- S[obsmask, isindep]
  Sred1 <- S[adjustable, isindep]
  invPostU_red <- forceSymmetric(crossprod(Sred1, solve(Ured, Sred1, sparse=TRUE)))

  s1 <- solve(Ured, v - vref)
  s2 <- t(Sred1) %*% s1
  s3 <- solve(invPostU_red, s2)

  zpost <- zprior
  zpost[isindep] <- zref[isindep] + s3
  # calculate posterior of experimental values
  zpost[obsmask] <- obs[obsmask] - (yref[obsmask] + Sred %*% (zpost[isindep] - zref[isindep]))

  if (ret.list) {
    return(list(
      zpost = zpost,
      ypost = yref + as.vector(S %*% (zpost - zref))
    ))
  } else {
    return(zpost)
  }
}
