gls <- function(map, zprior, U, obs, zref=zprior, damp=0, ret.list=FALSE) {

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
  stopifnot(all(zref[!adjustable] == zprior[!adjustable]))

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
  # if run in the loop of Levenberg-Marquart optimization
  diag(invPostU_red) <- diag(invPostU_red) + damp

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


LMalgo <- function(map, zprior, U, obs, zref=zprior, print.info=FALSE, ret.invcov=FALSE) {

  # LM parameters
  reltol <- 1e-5
  maxcount <- 200
  mincount <- 10
  tau <- 1e-8

  adjustable <- diag(U) > 0
  observed <- !is.na(obs)
  isindep <- adjustable & !observed
  adjobs <- observed[adjustable]
  Ured <- U[adjustable, adjustable]

  # determine the apriori scale for the damping term
  S <- map$jacobian(zref, with.id=TRUE)
  S <- S[adjustable, isindep]
  invUred_post <- crossprod(S, solve(Ured, S))
  lambda <- tau * max(diag(invUred_post))
  remove(S, invUred_post)

  # run the LM optimization loop
  yref <- map$propagate(zref)
  zref[observed] <- obs[observed] - yref[observed]
  last_d <- zref[adjustable] - zprior[adjustable]
  init_fex <- as.vector(crossprod(last_d, solve(Ured, last_d)))
  last_fapx <- last_fex <- init_fex
  cnt <- 0
  last_reject <- TRUE
  rel_gain <- Inf
  while (cnt < mincount ||
         (cnt < maxcount &&
          abs(relgain) > 1e-14 &&
          !(last_reject == FALSE && relgain <= reltol))) {
    cnt <- cnt + 1
    zprop <- gls(map, zprior, U, obs, zref=zref, damp=lambda)
    # calculate improvement according to linear approximation
    dapx <- zprop[adjustable] - zprior[adjustable]
    fapx <- as.vector(crossprod(dapx, solve(Ured, dapx)))
    # calculate improvement according to exact nonlinear map
    excess <- map$propagate(zprop)[observed] - obs[observed]
    dex <- dapx
    dex[adjobs] <- dex[adjobs] - excess
    fex <- as.vector(crossprod(dex, solve(Ured, dex)))
    # calculate the gain
    gain <- (last_fex - fex) / (last_fapx - fapx)
    relgain <- (last_fex - fex) / fex
    # adjust damping
    if (gain < 0.25) {
      lambda <- lambda * 2
    } else if (gain > 0.75) {
      lambda <- lambda / 3
    }
    # accept proposal if better (smaller) than current one
    if (fex - last_fex < 0) {
      zref <- zprop
      last_fapx <- fapx
      last_fex <- fex
      last_reject <- FALSE
    } else {
      last_reject <- TRUE
    }
    # print information
    if (print.info) {
      cat(paste0("iter ", cnt, " - gain: ", gain, " - relgain: ", relgain, " - fex: ", last_fex, "\n"))
      cat(paste0("last_reject: ", last_reject, "\n"))
    }
  }

  res <- list(
    zpost = zref,
    init_val = init_fex,
    final_val = last_fex,
    numiter = cnt
  )
  # determine the a posteriori covariance matrix
  if (ret.invcov) {
    S <- map$jacobian(zref, with.id=TRUE)
    S <- S[adjustable, isindep]
    invU <- crossprod(S, solve(Ured, S))
    mask <- isindep
    res[["invcov"]] <- invU
    res[["covmask"]] <- mask
  }
  return(res)
}

