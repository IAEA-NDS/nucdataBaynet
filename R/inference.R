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


get_posterior_cov <- function(map, zpost, U, obs, row_idcs, col_idcs) {

  stopifnot(length(zpost) == nrow(U))
  stopifnot(length(zpost) == ncol(U))

  if (is.logical(row_idcs)) {
    row_idcs <- which(row_idcs)
  }
  if (is.logical(col_idcs)) {
    col_idcs <- which(col_idcs)
  }

  adjustable <- diag(U) > 0
  observed <- !is.na(obs)
  isindep <- adjustable & !observed
  isfixed <- !adjustable & !observed
  adjobs <- observed[adjustable]

  is_indep_row <- isindep[row_idcs]
  is_indep_col <- isindep[col_idcs]
  is_obs_row <- observed[row_idcs]
  is_obs_col <- observed[col_idcs]

  is_fixed_row <- isfixed[row_idcs]
  is_fixed_col <- isfixed[col_idcs]

  map_idcs_to_indep <- rep(NA, length(zpost))
  map_idcs_to_indep[isindep] <- seq_len(sum(isindep))

  # determine uncertainties of independent variables
  S <- drop0(map$jacobian(zpost, with.id=TRUE))
  Sred <- S[adjustable, isindep]
  Ured <- U[adjustable, adjustable]
  invPostU_red <- forceSymmetric(crossprod(Sred, solve(Ured, Sred, sparse=TRUE)))

  S_left <- sparseMatrix(i=seq_along(row_idcs)[is_indep_row],
                         j=map_idcs_to_indep[row_idcs[is_indep_row]],
                         x=1,
                         dims=c(length(row_idcs), ncol(Sred)))
  S_left[is_obs_row,] <- S[row_idcs[is_obs_row], isindep]
  S_left <- drop0(S_left)

  S_right <- sparseMatrix(i=seq_along(col_idcs)[is_indep_col],
                          j=map_idcs_to_indep[col_idcs[is_indep_col]],
                          x=1,
                          dims=c(length(col_idcs), ncol(Sred)))
  S_right[is_obs_col,] <- S[col_idcs[is_obs_col], isindep]
  S_right <- drop0(S_right)

  Upost <- S_left %*% solve(invPostU_red, t(S_right))
  # fix the covariance matrix of error variables attached to observed nodes
  Upost[is_obs_row,] <- (-Upost[is_obs_row,])
  Upost[,is_obs_col] <- (-Upost[,is_obs_col])
  return(Upost)
}


get_posterior_sample <- function(map, zpost, U, obs, num) {

  adjustable <- diag(U) > 0
  observed <- !is.na(obs)
  isindep <- adjustable & !observed
  adjobs <- observed[adjustable]

  numindep <- sum(isindep)
  numobs <- sum(observed)

  S <- drop0(map$jacobian(zpost, with.id=TRUE))
  Sred <- S[adjustable, isindep]
  Ured <- U[adjustable, adjustable]
  invPostU_red <- forceSymmetric(crossprod(Sred, solve(Ured, Sred, sparse=TRUE)))

  Lfact <- chol(invPostU_red)
  rvec <- matrix(rnorm(numindep*num), nrow=numindep)
  zpost[observed] <- 0
  smpl <- matrix(rep(zpost, num), ncol=num, nrow=length(zpost))
  smpl[isindep,] <- smpl[isindep,] + as.matrix(solve(Lfact, rvec, system="L"))
  # calculate the errors associated with observed nodes
  for (i in seq_len(num)) {
    smpl[observed,i] <- obs[observed] - map$propagate(zpost, with.id=TRUE)[observed]
  }
  return(smpl)

  Sred <- S[adjustable, isindep]
  invPostU_red <- forceSymmetric(crossprod(Sred, solve(Ured, Sred, sparse=TRUE)))

  A <- matrix(runif(9), 3, 3)
  A <- A %*% t(A)
  x <- runif(3)


  L <- chol(A)
  t(L) %*% L

  A2 <- as(A, "sparseMatrix")
  L2 <- Cholesky(A2)

  solve(L, x)
  solve(L2, x, system="L")

  invA <- solve(A)
  res1 <- rmvnorm(1e6, rep(0, 3), invA)
  cov.wt(res1)

  rvec <- matrix(rnorm(3*10000), nrow=3)
  rvec2 <- solve(L2, rvec, system="L")
  cov.wt(as.matrix(t(rvec2)))

  S <- matrix(runif(9), 3, 3)
}
