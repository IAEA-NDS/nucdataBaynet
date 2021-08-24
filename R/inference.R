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


LMalgo <- function(map, zprior, U, obs, zref=zprior, print.info=FALSE, adjust_idcs = NULL,
                   ret.invcov=FALSE, control=list()) {

  # LM parameters
  defcontrol = list(
    reltol = 1e-6, maxcount = 100, mincount = 10, tau = 1e-10,
    reltol_steps = 3, reltol2 = 1e-12)
  control <- modifyList(defcontrol, control)

  reltol <- control$reltol
  reltol2 <- control$reltol2
  maxcount <- control$maxcount
  mincount <- control$mincount
  tau <- control$tau

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
  relgain_hist <- rep(Inf, control$reltol_steps)
  while (cnt < mincount ||
         (cnt < maxcount &&
          max(relgain_hist) > abs(reltol) &&
          abs(relgain) > reltol2)) {
    cnt <- cnt + 1
    zprop <- gls(map, zprior, U, obs, zref=zref, damp=lambda)
    if (!is.null(adjust_idcs)) {
      tmp <- zprop
      zprop <- zref
      zprop[adjust_idcs] <- tmp[adjust_idcs]
      zprop[observed] <- tmp[observed]
    }
    # calculate improvement according to linear approximation
    dapx <- zprop[adjustable] - zprior[adjustable]
    fapx <- as.vector(crossprod(dapx, solve(Ured, dapx)))
    # calculate improvement according to exact nonlinear map
    fun_error <- FALSE
    tryCatch({
      excess <- map$propagate(zprop)[observed] - obs[observed]
    }, error = function(e) {
      fun_error <<- TRUE
    })
    if (fun_error) {
      last_reject <- TRUE
      gain <- -Inf
      lambda <- lambda * 10
      if (print.info) {
        cat(paste0("iter ", cnt, " - objective function crashed with proposed parameters\n"))
        cat(paste0("last_reject: ", last_reject, " - lambda: ", lambda, "\n"))
      }
      next
    }
    dex <- dapx
    dex[adjobs] <- dex[adjobs] - excess
    fex <- as.vector(crossprod(dex, solve(Ured, dex)))
    # calculate the gain
    gain <- (last_fex - fex) / (abs(last_fapx - fapx) + 1e-15)
    relgain <- (last_fex - fex) / fex
    if (gain > 0) {
      relgain_hist <- c(tail(relgain_hist, n=9), relgain)
    }
    # adjust damping
    if (gain < -2) {
      lambda <- lambda * 2^5
    } else if (gain < 0.25) {
      lambda <- lambda * 2
    } else if (gain > 0.75) {
      lambda <- lambda / 3
    }
    # print information
    if (print.info) {
      cat(paste0("### iter ", cnt, " ###\n",
                 "gain: ", gain, " - relgain: ", relgain, "\n",
                 "fex: ", fex, " - last fex: ", last_fex, "\n",
                 "accept: ", fex < last_fex, " - new lambda: ", lambda, "\n"))
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


get_posterior_cov <- function(map, zpost, U, obs, row_idcs, col_idcs, ret.dep=FALSE) {

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

  S_right <- sparseMatrix(i=seq_along(col_idcs)[is_indep_col],
                          j=map_idcs_to_indep[col_idcs[is_indep_col]],
                          x=1,
                          dims=c(length(col_idcs), ncol(Sred)))

  if (!ret.dep) {
    # fixed (independent nodes) have zero prior uncertainty
    # therefore we can exclude them from the sensitivity matrix if ret.dep=FALSE
    # for independent nodes attached to observed nodes, we need to determine
    # their posterior uncertainty by forward propagating the independent nodes
    S_left[is_obs_row,] <- S[row_idcs[is_obs_row], isindep]
    S_right[is_obs_col,] <- S[col_idcs[is_obs_col], isindep]
  } else {
    # observed dependent nodes have zero posterior uncertainty
    # so we can exclude them from the sensitivity matrix if ret.dep = TRUE
    # we need to include however the fixed nodes because their
    # uncertainties are given by forward propagating independent nodes
    S_left[is_fixed_row,] <- S[row_idcs[is_fixed_row], isindep]
    S_right[is_fixed_col,] <- S[col_idcs[is_fixed_col], isindep]
  }

  S_left <- drop0(S_left)
  S_right <- drop0(S_right)

  Upost <- S_left %*% solve(invPostU_red, t(S_right))
  # fix the covariance matrix of error variables attached to observed nodes
  if (!ret.dep) {
    Upost[is_obs_row,] <- (-Upost[is_obs_row,])
    Upost[,is_obs_col] <- (-Upost[,is_obs_col])
  }
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
  Lfact <- Cholesky(invPostU_red, perm=TRUE, LDL=FALSE)
  Pmat <- as(Lfact, "pMatrix")
  rvec <- matrix(rnorm(numindep*num), nrow=numindep)
  zpost[observed] <- 0
  smpl <- matrix(rep(zpost, num), ncol=num, nrow=length(zpost))
  smpl[isindep,] <- smpl[isindep,] + as.matrix(crossprod(Pmat, solve(Lfact, rvec, system="Lt")))
  # calculate the errors associated with observed nodes
  for (i in seq_len(num)) {
    smpl[observed,i] <- obs[observed] - map$propagate(smpl[,i], with.id=TRUE)[observed]
  }
  return(smpl)
}
