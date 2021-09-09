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
