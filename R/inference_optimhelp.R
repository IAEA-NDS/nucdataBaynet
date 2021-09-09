evaluate_logpostpdf <- function(map, zprior, U, obs, zref) {

  # error checking of arguments
  stopifnot(length(obs) == length(zprior))
  stopifnot(nrow(U) == length(zprior))
  stopifnot(ncol(U) == length(zprior))
  stopifnot(length(zref) == length(zprior))

  # compute helpful quantities
  adjustable <- diag(U) > 0
  observed <- !is.na(obs)
  isindep <- adjustable & !observed
  adjobs <- observed[adjustable]
  Ured <- U[adjustable, adjustable]

  # we ignore the values in the attached noise nodes provided by the user
  # in zref, and compute their value here for consistency
  zref[observed] <- 0
  yref <- map$propagate(zref)
  zref[observed] <- obs[observed] - yref[observed]

  # NOTE: indices referred to by adjustable also include the observations
  d <- zref[adjustable] - zprior[adjustable]
  fex <- (-0.5) * as.vector(crossprod(d, solve(Ured, d)))

  # calculate the gradient
  # define the quantities as described in evaluation with Bayesian network paper
  v <- zprior
  v[observed] <- obs[observed] - zprior[observed]
  vref <- zref
  vref[observed] <- yref[observed]
  # select only adjustable (including observed)
  v <- v[adjustable]
  vref <- vref[adjustable]

  # compute inverse posterior matrix
  S <- map$jacobian(zref, with.id=TRUE)
  Sred <- S[observed, isindep]
  S <- S[adjustable, isindep]
  s1 <- solve(Ured, v - vref)
  s2 <- as.vector(t(S) %*% s1)

  return(list(
    val = fex,
    jac = s2,
    indepidx = isindep
  ))
}


generate_logpostpdf_funs <- function(map, zprior, U, obs, zref) {

  last_val <- NULL
  last_jac <- NULL
  isindep <- which(diag(U) > 0 & is.na(obs))

  update <- function(x) {
    stopifnot(length(isindep) == length(x))
    if (is.null(last_val) || any(x != last_val)) {
      thisz <- zref
      thisz[isindep] <- x
      ret <- evaluate_logpostpdf(map, zprior, U, obs, thisz)
      last_val <<- ret$val
      last_jac <<- ret$jac
    }
  }

  fun <- function(x) {
    update(x)
    return(last_val)
  }

  jac <- function(x) {
    update(x)
    return(last_jac)
  }

  return(list(
    fun = fun,
    jac = jac,
    indep_idx = isindep
  ))
}
