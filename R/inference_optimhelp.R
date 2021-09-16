#' Evaluate posterior density function
#'
#' Evaluate the posterior density function for a given vector of
#' independent variables.
#'
#' @param map Mapping object. Usually a compound map, see \code{\link{create_compound_map}}
#' @param zprior Vector of prior estimates of the independent variables (i.e., associated with nodes without parent nodes)
#' @param U Prior covariance matrix of the independent variables
#' @param obs Vector with observed values of dependent nodes. Must be of same
#'            size as \code{zprior}. An \code{NA} value in this vector means that the
#'            corresponding variable was not observed.
#' @param zref Vector with values of independent variables used as reference vector
#'             in Taylor expansion. The posterior pdf is also evaluated exactly
#'             at this point.
#'
#' @return
#' Return a list with the following elements:
#' \tabular{ll}{
#'   \code{val} \tab Value of the posterior probability density function evaluated at \code{zref} \cr
#'   \code{jac} \tab Gradient of the posteroir probability density function evaluated at \code{zref}
#' }
#'
#' @export
#'
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


#' Generate helper functions for optimization routines
#'
#' Generate a list with two function to evaluate the posterior probability
#' density function and its gradient. These functions can be used in
#' optimization to find the posterior maximum.
#'
#' @param map Mapping object. Usually a compound map, see \code{\link{create_compound_map}}
#' @param zprior Vector of prior estimates of the independent variables (i.e., associated with nodes without parent nodes)
#' @param U Prior covariance matrix of the independent variables
#' @param obs Vector with observed values of dependent nodes. Must be of same
#'            size as \code{zprior}. An \code{NA} value in this vector means that the
#'            corresponding variable was not observed.
#' @param zref Vector with values of independent variables used as reference vector
#'             in Taylor expansion. The posterior pdf is also evaluated exactly
#'             at this point.
#'
#' @return
#' Return a list with two functions:
#' \tabular{ll}{
#'   \code{fun(x)} \tab Function that evaluates the posterior pdf at \code{x} \cr
#'   \code{jac(x)} \tab Function that evaluates the gradient of the posterior pdf at \code{x}
#' }
#' @export
#'
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
