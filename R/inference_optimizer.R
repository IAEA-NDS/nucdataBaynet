#' Generalized Least Squares fitting
#'
#' Estimate the variables in a Bayesian network via the
#' Generalized Least Squares method.
#'
#' @param map Mapping object. Usually a compound map, see \code{\link{create_compound_map}}
#' @param zprior Vector of prior estimates of the independent variables (i.e., associated with nodes without parent nodes)
#' @param U Prior covariance matrix of the independent variables
#' @param obs Vector with observed values of dependent nodes. Must be of same
#'            size as \code{zprior}. An \code{NA} value in this vector means that the
#'            corresponding variable was not observed.
#' @param zref Values of the independent variable that should be used as reference point
#'             to calculate the Taylor approximation (using the \code{jacobian} function of \code{map})
#' @param adjust_idcs Indices of variables that should be adjusted, all other independent
#'                    variables are assumed to be fixed.
#' @param damp A value indicating the magnitude of the damping term applied to the
#'             inverse posterior covariance matrix. For a pure GLS update, leave it
#'             zero. This parameter is used by the
#'             in the case of a non-linear mapping.
#' @param ret.list If \code{FALSE}, return a vector with the estimated variables of the
#'                 independent variables. Otherwise, return a list with more information.
#'                 Default is \code{FALSE}
#'
#' @return
#' Return a vector with the posterior estimates of the independent variables if
#' \code{ret.list=FALSE}. Otherwise return a list with the following fields:
#' \tabular{ll}{
#' \code{zpost} \tab Vector of posterior estimates of the independent variables\cr
#' \code{ypost} \tab Vector of posterior estimates of the dependent variables
#' }
#'
#' @export
#'
#' @seealso \code{\link{LMalgo}}
#' @example man/examples/example_inference_01.R
#'
glsalgo <- function(map, zprior, U, obs, zref=zprior,
                adjust_idcs=NULL, damp=0, ret.list=FALSE) {

  stopifnot(length(obs) == length(zprior))
  stopifnot(nrow(U) == length(zprior))
  stopifnot(ncol(U) == length(zprior))
  stopifnot(length(zref) == length(zprior))

  obsmask <- !is.na(obs)

  # prepare the Taylor expansion
  zref[obsmask] <- 0
  S <- map$jacobian(zref, with.id=TRUE)
  S <- drop0(S)
  yref <- map$propagate(zref, with.id=TRUE)

  # prepare the statistical model description
  # and auxiliary quantities
  adjustable <- diag(U) > 0
  # tweak zprior and U if user requested to return the
  # conditional estimate given the elements referenced by adjust_idcs
  if (!is.null(adjust_idcs)) {
    adjust_mask <- rep(FALSE, length(zprior))
    adjust_mask[adjust_idcs] <- TRUE
    sel <- !obsmask & adjustable & !adjust_mask
    if (any(sel)) {
      zprior[sel] <- zref[sel]
      U[sel,] <- 0
      U[,sel] <- 0
      adjustable <- diag(U) > 0
    }
  }
  isindep <- !obsmask & adjustable

  stopifnot(all(diag(U)[obsmask] > 0))
  stopifnot(all(zref[!adjustable] == zprior[!adjustable]))

  # define the quantities as described in evaluation with Bayesian network paper
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
      ypost = yref + as.vector(S %*% (zpost - zref)),
      post_pdf_val = 1,
      post_pdf_grad = s2
    ))
  } else {
    return(zpost)
  }
}


#' Maximum a posterior estimate via LM algorithm
#'
#' Finds the values of the independent variables that maximize
#' the value of the posterior probability density function
#' using a customized Levenberg-Marquardt algorithm.
#'
#' The convergence criteria of the Levenberg-Marquardt algorithm can be
#' tweaked by several parameters in the \code{control} list passed as
#' argument. The following fields are available:
#' \tabular{ll}{
#'   \code{reltol} \tab If the relative improvement of posterior pdf value falls
#'                      below this value for a number of subsequent iterations
#'                      given by \code{reltol_steps}, terminate optimization.
#'                      (Default 1e-6) \cr
#'   \code{reltol_steps} \tab Number of subsequent iteratios with a relative
#'                            improvement of less than \code{reltol} necessary
#'                            to terminate optimization procedure.
#'                            (Default 3) \cr
#'   \code{reltol2} \tab If relative improvement falls below this value,
#'                       immediately terminate the optimization procedure.
#'                       (Default 1e-12) \cr
#'   \code{mincount} \tab Minimum number of iterations that should be performed
#'                        regardless of relative improvement.
#'                        (Default 10) \cr
#'   \code{maxcount} \tab Maximum number of iterations that should be performed
#'                        regardless of relative improvement.
#'                        (Default 100) \cr
#'   \code{maxreject} \tab Number of subsequent rejections of proposal vectors
#'                         before the LM algorithm gives up.
#'                         (Default 10) \cr
#'   \code{tau} \tab Parameter that influences the initial choice of the damping term.
#'                   The larger this parameter, the larger the initial damping applied
#'                   to the inverse posterior covariance matrix.
#'                   (Default 1e-10) \cr
#' }
#'
#' @param map Mapping object. Usually a compound map, see \code{\link{create_compound_map}}.
#' @param zprior Vector of prior estimates of the independent variables (i.e., associated with nodes without parent nodes)
#' @param U Prior covariance matrix of the independent variables
#' @param obs Vector with observed values of dependent nodes. Must be of same
#'            size as \code{zprior}. An \code{NA} value in this vector means that the
#'            corresponding variable was not observed.
#' @param zref Values of the independent variable that should be used as the initial reference point
#'             to calculate the Taylor approximation (using the \code{jacobian} function of \code{map})
#' @param print.info Display output regarding the progress of the optimization procedure
#' @param adjust_idcs Indices of variables that should be adjusted, all other independent
#'                    variables are assumed to be fixed. The default value \code{NULL} means
#'                    that all independent variables are adjusted.
#' @param ret.invcov If \code{TRUE}, also return the inverse posterior covariance matrix
#' @param control A list with control parameters to tweak the convergence criteria, see \strong{Details}.
#'
#' @return
#' Return a list with the following fields:
#' \tabular{ll}{
#'   \code{zpost} \tab Variable assignment corresponding to (potentially local) maximum of
#'   posterior density function. \cr
#'   \code{init_val} \tab Initial value of objective function (proportional to posterior pdf value) \cr
#'   \code{final_val} \tab Final value of objective function (proportional to posterior pdf value) \cr
#'   \code{numiter} \tab Number of iterations
#' }
#' @export
#'
#' @example man/examples/example_inference_01.R
#'
LMalgo <- function(map, zprior, U, obs, zref=zprior, print.info=FALSE, adjust_idcs = NULL,
                   ret.invcov=FALSE, control=list()) {

  # LM parameters
  defcontrol = list(
    reltol = 1e-6, maxcount = 100, mincount = 10, tau = 1e-10,
    reltol_steps = 3, reltol2 = 1e-12, maxreject=10)
  control <- modifyList(defcontrol, control)

  reltol <- control$reltol
  reltol2 <- control$reltol2
  maxcount <- control$maxcount
  mincount <- control$mincount
  tau <- control$tau
  maxreject <- control$maxreject

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
  # we ignore the values in the attached noise nodes provided by the user
  # in zref, and compute their value here for consistency
  zref[observed] <- 0
  yref <- map$propagate(zref)
  zref[observed] <- obs[observed] - yref[observed]
  # NOTE: indices referred to by adjustable also include the observations
  last_d <- zref[adjustable] - zprior[adjustable]
  init_fex <- as.vector(crossprod(last_d, solve(Ured, last_d)))
  last_fapx <- last_fex <- init_fex
  cnt <- 0
  cnt2 <- 0
  reject_cnt <- 0
  last_reject <- TRUE
  rel_gain <- Inf
  relgain_hist <- rep(Inf, control$reltol_steps)
  while ((reject_cnt < maxreject) &&
         ((mincount >= 0 && cnt < mincount && cnt < maxcount) ||
          (mincount < 0 && cnt2 < (-mincount) && cnt < maxcount) ||
          (cnt < maxcount &&
           max(relgain_hist) > abs(reltol) &&
           abs(relgain) > reltol2))) {
    cnt <- cnt + 1
    zprop <- glsalgo(map, zprior, U, obs, zref=zref,
                 adjust_idcs=adjust_idcs, damp=lambda)
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
      lambda <- lambda * 2^5
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
      relgain_hist <- c(tail(relgain_hist, n=(control$reltol_steps-1)), relgain)
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
      cnt2 <- cnt2 + 1
      reject_cnt <- 0
    } else {
      last_reject <- TRUE
      reject_cnt <- reject_cnt + 1
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

