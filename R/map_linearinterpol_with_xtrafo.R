#' Create a linear interpolation mapping including energy calibration
#'
#' Creates a map to linearly interpolate the values at the source indices
#' given on a one-dimensional mesh to the one-dimensional mesh associated with
#' the variables at the target indices.
#' It is possible to apply a shift and scaling to the target mesh to account
#' for, e.g., an energy calibration error of an experiment.
#' The transformation is given by \eqn{x = \alpha + \beta x'} where
#' \eqn{x'} is an x-value of the target mesh as stated by the user and the
#' resulting \eqn{x} is the \emph{correct} x-value that should be used for
#' the linear interpolation.
#'
#' The following fields are required in the parameter list to initialize the mapping:
#' \tabular{ll}{
#' \code{mapname} \tab Name of the mapping \cr
#' \code{maptype} \tab Must be \code{"linearinterpol_with_xtrafo_map"} \cr
#' \code{src_idx} \tab Vector of source indices \cr
#' \code{tar_idx} \tab Vector of target indices \cr
#' \code{src_x} \tab Vector with the mesh associated with the source indices \cr
#' \code{tar_x} \tab Vector with the mesh associated with the target indices \cr
#' \code{zero_outside} \tab Default is \code{FALSE}. If TRUE, y-values of target x-values outside
#'                          the limits of the source mesh will be zero, otherwise this situation
#'                          is not allowed. \cr
#' \code{shiftx_idx} \tab Index associated with the variable that contains \eqn{\alpha}. \cr
#' \code{scalex_idx} \tab Index associated with the variable that contains \eqn{\beta}
#' }
#'
#' @return
#' Returns a list of functions to operate with the mapping, see \code{\link{create_maptype_map}}.
#' @export
#'
#' @family mappings
#' @examples
#' params <- list(
#'   mapname = "mylinearintmap",
#'   maptype = "linearinterpol_with_xtrafo_map",
#'   src_idx = 1:3,
#'   tar_idx = 4:6,
#'   src_x = c(1,5,10),
#'   tar_x = c(4,5,6),
#'   shiftx_idx = 7,
#'   scalex_idx = 8
#' )
#' mymap <- create_linearinterpol_with_xtrafo_map()
#' mymap$setup(params)
#' x <- c(1,2,3,0,0,0,0.5,1.1)
#' mymap$propagate(x)
#' mymap$jacobian(x)
#'
create_linearinterpol_with_xtrafo_map <- function() {

  linmap <- NULL
  xtrafo_params <- NULL
  last_shiftx <- NULL
  last_scalex <- NULL
  last_src_x <- NULL
  last_with.id <- NULL
  energyderiv_coeffs <- NULL


  setup <- function(params) {
    stopifnot(params[["maptype"]] == getType())
    stopifnot(c("shiftx_idx", "scalex_idx") %in% names(params))
    xtrafo_params <<- list(
      shiftx_idx = params[["shiftx_idx"]],
      scalex_idx = params[["scalex_idx"]],
      claimed_tar_x = params[["tar_x"]]
    )
    # create the basic linear interpolation map
    params[["maptype"]] <- "linearinterpol_map"
    params[["shiftx_idx"]] <- NULL
    params[["scalex_idx"]] <- NULL
    linmap <<- create_map(params)
  }


  getType <- function() {
    return("linearinterpol_with_xtrafo_map")
  }


  getName <- function() {
    return(linmap$getName())
  }


  getDescription <- function() {
    return(linmap$getDescription())
  }


  is_linear <- function() {
    return(FALSE)
  }


  get_src_idx <- function() {
    return(c(linmap$get_src_idx(),
             xtrafo_params$shiftx_idx,
             xtrafo_params$scalex_idx))
  }


  get_tar_idx <- function() {
    return(linmap$get_tar_idx())
  }


  propagate <- function(x, with.id=TRUE) {
    update_linmap(x)
    return(linmap$propagate(x, with.id))
  }


  jacobian <- function(x, with.id=TRUE) {
    update_linmap(x)
    tar_idx <- linmap$get_tar_idx()
    S <- linmap$jacobian(x, with.id)
    matidcs_sel <- cbind(rep(tar_idx, 2),
                         c(rep(xtrafo_params$shiftx_idx, length(tar_idx)),
                           rep(xtrafo_params$scalex_idx, length(tar_idx))))
    S[matidcs_sel] <- c(energyderiv_coeffs,
                        energyderiv_coeffs * xtrafo_params$claimed_tar_x)
    return(S)
  }

  # internal utility functions

  update_linmap <- function(x) {
    cur_shiftx <- x[xtrafo_params$shiftx_idx]
    cur_scalex <- x[xtrafo_params$scalex_idx]
    src_idx <- linmap$get_src_idx()
    if (!isTRUE(last_shiftx == cur_shiftx) ||
        !isTRUE(last_scalex == cur_scalex) ||
        !isTRUE(all(x[src_idx] == last_src_x))) {
      # compute target xs derivatives with respect to target energy
      src_x <- linmap$get_src_x()
      tar_x <- linmap$get_tar_x()
      low_idx <- findInterval(tar_x, src_x, rightmost.closed=TRUE)
      high_idx <- low_idx + 1
      stopifnot(all(low_idx >= 1) && all(high_idx <= length(src_x)))
      xdiff <- src_x[high_idx] - src_x[low_idx]
      energyderiv_coeffs <<- (-x[src_idx[low_idx]] + x[src_idx[high_idx]]) / xdiff
      # update the linearinterpolation map
      claimed_tar_x <- xtrafo_params$claimed_tar_x
      new_tar_x <- cur_shiftx + cur_scalex * claimed_tar_x
      linmap$set_tar_x(new_tar_x)
      # keep track of current transformation parameters
      last_shiftx <<- cur_shiftx
      last_scalex <<- cur_scalex
      last_src_x <<- x[src_idx]
      return(TRUE)
    }
    return(FALSE)
  }


  return(list(
    setup = setup,
    getType = getType,
    getName = getName,
    getDescription = getDescription,
    is_linear = is_linear,
    get_src_idx = get_src_idx,
    get_tar_idx = get_tar_idx,
    propagate = propagate,
    jacobian = jacobian
  ))
}
