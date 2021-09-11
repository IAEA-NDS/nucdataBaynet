#' Create a mapping to effect a non-linear transformation
#'
#' Creates a mapping that transforms the values at the source indices by
#' a predefined non-linear transformation and outpts them at the target
#' indices. Supported transformations at present are \emph{exp}, \emph{relu}
#' and \emph{limiter}.
#'
#' The following fields are required in the parameter list to initialize the mapping:
#' \tabular{ll}{
#' \code{mapname} \tab Name of the mapping \cr
#' \code{maptype} \tab Must be \code{"nonlinear_map"} \cr
#' \code{src_idx} \tab Vector of source indices \cr
#' \code{tar_idx} \tab Vector of target indices.
#'                     Must be of same length as \code{src_idx} \cr
#' \code{funname} \tab A string specifying the non-linear transformation.
#' Can be \code{"exp"}, \code{"relu"} or \code{"limiter"}. \cr
#' \code{minvalue} \tab If \code{funname = "limiter"}, the lower cut-off limit. \cr
#' \code{maxvalue} \tab If \code{funname = "limiter"}, the upper cut-off limit.
#' }
#'
#' @return
#' Returns a list of functions to operate with the mapping, see \code{\link{create_maptype_map}}.
#' @export
#'
#' @family mappings
#' @examples
#' params <- list(
#'   mapname = "mynonlinearmap",
#'   maptype = "nonlinear_map",
#'   src_idx = 1:3,
#'   tar_idx = 4:6,
#'   funname = "relu"
#' )
#' mymap <- create_nonlinear_map()
#' mymap$setup(params)
#' x <- c(1,-2,3,0,0,0)
#' mymap$propagate(x)
#' mymap$jacobian(x)
#'
create_nonlinear_map <- function() {

  map <- NULL
  fun <- NULL
  dfun <- NULL

  setup <- function(params) {
    stopifnot(params$maptype == getType())
    stopifnot(length(params$src_idx) == length(params$tar_idx))
    stopifnot(is.vector(params$src_idx))
    stopifnot(is.vector(params$tar_idx))
    stopifnot(params$funname %in% c("exp", "relu", "limiter"))
    map <<- list(maptype = params$maptype,
                 mapname = params$mapname,
                 description = params$description,
                 src_idx = params$src_idx,
                 tar_idx = params$tar_idx,
                 funname = params$funname,
                 minvalue = params$minvalue, maxvalue = params$maxvalue)
    if (params$funname == "exp") {
      fun <<- exp
      dfun <<- exp
    } else if (params$funname == "relu") {
      fun <<- function(x) { pmax(x, 0) }
      dfun <<- function(x) { as.numeric(x >= 0) }
    } else if (params$funname == "limiter") {
      fun <<- function(x) { pmin(pmax(x, map$minvalue), map$maxvalue) }
      dfun <<- function(x) { as.numeric(x >= map$minvalue & x <= map$maxvalue) }
    } else {
      stop(paste0("function of name ", params$funname, " not supported"))
    }
  }


  getType <- function() {
    return("nonlinear_map")
  }


  getName <- function() {
    return(map[["mapname"]])
  }


  getDescription <- function() {
    return(map[["description"]])
  }


  is_linear <- function() {
    return(FALSE)
  }


  get_src_idx <- function() {
    return(map$src_idx)
  }


  get_tar_idx <- function() {
    return(map$tar_idx)
  }


  propagate <- function(x, with.id=TRUE) {
    res <- if (with.id) x else rep(0, length(x))
    res[map$tar_idx] <- res[map$tar_idx] + fun(x[map$src_idx])
    return(res)
  }


  jacobian <- function(x, with.id=TRUE) {
    S <- sparseMatrix(i = map$tar_idx, j = map$src_idx,
                      x = dfun(x[map$src_idx]), dims = rep(length(x), 2))
    if (with.id) {
      diag(S) <- diag(S) + 1
    }
    return(S)
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
