#' Create a generic linear mapping
#'
#' Create a linear mapping between the variables at the source indices
#' to the variables at the target indices, which is of the form
#' \deqn{
#'   \vec{y} = \vec{y}_{\textrm{ref}} + S (\vec{p} - \vec{p}_{\textrm{ref}})
#' }
#'
#' The following fields are required in the parameter list to initialize the mapping:
#' \tabular{ll}{
#' \code{mapname} \tab Name of the mapping \cr
#' \code{maptype} \tab Must be \code{"linearmap_map"} \cr
#' \code{src_idx} \tab Vector of source indices \cr
#' \code{tar_idx} \tab Vector of target indices. \cr
#' \code{yref} \tab The vector \eqn{\vec{y}_{\textrm{ref}}} given in the formula. \cr
#' \code{pref} \tab The vector \eqn{\vec{p}_{\textrm{ref}}} given in the formula. \cr
#' \code{S} \tab The matrix \eqn{S} given in the formula.
#' }
#'
#' @return
#' Returns a list of functions to operate with the mapping, see \code{\link{create_maptype_map}}.
#' @export
#'
#' @family mappings
#' @examples
#' params <- list(
#'   mapname = "mylinearmap",
#'   maptype = "linearmap_map",
#'   src_idx = 1:2,
#'   tar_idx = 3:5,
#'   yref = rep(10, 3),
#'   pref = c(2, 4),
#'   S = matrix(c(1, 2, 3, 4, 5, 6), nrow=3, ncol=2)
#' )
#' mymap <- create_linearmap_map()
#' mymap$setup(params)
#' x <- c(1,-2,0,0,0)
#' mymap$propagate(x)
#' mymap$jacobian(x)
#'
create_linearmap_map <- function() {

  mapinfo <- NULL
  Sfull <- NULL

  setup <- function(params) {
    stopifnot(!is.null(params[["maptype"]]))
    stopifnot(all(c("src_idx","tar_idx") %in% names(params)))
    stopifnot(all(c("S","yref","pref") %in% names(params)))
    stopifnot(length(params[["tar_idx"]]) == nrow(params[["S"]]))
    stopifnot(length(params[["src_idx"]]) == ncol(params[["S"]]))
    stopifnot(length(params[["yref"]]) == length(params[["tar_idx"]]))
    stopifnot(length(params[["pref"]]) == length(params[["src_idx"]]))
    mapinfo <<- list(
      mapname = params[["mapname"]],
      description = params[["description"]],
      src_idx = params[["src_idx"]],
      tar_idx = params[["tar_idx"]],
      S = params[["S"]],
      yref = params[["yref"]],
      pref = params[["pref"]]
    )
  }


  getType <- function() {
    return("linearmap_map")
  }


  getName <- function() {
    return(mapinfo[["mapname"]])
  }


  getDescription <- function() {
    return(mapinfo[["description"]])
  }


  is_linear <- function() {
    return(TRUE)
  }


  get_src_idx <- function() {
    return(mapinfo[["src_idx"]])
  }


  get_tar_idx <- function() {
    return(mapinfo[["tar_idx"]])
  }


  propagate <- function(x, with.id=TRUE) {
    m <- mapinfo
    if (with.id) res <- x else res <- rep(0, length(x))
    mres <- m$yref + m$S %*% (x[m$src_idx] - mapinfo$pref)
    res[m$tar_idx] <- res[m$tar_idx] + mres
    return(as.vector(res))
  }


  jacobian <- function(x, with.id=TRUE) {
    if (is.null(Sfull)) {
      m <- mapinfo
      Sfull <<- sparseMatrix(i=1, j=1, x=0, dims=rep(length(x),2))
      Sfull[m$tar_idx, m$src_idx] <<- m$S
      diag(Sfull) <- diag(Sfull) + 1
      Sfull <<- drop0(Sfull)
    }
    Sres <- Sfull
    if (!with.id) {
      diag(Sres) <- diag(Sres) - 1
    }
    return(Sres)
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
