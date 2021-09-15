#' Create a mapping from function values to the first derivative
#'
#' Create a mapping from a discretized version of a function with values
#' given on an x-mesh to the first derivatives for a subset of the x-values
#' of the source mesh.
#'
#' The following fields are required in the parameter list to initialize the mapping:
#' \tabular{ll}{
#' \code{mapname} \tab Name of the mapping \cr
#' \code{maptype} \tab Must be \code{"derivative_map"} \cr
#' \code{src_idx} \tab Vector of source indices \cr
#' \code{tar_idx} \tab Vector of target indices. \cr
#' \code{src_x} \tab Vector of x-values of the source mesh \cr
#' \code{tar_x} \tab Vector of x-values of the target mesh.
#'                   Must be a subset of the x-values of the source mesh
#'                   and must not include the largest x-value
#'                   of the source mesh.
#' }
#'
#' \loadmathjax
#' Let \mjseqn{x_i} denote the x-values of the mesh and \mjseqn{y_i} the
#' associated y-vales of the function. A finite-difference approximation of
#' the first derivative is given by
#' \mjsdeqn{
#'   \Delta_i = \frac{\vec{y}_{i+1} - \vec{y}_i}{\vec{x}_{i+1}-\vec{x}_i}
#' }
#'
#'
#' @return
#' Returns a list of functions to operate with the mapping, see \code{\link{create_maptype_map}}.
#' @export
#' @family mappings
#'
#' @examples
#' params <- list(
#'   mapname = "myderivmap",
#'   maptype = "derivative_map",
#'   src_idx = 1:4,
#'   tar_idx = 5:6,
#'   src_x = 1:4,
#'   tar_x = 2:3
#' )
#' mymap <- create_derivative_map()
#' mymap$setup(params)
#' x <- c(1,2,4,7,0,0)
#' mymap$propagate(x)
#' mymap$jacobian(x)
#'
create_derivative_map <- function() {

  map <- NULL
  S <- NULL
  is_with_id <- FALSE

  setup <- function(params) {
    stopifnot(all(c('maptype','src_idx','tar_idx','src_x','tar_x') %in% names(params)))
    stopifnot(params$maptype == getType())
    src_ord <- order(params$src_x)
    tar_ord <- order(params$tar_x)
    map <<- list(
      maptype = params$maptype,
      mapname = params$mapname,
      description = params$description,
      src_idx = params$src_idx[src_ord],
      tar_idx = params$tar_idx[tar_ord],
      src_x = params$src_x[src_ord],
      tar_x = params$tar_x[tar_ord]
    )
    idx <- findInterval(map$tar_x, map$src_x)
    stopifnot(all(map$src_x[idx]==map$tar_x))
    stopifnot(all(map$src_idx[idx] < max(map$src_idx)))
  }


  getType <- function() {
    return("derivative_map")
  }


  getName <- function() {
    return(map[["mapname"]])
  }


  getDescription <- function() {
    return(map[["description"]])
  }


  is_linear <- function() {
    return(TRUE)
  }


  get_src_idx <- function() {
    return(map$src_idx)
  }

  get_tar_idx <- function() {
    return(map$tar_idx)
  }

  propagate <- function(x, with.id=TRUE) {
    S <- jacobian(x, with.id)
    return(as.vector(S %*% x))
  }

  jacobian <- function(x, with.id=TRUE) {
    if (is.null(S) || with.id != is_with_id) {
      idx1 <- findInterval(map$tar_x, map$src_x)
      idx2 <- idx1 + 1
      delta <- map$src_x[idx2] - map$src_x[idx1]
      coeff1 <- (-1/delta)
      coeff2 <- 1/delta
      S <<- sparseMatrix(
        i = rep(map$tar_idx, 2),
        j = map$src_idx[c(idx1,idx2)],
        x = c(coeff1, coeff2),
        dims = rep(length(x), 2)
      )
      if (with.id) {
        diag(S) <- diag(S) + 1
      }
    }
    return(S)
  }

  list(
    setup = setup,
    getType = getType,
    getName = getName,
    getDescription = getDescription,
    is_linear = is_linear,
    get_src_idx = get_src_idx,
    get_tar_idx = get_tar_idx,
    propagate = propagate,
    jacobian = jacobian
  )
}
