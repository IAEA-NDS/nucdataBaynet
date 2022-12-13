#' Create a linear interpolation mapping
#'
#' Creates a map to linearly interpolate the values at the source indices
#' given on a one-dimensional mesh to the one-dimensional mesh associated with
#' the variables at the target indices.
#'
#' The following fields are required in the parameter list to initialize the mapping:
#' \tabular{ll}{
#' \code{mapname} \tab Name of the mapping \cr
#' \code{maptype} \tab Must be \code{"linearinterpol_map"} \cr
#' \code{src_idx} \tab Vector of source indices \cr
#' \code{tar_idx} \tab Vector of target indices \cr
#' \code{src_x} \tab Vector with the mesh associated with the source indices \cr
#' \code{tar_x} \tab Vector with the mesh associated with the target indices \cr
#' \code{zero_outside} \tab Default is \code{FALSE}. If TRUE, y-values of target x-values outside
#'                          the limits of the source mesh will be zero, otherwise this situation
#'                          is not allowed. \cr
#' \code{scalefact} \tab A scalar that is multiplied with the resulting numbers from linear interpolation
#'                       to obtain the final result of the mapping. Default is 1.
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
#'   maptype = "linearinterpol_map",
#'   src_idx = 1:3,
#'   tar_idx = 4:6,
#'   src_x = c(1,5,10),
#'   tar_x = c(4,5,6)
#' )
#' mymap <- create_linearinterpol_map()
#' mymap$setup(params)
#' x <- c(1,2,3,0,0,0)
#' mymap$propagate(x)
#' mymap$jacobian(x)
#'
create_linearinterpol_map <- function() {

  map <- NULL
  S <- NULL

  setup <- function(params) {
    stopifnot(params$maptype == getType())
    stopifnot(all(c('maptype','src_idx','tar_idx','src_x','tar_x') %in% names(params)))
    stopifnot(length(params$src_idx) > 0)
    stopifnot(length(params$tar_idx) > 0)
    stopifnot(length(params$src_idx) == length(params$src_x))
    stopifnot(length(params$tar_idx) == length(params$tar_x))
    stopifnot(is.null(params$scalefact) || length(params$scalefact)==1)
    src_ord <- order(params$src_x)
    map <<- list(
      maptype = params$maptype,
      mapname = params$mapname,
      description = params$description,
      src_idx = params$src_idx[src_ord],
      tar_idx = params$tar_idx,
      src_x = params$src_x[src_ord],
      tar_x = params$tar_x,
      zero_outside = params$zero_outside,
      scalefact = if (!is.null(params$scalefact)) params$scalefact else 1
    )
    if (!isTRUE(map$zero_outside)) {
      check_mesh(map$tar_x)
    }
  }

  getType <- function() {
    return("linearinterpol_map")
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

  propagate <- function(x, with.id=TRUE) {
    if (is.null(S)) {
      update_jacobian(length(x))
    }
    if (isTRUE(with.id)) {
      return(x + as.vector(S %*% x))
    } else {
      return(as.vector(S %*% x))
    }
  }

  jacobian <- function(x, with.id=TRUE) {
    if (is.null(S)) {
      update_jacobian(length(x))
    }
    if (isTRUE(with.id)) {
      return(S + Diagonal(n=length(x), x=1))
    } else {
      return(S)
    }
  }

  # internal functions

  check_mesh <- function(tar_x) {
    if (min(tar_x) < map$src_x[1]) {
      stop(paste0("some target x are smaller than the lowest x of the mesh",
                  " (e.g., x = ", min(tar_x), ")"))
    }
    if (max(tar_x) > tail(map$src_x, n=1)) {
      stop(paste0("some target x are larger than the largest x of the mesh",
                  " (e.g., x = ", max(tar_x), ")"))
    }
  }


  update_jacobian <- function(dim) {

    idx1 <- findInterval(map$tar_x, map$src_x, rightmost.closed=TRUE)
    idx2 <- idx1 + 1
    # special treatment of mesh with one point
    idrowsel <- idx1==0 & map$tar_x == map$src_x[idx2]
    gidrowidx <- map$tar_idx[idrowsel]
    gidcolidx <- map$src_idx[idx2[idrowsel]]
    idx1 <- idx1[!idrowsel]
    idx2 <- idx2[!idrowsel]

    i <- integer(0)
    j <- integer(0)
    coeff <- numeric(0)
    if (length(idx1) > 0) {
      if (isTRUE(map$zero_outside)) {
        sel <- idx1 >= 1 & idx2 <= length(map$src_x)
        idx1 <- idx1[sel]
        idx2 <- idx2[sel]
        tar_idx <- map$tar_idx[sel]
        xp <- map$tar_x[sel]
      } else {
        tar_idx <- map$tar_idx
        xp <- map$tar_x
      }
      gidx1 <- map$src_idx[idx1]
      gidx2 <- map$src_idx[idx2]
      x1 <- map$src_x[idx1]
      x2 <- map$src_x[idx2]
      delta <- x2 - x1
      i <- rep(tar_idx, 2)
      j <- c(gidx1, gidx2)
      coeff <- c((x2-xp)/delta, (xp-x1)/delta)
    }
    # add the cases with one-point meshes
    i <- c(i, gidrowidx)
    j <- c(j, gidcolidx)
    coeff <- c(coeff, rep(1, length(gidrowidx)))
    # apply scaling factor
    coeff <- coeff * map$scalefact
    # create the jacobian
    S <<- sparseMatrix(i = i, j = j, x = coeff,
                       dims = rep(dim,2))
  }

  # additional functions not part of standard interface

  get_src_idx <- function() {
    return(map$src_idx)
  }

  get_tar_idx <- function() {
    return(map$tar_idx)
  }

  get_src_x <- function() {
    return(map$src_x)
  }

  get_tar_x <- function() {
    return(map$tar_x)
  }

  set_tar_x <- function(tar_x) {
    if (!isTRUE(map$zero_outside)) {
      check_mesh(tar_x)
    }
    map[["tar_x"]] <<- tar_x
    S <<- NULL
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
    jacobian = jacobian,
    # additional functions not part of standard interface
    get_src_idx = get_src_idx,
    get_tar_idx = get_tar_idx,
    get_src_x = get_src_x,
    get_tar_x = get_tar_x,
    set_tar_x = set_tar_x
  )
}
