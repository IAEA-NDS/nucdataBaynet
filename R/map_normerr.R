#' Create a mapping for absolute energy-independent normalization errors
#'
#' Creates a map to distribute the absolute normalization errors given
#' at the source indices to the experimental data at the target indices.
#'
#' The following fields are required in the parameter list to initialize the mapping:
#' \tabular{ll}{
#' \code{mapname} \tab Name of the mapping \cr
#' \code{maptype} \tab Must be \code{"normerr_map"} \cr
#' \code{src_idx} \tab Vector of source indices \cr
#' \code{tar_idx} \tab Vector of target indices \cr
#' \code{src_feat} \tab Vector with the names of the normalization errors at the
#' source indices. This vector must not contain duplicated strings. \cr
#' \code{tar_feat} \tab Vector with the names of the normalization error associated
#' with the variables at the target indices.
#' }
#'
#' @note
#' This mapping may be used in a product mapping
#' to define \emph{relative} normalization errors,
#' see \code{\link{create_product_mapping}}
#'
#' @return
#' Returns a list of functions to operate with the mapping, see \code{\link{create_maptype_map}}.
#' @export
#'
#' @family mappings
#' @examples
#' params <- list(
#'   mapname = "mymap",
#'   maptype = "normerr_map",
#'   src_idx = 1:3,
#'   tar_idx = 4:9,
#'   src_feat = c("expnorm1", "expnorm2", "expnorm3"),
#'   tar_feat = c(rep("expnorm1", 2), rep("expnorm2", 3),
#'                rep("expnorm3", 1))
#' )
#' mymap <- create_normerr_map()
#' mymap$setup(params)
#' x <- c(1, 5, 7, rep(0, 6))
#' mymap$propagate(x)
#' mymap$jacobian(x)
#'
create_normerr_map <- function() {

  map <- NULL
  S <- NULL
  is_with_id <- FALSE

  setup <- function(params) {
    stopifnot(params$maptype == getType())
    stopifnot(!anyDuplicated(params$src_feat))
    S <<- NULL
    map <<- list(
      maptype = params$maptype,
      mapname = params$mapname,
      description = params$description,
      src_idx = params$src_idx,
      tar_idx = params$tar_idx,
      src_feat = params$src_feat,
      tar_feat = params$tar_feat
    )
  }

  getType <- function() {
    "normerr_map"
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
    map$src_idx
  }

  get_tar_idx <- function() {
    map$tar_idx
  }

  propagate <- function(x, with.id=TRUE) {
    S <- jacobian(x, with.id)
    as.vector(S %*% x)
  }

  jacobian <- function(x, with.id=TRUE) {
    if (is.null(S) || is_with_id != with.id) {
      idxassoc <- match(map$tar_feat, map$src_feat)
      S <<- sparseMatrix(i=map$tar_idx,
                        j=map$src_idx[idxassoc], x=1,
                        dims=rep(length(x),2))
      if (with.id) {
        diag(S) <<- diag(S) + 1
      }
      is_with_id <<- with.id
    }
    S
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
