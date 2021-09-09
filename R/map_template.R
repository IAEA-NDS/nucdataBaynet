#' Template for the creation of a mapping function
#'
#'  \code{create_maptype_map()} is a prototype function to construct a mapping skeleton.
#'  Its source code can be used as a template to create new mapping types.
#'
#' @return
#' Returns a named list with the following functions:
#'  \tabular{ll}{
#'  \code{setup(params)} \tab Takes a named list with the parameters required to
#'  set up the mapping. The required fields in the list depend on the mapping,
#'  but most mappings take at least a vector with source indices and a vector
#'  with target indices. \cr
#'  \code{getType()} \tab Returns a string that identifies the type of mapping. \cr
#'  \code{getName()} \tab Returns a string with the name of this specific mapping instance. \cr
#'  \code{getDescription()} \tab Returns a string with a more details about this specific mapping instance. \cr
#'  \code{is_linear()} \tab Must return TRUE if the mapping is linear, otherwise FALSE. \cr
#'  \code{get_src_idx()} \tab Returns the source indices of the mapping. \cr
#'  \code{get_tar_idx()} \tab Returns the target indices of the mapping. \cr
#'  \code{propagate(x, with.id=TRUE)} \tab Propagates the values in the vector \code{x}
#'  according to the mapping and returns the resulting vector. If \code{with.id=FALSE} only
#'  returns the changes caused by the mapping, otherwise it returns a vector with
#'  the changes added to \code{x}. \cr
#'  \code{jacobian(x, with.id=TRUE)} \tab Returns the Jacobian matrix of the mapping
#'  evaluated at \code{x}. If \code{with.id=FALSE}, the Jacobian only reflects the
#'  change of the mapping, otherwise it is augmented by the identity matrix.
#'  }
#'
#' @export
#' @family mappings
#'
#' @examples
#' \dontrun{
#'     # create the list with essential mapping functions
#'     mapping <- create_maptype_map()
#'     # some prototype parameter list
#'     params <- list(mapname = "somename",
#'                    src_idx = 1:5, tar_idx = 6:10)
#'     # initialize the mapping
#'     mapping$setup(params)
#'     # call some mapping functions
#'     mapping$get_src_idx()
#'     mapping$get_tar_idx()
#' }
create_maptype_map <- function() {

  mapinfo <- NULL

  setup <- function(params) {
    # check parameter list
    stopifnot(is.character(params$name))
    stopifnot(is.numeric(params$src_idx))
    stopifnot(is.numeric(params$tar_idx))
    # store parameters
    mapinfo <<- list(
      mapname = params$name,
      description = params$description,
      src_idx = params$src_idx,
      tar_idx = params$tar_idx,
      maptype = "example_maptype",
      is_linear = TRUE
    )
  }


  getType <- function() {
    return(mapinfo$maptype)
  }


  getName <- function() {
    return(mapinfo$mapname)
  }


  getDescription <- function() {
    return(mapinfo$description)
  }


  is_linear <- function() {
    return(mapinfo$is_linear)
  }


  get_src_idx <- function() {
    return(mapinfo$src_idx)
  }


  get_tar_idx <- function() {
    return(mapinfo$tar_idx)
  }


  propagate <- function(x, with.id=TRUE) {
    # TODO
  }


  jacobian <- function(x, with.id=TRUE) {
    # TODO
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
