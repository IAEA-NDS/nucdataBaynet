#' Create a product mapping
#'
#' Create a mapping that is formed by combining two mappings so that the output
#' of the combined mapping is the product of the outputs of the individual mappings.
#' For example, a product mapping can be used to define relative normalization errors.
#'
#' The following fields are required in the parameter list to initialize the mapping:
#' \tabular{ll}{
#' \code{mapname} \tab Name of the mapping \cr
#' \code{maptype} \tab Must be \code{"product_map"} \cr
#' \code{maps} \tab A list containing the parameter specifications of the
#' individual mappings.
#' }
#'
#' @family mappings
#' @return
#' Returns a list of functions to operate with the mapping, see \code{\link{create_maptype_map}}.
#' @export
#'
#' @examples
#' params1 <- list(
#'   mapname = "mymap1",
#'   maptype = "linearinterpol_map",
#'   src_idx = 1:3,
#'   tar_idx = 4:6,
#'   src_x = c(1,5,10),
#'   tar_x = c(4,5,6)
#' )
#' params2 <- list(
#'   mapname = "mymap1",
#'   maptype = "linearinterpol_map",
#'   src_idx = 7:9,
#'   tar_idx = 4:6,
#'   src_x = c(1,5,10),
#'   tar_x = c(4,5,6)
#' )
#' product_params <- list(
#'   mapname = "myprodmap",
#'   maptype = "product_map",
#'   maps = list(params1, params2)
#' )
#' mymap <- create_product_map()
#' mymap$setup(product_params)
#' x <- c(1,2,3,0,0,0,6,7,8)
#' mymap$propagate(x)
#' mymap$jacobian(x)
#'
create_product_map <- function() {

  map <- NULL

  setup <- function(params) {
    stopifnot(params$maptype == getType())
    stopifnot(!is.null(params$maps))
    map <<- list(
      maptype = params$maptype,
      mapname = params$mapname,
      description = params$description,
      map_list = lapply(params$maps, create_map)
    )
    map[["src_idx"]] <<- sort(unique(unlist(lapply(map$map_list, function(x) x$get_src_idx()))))
    map[["tar_idx"]] <<- sort(unique(unlist(lapply(map$map_list, function(x) x$get_tar_idx()))))
    stopifnot(anyDuplicated(c(map[["src_idx"]], map[["tar_idx"]])) == 0)
  }


  getType <- function() {
    return("product_map")
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
    return(map[["src_idx"]])
  }


  get_tar_idx <- function() {
    return(map[["tar_idx"]])
  }


  propagate <- function(x, with.id=TRUE) {
    res <- 1
    for (curmap in map[["map_list"]]) {
      res <- res * curmap$propagate(x, with.id=FALSE)
    }
    if (with.id) {
      res <- res + x
    }
    return(res)
  }


  jacobian <- function(x, with.id=TRUE) {
    t <- lapply(map$map_list, function(curmap) {
      tmpres <- curmap$propagate(x, with.id=FALSE)
      return(as(tmpres, "sparseVector"))
    })
    compS <- 0
    for (idx in seq_along(map$map_list)) {
      curmap <- map$map_list[[idx]]
      curS <- curmap$jacobian(x, with.id=FALSE)
      tprod <- 1
      for (tc in t[-idx]) {
        tprod <- tprod * tc
      }
      tprod <- as.vector(tprod)
      compS <- compS + curS * tprod
    }
    if (with.id) {
      compS <- compS + Diagonal(n = length(x), x = 1)
    }
    return(compS)
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
