create_relativemap_map <- function() {

  mapinfo <- list()
  basemap <- NULL

  setup <- function(params) {
    stopifnot(params$maptype == "relativemap_map")
    stopifnot(!is.null(params$basemap))
    basemap <<- create_map(params$basemap)
    stopifnot(!anyDuplicated(c(basemap$get_src_idx(), basemap$get_tar_idx())))

    mapinfo <<- list(
      mapname = params[["mapname"]],
      description = params[["description"]],
      src_idx = sort(c(basemap$get_src_idx(), basemap$get_tar_idx())),
      tar_idx = basemap$get_tar_idx()
    )
  }


  getType <- function() {
    return("relativemap_map")
  }


  getName <- function() {
    return(mapinfo[["mapname"]])
  }


  getDescription <- function() {
    return(mapinfo[["description"]])
  }


  is_linear <- function() {
    return(FALSE)
  }


  get_src_idx <- function() {
    return(mapinfo[["src_idx"]])
  }


  get_tar_idx <- function() {
    return(mapinfo[["tar_idx"]])
  }


  propagate <- function(x, with.id=TRUE) {
    src_idx <- basemap$get_src_idx()
    tar_idx <- basemap$get_tar_idx()
    if (!with.id) {
      res <- rep(0, length(x))
    } else {
      res <- x
    }
    subres <- x[tar_idx] * (basemap$propagate(x, with.id)[tar_idx])
    res[tar_idx] <- res[tar_idx] + subres
    return(res)
  }


  jacobian <- function(x, with.id=TRUE) {
    src_idx <- basemap$get_src_idx()
    tar_idx <- basemap$get_tar_idx()

    propx <- basemap$propagate(x, with.id=FALSE)[tar_idx]
    S <- basemap$jacobian(x, with.id=FALSE)

    S[tar_idx, src_idx] <- S[tar_idx, src_idx] * x[tar_idx]
    S[tar_idx, tar_idx] <- Diagonal(n = length(tar_idx), x = propx)
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
