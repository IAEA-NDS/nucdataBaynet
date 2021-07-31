create_product_map <- function() {

  map <- NULL

  setup <- function(params) {
    stopifnot(params$maptype == getType())
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
      tprod <- as(tprod, "sparseVector")
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
