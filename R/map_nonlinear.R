create_nonlinear_map <- function() {

  map <- NULL

  setup <- function(params) {
    stopifnot(params$maptype == getType())
    stopifnot(length(params$src_idx) == length(params$tar_idx))
    stopifnot(!is.unsorted(params$src_idx))
    stopifnot(!is.unsorted(params$tar_idx))
    stopifnot(params$funname == 'exp')
    map <<- list(maptype = params$maptype,
                 mapname = params$mapname,
                 description = params$description,
                 src_idx = params$src_idx,
                 tar_idx = params$tar_idx)
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


  get_src_idx <- function() {
    return(map$src_idx)
  }


  get_tar_idx <- function() {
    return(map$tar_idx)
  }


  propagate <- function(x, with.id=TRUE) {
    res <- if (with.id) x else rep(0, length(x))
    res[map$tar_idx] <- res[map$tar_idx] + exp(x[map$src_idx])
    return(res)
  }


  jacobian <- function(x, with.id=TRUE) {
    S <- sparseMatrix(i = map$tar_idx, j = map$src_idx,
                      x = exp(x[map$src_idx]), dims = rep(length(x), 2))
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
    get_src_idx = get_src_idx,
    get_tar_idx = get_tar_idx,
    propagate = propagate,
    jacobian = jacobian
  ))
}
