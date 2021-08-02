create_linearinterpol_map <- function() {

  map <- NULL
  S <- NULL
  is_with_id <- FALSE

  setup <- function(params) {
    stopifnot(params$maptype == getType())
    stopifnot(!is.unsorted(params$src_x))
    stopifnot(min(params$src_x) <= min(params$tar_x))
    stopifnot(max(params$src_x) >= max(params$tar_x))
    stopifnot(all(c('maptype','src_idx','tar_idx','src_x','tar_x') %in% names(params)))
    S <<- NULL
    map <<- list(
      maptype = params$maptype,
      mapname = params$mapname,
      description = params$description,
      src_idx = params$src_idx,
      tar_idx = params$tar_idx,
      src_x = params$src_x,
      tar_x = params$tar_x
    )
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
    S <- jacobian(x, with.id)
    return(as.vector(S %*% x))
  }

  jacobian <- function(x, with.id=TRUE) {
    if (is.null(S) || is_with_id != with.id) {
      idx1 <- findInterval(map$tar_x, map$src_x, rightmost.closed=TRUE)
      idx2 <- idx1 + 1
      gidx1 <- map$src_idx[idx1]
      gidx2 <- map$src_idx[idx2]
      x1 <- map$src_x[idx1]
      x2 <- map$src_x[idx2]
      xp <- map$tar_x
      delta <- x2 - x1
      i <- rep(map$tar_idx, 2)
      j <- c(gidx1, gidx2)
      coeff <- c((x2-xp)/delta, (xp-x1)/delta)
      S <<- sparseMatrix(i = i, j = j, x = coeff,
                        dims = rep(length(x),2))
      if (with.id) {
        diag(S) <<- diag(S) + 1
      }
      is_with_id <<- with.id
    }
    return(S)
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
    min_src_x <- map$src_x[1]
    max_src_x <- tail(map$src_x, n=1)
    stopifnot(all(tar_x >= min_src_x) && all(tar_x <= max_src_x))
    S <<- NULL
    map[["tar_x"]] <<- tar_x
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
