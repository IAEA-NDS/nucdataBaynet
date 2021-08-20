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
    src_ord <- order(params$src_x)
    map <<- list(
      maptype = params$maptype,
      mapname = params$mapname,
      description = params$description,
      src_idx = params$src_idx[src_ord],
      tar_idx = params$tar_idx,
      src_x = params$src_x[src_ord],
      tar_x = params$tar_x,
      zero_outside = params$zero_outside
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
    if (isTRUE(map$zero_outside)) {
      sel <- idx1 >= 1 & idx2 < length(map$tar_x)
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
