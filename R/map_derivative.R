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
