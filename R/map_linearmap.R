create_linearmap_map <- function() {

  mapinfo <- NULL
  Sfull <- NULL

  setup <- function(params) {
    stopifnot(!is.null(params[["maptype"]]))
    stopifnot(all(c("src_idx","tar_idx") %in% names(params)))
    stopifnot(all(c("S","yref","pref") %in% names(params)))
    stopifnot(length(params[["tar_idx"]]) == nrow(params[["S"]]))
    stopifnot(length(params[["src_idx"]]) == ncol(params[["S"]]))
    stopifnot(length(params[["yref"]]) == length(params[["tar_idx"]]))
    stopifnot(length(params[["pref"]]) == length(params[["src_idx"]]))
    mapinfo <<- list(
      mapname = params[["mapname"]],
      description = params[["description"]],
      src_idx = params[["src_idx"]],
      tar_idx = params[["tar_idx"]],
      S = params[["S"]],
      yref = params[["yref"]],
      pref = params[["pref"]]
    )
  }


  getType <- function() {
    return("linearmap_map")
  }


  getName <- function() {
    return(mapinfo[["mapname"]])
  }


  getDescription <- function() {
    return(mapinfo[["description"]])
  }


  is_linear <- function() {
    return(TRUE)
  }


  get_src_idx <- function() {
    return(mapinfo[["src_idx"]])
  }


  get_tar_idx <- function() {
    return(mapinfo[["tar_idx"]])
  }


  propagate <- function(x, with.id=TRUE) {
    m <- mapinfo
    if (with.id) res <- x else res <- rep(0, length(x))
    mres <- m$yref + m$S %*% (x[m$src_idx] - mapinfo$pref)
    res[m$tar_idx] <- res[m$tar_idx] + mres
    return(as.vector(res))
  }


  jacobian <- function(x, with.id=TRUE) {
    if (is.null(Sfull)) {
      m <- mapinfo
      Sfull <<- sparseMatrix(i=1, j=1, x=0, dims=rep(length(x),2))
      Sfull[m$tar_idx, m$src_idx] <<- m$S
      diag(Sfull) <- diag(Sfull) + 1
      Sfull <<- drop0(Sfull)
    }
    Sres <- Sfull
    if (!with.id) {
      diag(Sres) <- diag(Sres) - 1
    }
    return(Sres)
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
