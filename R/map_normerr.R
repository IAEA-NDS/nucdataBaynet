create_normerr_map <- function() {

  map <- NULL
  S <- NULL
  is_with_id <- FALSE

  setup <- function(params) {
    stopifnot(params$mapname == getName())
    stopifnot(!anyDuplicated(params$src_feat))
    S <<- NULL
    map <<- list(
      mapname = params$mapname,
      src_idx = params$src_idx,
      tar_idx = params$tar_idx,
      src_feat = params$src_feat,
      tar_feat = params$tar_feat
    )
  }

  getName <- function() {
    "normerr_map"
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
    getName = getName,
    get_src_idx = get_src_idx,
    get_tar_idx = get_tar_idx,
    propagate = propagate,
    jacobian = jacobian
  )
}
