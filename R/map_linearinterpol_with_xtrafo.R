create_linearinterpol_with_xtrafo_map <- function() {

  linmap <- NULL
  xtrafo_params <- NULL
  last_shiftx <- NULL
  last_scalex <- NULL
  last_with.id <- NULL
  energyderiv_coeffs <- NULL


  setup <- function(params) {
    stopifnot(params[["maptype"]] == getType())
    stopifnot(c("shiftx_idx", "scalex_idx") %in% names(params))
    xtrafo_params <<- list(
      shiftx_idx = params[["shiftx_idx"]],
      scalex_idx = params[["scalex_idx"]],
      claimed_tar_x = params[["tar_x"]]
    )
    # create the basic linear interpolation map
    params[["maptype"]] <- "linearinterpol_map"
    params[["shiftx_idx"]] <- NULL
    params[["scalex_idx"]] <- NULL
    linmap <<- create_map(params)
  }


  getType <- function() {
    return("linearinterpol_with_xtrafo_map")
  }


  getName <- function() {
    return(linmap$getName())
  }


  getDescription <- function() {
    return(linmap$getDescription())
  }


  is_linear <- function() {
    return(FALSE)
  }


  get_src_idx <- function() {
    return(c(linmap$get_src_idx(),
             xtrafo_params$shiftx_idx,
             xtrafo_params$scalex_idx))
  }


  get_tar_idx <- function() {
    return(linmap$get_tar_idx())
  }


  propagate <- function(x, with.id=TRUE) {
    update_linmap(x)
    return(linmap$propagate(x, with.id))
  }


  jacobian <- function(x, with.id=TRUE) {
    is_updated <- update_linmap(x)
    if (is_updated) {
      src_x <- linmap$get_src_x()
      tar_x <- linmap$get_tar_x()
      low_idx <- findInterval(tar_x, src_x, rightmost.closed=TRUE)
      high_idx <- low_idx + 1
      stopifnot(all(low_idx >= 1) && all(high_idx <= length(src_x)))
      xdiff <- src_x[high_idx] - src_x[low_idx]
      energyderiv_coeffs <<- (-x[low_idx] + x[high_idx]) / xdiff
    }
    tar_idx <- linmap$get_tar_idx()
    S <- linmap$jacobian(x, with.id)
    matidcs_sel <- cbind(rep(tar_idx, 2),
                         c(rep(xtrafo_params$shiftx_idx, length(tar_idx)),
                           rep(xtrafo_params$scalex_idx, length(tar_idx))))
    S[matidcs_sel] <- c(energyderiv_coeffs,
                        energyderiv_coeffs * xtrafo_params$claimed_tar_x)
    return(S)
  }

  # internal utility functions

  update_linmap <- function(x) {
    cur_shiftx <- x[xtrafo_params$shiftx_idx]
    cur_scalex <- x[xtrafo_params$scalex_idx]
    if (!isTRUE(last_shiftx == cur_shiftx) ||
        !isTRUE(last_scalex == cur_scalex)) {
      claimed_tar_x <- xtrafo_params$claimed_tar_x
      linmap$set_tar_x(cur_shiftx + cur_scalex * claimed_tar_x)
      last_shiftx <<- cur_shiftx
      last_scalex <<- cur_scalex
      return(TRUE)
    }
    return(FALSE)
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
