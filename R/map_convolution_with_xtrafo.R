create_convolution_with_xtrafo_map <- function() {

  map_params <- NULL

  S <- NULL
  flo <- NULL
  fhi <- NULL
  last_winsize <- NULL
  last_shiftx <- NULL
  last_scalex <- NULL
  last_src_y <- NULL
  need_energy_derivatives <- FALSE
  energy_deriv_info <- NULL

  setup <- function(params) {
    stopifnot(params$maptype == getType())
    stopifnot(c("src_x", "tar_x", "src_idx", "tar_idx") %in% names(params))
    stopifnot(xor(is.null(params[["winsize"]]), is.null(params[["winsize_idx"]])))
    stopifnot(xor(is.null(params[["shiftx"]]), is.null(params[["shiftx_idx"]])))
    stopifnot(xor(is.null(params[["scalex"]]), is.null(params[["scalex_idx"]])))
    stopifnot(anyDuplicated(params[["src_x"]]) == 0)
    src_order <- order(params$src_x)
    map_params <<- list(
      maptype = params$maptype,
      mapname = params$mapname,
      description = params$description,
      src_idx = params$src_idx[src_order],
      tar_idx = params$tar_idx,
      src_x = params$src_x[src_order],
      orig_tar_x = params$tar_x
    )
    map_params[["shiftx_idx"]] <<- params[["shiftx_idx"]]
    map_params[["shiftx"]] <<- params[["shiftx"]]
    map_params[["scalex_idx"]] <<- params[["scalex_idx"]]
    map_params[["scalex"]] <<- params[["scalex"]]
    map_params[["winsize_idx"]] <<- params[["winsize_idx"]]
    map_params[["winsize"]] <<- params[["winsize"]]
    need_energy_derivatives <<- !is.null(map_params[["shiftx_idx"]]) ||
      !is.null(map_params[["scalex_idx"]]) || !is.null(map_params[["winsize_idx"]])
  }


  getType <- function() {
    return("convolution_with_xtrafo_map")
  }


  getName <- function() {
    return(map_params[["mapname"]])
  }


  getDescription <- function() {
    return(map_params[["description"]])
  }


  is_linear <- function() {
    # only linear if energy related parameters winsize, shiftx, scalex are fixed
    return(!need_energy_derivatives)
  }


  get_src_idx <- function() {
    return(c(map_params[["src_idx"]],
             map_params[["shiftx_idx"]],
             map_params[["scalex_idx"]],
             map_params[["winsize_idx"]]))
  }


  get_tar_idx <- function() {
    return(map_params[["tar_idx"]])
  }


  propagate <- function(x, with.id=TRUE) {
    update_map(x)
    if (with.id) {
      return(as.vector(x + S %*% x))
    } else {
      return(as.vector(S %*% x))
    }
  }


  jacobian <- function(x, with.id=TRUE) {
    update_map(x)
    row_idcs <- NULL
    col_idcs <- NULL
    mat_vals <- NULL
    tar_idx <- map_params[["tar_idx"]]
    dinfo <- energy_deriv_info
    winsize_idx <- map_params[["winsize_idx"]]
    if (!is.null(winsize_idx)) {
      col_idcs <- c(col_idcs, winsize_idx)
      mat_vals <- cbind(mat_vals, (dinfo$flo + dinfo$fhi) / 2)
    }
    shiftx_idx <- map_params[["shiftx_idx"]]
    if (!is.null(shiftx_idx)) {
      col_idcs <- c(col_idcs, shiftx_idx)
      mat_vals <- cbind(mat_vals, (-dinfo$flo + dinfo$fhi))
    }
    scalex_idx <- map_params[["scalex_idx"]]
    if (!is.null(scalex_idx)) {
      col_idcs <- c(col_idcs, scalex_idx)
      orig_tar_x <- map_params[["orig_tar_x"]]
      mat_vals <- cbind(mat_vals, (-dinfo$flo*orig_tar_x + dinfo$fhi*orig_tar_x))
    }
    Scur <- S
    if (!is.null(mat_vals)) {
      Scur[tar_idx, col_idcs] <- mat_vals
    }
    if (with.id) {
      diag(Scur) <- diag(Scur) + 1
    }
    return(Scur)
  }

  # internal functions

  update_map <- function(x) {
    # take either hard-coded value during setup or value in x specified by *_idx params
    cur_winsize <- if (!is.null(map_params[["winsize"]]))
      map_params[["winsize"]] else x[map_params[["winsize_idx"]]]
    cur_shiftx <- if (!is.null(map_params[["shiftx"]]))
      map_params[["shiftx"]] else x[map_params[["shiftx_idx"]]]
    cur_scalex <- if (!is.null(map_params[["scalex"]]))
      map_params[["scalex"]] else x[map_params[["scalex_idx"]]]
    tar_idx <- get_tar_idx()
    cur_src_y <- x[map_params[["src_idx"]]]

    # any change since last call
    if (is.null(S) || (need_energy_derivatives && (
      !isTRUE(last_winsize == cur_winsize) ||
      !isTRUE(last_shiftx == cur_shiftx) ||
      !isTRUE(last_scalex == cur_scalex) ||
      !isTRUE(all(last_src_y == cur_src_y))))) {

      src_x <- map_params[["src_x"]]
      src_idx <- map_params[["src_idx"]]
      new_tar_x <- cur_shiftx + (1+cur_scalex) * map_params[["orig_tar_x"]]
      S <<- get_convolute_rect_matrix(src_x, new_tar_x, cur_winsize,
                                      src_idx = map_params[["src_idx"]],
                                      tar_idx = get_tar_idx(),
                                      dims = rep(length(x), 2))
      if (need_energy_derivatives) {
        energy_deriv_info <<- get_convolute_rect_derivative(src_x, x[src_idx],
                                                            new_tar_x, cur_winsize)
      }
      last_winsize <<- cur_winsize
      last_shiftx <<- cur_shiftx
      last_scalex <<- cur_scalex
      last_src_y <<- cur_src_y
    }
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
