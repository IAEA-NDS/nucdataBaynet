#' Create a mapping for windowed averages including energy calibration
#'
#' Creates a map to calulate averages of a linear piecewise function whose
#' functions values are given at the source indices.
#' The averages are evaluated at the x-values of the mesh associated with the
#' target indices.
#' It is also possible to apply the transformation \eqn{x = \alpha + (1-\beta)x'}
#' to the x'-values of the target mesh before the averaging is performed.
#'
#' The following fields are required in the parameter list to initialize the mapping:
#' \tabular{ll}{
#' \code{mapname} \tab Name of the mapping \cr
#' \code{maptype} \tab Must be \code{"convolution_with_xtrafo_map"} \cr
#' \code{src_idx} \tab Vector of source indices \cr
#' \code{tar_idx} \tab Vector of target indices \cr
#' \code{src_x} \tab Vector with the mesh associated with the source indices \cr
#' \code{tar_x} \tab Vector with the mesh associated with the target indices \cr
#' \code{winsize} \tab The average of the function at the source mesh is taken
#' between x-winsize/2 and x+winsize/2 \cr
#' \code{winsize_idx} \tab The index of the variable containing the window size.
#' Only one of \code{winsize} and \code{winsize_idx} must be present. \cr
#' \code{shiftx} \tab Value of \eqn{\alpha} \cr
#' \code{shiftx_idx} \tab Index associated with the variable that contains \eqn{\alpha}.
#' Only one of \code{shiftx} and \code{shiftx_idx} must be present. \cr
#' \code{scalex} \tab Value of \eqn{\beta}. \cr
#' \code{scalex_idx} \tab Index associated with the variable that contains \eqn{\beta}.
#' Only one of \code{scalex} and \code{scalex_idx} must be present.
#' }
#'
#' @return
#' Returns a list of functions to operate with the mapping, see \code{\link{create_maptype_map}}.
#' @export
#'
#' @family mappings
#' @examples
#' params <- list(
#'   mapname = "mylinearintmap",
#'   maptype = "convolution_with_xtrafo_map",
#'   src_idx = 1:3,
#'   tar_idx = 4:6,
#'   src_x = c(1,5,10),
#'   tar_x = c(4,5,6),
#'   winsize_idx = 7,
#'   shiftx_idx = 8,
#'   scalex_idx = 9
#' )
#' mymap <- create_convolution_with_xtrafo_map()
#' mymap$setup(params)
#' x <- c(1,2,3,0,0,0,2,1,0.1)
#' mymap$propagate(x)
#' mymap$jacobian(x)
#'
create_convolution_with_xtrafo_map <- function() {

  map_params <- NULL

  S <- NULL
  last_winsize <- NA
  last_shiftx <- NA
  last_scalex <- NA
  last_src_y <- NA
  last_tar_x <- NA
  last_prop_y <- NA

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
      return(as.vector(x + last_prop_y))
    } else {
      return(last_prop_y)
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
    winsize <- last_winsize
    tar_x <- last_tar_x
    src_x <- map_params[["src_x"]]
    tar_x_min <- tar_x - last_winsize / 2
    tar_x_max <- tar_x + last_winsize / 2
    Fval <- last_prop_y[tar_idx]
    if (!is.null(winsize_idx)) {
      col_idcs <- c(col_idcs, winsize_idx)
      deriv <- (dinfo$flo + dinfo$fhi) / 2 / winsize
      deriv <- deriv - Fval / winsize
      mat_vals <- cbind(mat_vals, deriv)
    }
    shiftx_idx <- map_params[["shiftx_idx"]]
    if (!is.null(shiftx_idx)) {
      col_idcs <- c(col_idcs, shiftx_idx)
      deriv <- (-dinfo$flo + dinfo$fhi) / winsize
      mat_vals <- cbind(mat_vals, deriv)
    }
    scalex_idx <- map_params[["scalex_idx"]]
    if (!is.null(scalex_idx)) {
      col_idcs <- c(col_idcs, scalex_idx)
      orig_tar_x <- map_params[["orig_tar_x"]]
      deriv <- (-dinfo$flo*orig_tar_x + dinfo$fhi*orig_tar_x) / winsize
      mat_vals <- cbind(mat_vals, deriv)
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
    src_y_changed <- !isTRUE(all(last_src_y == cur_src_y))
    if (is.null(S) || (need_energy_derivatives && (
      !isTRUE(last_winsize == cur_winsize) ||
      !isTRUE(last_shiftx == cur_shiftx) ||
      !isTRUE(last_scalex == cur_scalex) ||
      src_y_changed))) {

      src_x <- map_params[["src_x"]]
      src_idx <- map_params[["src_idx"]]
      new_tar_x <- cur_shiftx + (1+cur_scalex) * map_params[["orig_tar_x"]]
      tryCatch({
      S <<- get_convolute_rect_matrix(src_x, new_tar_x, cur_winsize,
                                      src_idx = map_params[["src_idx"]],
                                      tar_idx = get_tar_idx(),
                                      dims = rep(length(x), 2))
      }, error = function(e) {
        e$message <- paste0(e$message, " (in map with name ", getName(),
                            " and type ", getType(), ")")
        stop(e)
      })
      if (need_energy_derivatives) {
        energy_deriv_info <<- get_convolute_rect_derivative(src_x, x[src_idx],
                                                            new_tar_x, cur_winsize)
      }
      last_winsize <<- cur_winsize
      last_shiftx <<- cur_shiftx
      last_scalex <<- cur_scalex
      last_src_y <<- cur_src_y
      last_tar_x <<- new_tar_x
      last_prop_y <<- as.vector(S %*% x)
    } else if (src_y_changed) {
      last_prop_y <<- as.vector(S %*% x)
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
