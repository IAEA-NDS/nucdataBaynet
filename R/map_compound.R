#' Create a compound mapping
#'
#' Create a compound mapping by merging the effect of individual mappings.
#'
#'
#' The following fields are required in the parameter list to initialize the mapping:
#' \tabular{ll}{
#' \code{mapname} \tab Name of the mapping \cr
#' \code{maptype} \tab Must be \code{"compound_map"} \cr
#' \code{maps} \tab a list with the parameter lists of the individual maps
#' }
#'
#' This mapping will also take care of the correct ordering of the individual
#' mappings.
#' Bundling together individual mapping specifications to a compound map is
#' the last step before doing Bayesian inference.
#'
#' @return
#' Returns a list of functions to operate with the mapping, see \code{\link{create_maptype_map}}.
#' @export
#' @family mappings
#'
#' @examples
#' params1 <- list(
#'   mapname = "mymap1",
#'   maptype = "linearinterpol_map",
#'   src_idx = 1:3,
#'   tar_idx = 7:8,
#'   src_x = 1:3,
#'   tar_x = 2:3
#' )
#' params2 <- list(
#'   mapname = "mymap2",
#'   maptype = "linearinterpol_map",
#'   src_idx = 4:6,
#'   tar_idx = 7:8,
#'   src_x = 1:3,
#'   tar_x = 2:3
#' )
#' compmap_params <- list(
#'   mapname = "mycompmap",
#'   maptype = "compound_map",
#'   maps = list(params1, params2)
#' )
#' mycompmap <- create_compound_map()
#' mycompmap$setup(compmap_params)
#' x <- c(1,2,3,5,5,5,0,0)
#' mycompmap$propagate(x)
#' mycompmap$jacobian(x)
#'
#'
create_compound_map <- function() {

  mapinfo <- list()
  map_list <- list()
  pure_sources <- NULL
  pure_targets <- NULL

  setup <- function(params) {

    stopifnot(!is.null(params$maps))
    maps <- lapply(params$maps, function(curparams) {
      tryCatch(create_map(curparams), error = function(e) {
        e$message <- paste0(e$message, " (during creation of map with name ",
                            curparams$mapname, " of type ", curparams$maptype)
        stop(e)
      })
    })
    max_idx <- 0
    # ensure no duplicates in indices
    for (i in seq_along(maps)) {
      if (anyDuplicated(maps[[i]]$get_src_idx()))
        stop(paste("Duplicate source indices in map",i))
      if (anyDuplicated(maps[[i]]$get_tar_idx()))
        stop(paste("Duplicate target indices in map", j))
      tryCatch({
        max_idx <- max(max_idx, maps[[i]]$get_src_idx(), maps[[i]]$get_tar_idx())
      }, error = function(e) {
        e$message <- paste0(e$message, " (in map with name ", maps[[i]]$getName(),
                            " of type ", maps[[i]]$getType())
        stop(e)
      })
    }
    # determine indices which are pure sources
    is_dest <- rep(FALSE, max_idx)
    is_src <- rep(FALSE, max_idx)
    for (curmap in maps) {
      is_dest[curmap$get_tar_idx()] <- TRUE
      is_src[curmap$get_src_idx()] <- TRUE
    }
    pure_sources <<- which(is_src & !is_dest)
    pure_targets <<- which(is_dest)

    map_list <<- order_maps(maps)
    mapinfo[["mapname"]] <- params[["mapname"]]
    mapinfo[["description"]] <- params[["description"]]
  }


  getType <- function() {
    return("compound_map")
  }


  getName <- function() {
    return(mapinfo[["mapname"]])
  }


  getDescription <- function() {
    return(mapinfo[["description"]])
  }


  getMaps <- function() {
    return(map_list)
  }


  get_src_idx <- function() {
    return(pure_sources)
  }


  get_tar_idx <- function() {
    return(pure_targets)
  }


  get_map_order <- function() {
    return(unlist(lapply(map_list, function(x) x$getName())))
  }


  propagate <- function(x, with.id = TRUE) {
    res <- x
    if (!with.id) {
      res[-pure_sources] <- 0
    }
    initialres <- res
    for (curmap in map_list) {
      cur_src_idx <- curmap$get_src_idx()
      cur_tar_idx <- curmap$get_tar_idx()
      # special treatment of maps with overlapping src and tar idcs
      self_map_mask <- get_selfmap_mask(cur_src_idx, cur_tar_idx, n=length(x))

      res[self_map_mask] <- res[self_map_mask] - initialres[self_map_mask]
      tmpres <- res
      res <- curmap$propagate(res, FALSE)
      res[!self_map_mask] <- res[!self_map_mask] + tmpres[!self_map_mask]
      res[self_map_mask] <- res[self_map_mask] + initialres[self_map_mask]
    }
    if (!with.id) {
      res[pure_sources] <- 0
    }
    return(res)
  }


  jacobian <- function(x, with.id = TRUE) {
    S <- NULL
    if (!with.id) {
      x[-pure_sources] <- 0
    }
    orig_x <- x
    for (curmap in map_list) {
      if (is.null(S))
      {
        if (with.id) {
          S <- curmap$jacobian(x, TRUE)
        } else {
          S <- curmap$jacobian(x, FALSE)
          diag_els <- rep(0, length(x))
          diag_els[pure_sources] <- 1
          S <- S + Diagonal(x = diag_els)
        }
      }
      else
      {
        # self_map_flag <- is_self_map(curmap)
        cur_src_idx <- curmap$get_src_idx()
        cur_tar_idx <- curmap$get_tar_idx()
        selfmap_mask <- get_selfmap_mask(cur_src_idx, cur_tar_idx, n=length(x))
        diag(S)[selfmap_mask] <- 0
        x[selfmap_mask] <- x[selfmap_mask] - orig_x[selfmap_mask]
        curS <- curmap$jacobian(x, TRUE)
        diag(curS)[selfmap_mask] <- diag(curS)[selfmap_mask] - 1
        S <- curS %*% S
        if (with.id) {
          diag(S)[selfmap_mask] <- 1
          x[selfmap_mask] <- x[selfmap_mask] + orig_x[selfmap_mask]
        }
      }
      x <- curmap$propagate(x, with.id=TRUE)
    }
    if (!with.id) {
      diag(S)[pure_sources] <- 0
    }
    return(S)
  }


  list(setup = setup,
       getType = getType,
       getName = getName,
       getMaps = getMaps,
       getDescription = getDescription,
       get_map_order = get_map_order,
       get_src_idx = get_src_idx,
       get_tar_idx = get_tar_idx,
       propagate = propagate,
       jacobian = jacobian)
}
