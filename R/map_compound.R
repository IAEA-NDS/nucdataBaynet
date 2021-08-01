create_compound_map <- function() {

  mapinfo <- list()
  map_list <- list()
  pure_sources <- NULL
  pure_targets <- NULL

  setup <- function(params) {

    maps <- lapply(params$maps, function(curparams) {
        return(create_map(curparams))
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
    return(unlist(lapply(map_list, function(x) x$getType())))
  }


  propagate <- function(x, with.id = TRUE) {
    res <- x
    if (!with.id) {
      res[-pure_sources] <- 0
    }
    initialres <- res
    for (curmap in map_list) {
      self_map_flag <- is_self_map(curmap)
      cur_tar_idx <- curmap$get_tar_idx()
      if (self_map_flag) {
        res[cur_tar_idx] <- res[cur_tar_idx] - initialres[cur_tar_idx]
      }
      if (with.id) {
        if (!self_map_flag) {
          res <- curmap$propagate(res, TRUE)
        } else {
          res[cur_tar_idx] <- curmap$propagate(res, FALSE)[cur_tar_idx]
        }
      } else {
        if (!self_map_flag) {
          res <- res + curmap$propagate(res, FALSE)
        } else {
          res[cur_tar_idx] <- curmap$propagate(res, FALSE)[cur_tar_idx]
        }
      }
      if (self_map_flag) {
        res[cur_tar_idx] <- res[cur_tar_idx] + initialres[cur_tar_idx]
      }
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
        self_map_flag <- is_self_map(curmap)
        cur_tar_idx <- curmap$get_tar_idx()
        if (self_map_flag) {
          diag(S)[cur_tar_idx] <- 0
          x[cur_tar_idx] <- x[cur_tar_idx] - orig_x[cur_tar_idx]
        }
        curS <- curmap$jacobian(x, TRUE)
        if (self_map_flag) {
          diag(curS)[cur_tar_idx] <- diag(curS)[cur_tar_idx] - 1
        }
        S <- curS %*% S
        if (self_map_flag) {
          if (with.id) {
            diag(S)[cur_tar_idx] <- 1
          }
          x[cur_tar_idx] <- x[cur_tar_idx] + orig_x[cur_tar_idx]
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
