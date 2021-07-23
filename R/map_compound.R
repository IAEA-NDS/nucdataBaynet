create_compound_map <- function() {

  registered_map_creators <- list(
    "normerr_map" = create_normerr_map,
    "linmod_map" = create_linmod_map,
    "nonlinear_map" = create_nonlinear_map
  )

  map_list <- list()
  pure_sources <- NULL
  pure_targets <- NULL

  order_maps <- function(maps) {
    ordmaps <- maps
    if (length(ordmaps) > 1)
    {
      i <- 1
      j <- 2
      while (TRUE) {
        if (anyDuplicated(c(ordmaps[[i]]$get_src_idx(),
                            ordmaps[[j]]$get_tar_idx()))) {
          tmp <- ordmaps[[j]]
          ordmaps[[j]] <- ordmaps[[i]]
          ordmaps[[i]] <- tmp
          j <- i + 1
        } else {
          j <- j + 1
        }
        if (j > length(ordmaps)) {
          i <- i + 1
          j <- i + 1
        }
        if (j > length(ordmaps)) {
          break
        }
      }
    }
    ordmaps
  }

  setup <- function(params) {

    maps <- lapply(params$maps, function(curparams) {
      curmap <- registered_map_creators[[curparams$mapname]]()
      curmap$setup(curparams)
      curmap
    })

    max_idx <- 0
    # ensure no duplicates in indices
    for (i in seq_along(maps)) {
      if (anyDuplicated(maps[[i]]$get_src_idx()))
        stop(paste("Duplicate source indices in map",i))
      if (anyDuplicated(maps[[i]]$get_tar_idx()))
        stop(paste("Duplicate target indices in map", j))
      max_idx <- max(max_idx, maps[[i]]$get_src_idx(), maps[[i]]$get_tar_idx())
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
  }


  getName <- function() {
    return("compound_map")
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

  is_self_map <- function(map) {
    numDups <- anyDuplicated(c(map$get_src_idx(),
                               map$get_tar_idx()))
    return(numDups > 0)
  }


  propagate <- function(x, with.id = TRUE) {
    res <- x
    if (!with.id) {
      res[pure_targets] <- 0
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
          save_diag <- diag(S)[cur_tar_idx]
          diag(S)[cur_tar_idx] <- 0
        }
        S <- curmap$jacobian(x, TRUE) %*% S
        if (self_map_flag) {
          diag(S)[cur_tar_idx] <- save_diag
        }
      }
    }
    if (!with.id) {
      diag(S)[pure_sources] <- 0
    }
    return(S)
  }


  list(setup = setup,
       getName = getName,
       get_map_order = get_map_order,
       get_src_idx = get_src_idx,
       get_tar_idx = get_tar_idx,
       propagate = propagate,
       jacobian = jacobian)
}
