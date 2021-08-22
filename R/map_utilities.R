create_map <- function(params) {
  if (! "maptype" %in% names(params)) {
    print(params)
    stop(paste0("maptype missing in parameters of map named '", params$mapname, "'"))
  }
  tryCatch({
    map_generator <- get_map_generator(params$maptype)
  }, error = function(e) {
    e$message <- paste0(e$message, " (during creation of map with name ",
                        params$mapname, " of type ", params$maptype)
    stop(e)
  })
  map <- map_generator()
  map$setup(params)
  return(map)
}


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


is_self_map <- function(map) {
  numDups <- anyDuplicated(c(map$get_src_idx(),
                             map$get_tar_idx()))
  return(numDups > 0)
}

get_selfmap_mask <- function(x, src_idx, tar_idx) {
  self_map_mask <- rep(FALSE, length(x))
  self_map_mask[src_idx] <- TRUE
  self_map_mask[tar_idx] <- self_map_mask[tar_idx] & TRUE
  self_map_mask[-tar_idx] <- self_map_mask[-tar_idx] & FALSE
  return(self_map_mask)
}


get_network_structure <- function(maplist, nodes, obs, nonlinear_ind=TRUE) {
  unique_nodes <- unique(nodes)
  is_observed <- rep(FALSE, length(unique_nodes))

  node_idcs_map <- match(nodes, unique_nodes)
  is_observed[node_idcs_map[!is.na(obs)]] <- TRUE
  adjmat <- matrix(0, nrow=length(unique_nodes),
                   ncol=length(unique_nodes))
  for (curmap in maplist) {
    if (is_self_map(curmap)) { next }
    src_idx <- curmap$get_src_idx()
    tar_idx <- curmap$get_tar_idx()
    src_nodes <- unique(node_idcs_map[src_idx])
    tar_nodes <- unique(node_idcs_map[tar_idx])
    adjmat[node_idcs_map[src_idx], node_idcs_map[tar_idx]] <- 1
    if (length(src_nodes) > 1 && !curmap$is_linear() && nonlinear_ind) {
      selmat1 <- upper.tri(adjmat, diag=FALSE)
      selmat2 <- matrix(FALSE, nrow=length(unique_nodes),
                        ncol=length(unique_nodes))
      selmat2[src_nodes, src_nodes] <- TRUE
      selmat <- selmat1 & selmat2
      adjmat[selmat] <- 2
    }
  }
  colnames(adjmat) <- rownames(adjmat) <- unique_nodes
  grph <- graph_from_adjacency_matrix(adjmat, weighted="tmp")
  E(grph)$lty <- ifelse(E(grph)$tmp == 1, "solid", "dashed")
  E(grph)$arrow.mode <- ifelse(E(grph)$tmp == 1, ">", "-")
  V(grph)$color <- ifelse(is_observed, "lightgray", "white")
  V(grph)$observed <- is_observed
  return(grph)
}
