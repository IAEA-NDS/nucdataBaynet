create_map <- function(params) {
  map_generator <- get_map_generator(params$maptype)
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


get_network_structure <- function(maplist, node_names) {
  unique_nodes <- unique(node_names)
  node_idcs_map <- match(node_names, unique_nodes)
  adjmat <- matrix(FALSE, nrow=length(unique_nodes),
                   ncol=length(unique_nodes))
  for (curmap in maplist) {
    if (is_self_map(curmap)) { next }
    src_idx <- curmap$get_src_idx()
    tar_idx <- curmap$get_tar_idx()
    adjmat[node_idcs_map[src_idx], node_idcs_map[tar_idx]] <- TRUE
  }
  colnames(adjmat) <- rownames(adjmat) <- unique_nodes
  return(adjmat)
}
