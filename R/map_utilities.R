create_map <- function(params) {
  map_generator <- get_map_generator(params$mapname)
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
