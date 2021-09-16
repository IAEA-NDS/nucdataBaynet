#' Create a mapping
#'
#' This function creates a mapping from a parameter list.
#'
#' It calls the mapping creation function responsible for mappings
#' of \code{maptype} stated in the parameter list. It also
#' invokes the \code{setup} function to initialize the mapping,
#' provided by all mapping objects, see the \strong{Value} section
#' of \code{\link{create_maptype_map}}.
#' For all maps the fields \code{mapname} and \code{maptype} must be present
#' in the parameter list.
#' The other parameters depend on the specifics
#' of the mapping, see e.g., \code{\link{create_linearinterpol_map}}.
#'
#' @note
#' This function should be preferred over the specialized mapping creation
#' functions, such as \code{\link{create_linearinerpol_mapping}}.
#'
#' @param params Parameter list with the mapping definition, see \strong{Details}.
#'
#' @return
#' Returns a list of functions to operate with the mapping, see \code{\link{create_maptype_map}}.
#'
#' @export
#'
#' @family mappings
#' @examples
#' params <- list(
#'   mapname = "mylinearintmap",
#'   maptype = "linearinterpol_map",
#'   src_idx = 1:3,
#'   tar_idx = 4:6,
#'   src_x = c(1,5,10),
#'   tar_x = c(4,5,6)
#' )
#' mymap <- create_map(params)
#' x <- c(1,2,3,0,0,0)
#' mymap$propagate(x)
#' mymap$jacobian(x)
#'
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


#' Sort a list of mappings
#'
#' Sort a list of mappings so that the sequential application of
#' the mappings yields the correct result.
#'
#' Mappings whose source indices contain target indices of another mapping
#' must not be applied before that mapping. This function puts the mappings
#' into an order to ensure this criterion is met.
#'
#' @param maps List of mapping objects, e.g., as obtained by \code{\link{create_map}}.
#'
#' @return
#' Return a list with the ordered maps
#' @export
#'
#' @examples
#' params1 <- list(
#'   mapname = "mylinearmap1",
#'   maptype = "linearinterpol_map",
#'   src_idx = 1:3,
#'   tar_idx = 4:6,
#'   src_x = c(1,5,10),
#'   tar_x = c(4,5,6)
#' )
#' params2 <- list(
#'   mapname = "mylinearmap2",
#'   maptype = "linearinterpol_map",
#'   src_idx = 4:6,
#'   tar_idx = 7:9,
#'   src_x = c(4,5,7),
#'   tar_x = c(4.5, 5, 5.5)
#' )
#' mymap1 <- create_map(params1)
#' mymap2 <- create_map(params2)
#' ordmaps <- order_maps(list(mymap2, mymap1))
#' lapply(ordmaps, function(x) x$getName())
#'
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


#' Get network structre
#'
#' Create a Bayesian network graph
#'
#' @param maplist A list of mapping objects
#' @param nodes A vector with the node name associated to the variables.
#'              Several variables can belong to the same node.
#' @param obs A vector of observed values. Unobserved values are represented by
#'            \code{NA}.
#' @param nonlinear_ind Should non-linear interactions between nodes be
#'                      indicated by dashed lines.
#'
#' @return
#' Return an igraph graph object
#' @export
#'
#' @examples
#' \dontrun{
#'   library(igraph)
#'   params1 <- list(
#'     mapname = "mylinearmap1",
#'     maptype = "linearinterpol_map",
#'     src_idx = 1:3,
#'     tar_idx = 4:6,
#'     src_x = c(1,5,10),
#'     tar_x = c(4,5,6)
#'   )
#'   params2 <- list(
#'     mapname = "mylinearmap2",
#'     maptype = "linearinterpol_map",
#'     src_idx = 4:6,
#'     tar_idx = 7:9,
#'     src_x = c(4,5,7),
#'     tar_x = c(4.5, 5, 5.5)
#'   )
#'   compparams <- list(
#'     mapname = "mycompmap",
#'     maptype = "compound_map",
#'     maps = list(params1, params2)
#'   )
#'
#'   mymap1 <- create_map(params1)
#'   mymap2 <- create_map(params2)
#'
#'   mymap <- create_map(compparams)
#'   maplist <- mymap$getMaps()
#'   nodes <- c(rep("mynode1", 3),
#'              rep("mynode2", 3),
#'              rep("obsnode", 3))
#'   obs <- rep(NA, 9)
#'
#'   grph <- get_network_structure(maplist)
#' }
#'
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


is_self_map <- function(map) {
  numDups <- anyDuplicated(c(map$get_src_idx(),
                             map$get_tar_idx()))
  return(numDups > 0)
}

get_selfmap_mask <- function(src_idx, tar_idx, n=max(src_idx, tar_idx)) {
  self_map_mask <- rep(FALSE, n)
  self_map_mask[src_idx] <- TRUE
  self_map_mask[tar_idx] <- self_map_mask[tar_idx] & TRUE
  self_map_mask[-tar_idx] <- self_map_mask[-tar_idx] & FALSE
  return(self_map_mask)
}
