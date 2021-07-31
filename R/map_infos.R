get_map_generators <- function(maptypes=character(0)) {
  registered_map_creators <- list(
    "compound_map" = create_compound_map,
    "derivative_map" = create_derivative_map,
    "derivative2nd_map" = create_derivative2nd_map,
    "linmod_map" = create_linmod_map,
    "nonlinear_map" = create_nonlinear_map,
    "normerr_map" = create_normerr_map,
    "product_map" = create_product_map
  )
  if (length(maptypes) == 0) {
    maptypes <- names(registered_map_creators)
  }
  lapply(maptypes, function(m) {
    if (! m %in% names(registered_map_creators)) {
      stop(paste0("Maptype ", m, " is not known"))
    }
  })
  return(registered_map_creators[maptypes])
}


get_map_generator <- function(maptype) {
  return(get_map_generators(maptype)[[1]])
}


get_map_params_example <- function(maptype) {
  if (! maptype %in% names(get_map_generators())) {
    stop(paste0("No map of name ", maptype, " known"))
  }
  examples <- list(
    "linmod_map" = list(
      maptype = "linmod_map",
      src_idx = 1:4,
      tar_idx = 6:10,
      src_x = 11:14,
      tar_x = c(12,12.5,13,13.5,14)
    ),
    "nonlinear_map" = list(
      maptype = "nonlinear_map",
      src_idx = 1:5,
      tar_idx = 6:10,
      funname = "exp"
    )
  )
  if (! maptype %in% names(examples)) {
    stop(paste0("No example for map ", maptype, " available"))
  }
  return(examples[[maptype]])
}
