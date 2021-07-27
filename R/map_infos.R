get_map_generators <- function(mapnames=character(0)) {
  registered_map_creators <- list(
    "compound_map" = create_compound_map,
    "derivative_map" = create_derivative_map,
    "derivative2nd_map" = create_derivative2nd_map,
    "linmod_map" = create_linmod_map,
    "nonlinear_map" = create_nonlinear_map,
    "normerr_map" = create_normerr_map,
    "product_map" = create_product_map
  )
  if (length(mapnames) == 0) {
    mapnames <- names(registered_map_creators)
  }
  lapply(mapnames, function(m) {
    if (! m %in% names(registered_map_creators)) {
      stop(paste0("Maptype ", m, " is not known"))
    }
  })
  return(registered_map_creators[mapnames])
}


get_map_generator <- function(mapname) {
  return(get_map_generators(mapname)[[1]])
}


get_map_params_example <- function(mapname) {
  if (! mapname %in% names(get_map_generators())) {
    stop(paste0("No map of name ", mapname, " known"))
  }
  examples <- list(
    "linmod_map" = list(
      mapname = "linmod_map",
      src_idx = 1:4,
      tar_idx = 6:10,
      src_x = 11:14,
      tar_x = c(12,12.5,13,13.5,14)
    ),
    "nonlinear_map" = list(
      mapname = "nonlinear_map",
      src_idx = 1:5,
      tar_idx = 6:10,
      funname = "exp"
    )
  )
  if (! mapname %in% names(examples)) {
    stop(paste0("No example for map ", mapname, " available"))
  }
  return(examples[[mapname]])
}
