################################################################################
#                  Define the mappings
################################################################################

truexs_to_expdata_mapping <- list(
  maptype = "convolution_with_xtrafo_map",
  mapname = "aconvolution",
  src_idx = node_dt[NODE=="truexs", IDX],
  tar_idx = node_dt[NODE=="exp", IDX],
  src_x = node_dt[NODE=="truexs", ENERGY],
  tar_x = node_dt[NODE=="exp", ENERGY],
  winsize = enres, scalex = 0, shiftx = 0
)

truexs_to_truexs2nd_mapping <- list(
  maptype = "derivative2nd_map",
  mapname = "aderivativemap",
  src_idx = node_dt[NODE=="truexs", IDX],
  tar_idx = node_dt[NODE=="truexs2nd", IDX],
  src_x = node_dt[NODE=="truexs", ENERGY],
  tar_x = node_dt[NODE=="truexs2nd", ENERGY]
)

compound_mapping <- list(
  maptype = "compound_map",
  mapname = "acompoundmap",
  maps = list(
    truexs_to_expdata_mapping,
    truexs_to_truexs2nd_mapping
  )
)

