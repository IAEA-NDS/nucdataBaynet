################################################################################
#                  Define the mappings
################################################################################

truexs_avg_to_truexs2nd_avg_mappings <- lapply(indep_channels, function(curreac) {
  src_row_sel <- node_dt[, NODE==paste0("truexs_avg_", curreac) & REAC == curreac]
  tar_row_sel <- node_dt[, NODE==paste0("truexs2nd_avg_", curreac) & REAC == curreac]
  list(
    maptype = "derivative2nd_map",
    mapname = paste0("truexs_avg_to_truexs2nd_avg_for_", curreac),
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    src_x = node_dt[src_row_sel, ENERGY],
    tar_x = node_dt[tar_row_sel, ENERGY]
  )
})

truexs_hires_to_truexs2nd_hires_mappings <- lapply(indep_channels, function(curreac) {
  src_row_sel <- node_dt[, NODE==paste0("truexs_hires_", curreac) & REAC == curreac]
  tar_row_sel <- node_dt[, NODE==paste0("truexs2nd_hires_", curreac) & REAC == curreac]
  list(
    maptype = "derivative2nd_map",
    mapname = paste0("truexs_hires_to_truexs2nd_hires_for_", curreac),
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    src_x = node_dt[src_row_sel, ENERGY],
    tar_x = node_dt[tar_row_sel, ENERGY]
  )
})

truexs_hires_to_inttruexs_hires_mappings <- lapply(all_channels, function(curreac) {
  src_row_sel <- node_dt[, NODE==paste0("truexs_hires_", curreac) & REAC == curreac]
  tar_row_sel <- node_dt[, NODE==paste0("inttruexs_hires_", curreac) & REAC == curreac]
  list(
    maptype = "convolution_with_xtrafo_map",
    mapname = paste0("truexs_hires_to_inttruexs_hires_for_", curreac),
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    src_x = node_dt[src_row_sel, ENERGY],
    tar_x = node_dt[tar_row_sel, ENERGY],
    winsize = node_dt[tar_row_sel, WINSIZE[1]],
    scalex = 0, shiftx = 0
  )
})

truexs_to_expdata_mappings <- lapply(all_channels, function(curreac) {
  src_row_sel <- node_dt[, NODE==paste0("truexs_", curreac) & REAC == curreac]
  tar_row_sel <- node_dt[, grepl("^expdata_", NODE) & REAC == curreac]
  list(
    maptype = "convolution_with_xtrafo_map",
    mapname = paste0("truexs_to_expdata_for_", curreac),
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    src_x = node_dt[src_row_sel, ENERGY],
    tar_x = node_dt[tar_row_sel, ENERGY],
    winsize = 0.003, scalex=0, shiftx=0
  )
})

comptruexs_to_truexs_mappings <- apply(expand.grid(indep_channels, c("avg", "hires")), 1, function(x) {
  curreac <- x[1]; comp <- x[2]
  src_row_sel <- node_dt[, NODE == paste0("truexs_", comp, "_", curreac) & REAC == curreac]
  tar_row_sel <- node_dt[, NODE == paste0("truexs_", curreac) & REAC == curreac]
  list(
    maptype = "linearinterpol_map",
    mapname = paste0("truexs_", curreac, "_", comp, "_to_", curreac),
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    src_x = node_dt[src_row_sel, ENERGY],
    tar_x = node_dt[tar_row_sel, ENERGY]
  )
}, simplify=FALSE)

indep_to_dep_channel_mappings <- lapply(indep_channels, function(curreac) {
  src_row_sel <- node_dt[, NODE==paste0("truexs_", curreac) & REAC == curreac]
  tar_row_sel <- node_dt[, NODE==paste0("truexs_", "TOT") & REAC == "TOT"]
  list(
    maptype = "linearinterpol_map",
    mapname = paste0("truexs_", curreac, "_to_TOT"),
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    src_x = node_dt[src_row_sel, ENERGY],
    tar_x = node_dt[tar_row_sel, ENERGY]
  )
})

indepcomp_to_depcomp_channel_mappings <- apply(expand.grid(indep_channels,c("avg","hires")), 1, function(x) {
  curreac <- x[1]
  comp <- x[2]
  src_row_sel <- node_dt[, NODE==paste0("truexs_", comp, "_", curreac) & REAC == curreac]
  tar_row_sel <- node_dt[, NODE==paste0("truexs_", comp, "_TOT") & REAC == "TOT"]
  list(
    maptype = "linearinterpol_map",
    mapname = paste0("truexs_", curreac, "_to_TOT"),
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    src_x = node_dt[src_row_sel, ENERGY],
    tar_x = node_dt[tar_row_sel, ENERGY]
  )
}, simplify=FALSE)

positive_mappings <- lapply(indep_channels, function(curreac) {
  src_row_sel <- node_dt[NODE==paste0("truexs_",curreac), IDX]
  tar_row_sel <- node_dt[NODE==paste0("truexs_",curreac), IDX]
  list(
    maptype = "nonlinear_map",
    mapname = paste0("positive_map_", curreac),
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    funname = "relu"
  )
})

normerr_to_expdata_mapping <- {
  src_row_sel <- node_dt[, NODE=="normerr"]
  tar_row_sel <- node_dt[, grepl("^expdata_", NODE)]
  list(
  maptype = "normerr_map",
  mapname = "normerr_map",
  src_idx = node_dt[src_row_sel, IDX],
  tar_idx = node_dt[tar_row_sel, IDX],
  src_feat = node_dt[src_row_sel, EXPREF],
  tar_feat = node_dt[tar_row_sel, EXPREF]
  )
}


compound_mapping <- list(
  maptype = "compound_map",
  mapname = "compoundmap",
  maps = c(
    truexs_to_expdata_mappings,
    comptruexs_to_truexs_mappings,
    indepcomp_to_depcomp_channel_mappings,
    indep_to_dep_channel_mappings,
    truexs_avg_to_truexs2nd_avg_mappings,
    truexs_hires_to_truexs2nd_hires_mappings,
    truexs_hires_to_inttruexs_hires_mappings,
    positive_mappings,
    list(normerr_to_expdata_mapping)
  )
)
