################################################################################
#                  Define the mappings
################################################################################

# the following two maps take care of the smoothness property of the true cross section

# map from average true cross section to second derivative of average true cross section
truexs_to_truexs2nd_avg_mapping <- {
  src_row_sel <- node_dt[, NODE == "truexs_avg"]
  tar_row_sel <- node_dt[, NODE == "truexs2nd_avg"]
  list(
    maptype = "derivative2nd_map",
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    src_x = node_dt[src_row_sel, ENERGY],
    tar_x = node_dt[tar_row_sel, ENERGY]
  )
}

# map from hi-resolution true cross section to second derivative of high-resolution true cross section
truexs_to_truexs2nd_hires_mapping <- {
  src_row_sel <- node_dt[, NODE == "truexs_hires"]
  tar_row_sel <- node_dt[, NODE == "truexs2nd_hires"]
  list(
    maptype = "derivative2nd_map",
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    src_x = node_dt[src_row_sel, ENERGY],
    tar_x = node_dt[tar_row_sel, ENERGY]
  )
}

# map from high-resolution true cross section to convoluted cross section with pseudo observation
truexs_to_inttruexs_hires_mapping <- {
  src_row_sel <- node_dt[, NODE == "truexs_hires"]
  tar_row_sel <- node_dt[, NODE == "inttruexs_hires"]
  list(
    maptype = "convolution_with_xtrafo_map",
    mapname = "truexs_hires_to_inttruexs_map",
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    src_x = node_dt[src_row_sel, ENERGY],
    tar_x = node_dt[tar_row_sel, ENERGY],
    winsize = node_dt[tar_row_sel, WINSIZE[1]],
    scalex = 0, shiftx = 0
  )
}

# compose the overall cross section from the hires and the avg xs component

# the following two mappings map the true cross section (avg + hires) to the measurements
truexs_avg_to_truexs_mapping <- {
  src_row_sel <- node_dt[, NODE == "truexs_avg"]
  tar_row_sel <- node_dt[, NODE == "truexs"]
  list(
    maptype = "linearinterpol_map",
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    src_x = node_dt[src_row_sel, ENERGY],
    tar_x = node_dt[tar_row_sel, ENERGY]
  )
}

truexs_hires_to_truexs_mapping <- {
  src_row_sel <- node_dt[, NODE == "truexs_hires"]
  tar_row_sel <- node_dt[, NODE == "truexs"]
  list(
    maptype = "linearinterpol_map",
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    src_x = node_dt[src_row_sel, ENERGY],
    tar_x = node_dt[tar_row_sel, ENERGY]
  )
}

truexs_positivity_mapping <- {
  src_row_sel <- node_dt[, NODE == "truexs"]
  tar_row_sel <- src_row_sel
  list(
    maptype = "nonlinear_map",
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    funname = "relu"
  )
}

# map the average cross section to the experimental node
truexs_to_expdata_mappings <- lapply(unique(expdata_dt$EXPREF), function(curexp) {
  tar_row_sel <- node_dt[, NODE =="expdata" & EXPREF == curexp]
  src_row_sel <- node_dt[, NODE == "truexs"]
  list(
    maptype = "convolution_with_xtrafo_map",
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    src_x = node_dt[src_row_sel, ENERGY],
    tar_x = node_dt[tar_row_sel, ENERGY],
    shiftx_idx = node_dt[NODE=="encalib" & EXPREF == curexp & PARNAME == "shiftx", IDX],
    scalex_idx = node_dt[NODE=="encalib" & EXPREF == curexp & PARNAME == "scalex", IDX],
    winsize_idx = node_dt[NODE=="encalib" & EXPREF == curexp & PARNAME == "winsize", IDX]
  )
})
names(truexs_to_expdata_mappings) <- unique(expdata_dt$EXPREF)

truexs_avg_to_expdata_mappings <- lapply(unique(expdata_dt$EXPREF), function(curexp) {
  tar_row_sel <- node_dt[, NODE =="expdata" & EXPREF == curexp]
  src_row_sel <- node_dt[, NODE == "truexs_avg"]
  list(
    maptype = "convolution_with_xtrafo_map",
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    src_x = node_dt[src_row_sel, ENERGY],
    tar_x = node_dt[tar_row_sel, ENERGY],
    shiftx_idx = node_dt[NODE=="encalib" & EXPREF == curexp & PARNAME == "shiftx", IDX],
    scalex_idx = node_dt[NODE=="encalib" & EXPREF == curexp & PARNAME == "scalex", IDX],
    winsize_idx = node_dt[NODE=="encalib" & EXPREF == curexp & PARNAME == "winsize", IDX]
  )
})
names(truexs_avg_to_expdata_mappings) <- unique(expdata_dt$EXPREF)

# the following mappings establish the contributions of the relative normalization errors
# in principle we could take the relative normalization error with respect to
# the true cross section (avg + hires), but here we define as relative to the
# average true cross section. The construction we use is to define a mapping
# directly from the true average cross section to the convoluted experimental data,
# and then use a product mapping.

syserr_to_expdata_mappings <- lapply(unique(expdata_dt$EXPREF), function(curexp) {
  src_row_sel <- node_dt[, NODE == "relsyserr" & EXPREF == curexp]
  tar_row_sel <- node_dt[, NODE == "expdata" & EXPREF == curexp]
  list(
    maptype = "convolution_with_xtrafo_map",
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    src_x = node_dt[src_row_sel, ENERGY],
    tar_x = node_dt[tar_row_sel, ENERGY],
    shiftx = 0, scalex = 0, winsize = 1
  )
})
names(syserr_to_expdata_mappings) <- unique(expdata_dt$EXPREF)

relsyserr_to_expdata_mappings <- lapply(unique(expdata_dt$EXPREF), function(curexp) {
  list(
    maptype = "product_map",
    maps = list(syserr_to_expdata_mappings[[curexp]],
                truexs_avg_to_expdata_mappings[[curexp]])
  )
})

relsyserr_to_relsyserr2nd_mappings <- lapply(unique(expdata_dt$EXPREF), function(curexp) {
  src_row_sel <- node_dt[, NODE == "relsyserr" & EXPREF == curexp]
  tar_row_sel <- node_dt[, NODE == "relsyserr2nd" & EXPREF == curexp]
  list(
    maptype = "derivative2nd_map",
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    src_x = node_dt[src_row_sel, ENERGY],
    tar_x = node_dt[tar_row_sel, ENERGY]
  )
})
names(relsyserr_to_relsyserr2nd_mappings) <- unique(expdata_dt$EXPREF)

staterr_to_expdata_mappings <- lapply(unique(expdata_dt$EXPREF), function(curexp) {
  src_row_sel <- node_dt[, NODE == "relstaterr" & EXPREF == curexp]
  tar_row_sel <- node_dt[, NODE == "expdata" & EXPREF == curexp]
  list(
    maptype = "linearinterpol_map",
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    src_x = node_dt[src_row_sel, ENERGY],
    tar_x = node_dt[tar_row_sel, ENERGY]
  )
})
names(staterr_to_expdata_mappings) <- unique(expdata_dt$EXPREF)

relstaterr_to_expdata_mappings <- lapply(unique(expdata_dt$EXPREF), function(curexp) {
  list(
    maptype = "product_map",
    maps = list(staterr_to_expdata_mappings[[curexp]],
                truexs_to_expdata_mappings[[curexp]])
  )
})

# assemble all the mappings to the compound map

compound_mapping <- list(
  maptype = "compound_map",
  maps = c(
    # smoothness
    list(
      truexs_to_truexs2nd_avg_mapping,
      truexs_to_truexs2nd_hires_mapping,
      truexs_to_inttruexs_hires_mapping,
      truexs_avg_to_truexs_mapping,
      truexs_hires_to_truexs_mapping,
      truexs_positivity_mapping
    ),
    # true xs to expdata mapping
    truexs_to_expdata_mappings,
    # experimental error mappings
    relsyserr_to_expdata_mappings,
    relsyserr_to_relsyserr2nd_mappings,
    relstaterr_to_expdata_mappings
  )
)
