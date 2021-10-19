################################################################################
#                  Define the mappings
################################################################################

auxmodpar_to_modpar_mapping <- list(
  maptype = "nonlinear_map",
  mapname = "auxmodpar_to_modpar",
  src_idx = node_dt[NODE == "auxmodpar", IDX],
  tar_idx = node_dt[NODE == "modpar", IDX],
  funname = "limiter",
  minvalue = 0.9, maxvalue=1.1
)

modpar_to_modpred_mapping <- list(
  maptype = "linearmap_map",
  mapname = "modpar_to_modpred",
  src_idx = node_dt[NODE=="modpar", IDX],
  tar_idx = node_dt[grepl("^pred_", NODE), IDX],
  yref = talysinfo$yref,
  pref = talysinfo$pref,
  S = talysinfo$S
)

modpred_to_truexs_mappings <- lapply(talysinfo$indep_reacs, function(curreac) {
  src_row_sel <- node_dt[, grepl("^pred_", NODE) & REAC == curreac]
  tar_row_sel <- node_dt[, grepl("^truexs_", NODE) & REAC == curreac]
  list(
    maptype = "linearinterpol_map",
    mapname = paste0("modpred_to_truexs_for_", curreac),
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    src_x = node_dt[src_row_sel, ENERGY],
    tar_x = node_dt[tar_row_sel, ENERGY]
  )
})
names(modpred_to_truexs_mappings) <- talysinfo$indep_reacs

indep_moddef_to_truexs_mappings <- lapply(talysinfo$indep_reacs, function(curreac) {
  src_row_sel <- node_dt[, grepl("^def_", NODE) & REAC == curreac]
  tar_row_sel <- node_dt[, grepl("^truexs_", NODE) & REAC == curreac]
  list(
    maptype = "linearinterpol_map",
    mapname = paste0("moddef_to_truexs_for_", curreac),
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    src_x = node_dt[src_row_sel, ENERGY],
    tar_x = node_dt[tar_row_sel, ENERGY],
    zero_outside = TRUE
  )
})
names(indep_moddef_to_truexs_mappings) <- talysinfo$indep_reacs

truexs_positivity_mapping <- {
  src_row_sel <- node_dt[, grepl("^truexs_", NODE) & REAC %in% talysinfo$indep_reacs]
  tar_row_sel <- src_row_sel
  list(
    maptype = "nonlinear_map",
    mapname = "truexs_positivity_map",
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    funname = "relu"
  )
}

indep_moddef_to_indep_moddef2nd_mappings <- lapply(talysinfo$indep_reacs, function(curreac) {
  src_row_sel <- node_dt[, grepl("^def_", NODE) & REAC == curreac]
  tar_row_sel <- node_dt[, grepl("^def2nd_", NODE) & REAC == curreac]
  list(
    maptype = "derivative2nd_map",
    mapname = paste0("moddef_to_moddef2nd_for_", curreac),
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    src_x = node_dt[src_row_sel, ENERGY],
    tar_x = node_dt[tar_row_sel, ENERGY]
  )
})

relmoddef_to_truexs_mappings <- lapply(talysinfo$indep_reacs, function(curreac) {
  list(
    maptype = "product_map",
    mapname = paste0("relmoddef_to_truexs_for_", curreac),
    maps = list(indep_moddef_to_truexs_mappings[[curreac]],
                modpred_to_truexs_mappings[[curreac]])
  )
})

truexs_to_expdata_mappings <- lapply(talysinfo$all_reacs, function(curreac) {
  src_row_sel <- node_dt[, grepl("^truexs_", NODE) & REAC == curreac]
  tar_row_sel <- node_dt[, grepl("^exp_", NODE) & REAC == curreac]
  list(
    maptype = "linearinterpol_map",
    mapname = paste0("truexs_to_expdata_for_", curreac),
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    src_x = node_dt[src_row_sel, ENERGY],
    tar_x = node_dt[tar_row_sel, ENERGY]
  )
})
names(truexs_to_expdata_mappings) <- talysinfo$all_reacs

staterr_to_expdata_mappings <- lapply(talysinfo$all_reacs, function(curreac) {
  src_row_sel <- node_dt[, grepl("^staterr_",NODE) & REAC == curreac]
  tar_row_sel <- node_dt[, grepl("^exp_",NODE) & REAC == curreac]
  list(
    maptype = "linearinterpol_map",
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    src_x = node_dt[src_row_sel, ENERGY],
    tar_x = node_dt[tar_row_sel, ENERGY]
  )
})
names(staterr_to_expdata_mappings) <- talysinfo$all_reacs

relstaterr_to_expdata_mappings <- lapply(talysinfo$all_reacs, function(curreac) {
  list(
    maptype = "product_map",
    mapname = paste0("relstaterr_for_", curreac),
    maps = list(staterr_to_expdata_mappings[[curreac]],
                truexs_to_expdata_mappings[[curreac]])
  )
})

normerr_to_expdata_mappings <- lapply(talysinfo$all_reacs, function(curreac) {
  src_row_sel <- node_dt[, grepl("^normerr_",NODE) & REAC == curreac]
  tar_row_sel <- node_dt[, grepl("^exp_",NODE) & REAC == curreac]
  list(
    maptype = "normerr_map",
    mapname = "normerr_map",
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    src_feat = node_dt[src_row_sel, EXPREF],
    tar_feat = node_dt[tar_row_sel, EXPREF]
  )
})
names(normerr_to_expdata_mappings) <- talysinfo$all_reacs

relnormerr_to_expdata_mappings <- lapply(talysinfo$all_reacs, function(curreac) {
  list(
    maptype = "product_map",
    mapname = paste0("relnormerr_for_", curreac),
    maps = list(normerr_to_expdata_mappings[[curreac]],
                truexs_to_expdata_mappings[[curreac]])
  )
})

truexs_to_truexs_tot_mappings <- lapply(talysinfo$indep_reacs, function(curreac) {
  src_row_sel <- node_dt[, grepl("^truexs_", NODE) & REAC == curreac]
  tar_row_sel <- node_dt[, grepl("^truexs_\\(N,TOT\\)", NODE) & REAC == "(N,TOT)"]
  list(
    maptype = "linearinterpol_map",
    mapname = paste0("truexs_to_expdata_for_", curreac),
    src_idx = node_dt[src_row_sel, IDX],
    tar_idx = node_dt[tar_row_sel, IDX],
    src_x = node_dt[src_row_sel, ENERGY],
    tar_x = node_dt[tar_row_sel, ENERGY]
  )
})


comp_mapping <- list(
  maptype = "compound_map",
  maps = c(
    # individual mappings must be enclosed by list
    list(
      modpar_to_modpred_mapping,
      auxmodpar_to_modpar_mapping
    ),
    # list of mapping definitions go here
    modpred_to_truexs_mappings,
    relmoddef_to_truexs_mappings,
    list(truexs_positivity_mapping),
    truexs_to_expdata_mappings,
    truexs_to_truexs_tot_mappings,
    indep_moddef_to_indep_moddef2nd_mappings,
    relstaterr_to_expdata_mappings,
    relnormerr_to_expdata_mappings
  )
)
