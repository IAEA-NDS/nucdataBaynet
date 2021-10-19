################################################################################
#                  Compose the node data table
################################################################################

energymesh <- unique(talysinfo$extNeedsDt[, L1])
energymesh <- sort(c(energymesh, unlist(talysinfo$energy_thresholds)))
names(energymesh) <- NULL

expdata_dt <- expDt[L1>5, list(
  NODE = paste0("exp_", REAC),
  PRIOR = 0,
  UNC = 1e-1,
  OBS = DATA,
  REAC = REAC,
  ENERGY = L1,
  EXPREF = EXPID
)]

relstaterr_dt <- expdata_dt[, list(
  NODE = paste0("staterr_", REAC),
  PRIOR = 0,
  UNC = 0.01,
  ENERGY = ENERGY,
  EXPREF = EXPREF
), by=REAC]

relnormerr_dt <- expdata_dt[, list(
  NODE = paste0("normerr_", REAC),
  PRIOR = 0,
  # we reduce the normalization uncertainty of some datapoints to make them
  # work as anchor datasets to determine the normalization errors.
  # anchor datasets were chosen by visual inspection of plots
  # UNC = ifelse(EXPREF %in% c("23171003","41614019","22316003","22976017"), 0.01, 0.10)
  UNC = 0.10
), by=c("REAC","EXPREF")]

auxmodpar_dt <- refParamDt[ADJUSTABLE==TRUE, list(
  NODE = "auxmodpar",
  PRIOR = 1,
  UNC = 0.05,
  PARNAME = PARNAME
)]

modpar_dt <- auxmodpar_dt[, list(
  NODE = "modpar",
  PRIOR = 0,
  UNC = 0,
  PARNAME = PARNAME
)]

modpred_dt <- extNeedsDt[REAC %in% talysinfo$indep_reacs, list(
  NODE = paste0("pred_", REAC),
  PRIOR = 0,
  UNC = 0,
  REAC = REAC,
  ENERGY = L1
)]

moddef_dt <- modpred_dt[, list(
  NODE = paste0("def_", REAC),
  PRIOR = 0,
  UNC = 0.2,
  ENERGY = energymesh[energymesh > talysinfo$energy_thresholds[[REAC]]]
), by=REAC]

moddef2nd_dt <- moddef_dt[, list(
  NODE = paste0("def2nd_", REAC),
  PRIOR = 0,
  UNC = 1e-2,
  OBS = 0,
  ENERGY = head(ENERGY[-1], n=-1)
), by=REAC]

truexs_dt <- extNeedsDt[REAC %in% talysinfo$all_reacs, list(
  NODE = paste0("truexs_", REAC),
  PRIOR = 0,
  UNC = 0,
  REAC = REAC,
  ENERGY = L1
)]

node_dt <- rbindlist(list(
  expdata_dt, relstaterr_dt, relnormerr_dt,
  auxmodpar_dt, modpar_dt, modpred_dt, moddef_dt, moddef2nd_dt,
  truexs_dt
), fill=TRUE)
node_dt[, IDX:=seq_len(.N)]
