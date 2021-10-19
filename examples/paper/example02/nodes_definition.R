################################################################################
#                  Compose the node data table
################################################################################

indep_channels <- c("EL", "INL")
dep_channels <- c("TOT")
all_channels <- c(indep_channels, dep_channels)

expdata_dt <- rawdata[, list(
  NODE = paste0("expdata_", REAC),
  PRIOR = 0,
  UNC = ifelse(REAC=="TOT", 3000*0.1, ifelse(REAC=="INL", 700*0.1, 2000*0.1)),
  OBS = DATA,
  REAC = REAC,
  EXPREF = EXPID,
  ENERGY = L1
)]

truexs_dt <- rawdata[, list(
  NODE = paste0("truexs_", REAC),
  PRIOR = 0,
  UNC = 0,
  ENERGY = seq(0.75, 2.25, by=0.001)
), by=REAC]

truexs_avg_dt <- rawdata[, list(
  NODE = paste0("truexs_avg_", REAC),
  PRIOR = 0,
  UNC = ifelse(REAC %in% indep_channels, 1e8, 0),
  ENERGY = seq(0.75, 2.25, by=0.05)
), by=REAC]

truexs2nd_avg_dt <- truexs_avg_dt[NODE %in% c("truexs_avg_EL", "truexs_avg_INL"), list(
  NODE = sub("truexs","truexs2nd",NODE[1]),
  PRIOR = 0,
  UNC = 1e4,
  OBS = 0,
  ENERGY = head(ENERGY[-1], n=-1)
), by=REAC]

truexs_hires_dt <- rawdata[, list(
  NODE = paste0("truexs_hires_", REAC),
  PRIOR = 0,
  #UNC = ifelse(REAC=="EL", 1e3, ifelse(REAC=="INL", 500, 0)),
  UNC = ifelse(REAC=="EL", 1e4, ifelse(REAC=="INL", 1e4, 0)),
  ENERGY = seq(0.75, 2.25, by=0.001)
), by=REAC]

truexs2nd_hires_dt <- truexs_hires_dt[NODE %in% c("truexs_hires_EL", "truexs_hires_INL"), list(
  NODE = sub("truexs","truexs2nd",NODE[1]),
  PRIOR = 0,
  UNC = 1e8,
  OBS = 0,
  ENERGY = head(ENERGY[-1], n=-1)
), by=REAC]

int_truexs_hires_dt <- truexs_hires_dt[, list(
  NODE = sub("truexs", "inttruexs", NODE[1]),
  PRIOR = 0,
  UNC = 5e1,
  OBS = 0,
  ENERGY = seq(1.1, 1.9, by=0.1),
  WINSIZE = 0.2
), by=REAC]

normerr_dt <- rawdata[, list(
  NODE = "normerr",
  PRIOR = 0,
  UNC = 100
), by=c("REAC","EXPID")]
setnames(normerr_dt, "EXPID", "EXPREF")

node_dt <- rbindlist(list(
  expdata_dt,
  truexs_dt,
  normerr_dt,
  truexs_hires_dt,
  int_truexs_hires_dt,
  truexs2nd_hires_dt,
  truexs_avg_dt,
  truexs2nd_avg_dt
), fill=TRUE)
node_dt[, IDX := seq_len(.N)]
