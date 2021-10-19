################################################################################
#                  Compose the node data table
################################################################################

expdata_dt <- rawdt[Energy > 7000 & Energy < 12000, list(
  NODE = "expdata",
  PRIOR = 0,
  UNC = 1e-2,
  OBS = Data,
  EXPREF = Reference,
  EXFOR = Entry,
  ENERGY = Energy
)]

encalib_dt <- expdata_dt[, list(
  NODE = "encalib",
  PRIOR = c(4, 0, 0),
  UNC = c(2, 0.05, 0.05),
  OBS = NA,
  PARNAME = c("winsize", "shiftx", "scalex")
), by=EXPREF]

# true average cross section mesh and second derivative mesh

truexs_avg_dt <- data.table(
  NODE = "truexs_avg",
  PRIOR = 0,
  UNC = 1e4,
  OBS = NA,
  ENERGY = seq(6000, 14000, by=50)
)

truexs2nd_avg_dt <- truexs_avg_dt[, list(
  NODE = "truexs2nd_avg",
  PRIOR = 0,
  UNC = 1e-5,
  OBS = 0,
  ENERGY = head(ENERGY[-1], n=-1)
)]

# true hi-resolution cross section component mesh and second derivative mesh
# our modeling assumption will be that this component needs to be added
# to the average cross section to get to the experimental cross section

truexs_hires_dt <- data.table(
  NODE = "truexs_hires",
  PRIOR = 0,
  UNC = 1000,
  OBS = NA,
  ENERGY = seq(6000, 14000, by=1)
)

truexs2nd_hires_dt <- truexs_hires_dt[, list(
  NODE = "truexs2nd_hires",
  PRIOR = 0,
  UNC = 1,
  OBS = 0,
  ENERGY = head(ENERGY[-1], n=-1)
)]

int_truexs_hires_dt <- truexs_hires_dt[, list(
  NODE = sub("truexs", "inttruexs", NODE[1]),
  PRIOR = 0,
  UNC = 1e-3,
  OBS = 0,
  ENERGY = seq(min(ENERGY)+101, max(ENERGY)-101, by=100),
  WINSIZE = 200
)]

truexs_dt <- truexs_hires_dt[, list(
  NODE = "truexs",
  PRIOR = 0,
  UNC = 0,
  OBS = NA,
  ENERGY = ENERGY
)]

# these errors explain the difference between
# true cross section and experimental values

relstaterr_dt <- expdata_dt[, list(
  NODE = "relstaterr",
  PRIOR = 0,
  UNC = 3e-2,
  OBS = NA,
  ENERGY = ENERGY
), by=EXPREF]

relsyserr_dt <- expdata_dt[, list(
  NODE = "relsyserr",
  PRIOR = 0,
  UNC = ifelse(grepl("Paradela", EXPREF), 0, 5e-2),
  OBS = NA,
  ENERGY = seq(min(ENERGY)-150, max(ENERGY)+150, by=50)
), by=EXPREF]

relsyserr2nd_dt <- relsyserr_dt[, list(
  NODE = "relsyserr2nd",
  PRIOR = 0,
  UNC = 1e-6,
  OBS = 0,
  ENERGY = head(ENERGY[-1], n=-1)
), by=EXPREF]

# assemble all of these data tables together to the node data table

node_dt <- rbindlist(list(
  expdata_dt,
  relstaterr_dt,
  relsyserr_dt,
  relsyserr2nd_dt,
  truexs_dt,
  truexs_avg_dt,
  truexs2nd_avg_dt,
  truexs_hires_dt,
  truexs2nd_hires_dt,
  int_truexs_hires_dt,
  encalib_dt
), fill=TRUE)
node_dt[, IDX:=seq_len(.N)]
