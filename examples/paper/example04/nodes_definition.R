################################################################################
#                  Compose the node data table
################################################################################

set.seed(27)

enres <- 500
mesh_energy <- seq(300, 25000, length=1000)
exp_energy <- seq(300+enres/2, 25000-enres/2, by=2000) 
# +1500 because of energy convolution
exp_measurement <- sin(exp_energy *2*pi / diff(range(exp_energy))*2) + rnorm(length(exp_energy), 0, 0.2)

expdata_dt <- data.table(
  NODE = "exp",
  PRIOR = 0,
  UNC = 1e-1,
  OBS = exp_measurement,
  ENERGY = exp_energy
)

truexs_dt <- data.table(
  NODE = "truexs",
  PRIOR = 0,
  UNC = 1e4,
  OBS = NA,
  ENERGY = mesh_energy
)

truexs2nd_dt <- truexs_dt[, list(
  NODE = "truexs2nd",
  PRIOR = 0,
  UNC = 1e-6,
  OBS = 0,
  ENERGY = head(mesh_energy[-1], n=-1)
)]

node_dt <- rbindlist(list(expdata_dt, truexs_dt, truexs2nd_dt), fill=TRUE)
node_dt[, IDX := seq_len(.N)]
