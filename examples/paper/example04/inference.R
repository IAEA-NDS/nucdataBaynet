################################################################################
#                  Perform the Bayesian inference
################################################################################

compmap <- create_map(compound_mapping)

zprior <- node_dt$PRIOR
obs <- node_dt$OBS
U <- Diagonal(n = nrow(node_dt), x = node_dt$UNC^2)

# three different smoothness scales

smoothval <- c(1e-6, 1e-7, 1e-8)

node_dt[NODE=="truexs2nd", UNC := smoothval[1]]
U <- Diagonal(n = nrow(node_dt), x = node_dt$UNC^2)
optres <- glsalgo(compmap, zprior, U, obs, ret.list=TRUE)
node_dt[, ZPOST1 := as.vector(optres$zpost)]

node_dt[NODE=="truexs2nd", UNC := smoothval[2]]
U <- Diagonal(n = nrow(node_dt), x = node_dt$UNC^2)
optres <- glsalgo(compmap, zprior, U, obs, ret.list=TRUE)
node_dt[, ZPOST2 := as.vector(optres$zpost)]

node_dt[NODE=="truexs2nd", UNC := smoothval[3]]
U <- Diagonal(n = nrow(node_dt), x = node_dt$UNC^2)
optres <- glsalgo(compmap, zprior, U, obs, ret.list=TRUE)
node_dt[, ZPOST3 := as.vector(optres$zpost)]

names(smoothval) <- c("ZPOST1", "ZPOST2", "ZPOST3")
