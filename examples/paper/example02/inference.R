################################################################################
#                  Do the inference
################################################################################

compmap <- create_map(compound_mapping)

zprior <- node_dt$PRIOR
obs <- node_dt$OBS
U <- Diagonal(n = nrow(node_dt), x = node_dt$UNC^2)

optres <- LMalgo(compmap, zprior, U, obs, print.info=TRUE)
zpost <- optres$zpost

node_dt[, ZPOST := zpost]
node_dt[, PRED := compmap$propagate(zpost)]

# compute selected uncertainties

# truexs avg uncertainties
ensel <- node_dt[, ENERGY > 0.9 & ENERGY < 2.1]
sel <- node_dt[, grepl("^truexs_avg_(TOT|EL|INL)", NODE) & ensel]
covmat <- get_posterior_cov(compmap, zpost, U, obs, sel, sel, ret.dep=TRUE)
node_dt[sel, PREDUNC := sqrt(diag(covmat))]

# truexs hires uncertainties
ensel <- node_dt[, ENERGY > 0.9 & ENERGY < 2.1]
sel <- node_dt[, grepl("^truexs_hires_(TOT|EL|INL)", NODE) & ensel]
covmat <- get_posterior_cov(compmap, zpost, U, obs, sel, sel, ret.dep=TRUE)
node_dt[sel, PREDUNC := sqrt(diag(covmat))]

# truexs uncertainties
ensel <- node_dt[, ENERGY > 0.9 & ENERGY < 2.1]
sel <- node_dt[, grepl("^truexs_(TOT|EL|INL)", NODE) & ensel]
covmat <- get_posterior_cov(compmap, zpost, U, obs, sel, sel, ret.dep=TRUE)
node_dt[sel, PREDUNC := sqrt(diag(covmat))]

# create a sample from posterior distribution
postsmpl <- get_posterior_sample(compmap, zpost, U, obs, 10)
postsmpl <- apply(postsmpl, 2, compmap$propagate)
node_dt <- cbind(node_dt, postsmpl)
