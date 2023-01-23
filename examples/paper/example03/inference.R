################################################################################
#                  Perform the Bayesian inference
################################################################################

compmap <- create_map(comp_mapping)

zprior <- node_dt$PRIOR
obs <- node_dt$OBS
U <- Diagonal(n = nrow(node_dt), x = node_dt$UNC^2)

# optimization schedule
nodes <- node_dt[, unique(NODE)]
optim_scheme <- list(
  list(c(grep("^def_", nodes, value=TRUE),
         grep("^staterr_", nodes, value=TRUE),
         grep("^normerr_", nodes, value=TRUE)), 30, 10),
  list(c(grep("^auxmodpar", nodes, value=TRUE),
         grep("^normerr_", nodes, value=TRUE),
         grep("^staterr_", nodes, value=TRUE),
         grep("^def_", nodes, value=TRUE)), 100, 10),
  list(nodes, 30, 10)
)

zpost <- zprior
for (cur_stage in optim_scheme) {
  cat(paste0("############################\n",
             "optimizing the following nodes:\n",
             paste0(cur_stage[[1]], collapse=", "), "\n",
             "############################\n"))
  zref <- zpost
  nodesel <- cur_stage[[1]]
  maxiter <- cur_stage[[2]]
  miniter <- cur_stage[[3]]
  adjust_idcs <- node_dt[NODE %in% nodesel, IDX]
  optres <- LMalgo(compmap, zprior, U, obs, print.info=TRUE, zref = zref,
                   adjust_idcs = adjust_idcs, control=list(maxcount=maxiter, mincount=miniter, tau=1e-6),
                   must.converge=FALSE)
  zpost <- optres$zpost
}

node_dt[, ZPOST := zpost]
node_dt[, PRED := compmap$propagate(zpost)]
node_dt[, PRIORPRED := compmap$propagate(PRIOR)]

S <- compmap$jacobian(optres$zpost)

# propagate prior experiment uncertainty
idxred <- node_dt[grepl("^staterr_",NODE) | grepl("^normerr_",NODE) | grepl("^exp_",NODE), IDX]
expsel <- node_dt[, grepl("^exp_",NODE)]
Ured <- U
Ured[-idxred,-idxred] <- 0
err <- sqrt(diag(S %*% Ured %*% t(S)))
err <- err[expsel]
node_dt[expsel, PRIORUNC := err]

# get relevant prior covariance matrices
idxred <- node_dt[grepl("^auxmodpar$",NODE) | grepl("^def_",NODE), IDX]
idxtar <- node_dt[grepl("^truexs_",NODE), IDX]
Ured <- U
Ured[-idxred,-idxred] <- 0
err <- sqrt(diag(S %*% Ured %*% t(S)))
err <- err[idxtar]
node_dt[idxtar, PRIORUNC := err]

# get relevant posterior covariance matrices
sel <- node_dt[, grepl("^truexs_", NODE)]
covmat <- get_posterior_cov(compmap, zpost, U, obs, sel, sel, ret.dep=TRUE)
node_dt[sel, PREDUNC := sqrt(diag(covmat))]

sel <- node_dt[, grepl("^def_", NODE)]
covmat <- get_posterior_cov(compmap, zpost, U, obs, sel, sel, ret.dep=TRUE)
node_dt[sel, PREDUNC := sqrt(diag(covmat))]

sel <- node_dt[, grepl("^auxmodpar", NODE)]
covmat <- get_posterior_cov(compmap, zpost, U, obs, sel, sel, ret.dep=TRUE)
node_dt[sel, PREDUNC := sqrt(diag(covmat))]

sel <- node_dt[, grepl("^normerr_", NODE)]
covmat <- get_posterior_cov(compmap, zpost, U, obs, sel, sel, ret.dep=TRUE)
node_dt[sel, PREDUNC := sqrt(diag(covmat))]
