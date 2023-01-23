################################################################################
#                  Do the inference
################################################################################

save_dir <- "results"
save_prefix <- "run02_test_"


compmap <- create_map(compound_mapping)
zprior <- node_dt$PRIOR
U <- Diagonal(x = node_dt$UNC^2)
obs <- node_dt$OBS

# optimization schedule
optim_scheme <- list(
  list("truexs_avg", 10),
  list(c("truexs_avg", "relnormerr", "relstaterr"), 30),
  list(c("relnormerr", "relstaterr", "truexs_hires"), 30),
  list(unique(node_dt$NODE), 300)
)

zpost <- zprior
# zpost[node_dt[,grepl("^truexs_",NODE)]] <- 1
for (cur_stage in optim_scheme) {
  cat(paste0("############################\n",
             "optimizing the following nodes:\n",
             paste0(cur_stage[[1]], collapse=", "), "\n",
             "############################\n"))
  zref <- zpost
  nodesel <- cur_stage[[1]]
  maxiter <- cur_stage[[2]]
  adjust_idcs <- node_dt[NODE %in% nodesel, IDX]
  optres <- LMalgo(compmap, zprior, U, obs, print.info=TRUE, zref = zref,
                   adjust_idcs = adjust_idcs,
                   control=list(maxcount=maxiter, mincount=-3, tau=1e-6),
                   must.converge=FALSE)
  zpost <- optres$zpost
}

# save_file <- paste0(save_prefix, "optres.rda")
# save_path <- file.path(save_dir, save_file)
# saveRDS(optres, file = save_path)

# calculate relevant posterior covariance blocks

node_dt[, PRED := as.vector(compmap$propagate(zpost))]

#sel <- node_dt[NODE=="truexs_hires" & ENERGY > 7500 & ENERGY < 7700,IDX]
sel <- node_dt[NODE %in% c("relsyserr","truexs_avg") & ENERGY > 6500 & ENERGY < 9500,IDX]
postcov <- get_posterior_cov(compmap, optres$zpost, U, obs, sel, sel)
node_dt[sel, PREDUNC := sqrt(diag(postcov))]
