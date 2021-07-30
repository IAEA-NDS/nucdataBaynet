set.seed(15)

truexs_dt <- data.table(node="true",
                        energy = 1:5,
                        data = 0,
                        obs = exp(1:5))

normerr_dt <- data.table(node="normerr",
                         energy = NA,
                         data = 0,
                         unc = 0.2,
                         obs = NA)

expxs_dt <- data.table(node="exp",
                       energy = truexs_dt$energy,
                       data = 0,
                       unc = 0.05,
                       obs = exp(1:5)*(1 + rnorm(1, sd=normerr_dt$unc[1])) +
                                         rnorm(5, sd=0.05))

modxs_dt <- data.table(node="mod",
                       energy = c(1,5),
                       data = 0,
                       unc = 1e8,
                       obs = NA)

modxs2_dt <- data.table(node="mod2",
                       energy = 1:5,
                       data = 0,
                       unc = 0,
                       obs = NA)


dt <- rbindlist(list(modxs_dt, modxs2_dt, normerr_dt, expxs_dt), use.names=TRUE)
dt[, IDX := seq_len(.N)]


linmod_params <- list(
  mapname = "linmod_map",
  src_idx = dt[node=="mod",IDX],
  tar_idx = dt[node=="mod2",IDX],
  src_x = dt[node=="mod",energy],
  tar_x = dt[node=="mod2",energy]
)


nonlinear_params <- list(
  mapname = "nonlinear_map",
  src_idx = dt[node=="mod2",IDX],
  tar_idx = dt[node=="exp",IDX],
  funname = "exp"
)


normerr_params <- normerr_map <- list(
  mapname = "normerr_map",
  src_idx = dt[node=="normerr",IDX],
  tar_idx = dt[node=="exp",IDX],
  src_feat = 1,
  tar_feat = dt[node=="exp", rep(1, .N)]
)


relnormerr_params <- list(
  mapname = "product_map",
  maps = list(nonlinear_params, normerr_params)
)

compound_params <- list(
  mapname = "compound_map",
  maps = list(linmod_params, nonlinear_params, relnormerr_params)
)

comp_map <- create_compound_map()
comp_map$setup(compound_params)


test_that("compound map jacobian corresponds to numeric jacobian", {
  myinp <- runif(nrow(dt))
  res1 <- jacobian(comp_map$propagate, myinp, with.id=TRUE)
  res2 <- as.matrix(comp_map$jacobian(myinp, with.id=TRUE))
  dimnames(res1) <- NULL
  dimnames(res2) <- NULL
  expect_equal(res1, res2)
})


test_that("LM algorithm converges to same solution from different starting points", {
  zprior <- dt$data
  U <- Diagonal(x = dt$unc^2)
  zstart1 <- zstart2 <- zprior
  sel <- dt[node %in% c("mod", "normerr"),IDX]
  zstart1[sel] <- c(4, 2, 1)
  zstart2[sel] <- c(1, 3, 3)
  curres1 <- LMalgo(comp_map, zprior, U, dt$obs, zref=zstart1)
  curres2 <- LMalgo(comp_map, zprior, U, dt$obs, zref=zstart2)
  expect_equal(curres1$final_val, curres2$final_val)
})
