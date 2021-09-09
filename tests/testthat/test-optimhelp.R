set.seed(23)

# same definition as in test-LMalgo.R

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
  maptype = "linearinterpol_map",
  src_idx = dt[node=="mod",IDX],
  tar_idx = dt[node=="mod2",IDX],
  src_x = dt[node=="mod",energy],
  tar_x = dt[node=="mod2",energy]
)


nonlinear_params <- list(
  maptype = "nonlinear_map",
  src_idx = dt[node=="mod2",IDX],
  tar_idx = dt[node=="exp",IDX],
  funname = "exp"
)


normerr_params <- normerr_map <- list(
  maptype = "normerr_map",
  src_idx = dt[node=="normerr",IDX],
  tar_idx = dt[node=="exp",IDX],
  src_feat = 1,
  tar_feat = dt[node=="exp", rep(1, .N)]
)


relnormerr_params <- list(
  maptype = "product_map",
  maps = list(nonlinear_params, normerr_params)
)

compound_params <- list(
  maptype = "compound_map",
  maps = list(linmod_params, nonlinear_params, relnormerr_params)
)

comp_map <- create_compound_map()
comp_map$setup(compound_params)


test_that("numeric jacobian corresponds to analytic jacobian", {
  zprior <- dt$data
  U <- Diagonal(x = dt$unc^2)
  optimfuns <- generate_logpostpdf_funs(comp_map, zprior, U, dt$obs, zref = zprior)
  # test 1
  ztest <- runif(length(optimfuns$indep_idx))
  expres <- as.vector(jacobian(optimfuns$fun, ztest))
  res <- optimfuns$jac(ztest)
  expect_equal(res, expres)
  # test 2
  ztest <- runif(length(optimfuns$indep_idx))
  expres <- as.vector(jacobian(optimfuns$fun, ztest))
  res <- optimfuns$jac(ztest)
  expect_equal(res, expres)
  # test 3
  ztest <- runif(length(optimfuns$indep_idx))
  expres <- as.vector(jacobian(optimfuns$fun, ztest))
  res <- optimfuns$jac(ztest)
  expect_equal(res, expres)
})
