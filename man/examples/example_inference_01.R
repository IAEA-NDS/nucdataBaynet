library(Matrix)
params <- list(
  mapname = "mymap",
  maptype = "linearinterpol_map",
  src_idx = 1:10,
  tar_idx = 11:15,
  src_x = 1:10,
  tar_x = 3:7
)
mymap <- create_linearinterpol_map()
mymap$setup(params)

U <- Diagonal(n=15, x=c(rep(1e3, 10), rep(1, 5)))
zprior <- rep(0, 15)
zref <- rep(0, 15)
obs <- c(rep(NA,10), 5:9)

# glsalgo only works for linear relationships
# LMalgo can also deal with non-linear relationships
zpost <- glsalgo(mymap, zprior, U, obs)
optres <- LMalgo(mymap, zprior, U, obs)
zpost2 <- optres$zpost

# posterior estimates of values on computational grid
zpost[1:10]
# posterior estimates of error variables associated with observations
zpost[11:15]
# get posterior covariance block of independent variables
get_posterior_cov(mymap, zpost, U, obs, 1:5, 5:10)
# draw samples of independent variables from posterior distribution
get_posterior_sample(mymap, zpost, U, obs, 10)
