linmod_params <- list(
  maptype = "linearinterpol_map",
  src_idx = 1:4,
  tar_idx = 6:10,
  src_x = 11:14,
  tar_x = c(12,12.5,13,13.5,14)
)


normerr_params <- list(maptype = "normerr_map",
               src_idx = 5,
               tar_idx = 6:10,
               src_feat = 1,
               tar_feat = c(1,1,1,1,1))


prod_params <- list(maptype = "product_map",
                    maps = list(linmod_params, normerr_params))


test_that("product map propagates properly with with.id=FALSE", {
  prodmap <- create_product_map()
  prodmap$setup(prod_params)
  linmod_map <- create_linearinterpol_map()
  linmod_map$setup(linmod_params)
  normerr_map <- create_normerr_map()
  normerr_map$setup(normerr_params)
  inp <- 1:10
  res1 <- normerr_map$propagate(inp, with.id=FALSE)
  res2 <- linmod_map$propagate(inp, with.id=FALSE)
  expres <- res1 * res2
  res <- prodmap$propagate(inp, with.id=FALSE)
  expect_equal(res, expres)
})


test_that("product map propagates properly with with.id=TRUE", {
  prodmap <- create_product_map()
  prodmap$setup(prod_params)
  linmod_map <- create_linearinterpol_map()
  linmod_map$setup(linmod_params)
  normerr_map <- create_normerr_map()
  normerr_map$setup(normerr_params)
  inp <- 1:10
  res1 <- normerr_map$propagate(inp, with.id=FALSE)
  res2 <- linmod_map$propagate(inp, with.id=FALSE)
  expres <- inp + res1 * res2
  res <- prodmap$propagate(inp, with.id=TRUE)
  expect_equal(res, expres)
})


test_that("product map propagates properly if one map yields zero values with with.id=FALSE", {
  prodmap <- create_product_map()
  prodmap$setup(prod_params)
  linmod_map <- create_linearinterpol_map()
  linmod_map$setup(linmod_params)
  normerr_map <- create_normerr_map()
  normerr_map$setup(normerr_params)
  inp <- c(0,0, 3:10)
  res1 <- normerr_map$propagate(inp, with.id=FALSE)
  res2 <- linmod_map$propagate(inp, with.id=FALSE)
  expres <- res1 * res2
  res <- prodmap$propagate(inp, with.id=FALSE)
  expect_equal(res, expres)
})


test_that("product map evaluates jacobian correctly with with.id=FALSE", {
  prodmap <- create_product_map()
  prodmap$setup(prod_params)
  jacfun <- function(x) {
    inp <- rep(0, 10)
    inp[1:5] <- x
    res <- prodmap$propagate(inp, with.id=FALSE)[6:10]
    return(res)
  }
  expres <- jacobian(jacfun, 1:5)
  inp <- 1:10
  res <- as.matrix(prodmap$jacobian(inp, with.id=FALSE))[6:10,1:5]
  dimnames(expres) <- NULL
  dimnames(res) <- NULL
  expect_equal(res, expres)
})


test_that("product map evaluates jacobian correctly with with.id=TRUE", {
  prodmap <- create_product_map()
  prodmap$setup(prod_params)
  jacfun <- function(x) {
    res <- prodmap$propagate(x, with.id=TRUE)
    return(res)
  }
  inp <- 1:10
  expres <- jacobian(jacfun, inp)
  res <- as.matrix(prodmap$jacobian(inp, with.id=TRUE))
  dimnames(expres) <- NULL
  dimnames(res) <- NULL
  expect_equal(res, expres)
})


test_that("product map evaluates jacobian correctly if zero elements present with with.id=FALSE", {
  prodmap <- create_product_map()
  prodmap$setup(prod_params)
  jacfun <- function(x) {
    inp <- rep(0, 10)
    inp[1:5] <- x
    res <- prodmap$propagate(inp, with.id=FALSE)[6:10]
    return(res)
  }
  inp <- c(0,0,3:10)
  expres <- jacobian(jacfun, inp[1:5])
  res <- as.matrix(prodmap$jacobian(inp, with.id=FALSE))[6:10,1:5]
  dimnames(expres) <- NULL
  dimnames(res) <- NULL
  expect_equal(res, expres)
})
