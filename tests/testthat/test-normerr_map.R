params <- list(mapname = "normerr_map",
               src_idx = 1:3,
               tar_idx = 6:10,
               src_feat = c(1,2,3),
               tar_feat = c(1,1,3,2,2))


test_that("normerr_map propagation works correctly without identity matrix", {
  map <- create_normerr_map()
  map$setup(params)
  res <- map$propagate(c(5,9,2,0,0,0,0,0,0,0), with.id=FALSE)
  expect_equal(res[6:10], c(5,5,2,9,9))
  expect_equal(res[1:5], rep(0,5))
})


test_that("normerr_map propagation works correctly with identity matrix", {

  map <- create_normerr_map()
  map$setup(params)
  res <- map$propagate(c(5,9,2,0,0,0,0,0,0,0), with.id=TRUE)
  expect_equal(res[6:10], c(5,5,2,9,9))
  expect_equal(res[1:5], c(5,9,2,0,0))
})
