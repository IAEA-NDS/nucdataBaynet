map_creators <- get_map_generators()
map_creator_names <- names(map_creators)


test_that("all mappings have required methods", {
  for (idx in seq_along(map_creators)) {
    curname <- map_creator_names[idx]
    curmap <- map_creators[[idx]]()
    expect_true(all(c("setup", "getName", "getDescription", "getType", "get_src_idx", "get_tar_idx",
                  "propagate", "jacobian") %in% names(curmap)),
                label = paste0("mapping ", curname, " has all required methods"))
  }
})


test_that("all mappings have correct formals in propagate method", {
  for (idx in seq_along(map_creators)) {
    curname <- map_creator_names[idx]
    curmap <- map_creators[[idx]]()
    expect_true(all(c("x", "with.id") %in% names(formals(curmap$propagate))),
                label = paste0("propagate method of ", curname, " has required arguments"))
    expect_true(formals(curmap$propagate)[["with.id"]],
                label = paste0("with.id argument of propagate method of ", curname, " defaults to TRUE"))
  }
})


test_that("all mappings have correct formals in jacobian method", {
  for (idx in seq_along(map_creators)) {
    curname <- map_creator_names[idx]
    curmap <- map_creators[[idx]]()
    expect_true(all(c("x", "with.id") %in% names(formals(curmap$jacobian))),
                label = paste0("jacobian method of ", curname, " has required arguments"))
    expect_true(formals(curmap$propagate)[["with.id"]],
                label = paste0("with.id argument of jacobian method of ", curname, " defaults to TRUE"))
  }
})
