test_that("carmon works", {
  old_seed <- get0(".Random.seed", envir = .GlobalEnv)
  on.exit({
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  })
  set.seed(42)
  c_obj <- carmon(multi_omics_small, nlambda_w = 3, nlambda_b = 3, nc = 3, verbose = FALSE)
  expect_equal(c(c_obj$marginals[1], c_obj$sel_icov[15,11], c_obj$sel_index_lw, length(c_obj$report)),
               c("negative binomial", 0.0414107035728185, 2, 6))
})

test_that("verbosity 2 and 0 work", {
  expect_snapshot(c_obj <- carmon(multi_omics_micro, net_method = "correlation", cor_cutoff = 0.7, verbose = 2, plot = FALSE))
  expect_silent(c_obj <- carmon(multi_omics, net_method = "correlation", cor_quant = 0.05, analyse = FALSE, verbose = 0, plot = FALSE))
})

test_that("verbosity 1 works", {
  expect_snapshot(c_obj <- carmon(multi_omics_micro, net_method = "correlation", cor_cutoff = 0.7, verbose = 1, plot = FALSE))
})

test_that("print.carmon without centrality works", {
  c_obj <- carmon(multi_omics, net_method = "correlation", cor_quant = 0.05, analyse = FALSE, verbose = 0, plot = FALSE)
  expect_snapshot(print(c_obj))
})

test_that("print.carmon with centrality works", {
  c_obj <- carmon(multi_omics, net_method = "correlation", cor_quant = 0.05, verbose = 0, plot = FALSE)
  expect_snapshot(print(c_obj))
})
