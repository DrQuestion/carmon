test_that("centrality analysis and report making works", {
  expect_output(c_obj <- carmon(multi_omics, marginals = c("e", 0),
                         net_method = "correlation", cor_quant = 0.05,
                         plot = FALSE, verbose = 1))
  expect_no_error(rep <- centrality_report(c_obj))
})

test_that("error and warning triggering works", {
  # Unexisting measure requested
  expect_error(carmon(multi_omics, marginals = c("e", 0),
                         net_method = "correlation", cor_quant = 0.05,
                         plot = FALSE, verbose = 0, c_measures = "r"))
  expect_error(carmon(multi_omics, marginals = c("e", 0),
                      net_method = "correlation", cor_quant = 0.05,
                      plot = FALSE, verbose = 0, c_measures = ""))
})
