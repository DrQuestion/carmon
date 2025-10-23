test_that("centrality analysis and report making works", {
    data(multi_omics)
    suppressMessages(expect_message(c_obj <- carmon(multi_omics,
        marginals = c("e", 0),
        net_method = "correlation", cor_quant = 0.05, plot = FALSE,
        verbose = 1
    )))
    expect_no_error(rep <- centrality_report(c_obj))
})

test_that("error and warning triggering works", {
    # Unexisting measure requested
    data(multi_omics)
    expect_error(carmon(multi_omics,
        marginals = c("e", 0),
        net_method = "correlation", cor_quant = 0.05, plot = FALSE,
        verbose = 0, c_measures = "r"
    ))
    expect_error(carmon(multi_omics,
        marginals = c("e", 0),
        net_method = "correlation", cor_quant = 0.05, plot = FALSE,
        verbose = 0, c_measures = ""
    ))
})
