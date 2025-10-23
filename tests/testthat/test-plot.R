test_that("highlighting multiple central nodes works", {
    data(multi_omics)
    expect_no_error(carmon(multi_omics,
        net_method = "correlation",
        cor_quant = 0.05, verbose = FALSE, plot_node_labels = FALSE
    ))
})

test_that("empty network warnings work for both plot functions", {
    data(multi_omics)
    expect_warning(c_obj <- carmon(multi_omics,
        net_method = "correlation", cor_cutoff = 0.9999,
        verbose = FALSE, analyse = FALSE
    ), "The reconstructed network is empty, it cannot be plotted.")
    expect_warning(plot_report(c_obj),
                "The reconstructed network is empty, cannot plot the results of
    the centrality analysis.")
})
