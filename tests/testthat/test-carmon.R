test_that("carmon works", {
    old_seed <- get0(".Random.seed", envir = .GlobalEnv)
    on.exit({
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
    })
    data(multi_omics_small)
    set.seed(42)
    c_obj <- carmon(multi_omics_small,
        nlambda_w = 3, nlambda_b = 3,
        nc = 3, verbose = FALSE
    )
    expect_equal(c(
        c_obj$marginals[1], round(c_obj$sel_icov[15, 11], 6),
        c_obj$sel_index_lw, unname(c_obj$report["central for"])
    ), c(
        "negative binomial",
        0.041411, 2, "be"
    ))
})

test_that("verbosity 2 and 0 work", {
    data(multi_omics_micro)
    data(multi_omics)
    expect_snapshot(c_obj <- carmon(multi_omics_micro,
        net_method = "correlation", cor_cutoff = 0.7, verbose = 2,
        plot = FALSE
    ))
    expect_silent(c_obj <- carmon(multi_omics,
        net_method = "correlation",
        cor_quant = 0.05, analyse = FALSE, verbose = 0,
        plot = FALSE
    ))
})

test_that("verbosity 1 works", {
    data(multi_omics_micro)
    expect_snapshot(c_obj <- carmon(multi_omics_micro,
        net_method = "correlation", cor_cutoff = 0.7, verbose = 1,
        plot = FALSE
    ))
})

test_that("print.carmon without centrality works", {
    data(multi_omics)
    c_obj <- carmon(multi_omics,
        net_method = "correlation", cor_quant = 0.05,
        analyse = FALSE, verbose = 0, plot = FALSE
    )
    expect_snapshot(print(c_obj))
})

test_that("print.carmon with centrality works", {
    data(multi_omics)
    c_obj <- carmon(multi_omics,
        net_method = "correlation", cor_quant = 0.05,
        verbose = 0, plot = FALSE
    )
    expect_snapshot(print(c_obj))
})
