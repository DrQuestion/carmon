test_that("layers check function returns layers unchanged if rownames are
          identical and nrow is equal",
    {
        m1 <- matrix(seq(1, 4), nrow = 2, dimnames = list(c("A", "B"), NULL))
        m2 <- matrix(seq(5, 8), nrow = 2, dimnames = list(c("A", "B"), NULL))
        layers <- list(m1, m2)
        result <- check_layers_dims(layers)
        expect_equal(result, layers)
    })

test_that("layers check function reorders layers if rownames match but are in
          different order",
    {
        m1 <- matrix(seq(1, 4), nrow = 2, dimnames = list(c("A", "B"), NULL))
        m2 <- matrix(seq(5, 8), nrow = 2, dimnames = list(c("B", "A"), NULL))
        layers <- list(m1, m2)
        result <- check_layers_dims(layers)
        expect_equal(rownames(result[[2]]), c("A", "B"))
    })

test_that("layers check function stops if rownames differ across layers", {
    m1 <- matrix(seq(1, 4), nrow = 2, dimnames = list(c("A", "B"), NULL))
    m2 <- matrix(seq(5, 8), nrow = 2, dimnames = list(c("C", "D"), NULL))
    layers <- list(m1, m2)
    expect_error(check_layers_dims(layers), "Samples should be along rows")
})

test_that("layers check function uses colnames if ncol matches but not nrow,
          returns just transposed if order same",
    {
        m1 <- matrix(seq(1, 4), nrow = 2, dimnames = list(NULL, c("A", "B")))
        m2 <- matrix(seq(5, 10), nrow = 3, dimnames = list(NULL, c("A", "B")))
        layers <- list(m1, m2)
        expect_warning(result <- check_layers_dims(layers),
        "samples will be assumed to be distributed along the
    columns")
        expect_equal(result, lapply(layers, t))
    })

test_that("layers check function uses colnames, reorders if same colnames in
          different order and transposes them",
    {
        m1 <- matrix(seq(1, 4), nrow = 2, dimnames = list(NULL, c("A", "B")))
        m2 <- matrix(seq(5, 10), nrow = 3, dimnames = list(NULL, c("B", "A")))
        layers <- list(m1, m2)
        expect_warning(result <- check_layers_dims(layers),
        "samples will be assumed to be distributed along the
    columns")
        expect_equal(rownames(result[[2]]), c("A", "B"))
    })

test_that("layers check function stops if neither nrow nor ncol matches",
    {
        m1 <- matrix(seq(1, 6), nrow = 2)
        m2 <- matrix(seq(7, 12), nrow = 3)
        layers <- list(m1, m2)
        expect_error(check_layers_dims(layers),
            "Dimensions do not match among layers")
    })

test_that("layers number check function stops if number is inconsistent", {
    layers <- list("a", "b")
    omics <- c("a", "b")
    marginals <- c("a", "b", "c")
    expect_error(suppressMessages(check_layers_num(layers, omics = omics,
                                                 marginals = marginals)),
        "Number of layers is inconsistent across input parameters. Please
    check your input.")
})

test_that("omics check function raises warnings for non implemented omics and
          returns FALSE",
    {
        expect_warning(check <- check_omics(omics = c("foo-seq", "baromics"),
                                            marginals = c(0, 0)))
        expect_false(check)
    })

test_that("layers name check function works when renaming duplicates", {
    layers <- list()
    layers[[1]] <- "foo"
    layers[[2]] <- "bar"
    omics <- c("a", "a")
    expect_equal(names(check_layers_names(layers, omics)), c("a_A", "a_B"))
})

test_that("unnamed features are renamed after their omics type",
    {
        data(multi_omics_micro)
        dataset <- multi_omics_micro
        colnames(dataset[[1]]) <- NULL
        dataset <- check_gen_colnames(dataset)
        actual_names <- unlist(lapply(dataset, colnames))
        expected_names <- c("rnaseq_1", "rnaseq_2",
            colnames(multi_omics_micro[[2]]))
        names(expected_names) <- names(actual_names)
        expect_equal(actual_names, expected_names)
    })

test_that("which_omics works", {
    expect_snapshot(which_omics())
})

test_that("which_marginals works", {
    expect_snapshot(which_marginals())
})

test_that("get_edge_weights works with correlation and mb", {
    data(multi_omics_micro)
    with_cor <- carmon(multi_omics_micro, net_method = "correlation",
        cor_quant = 0.75, plot = FALSE, analyse = FALSE, verbose = FALSE)
    with_mb <- carmon(multi_omics_micro, net_method = "mb", nlambda = 1,
        rep.num = 2, plot = FALSE, analyse = FALSE, verbose = FALSE)
    expect_no_error(get_edge_weights(with_cor))
    expect_no_error(get_edge_weights(with_mb))
})

test_that("assembling a carmon_obj from carmon_cop and carmon_rec works",
    {
        data(multi_omics_micro)
        carmon_cop_obj <- copulize(multi_omics_micro, verbose = FALSE)
        carmon_rec_obj <- reconstruct(carmon_cop_obj$layers,
            net_method = "correlation", cor_quant = 0.5, verbose = FALSE)
        carmon_obj <- assemble_carmon_obj(carmon_cop_obj, carmon_rec_obj)
        expect_true(inherits(carmon_obj, "carmon"))
    })
