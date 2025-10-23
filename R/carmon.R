#' Perform copula-aided reconstruction, analysis and plot of a multi-omics
#' network
#'
#' Single wrapper function encompassing all the functional units of the
#' \pkg{carmon} package: copula-based transition of non-normal data to the
#' normal realm, network reconstruction and selection, consensus centrality
#' analysis of the network to identify key omics features, and enriched plot of
#' the network and of the analysis results.
#'
#' @param layers The omics layers to analyze. Preferably provided as a named R
#'   list of \strong{non-normalized} omics data sets. If possible the names of
#'   the list
#'   should correspond to the respective omics type. If not possible, for
#'   example because of two layers from the same technology, please provide the
#'   omics types with the parameter `omics`.  To see a list of available
#'   omics types use the function `which_omics()`.\cr
#'   Each data set should be source-matched (same amount of matched samples or
#'   individuals across each data set). Placing of the samples (or individuals)
#'   should also be consistent: either along the rows for \emph{all} the data
#'   sets, or along the columns for \emph{all} the data sets, nothing in
#'   between. All the samples (or individuals) should also have consistent
#'   naming across the data sets.\cr
#'   `layers` can also be a single unified data set, but then it is necessary to
#'   specify the argument `p`.
#' @param p  Necessary \strong{only} in case `layers` is a single data set.
#'   A vector with with the number of variables for each omic layer of the
#'   data set (e.g. the number of transcripts, metabolites etc.), in the same
#'   order the layers have in the data set. If the sum of the elements in `p`
#'   is less than the real number of features in the data set, `carmon()`
#'   assumes that there is an omics layer more, with as many features as the
#'   real number of features minus `sum(p)`.
#' @param omics Highly recommended. A character vector of as many elements as
#'   the number of omics layers, naming which omics types each layer contains,
#'   in the same order as provided in the input `layers`
#'   (\emph{e.g.} `omics = c('RNA-seq', 'proteomics', 'metabolomics')`). To
#'   see a list of terms and omics technologies for which \pkg{carmon} is
#'   specifically tailored, use the function `which_omics()`.
#' @param marginals Optional, to be specified when the user prefers to use
#'   different marginal distributions than the default distribution \pkg{carmon}
#'   tailored for each omics layers. A vector of as many elements as the number
#'   of layers, specifying which marginal distribution should be used for each
#'   omics layer, in the same order as provided in the input `layers`. To see a
#'   list of available marginal distributions, use the function
#'   `which_marginals()`. For a mixed custom setting, place a `0` in the vector
#'   in the position corresponding to the omics layers for which the default
#'   distribution is desired. Otherwise, specify the desired marginal
#'   distribution (\emph{e.g.} `marginals = c(0, 'lognormal', 0)` means
#'   'default, lognormal, default').
#' @param noninv_method A placeholder for future functionalities of
#'   \pkg{carmon}, do not use.
#' @param copula A placeholder for future functionalities of \pkg{carmon}, do
#'   not use.
#' @param net_method The network reconstruction method to use. The four methods
#'   currently available are: `'coglasso'` for
#'   \emph{collaborative graphical lasso},
#'   `'glasso'` for \emph{graphical lasso}, `'mb'` for the Meinshausen-Buhlmann
#'   neighborhood selection, and `'correlation'` for a thresholded Pearson's
#'   correlation network. The default method is `'coglasso'`. See References for
#'   more information on the network reconstruction methods.
#' @param sel_method The network selection method. Each reconstruction procedure
#'   has its own set of model selection procedures available. See References for
#'   more information on the network selection methods. Here is the list of
#'   them:
#'   * for `net_method = 'coglasso'`: `'xstars'` for \emph{eXtended StARS},
#'   `'xestars'` for \emph{eXtended  Efficient StARS}, and `'ebic'` for Extended
#'   Bayesian Information Criterion. Default is `'xestars'`.
#'   * for `net_method = 'glasso'`: `'stars'` for \emph{StARS}, `'ric'` for
#'   Rotational Information Criterion, and `'ebic'` for Extended Bayesian
#'   Information Criterion. Default is `'stars'`.
#'   * for `net_method = 'mb'` or `'correlation'`: `'stars'` for \emph{StARS},
#'   and `'ric'` for Rotational Information Criterion. Default is `'stars'`.
#' @param ... The additional optional arguments to be given for the network
#'   reconstruction and selection procedure. The available arguments depend on
#'   the chosen network reconstruction and network selection methods. If using
#'   `'coglasso'` (default option), see [coglasso::bs()] for both reconstruction
#'   and selection arguments. If reconstructing with `'glasso'`, `'mb'`, or
#'   `'correlation'` without setting `cor_cutoff` or `cor_quant` see,
#'   respectively, [huge::huge.glasso()], [huge::huge.mb()], or
#'   [huge::huge.ct()]. For additional options for network selection from these
#'   last three reconstruction methods, see [huge::huge.select()].
#' @param cor_cutoff Optional for `net_method = 'correlation'`. The cutoff
#'   value for the absolute Pearson's correlation network. Any edge with an
#'   absolute correlation below the cutoff is excluded from the final network.
#'   Not used by default, as the cutoff is generated by internal calculations.
#' @param cor_quant Optional for `net_method = 'correlation'`, to set as an
#'   alternative to `cor_cutoff`. Determine the cutoff of correlation based on
#'   the top percentile indicated by the user. For example, `cor_quant = 0.2`
#'   would set as a correlation cutoff the 20th percentile of the absolute
#'   Pearson's correlation values, ordered from highest to lowest. Not used by
#'   default, as the cutoff is generated by internal calculations.
#' @param minimal_output Logical. Set to `TRUE` to get only a minimal output
#' from
#'   the network reconstruction and selection module of `carmon()`, which mainly
#'   differs on the basis of on the network reconstruction method set with
#'   `net_method.` Defaults to `FALSE`.
#' @param analyse Logical. Whether to perform the consensus centrality analysis
#'   to identify key omics features. The highlighted features are important
#'   based on a consensus of multiple centrality measures. Defaults to `TRUE`.
#' @param c_measures A string of characters, each one representing one of four
#'   possible centrality measures implemented in the consensus centrality
#'   analysis: `'d'` stands for degree centrality, `'b'` for betweenness
#'   centrality, `'c'` for closeness centrality, and `'e'` for eigenvector
#'   centrality. Default is `'dbce'`: all the four measures.
#' @param max_candidates_c_measures What is the highest amount of nodes that can
#'   be highlighted by each measure of centrality, \emph{before} finding the
#'   consensus? Default is 20. When given together with the
#'   `quantile_c_measures`
#'   argument, it overrides the `quantile_c_measures` argument when
#'   `max_candidates_c_measures` is smaller than  the chosen quantile, it is
#'   overrode by `quantile_c_measures` in the opposite case.
#' @param quantile_c_measures What is the top percentile of nodes that can be
#'   highlighted by each measure of centrality, \emph{before} finding the
#'   consensus? Default is the top 5%. When given together with the
#'   `max_candidates_c_measures`, it overrides the `max_candidates_c_measures`
#'   argument when the amount of nodes in the chosen top percentile is smaller
#'   than the chosen maximum amount of candidate nodes, it is overrode by
#'   `max_candidates_c_measures` in the opposite case.
#' @param scaled_c_measures Logical. Whether to compute centrality measures as
#'   0-1 scaled values. Defaults to `TRUE`.
#' @param plot Logical. Whether to plot the multi-omics network (enriched by the
#'   results of the centrality analysis when that is performed), and the results
#'   of the centrality analysis (when performed). Defaults to `TRUE`.
#' @param plot_node_labels Show node names in the plot of the network. Defaults
#'   to `TRUE`.
#' @param hide_isolated Hide from the plot nodes that are not connected to any
#'   other node. Defaults to `TRUE`.
#' @param plot_hot_nodes Highlight in the plot nodes that are selected as
#'   consensus-central nodes by the node centrality analysis, by drawing them
#'   bigger and adding the criteria for which they result to be central. The
#'   larger the consensus, the bigger the nodes are plotted. Defaults to `TRUE`.
#' @param verbose The level of verbosity of the `carmon()` process. `0`
#'   suppresses the information output, while `1` and `2` give progressively
#'   increasing amounts of information about the inner computations happening
#'   inside the package's modules. Defaults to `1`
#'
#' @returns `carmon()` returns an object of `S3` class `carmon`. The elements of
#' this object depend on the chosen network reconstruction (`'coglasso'`, by
#' default) and selection (`'xestars'`, by default) strategies, and on whether
#' or not the consensus centrality analysis is carried out (yes, by default).
#' \cr Every possible running mode of `carmon()` produces a `carmon` object
#' with the following elements:
#' \itemize{
#'   \item `layers` is an R list, each element being a data set of the
#'     corresponding omics layer, already copulized and transferred to the
#'     normal realm.
#'   \item `omics` is a vector containing the omics type assigned to each omics
#'     layer.
#'   \item `marginals` is a vector containing the marginal distributions used to
#'     transfer each omics layer to the normal realm.
#'   \item `sel_adj` is the adjacency matrix of the final selected \pkg{carmon}
#'     network.
#'   \item `net_method` is the chosen network reconstruction method.
#'   \item `sel_method` is the chosen model selection method.
#'   \item `network` is the `igraph` network object of the final selected
#'     \pkg{carmon} network.
#'   \item `call` is the matched call.
#' }
#'
#' Depending on the different network reconstruction and selection strategies
#' chosen, some of the other elements may change. When choosing `'coglasso'` or
#' `'glasso'` as network reconstruction method, an additional returned element
#' is:
#' \itemize{
#'   \item `sel_icov` is the inverse covariance matrix of the final selected
#'     network.
#' }
#' Otherwise, if choosing 'correlation', `carmon()` also returns:
#' \itemize{
#'   \item `cor` is the Pearson's correlation matrix of the final selected
#'     network.
#' }
#' When `minimal_output` is set to `TRUE`, these are also the only elements
#' returned in the object of `S3` class `carmon` given by `carmon()`. For
#' the non minimal output elements resulting from network reconstruction and
#' selection, please look at [coglasso::bs()] for `net_method = 'coglasso'`.
#' When using `'glasso'`, `'mb'`, and when using `'correlation'` without setting
#' any value for `cor_cutoff` or `cor_quant`, look at [huge::huge()] and at
#' [huge::huge.select()].\cr
#'
#' Moreover, when `analyse = TRUE` (default behaviour), the `carmon` object
#' contains two additional elements resulting from the consensus centrality
#' analysis:
#' \itemize{
#' \item `report` an R data frame. The rows correspond to the nodes identified
#'   to be central by the analysis, and they are ordered based on how large is
#'   the consensus among the different measures. The data frame has 6
#'   columns: \cr
#'   \emph{candidate}, the name of the central node; \cr
#'   \emph{degree}, the degree centrality of the node; \cr
#'   \emph{betweenness}, the betweenness centrality of the node; \cr
#'   \emph{closeness}, the closeness centrality of the node; \cr
#'   \emph{eigenvector}, the eigenvector centrality of the node; and \cr
#'   \emph{central for}, the a string reporting the first letter of all the
#'   measures according to which the node is central for.
#' \item `measures_list` an R named list of as many elements as the number of
#'   chosen centrality measures, the name of each element being the associated
#'   centrality measure. Each element is a named numerical vector, containing
#'   the measures of the top central nodes identified in the analysis. The
#'   name of each element of the vectors is the name of the node associated to
#'   the reported measure.
#' }
#'
#' @references For \emph{collaborative graphical lasso} and
#' \emph{eXtended StARS}, see
#' \href{https://doi.org/10.48550/arXiv.2403.18602
#' }{Albanese \emph{et al.} (2024)}
#' @references For \emph{graphical lasso}, see
#' \href{https://doi.org/10.1093/biostatistics/kxm045
#' }{Friedman \emph{et al.} (2007)}
#' @references For neighborhood selection, see
#' \href{https://doi.org/10.1214/009053606000000281
#' }{Meinshausen and Buhlmann (2006)}
#' @references For the \emph{StARS} selection procedure, see
#' \href{https://proceedings.neurips.cc/paper_files/paper/2010/file/
#' 301ad0e3bd5cb1627a2044908a42fdc2-Paper.pdf}{Liu \emph{et al.} (2010)}
#'
#'
#' @export
#'
#' @examples
#' # Suggested usage: provide input data as a named R list, each element being
#' # the data set of an omics layer, then let carmon do the rest. See the
#' # vignettes for a more custom usage!
#' data(multi_omics_small)
#' c_obj <- carmon(multi_omics_small, verbose = FALSE)
#'
carmon <- function(
    layers, ..., p = NULL, omics = NULL,
    marginals = NULL, noninv_method = NULL, copula = NULL,
    net_method = "coglasso", sel_method = NULL, cor_cutoff = NULL,
    cor_quant = NULL, minimal_output = FALSE, analyse = TRUE,
    c_measures = "dbce", max_candidates_c_measures = NULL,
    quantile_c_measures = NULL, scaled_c_measures = TRUE,
    plot = TRUE, plot_node_labels = TRUE, hide_isolated = TRUE,
    plot_hot_nodes = TRUE, verbose = 1) {

    carmon_obj <- copulize(layers, p = p, omics = omics, marginals = marginals,
            noninv_method = noninv_method, copula = copula, verbose = verbose)

    reconstruct_args <- list(layers = carmon_obj$layers,
            net_method = net_method, sel_method = sel_method,
            cor_cutoff = cor_cutoff, cor_quant = cor_quant,
            minimal_output = minimal_output, verbose = verbose)
    reconstruct_args <- c(reconstruct_args, list(...))

    carmon_rec <- do.call(reconstruct, reconstruct_args)
    carmon_obj <- c(carmon_obj, carmon_rec)
    carmon_obj$network <- get_network(carmon_obj)

    carmon_obj$call <- match.call()
    carmon_obj$reconstruct_call <- NULL
    carmon_obj$copulize_call <- NULL
    class(carmon_obj) <- "carmon"

    if (analyse) {
        carmon_obj <- compute_centrality(carmon_obj, c_measures,
            max_candidates_c_measures, quantile_c_measures,
            scaled_c_measures, verbose)
    }
    if (plot) {
        if (!is.null(carmon_obj$measures_list)) {
            plot_report(carmon_obj, scaled_c_measures)
        }
        plot(carmon_obj, plot_node_labels, hide_isolated, plot_hot_nodes)
    }
    return(carmon_obj)
}

#' Print function for the S3 class `carmon`
#'
#' Print information on the selected carmon network and the explored
#' hyperparameters and see next suggested step.
#'
#' @param x is the object of `S3` class `carmon`.
#' @param ... system required, not used.
#'
#' @noRd
#' @export
print.carmon <- function(x, ...) {
    carmon_name <- rlang::call_args(match.call())
    if (!is.null(x$sel_method)) {
        cat("Carmon network estimated with ", x$net_method,
            " and selected with ", x$sel_method, ".\n\n", sep = "")
    } else {
        cat("Carmon network estimated with ", x$net_method,
            " and selected with ", x$sel_method, ".\n\n", sep = "")
    }
    cat("The call was:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    cat("******************************************************\n\n")
    # cat('The model selection method was:\n', x$method, '\n', sep = '')
    # cat('The density of the selected network is:\n', x$sel_density, '\n\n',
    # sep = '')
    if (length(x$omics) > 2) {
        mes <- paste(x$omics[seq(1, length(x$omics) - 1)], collapse = ", ")
        mes <- paste(mes, ", and ", x$omics[length(x$omics)], sep = "")
    } else {
        mes <- paste(x$omics, collapse = " and ")
    }
    cat("The network is made of ", length(x$omics), " omics layers: ", mes,
        ".\n", sep = "")
    cat("It has a total of ", sum(vapply(x$layers,ncol, numeric(1))),
        " nodes.\n", sep = "")
    mes <- ""
    for (i in seq_len(length(x$omics))) {
        layer_name <- x$omics[i]
        mes_ <- paste("The ", layer_name, " layer has ", ncol(x$layers[[i]]),
            " nodes; \n    its chosen marginal distribution is the ",
            x$marginals[i], ".\n", sep = "")
        mes <- paste(mes, mes_, sep = "")
    }
    cat(mes)
    cat("\n")
    cat("******************************************************\n\n")
    if (is.null(x$report)) {
        cat("Find central nodes with:\n", carmon_name[[1]],
            " <- compute_centrality(", carmon_name[[1]], ")\n", sep = "")
        cat("Then p")
    } else {
        cat("P")
    }
    cat("rint a report of the central nodes with:\ncentrality_report(",
        carmon_name[[1]], ")\n\n", sep = "")
    cat("Plot a comparison of the central nodes with:\nplot_report(",
        carmon_name[[1]], ")\n\n", sep = "")
    cat("Plot the carmon network with:\nplot(", carmon_name[[1]], ")\n",
        sep = "")
}
