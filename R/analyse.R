#' Find key nodes in a `carmon` network
#'
#' Find the top central nodes according to multiple centrality measures of
#' choice Use a consensus of multiple centrality measures to identify nodes
#' that could be important in the network obtained through `carmon`
#'
#' @param carmon_obj The object of `S3` class `carmon`.
#' @param measures A string of characters, each one representing one of four
#'   possible centrality measures implemented: `'d'` stands for degree
#'   centrality, `'b'` for betweenness centrality, `'c'` for closeness
#'   centrality, and `'e'` for eigenvector centrality. Default is `'dbce'`, all
#'   the four measures.
#' @param max_candidates What is the highest amount of nodes that can be
#'   highlighted by each measure, \emph{before} finding the consensus? Default
#'   is 20. When given together with the `quant` argument, it overrides the
#'   `quant` argument when `max_candidates` is smaller than the chosen quantile,
#'   it is overrode by `quant` in the opposite case.
#' @param quant What is the top percentile of nodes that can be highlighted by
#'   each measure, \emph{before} finding the consensus? Default is the top 5%.
#'   When given together with the `max_candidates` argument, it overrides the
#'   `max_candidates` argument when the amount of nodes in the
#'   chosen top percentile is smaller than the chosen maximum amount of
#'   candidate nodes, it is overrode by `max_candidates` in the opposite case.
#' @param scaled Logical, whether to compute centrality measures as
#'   0-1 scaled values. Defaults to TRUE.
#' @param verbose The level of verbosity of the centrality analysis. `0`
#'   suppresses the information output, while `1` and `2` give progressively
#'   increasing amounts of information about the inner computations happening
#'   inside the analysis.
#'
#' @returns `compute_centrality()` returns an object of `S3` class `carmon`,
#'   consisting of the input `carmon_obj` enriched of two new elements:
#'   * `report` an R data frame. The rows correspond to the nodes identified to
#'     be central by the analysis, and they are ordered based on how large is
#'     the consensus among the different measures. The data frame has 6 or less
#'     columns, depending on whether the respective centrality measures have
#'     found any central node: \cr
#'     \emph{candidate}, the name of the central node; \cr
#'     \emph{degree}, the degree centrality of the node; \cr
#'     \emph{betweenness}, the betweenness centrality of the node; \cr
#'     \emph{closeness}, the closeness centrality of the node; \cr
#'     \emph{eigenvector}, the eigenvector centrality of the node; and \cr
#'     \emph{central for}, the a string reporting the first letter of all the
#'     measures according to which the node is central for.
#'   * `measures_list` an R named list of as many elements as the number of
#'     chosen centrality measures, the name of each element being the associated
#'     centrality measure. Each element is a named numerical vector, containing
#'     the measures of the top central nodes identified in the analysis. The
#'     name of each element of the vectors is the name of the node associated to
#'     the reported measure.
#' @export
#'
#' @examples
#' # compute_centrality() needs an object of S3 class 'carmon' as an input. The
#' # following line quickly obtains one.
#' data(multi_omics_small)
#' c_obj <- carmon(multi_omics_small,
#'     net_method = "correlation",
#'     cor_quant = 0.25, analysis = FALSE, verbose = FALSE
#' )
#' # Then to perform the consensus centrality analysis:
#' c_obj <- compute_centrality(c_obj)
#'
compute_centrality <- function(
    carmon_obj, measures = "dbce",
    max_candidates = NULL, quant = NULL, scaled = TRUE, verbose = FALSE) {
    if (verbose) {
        message("**************Beginning network analysis**************")
        utils::flush.console()
    }
    n_candidates <- determine_n_candidates(carmon_obj, max_candidates, quant)
    if (nchar(measures) > 0) {
        for (i in seq_len(nchar(measures))) {
            m <- substr(measures, i, i)
            if (!grepl(m, "dbce")) {
                stop(paste0("Unknown measure inserted: ", m, collapse = ""))
            }
        }
        if (verbose > 0 & verbose < 2) {
            message("Computing centrality measures....")
        }
        res <- compute_measures(carmon_obj, n_candidates, scaled, measures,
                                verbose)
        all_measures <- res$all_measures
        no_node_found <- which(res$no_node_found == 1)
        all_measures[no_node_found] <- NULL
        if (verbose > 0 & verbose < 2) {
            message("Centrality measures computed!")
        }
        if (length(no_node_found) == nchar(measures)) {
            warning("No node was found to be central. No centrality report
    will be generated.")
            if (verbose) {
                message(
                    "**************Network analysis complete***************")
                utils::flush.console()
            }
            return(carmon_obj)
        }
    } else {
        stop("Give at least a centrality measure for network analysis")
    }
    report <- compile_report(all_measures = all_measures, degrees = res$degrees,
                    betweenness = res$betweenness, closeness = res$closeness,
                    eigenvector = res$eigenvector, verbose = verbose)
    carmon_obj$report <- report
    carmon_obj$measures_list <- all_measures
    return(carmon_obj)
}

#' Compute the centrality measures
#'
#' This function works as a coordinator of the computation of each centrality
#' measure separately.
#'
#' @inheritParams compute_centrality
#' @param n_candidates The number of candidates defining the size of the top
#'     list for each centrality measure.
#'
#' @returns compute_measures() returns a named R list with the following
#'     elements:
#' \itemize {
#'     \item `all_measures` an R named list of as many elements as the number of
#'         chosen centrality measures, the name of each element being the
#'         associated centrality measure. Each element is a named numerical
#'         vector, containing the measures of the top central nodes identified
#'         in the analysis. The name of each element of the vectors is the name
#'         of the node associated to the reported measure.
#'     \item `no_node_found` Is a numerical vector of as many elements as the
#'         number of requested measures, in the same order as the input to the
#'         argument `measures` (or `c_measures` when calling directly as part of
#'         `carmon()`). It contains 1 if the corresponding measure did not find
#'         any node to be central, 0 otherwise.
#'     \item `degrees`, `betweenness`, `closeness`, and `eigenvector` are
#'         numerical vectors containing the corresponding centrality measure
#'         for each node of the network.
#' }
#'
#' @noRd
#'
compute_measures <- function(carmon_obj, n_candidates, scaled, measures,
                            verbose) {
    no_node_found <- rep(0, nchar(measures))
    all_measures <- vector(mode = "list", length = nchar(measures))
    degrees <- NULL; betweenness <- NULL; closeness <- NULL;
    eigenvector <- NULL
    for (i in seq_len(nchar(measures))) {
        all_measures[[i]] <- rep(0, n_candidates)
        m <- substr(measures, i, i)
        if (m == "d") {
            res <- compute_degree(carmon_obj = carmon_obj,
                all_measures = all_measures, n_candidates = n_candidates,
                scaled = scaled, i = i, measures = measures, verbose = verbose)
            all_measures <- res[[1]]
            degrees <- res[[2]]
        } else if (m == "b") {
            res <- compute_betweenness(carmon_obj = carmon_obj,
                all_measures = all_measures, n_candidates = n_candidates,
                scaled = scaled, i = i, measures = measures, verbose = verbose)
            all_measures <- res[[1]]
            betweenness <- res[[2]]
        } else if (m == "c") {
            res <- compute_closeness(carmon_obj = carmon_obj,
                all_measures = all_measures, n_candidates = n_candidates,
                scaled = scaled, i = i, measures = measures, verbose = verbose)
            all_measures <- res[[1]]
            closeness <- res[[2]]
        } else {
            res <- compute_eigenvector_c(carmon_obj = carmon_obj,
                all_measures = all_measures, n_candidates = n_candidates,
                scaled = scaled, i = i, measures = measures, verbose = verbose)
            all_measures <- res[[1]]
            eigenvector <- res[[2]]
        }
        if (length(all_measures[[i]]) == 0) {
            no_node_found[i] <- 1
            if (verbose) {
                message("No node was found to be central for ",
                                names(all_measures)[i], ".")
            }
        }
    }
    return(list("all_measures" = all_measures, "no_node_found" = no_node_found,
                "degrees" = degrees, "betweenness" = betweenness,
                "closeness" = closeness, "eigenvector" = eigenvector))
}

#' Compute the degree for each node of the input carmon network
#'
#' @inheritParams compute_centrality carmon_obj scaled measures verbose
#' @inheritParams compute_measures n_candidates
#' @param all_measures The R list of top central nodes, iteratively updated for
#'     each measure within the function `compute_measures()`.
#' @param i The index of the current centrality measure among those to compute.
#'
#' @returns It returns an R list of two elements: the updated `all_measures` R
#'     list (see arguments of the function); and `measure`, the vector of the
#'     degrees of each node in the network.
#'
#' @noRd
#'
compute_degree <- function(carmon_obj, all_measures, n_candidates, scaled, i,
                            measures, verbose) {
    if (verbose > 1) {
        mes <- paste(c(
            "Computing centrality measure ", i, " of ",
            nchar(measures), " (Degree)"
        ))
        message(mes)
    }
    names(all_measures)[i] <- "Degree"
    # maybe for degrees (or everything?) I should use the absolute
    # number, and row normalize the resulting table, it is hard to
    # compare conditions and fused otherwise
    degrees <- sort(igraph::degree(carmon_obj$network,
                                    normalized = scaled
    ), decreasing = TRUE)

    if(n_candidates >= length(degrees)) {
        n_candidates <- length(degrees)
    }

    d_limit_value <- degrees[n_candidates]

    if (sum(degrees >= d_limit_value) > n_candidates) {
        # Too many nodes to plot
        degrees_names <- names(degrees[degrees > d_limit_value])
        degrees_top <- degrees[degrees > d_limit_value]
        all_measures[[i]] <- rep(0, length(degrees_names))
    } else {
        degrees_names <- names(degrees[degrees >= d_limit_value])
        degrees_top <- degrees[degrees >= d_limit_value]
    }

    all_measures[[i]] <- degrees_top
    names(all_measures[[i]]) <- degrees_names

    return(list("all_measures" = all_measures, "measure" = degrees))
}

#' Compute the betweenness for each node of the input carmon network
#'
#' @inheritParams compute_centrality carmon_obj scaled measures verbose
#' @inheritParams compute_measures n_candidates
#' @param all_measures The R list of top central nodes, iteratively updated for
#'     each measure within the function `compute_measures()`.
#' @param i The index of the current centrality measure among those to compute.
#'
#' @returns It returns an R list of two elements: the updated `all_measures` R
#'     list (see arguments of the function); and `measure`, the vector of
#'     the betweenness value of each node in the network.
#'
#' @noRd
#'
compute_betweenness <- function(carmon_obj, all_measures, n_candidates, scaled,
                                i, measures, verbose) {
    if (verbose > 1) {
        mes <- paste(c(
            "Computing centrality measure ", i, " of ",
            nchar(measures), " (Betweenness) "
        ))
        message(mes)
    }
    names(all_measures)[i] <- "Betweenness"

    betweenness <- sort(igraph::betweenness(carmon_obj$network,
                                            normalized = scaled, weights = NA
    ), decreasing = TRUE)

    if(n_candidates >= length(betweenness)) {
        n_candidates <- length(betweenness)
    }

    b_limit_value <- betweenness[n_candidates]

    if (sum(betweenness >= b_limit_value) > n_candidates) {
        # Too many nodes to plot
        betweenness_names <- names(betweenness[betweenness > b_limit_value])
        betweenness_top <- betweenness[betweenness > b_limit_value]
        all_measures[[i]] <- rep(0, length(betweenness_names))
    } else {
        betweenness_names <- names(betweenness[betweenness >= b_limit_value])
        betweenness_top <- betweenness[betweenness >= b_limit_value]
    }

    all_measures[[i]] <- betweenness_top
    names(all_measures[[i]]) <- betweenness_names

    return(list("all_measures" = all_measures, "measure" = betweenness))
}

#' Compute the closeness for each node of the input carmon network
#'
#' @inheritParams compute_centrality carmon_obj scaled measures verbose
#' @inheritParams compute_measures n_candidates
#' @param all_measures The R list of top central nodes, iteratively updated for
#'     each measure within the function `compute_measures()`.
#' @param i The index of the current centrality measure among those to compute.
#'
#' @returns It returns an R list of two elements: the updated `all_measures` R
#'     list (see arguments of the function); and `measure`, the vector of
#'     the closeness value of each node in the network.
#'
#' @noRd
#'
compute_closeness <- function(carmon_obj, all_measures, n_candidates, scaled,
                                i, measures, verbose) {
    if (verbose > 1) {
        mes <- paste(c(
            "Computing centrality measure ", i, " of ",
            nchar(measures), " (Closeness) "
        ))
        message(mes)
    }
    names(all_measures)[i] <- "Closeness"

    closeness <- sort(igraph::closeness(carmon_obj$network,
                                        normalized = scaled, weights = NA
    ), decreasing = TRUE)

    if(n_candidates >= length(closeness)) {
        n_candidates <- length(closeness)
    }

    c_limit_value <- closeness[n_candidates]

    if (sum(closeness >= c_limit_value) > n_candidates) {
        # Too many nodes to plot
        closeness_names <- names(closeness[closeness > c_limit_value])
        closeness_top <- closeness[closeness > c_limit_value]
        all_measures[[i]] <- rep(0, length(closeness_names))
    } else {
        closeness_names <- names(closeness[closeness >= c_limit_value])
        closeness_top <- closeness[closeness >= c_limit_value]
    }

    all_measures[[i]] <- closeness_top
    names(all_measures[[i]]) <- closeness_names

    return(list("all_measures" = all_measures, "measure" = closeness))
}

#' Compute the eigenvector centrality for each node of the input carmon network
#'
#' @inheritParams compute_centrality carmon_obj scaled measures verbose
#' @inheritParams compute_measures n_candidates
#' @param all_measures The R list of top central nodes, iteratively updated for
#'     each measure within the function `compute_measures()`.
#' @param i The index of the current centrality measure among those to compute.
#'
#' @returns It returns an R list of two elements: the updated `all_measures` R
#'     list (see arguments of the function); and `measure`, the vector of
#'     the eigenvector centrality of each node in the network.
#'
#' @noRd
#'
compute_eigenvector_c <- function(carmon_obj, all_measures, n_candidates,
                                scaled, i, measures, verbose) {
    if (verbose > 1) {
        mes <- paste(c(
            "Computing centrality measure ", i, " of ", nchar(measures),
            " (Eigenvector) "))
        message(mes)
    }
    names(all_measures)[i] <- "Eigenvector Centrality"

    eigenvector <- sort(igraph::eigen_centrality(carmon_obj$network,
                                        weights = NA)$vector, decreasing = TRUE)

    if(n_candidates >= length(eigenvector)) {
        n_candidates <- length(eigenvector)
    }

    e_limit_value <- eigenvector[n_candidates]

    if (sum(eigenvector >= e_limit_value) > n_candidates) {
        # Too many nodes to plot
        eigenvector_names <- names(eigenvector[eigenvector > e_limit_value])
        eigenvector_top <- eigenvector[eigenvector > e_limit_value]
        all_measures[[i]] <- rep(0, length(eigenvector_names))
    } else {
        eigenvector_names <- names(eigenvector[eigenvector >= e_limit_value])
        eigenvector_top <- eigenvector[eigenvector >= e_limit_value]
    }

    all_measures[[i]] <- eigenvector_top
    names(all_measures[[i]]) <- eigenvector_names

    return(list("all_measures" = all_measures, "measure" = eigenvector))
}

#' Compile the centrality report from the computed centrality measures
#'
#' @param all_measures an R named list of as many elements as the number of
#'         chosen centrality measures, the name of each element being the
#'         associated centrality measure. Each element is a named numerical
#'         vector, containing the measures of the top central nodes identified
#'         in the analysis. The name of each element of the vectors is the name
#'         of the node associated to the reported measure.
#' @param degrees, betweenness, closeness, eigenvector The vectors containing
#'         the relative centrality measures for each node of the network.
#' @param verbose The verbosity of the report compiling.
#'
#' @returns It returns an R data frame. The rows correspond to the nodes
#'     identified to be central by the analysis, and they are ordered based on
#'     how large is the consensus among the different measures. The data frame
#'     has 6 or less columns, depending on whether the respective centrality
#'     measures have found any central node: \cr
#'     \emph{candidate}, the name of the central node; \cr
#'     \emph{degree}, the degree centrality of the node; \cr
#'     \emph{betweenness}, the betweenness centrality of the node; \cr
#'     \emph{closeness}, the closeness centrality of the node; \cr
#'     \emph{eigenvector}, the eigenvector centrality of the node; and \cr
#'     \emph{central for}, the a string reporting the first letter of all the
#'     measures according to which the node is central for.
#'
#' @noRd
#'
compile_report <- function(all_measures, degrees = NULL, betweenness = NULL,
                            closeness = NULL, eigenvector = NULL,
                            verbose = FALSE) {
    if (verbose) {
        message("Compiling the report....")
    }
    all_candidates <- unique(unlist(lapply(all_measures, names)))
    report <- matrix(0, nrow = length(all_candidates),
                    ncol = 4 + length(all_measures))
    for (i in seq_len(length(all_candidates))) {
        candidate <- all_candidates[i]
        report[i, 1] <- all_candidates[i]
        count <- 0
        string <- ""
        column <- 2
        which_measures <- rep(FALSE,4)
        measures_name <- c("Degree", "Betweenness", "Closeness",
                        "Eigenvector Centrality")
        measures <- list("d" = degrees, "b" = betweenness, "c" = closeness,
                            "e" = eigenvector)
        for (j in seq_along(measures_name)) {
            m <- measures_name[j]
            if(m %in% names(all_measures)) {
                if (candidate %in% names(all_measures[[m]])) {
                    string <- paste(string, tolower(substr(m, 1, 1)), sep = "")
                    count <- count + 1
                    report[i, column] <- paste(c(round(measures[[j]][candidate],
                                            digits = 4), "*"), collapse = "")
                } else {
                    report[i, column] <- as.character(
                        round(measures[[j]][candidate], digits = 4))
                }
                column <- column + 1
                which_measures[j] <- TRUE
            }
        }
        last_column_measure <- column - 1
        report[i, column] <- string
        report[i, column <- column + 1] <- count
        report[i, column <- column + 1] <- sum(as.numeric(strsplit(
                report[i, seq(2, last_column_measure)], split = "\\*")))
    }
    report <- format_report(report, column, which_measures)
    if (verbose) {
        message("Report compiled!")
        message("**************Network analysis complete***************")
        utils::flush.console()
    }
    return(report)
}

#' Format the centrality report as an R data frame or a named vector
#'
#' @param report The matrix containing the report generated inside
#'     compile_report().
#' @param column The index of the last column of the data frame or element of
#'     the vector.
#' @param which_measures A logical vector with TRUEs corresponding to those
#'     measures for which central nodes were found, FALSE otherwise.
#'
#' @returns It returns the formatted report, either as an R data.frame with a
#'     number of columns equal to 4 + the number of measures for which central
#'     nodes were found, or a named vector with the same number of elements.
#'
#' @noRd
#'
format_report <- function(report, column, which_measures) {
    report <- report[order(as.numeric(report[, column - 1]),
                            as.numeric(report[, column]), decreasing = TRUE), ]
    col_names <- c("degree", "betweenness", "closeness", "eigenvector")
    if (is.null(dim(report))) {
        # Only one candidate
        report <- report[-c(column - 1 , column)]
        names(report) <- c(
            "candidate", col_names[which(which_measures)], "central for")
    } else {
        report <- report[, -c(column - 1 , column)]
        colnames(report) <- c(
            "candidate", col_names[which(which_measures)], "central for")
        rownames(report) <- report[, 1]
        report <- as.data.frame(report)
    }
    return(report)
}

#' Display the table of results of the consensus centrality analysis of the
#' `carmon` network
#'
#' @param carmon_obj The object of `S3` class `carmon`, possibly already
#'   resulting from a \emph{full} `carmon()` run, or at least undergone
#'   centrality analysis through `compute_centrality()`.
#'
#' @returns `centrality_report()` returns the centrality report of the object of
#'   `S3` class `carmon` given as an input, in the form of an R data frame.
#'   The rows correspond to the nodes identified to be central by the analysis,
#'   and they are ordered based on how large is the consensus among the
#'   different measures. The data frame has 6 columns: \cr
#'   \emph{candidate}, the name of the central node; \cr
#'   \emph{degree}, the degree centrality of the node; \cr
#'   \emph{betweenness}, the betweenness centrality of the node; \cr
#'   \emph{closeness}, the closeness centrality of the node; \cr
#'   \emph{eigenvector}, the eigenvector centrality of the node; and \cr
#'   \emph{central for}, the a string reporting the first letter of all the
#'   measures according to which the node is central for.
#' @export
#'
#' @examples
#' # Let's build and analyse a carmon network:
#' data(multi_omics_small)
#' c_obj <- carmon(multi_omics_small,
#'     net_method = "correlation",
#'     cor_quant = 0.25, analysis = TRUE, plot = FALSE,
#'     # analysis is already TRUE by default
#'     verbose = FALSE
#' )
#' # To display the table of the results of the centrality analysis:
#' centrality_report(c_obj)
#'
centrality_report <- function(carmon_obj) {
    if (is.null(carmon_obj$report)) {
        warning("The centrality report was not compiled yet. \nCompiling it now,
    but it will not be saved in the carmon object.")
        carmon_obj <- compute_centrality(carmon_obj)
    }
    return(carmon_obj$report)
}

# Use stability of central measures as well. Perform subsampling to check how
# often that node is associated to a certain (for example) degree. Two ways:
# measure the degree in every subsample and draw a statistical distribution of
# the various subsamples to indicate how stable that is.  Otherwise, check in
# how many subsamples that node is in the top 5%, or 10 or whatever, nodes for
# that measure, and use that as a measure of consensus, somehow.
