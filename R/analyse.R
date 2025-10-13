#' Find key nodes in a `carmon` network
#'
#' Find the top central nodes according to multiple centrality measures of
#' choice Use a consensus of multiple centrality measures to identify nodes that could
#' be important in the network obtained through `carmon`
#'
#' @param carmon_obj The object of `S3` class `carmon`.
#' @param measures A string of characters, each one representing one of four
#'   possible centrality measures implemented: `"d"` stands for degree
#'   centrality, `"b"` for betweenness centrality, `"c"` for closeness
#'   centrality, and `"e"` for eigenvector centrality. Default is `"dbce"`, all
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
#'     the consensus among the different measures. The data frame has 6 columns: \cr
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
#' # compute_centrality() needs an object of S3 class "carmon" as an input. The
#' # following line quickly obtains one.
#' c_obj <- carmon(multi_omics_small, net_method = "correlation",
#'                 cor_quant = 0.25, analysis = FALSE, verbose = FALSE)
#' # Then to perform the consensus centrality analysis:
#' c_obj <- compute_centrality(c_obj)
#'
compute_centrality <- function(carmon_obj, measures = "dbce", max_candidates = NULL, quant = NULL, scaled = TRUE, verbose = FALSE) {
  if (verbose) {
    cat("**************Beginning network analysis**************")
    cat("\n")
    utils::flush.console()
  }

  if (is.null(max_candidates) & is.null(quant)) {
    n_candidates <- min(20, round(0.05 * sum(unlist(lapply(carmon_obj$layers, ncol)))))
  } else if (!is.null(max_candidates) & is.null(quant)) {
    n_candidates <- max_candidates
  } else if (is.null(max_candidates) & !is.null(quant)) {
    n_candidates <- round(quant * sum(unlist(lapply(carmon_obj$layers, ncol))))
  } else {
    n_candidates <- min(max_candidates, round(quant * sum(unlist(lapply(carmon_obj$layers, ncol)))))
  }

  if (nchar(measures) > 0) {
    no_node_found <- rep(0, nchar(measures))
    for (m in measures) {
      if (!grepl(m, "dbce")) {
        stop(paste0("Unknown measure inserted: ", m, collapse = ""))
      }
    }

    if (verbose > 0 & verbose < 2) {
      mes <- "Computing centrality measures...."
      cat(mes, "\r")
      cat("\n")
    }

    all_measures <- vector(mode = "list", length = nchar(measures))

    for (i in 1:nchar(measures)) {
      all_measures[[i]] <- rep(0, n_candidates)
      m <- substr(measures, i, i)
      if (m == "d") {
        if (verbose > 1) {
          mes <- paste(c("Computing centrality measure ", i, " of ", nchar(measures), " (Degree)"))
          cat(mes, "\r")
          cat("\n")
        }
        names(all_measures)[i] <- "Degree"
        # maybe for degrees (or everything?) I should use the absolute number,
        # and row normalize the resulting table, it is hard to compare
        # conditions and fused otherwise
        degrees <- sort(igraph::degree(carmon_obj$network, normalized = scaled),
          decreasing = TRUE
        )
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
      } else if (m == "b") {
        if (verbose > 1) {
          mes <- paste(c("Computing centrality measure ", i, " of ", nchar(measures), " (Betweenness)"))
          cat(mes, "\r")
          cat("\n")
        }
        names(all_measures)[i] <- "Betweenness"

        betweenness <- sort(igraph::betweenness(carmon_obj$network,
          normalized = scaled,
          weights = NA
        ), decreasing = TRUE)
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
      } else if (m == "c") {
        if (verbose > 1) {
          mes <- paste(c("Computing centrality measure ", i, " of ", nchar(measures), " (Closeness)"))
          cat(mes, "\r")
          cat("\n")
        }
        names(all_measures)[i] <- "Closeness"

        closeness <- sort(igraph::closeness(carmon_obj$network,
          normalized = scaled,
          weights = NA
        ), decreasing = TRUE)

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
      } else {
        if (verbose > 1) {
          mes <- paste(c("Computing centrality measure ", i, " of ", nchar(measures), " (Eigenvector)"))
          cat(mes, "\r")
          cat("\n")
        }
        names(all_measures)[i] <- "Eigenvector Centrality"

        eigenvector <- sort(
          igraph::eigen_centrality(carmon_obj$network,
            weights = NA
          )$vector,
          decreasing = TRUE
        )

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
      }
      if(length(all_measures[[i]]) == 0) {
        no_node_found[i] <- 1
        if(verbose) {
          cat("No node was found to be central for ", names(all_measures)[i], ".\n", sep = "")
        }
      }
    }

    no_node_found <- which(no_node_found == 1)
    all_measures[no_node_found] <- NULL

    if (verbose > 0 & verbose < 2) {
      mes <- "Computing centrality measures....done"
      cat(mes, "\r")
      cat("\n")
    }
    if (length(no_node_found) == nchar(measures)) {
      warning("No node was found to be central. No centrality report will be generated.")
      if (verbose) {
        cat("**************Network analysis complete***************")
        cat("\n")
        utils::flush.console()
      }
      return(carmon_obj)
    }
  } else {
    stop("Give at least a centrality measure for network analysis")
  }

  if (verbose) {
    mes <- "Compiling the report...."
    cat(mes, "\r")
  }

  all_candidates <- unique(unlist(lapply(all_measures, names)))
  report <- matrix(0, nrow = length(all_candidates), ncol = 8)
  for (i in 1:length(all_candidates)) {
    candidate <- all_candidates[i]
    report[i, 1] <- all_candidates[i]

    count <- 0
    string <- ""
    if (candidate %in% degrees_names) {
      string <- paste(string, "d", sep = "")
      count <- count + 1
      report[i, 2] <- paste(c(round(degrees[candidate], digits = 3), "*"), collapse = "")
    } else {
      report[i, 2] <- as.character(round(degrees[candidate], digits = 3))
    }
    if (candidate %in% betweenness_names) {
      string <- paste(string, "b", sep = "")
      count <- count + 1
      report[i, 3] <- paste(c(round(betweenness[candidate], digits = 3), "*"), collapse = "")
    } else {
      report[i, 3] <- as.character(round(betweenness[candidate], digits = 3))
    }
    if (candidate %in% closeness_names) {
      string <- paste(string, "c", sep = "")
      count <- count + 1
      report[i, 4] <- paste(c(round(closeness[candidate], digits = 3), "*"), collapse = "")
    } else {
      report[i, 4] <- as.character(round(closeness[candidate], digits = 3))
    }
    if (candidate %in% eigenvector_names) {
      string <- paste(string, "e", sep = "")
      count <- count + 1
      report[i, 5] <- paste(c(round(eigenvector[candidate], digits = 3), "*"), collapse = "")
    } else {
      report[i, 5] <- as.character(round(eigenvector[candidate], digits = 3))
    }
    report[i, 6] <- string
    report[i, 7] <- count
    report[i, 8] <- sum(as.numeric(strsplit(report[i, 2:5], split = "\\*")))
  }

  report <- report[order(as.numeric(report[, 7]), as.numeric(report[, 8]),
    decreasing = TRUE
  ), ]

  if (is.null(dim(report))) { # Only one gene
    report <- report[-c(7, 8)]
    names(report) <- c("candidate", "degree", "betweenness", "closeness", "eigenvector", "central for")
  } else {
    report <- report[, -c(7, 8)]

    colnames(report) <- c("candidate", "degree", "betweenness", "closeness", "eigenvector", "central for")
    rownames(report) <- report[, 1]

    report <- as.data.frame(report)
  }

  if (verbose) {
    mes <- "Compiling the report....done"
    cat(mes, "\r")
    cat("\n")
    cat("**************Network analysis complete***************")
    cat("\n")
    utils::flush.console()
  }

  carmon_obj$report <- report
  carmon_obj$measures_list <- all_measures
  return(carmon_obj)
}

#' Display the table of results of the consensus centrality analysis of the `carmon` network
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
#' c_obj <- carmon(multi_omics_small, net_method = "correlation",
#'                 cor_quant = 0.25, analysis = TRUE, plot = FALSE,
#'                 # analysis is already TRUE by default
#'                 verbose = FALSE)
#' # To display the table of the results of the centrality analysis:
#' centrality_report(c_obj)
#'
centrality_report <- function(carmon_obj) {
  if(is.null(carmon_obj$report)) {
    warning("The centrality report was not compiled yet. \nCompiling it now, but it will not be saved in the carmon object.")
    carmon_obj <- compute_centrality(carmon_obj)
  }
  return(carmon_obj$report)
}

# Use stability of central measures as well. Perform subsampling to check how
# often that node is associated to a certain (for example) degree. Two ways:
# measure the degree in every subsample and draw a statistical distribution of
# the various subsamples to indicate how stable that is. Otherwise, check in how
# many subsamples that node is in the top 5%, or 10 or whatever, nodes for that
# measure, and use that as a measure of consensus, somehow.
