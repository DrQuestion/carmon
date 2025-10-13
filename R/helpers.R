#' Check that layers are properly sample matched
#'
#' Check that samples/individuals are distributed coherently either along
#' columns or along rows among all layers, and that they share the same names.
#'
#' @param layers A named R list of the source-matched omics layers to analyze.
#'
#' @returns It raises an error if the samples/individuals are not
#'   source-matched, or if they don't share the same names across all layers.
#'
#' @noRd
check_layers_dims <- function(layers) {
  if (length(unique(sapply(layers, nrow))) == 1) {
    ref_names <- rownames(layers[[1]])
    all_same_names <- all(
      sapply(layers, function(l) setequal(rownames(l), ref_names))
    )
    if (!all_same_names) {
      stop("Samples should be along rows and they must share the same names across layers.")
    } else {
      same_order <- all(
        sapply(layers, function(l) identical(rownames(l), ref_names))
      )
      if (same_order) {
        return(layers)
      } else {
        return(lapply(layers, function(l) l[ref_names, , drop = FALSE]))
      }
    }
  } else if (length(unique(sapply(layers, ncol))) == 1) {
    ref_names <- colnames(layers[[1]])
    all_same_names <- all(
      sapply(layers, function(l) setequal(colnames(l), ref_names))
    )
    if (!all_same_names) {
      stop("Samples should be along rows and they must share the same names across layers.")
    } else {
      warning("The samples will be assumed to be distributed along the columns of the layers. Please, next time arrange them along rows.")
      same_order <- all(
        sapply(layers, function(l) identical(colnames(l), ref_names))
      )
      if (same_order) {
        return(lapply(layers, t))
      } else {
        return(lapply(layers, function(l) t(l[, ref_names, drop = FALSE])))
      }
    }
  } else {
    stop("Dimensions do not match among layers. Please check that all layers share the same samples.")
  }
}

#' Check that all input arguments have a consistent number of omics layers.
#'
#' @param layers A named R list of the source-matched omics layers to analyze.
#' @param p A vector with the number of omics features for each layer,
#'   separately.
#' @param omics A vector containing the omics type of each omics layer, in the
#'   same order as the list of layers given as an input to carmon.
#' @param marginals A vector of as many elements as the number
#'   of layers, specifying which marginal distribution should be used for each
#'   omics layer, in the same order as provided in the input `layers`.
#'
#' @return It returns nothing, but it raises an error if the check fails.
#'
#' @noRd
check_layers_num <- function (layers = NULL, p = NULL, omics = NULL, marginals = NULL) {
  # Function checking that number of layers is consistent across arguments
  input <- list(layers = layers, p = p, omics = omics, marginals = marginals)
  non_null <- input[!sapply(input, is.null)]
  if (length(unique(sapply(non_null, length))) != 1) {
    lengths <- lapply(non_null, length)
    print(lengths)
    stop("Number of layers is inconsistent across input parameters. Please check your input.")
  }
}

#' Split a single data set in a list of sub-data sets
#'
#' `split_layers()` splits a single data set in multiple ones, each one having
#' as many omics features as the corresponding entry of the vector `p`, and
#' collects them in a named R list. The names of the list are assigned based on
#' the elements of the vector of omics types given through `omics`.
#'
#' @param layers The input data set with multiple omics layers merged.
#' @param p A vector with the number of omics features of each layer, in the
#'   same order as the layers.
#' @param omics A vector containing the omics type of each omics layer, in the
#'   same order as the the layers in the data set given as an input to carmon.
#'
#' @returns `split_layers()` returns a named R list. Each element of the list
#'   the data set of the corresponding omics layer. The name of each element is
#'   based on the corresponding omics type.
#'
#' @noRd
split_layers <- function(layers, p, omics = NULL) {
  groups <- factor(rep(seq_along(p), times = p))
  if (!is.null(omics)) {
    levels(groups) <- omics
  }
  layers_idxs <- split(seq_len(ncol(layers)), groups)
  layers <- lapply(layers_idxs, function(idx) layers[, idx, drop = FALSE])
  return(layers)
}

#' Checks whether the omics types provided as input are compatible with carmon
#' or with its own omics nomenclature.
#'
#' @param omics A vector containing the omics type of each omics layer, in the
#'   same order as the list of layers given as an input to carmon.
#' @param marginals The vector of user-chosen marginal distributions for each
#'   omics layer.
#'
#' @returns `check_omics()` returns a logical. `TRUE` if all given omics are
#'   compatible with carmon's nomenclature, `FALSE` if any provided omics is
#'   not. Also raises a warning if for the non-compatible omics layer no custom
#'   marginal distribution is specified by the user, telling that the empirical
#'   marginal will be automatically chosen.
#'
#' @noRd
check_omics <- function(omics, marginals) {
  admitted_omics <- c("rna-seq", "rnaseq", "gene counts", "transcriptomics",
                      "proteomics", "protein fragments", "protein counts",
#                       "bs-seq", "bsseq", "wgbs", "methylomics",
                      "metabolomics", "lc-ms", "gc-ms", "ms")
  if(any(!(omics %in% admitted_omics))){
    if(all(marginals == 0)) {
      warning("Some omics provided are not supported by carmon.\n
              Check which_omics() for a list of supported omics technologies.\n
              In this run their marginal distributions will be modelled with\n
              the empirical marginal distribution.")
    }
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#' Name each unnamed element of the layers list with the corresponding omics type
#'
#' @param layers A named R list of the source-matched omics layers to analyze.
#' @param omics  A vector containing the omics type of each omics layer, in the
#'   same order as the list of layers given as an input to carmon.
#'
#' @returns `check_layers_names()` returns the R list of omics layers, with each
#'   element named after the omics type, if a name was not provided by the user
#'   already. In case of multiple omics layers of the same type, the given names
#'   follow the form "omics_A, omics_B, ...".
#'
#' @noRd
check_layers_names <- function(layers, omics) {
  layer_names <- omics
  if (length(unique(layer_names)) < length(layer_names)) {
    # There are multiple layers of the same omics type
    for (i in 1:length(layer_names)) {
      if (duplicated(layer_names)[i]) {
        omic <- layer_names[i]
        same <- layer_names[which(layer_names == omic)]
        same <- paste(same, LETTERS[1:length(same)], sep = "_")
        layer_names[which(layer_names == omic)] <- same
      }
    }
  }
  for (i in 1:length(layers)) {
    if(is.null(names(layers)[i])){
      names(layers)[i] <- layer_names[i]
    } else if (is.na(names(layers)[i])) {
      names(layers)[i] <- layer_names[i]
    }
  }
  return(layers)
}

#' Name the unnamed omics features in each layer after their omics type
#'
#' @param layers A named R list of the source-matched omics layers to analyze.
#'
#' @returns `check_gen_colnames()` returns the R list of omics layers, and if
#'   the features of the data set of an omics layer are not named (\emph{i.e.}
#'   there are no column names), the features will be named after the omics
#'   type, following the format "omics_colnumber".
#'
#' @noRd
check_gen_colnames <- function(layers) {
  for (i in 1:length(layers)) {
    if (is.null(colnames(layers[[i]]))) {
      colnames(layers[[i]]) <- paste(names(layers)[i], seq(1:ncol(layers[[i]])),
                                     sep = "_")
    }
  }
  return(layers)
}

#' Select from the ellipsis (...) only arguments valid for the target function
#'
#' @param fun The target function for whom the validity of arguments in the
#'   ellipsis must be checked.
#' @param dots The ellipsis, to be provided as `list(...)`.
#'
#' @returns `filter_args` returns the filtered list with only the arguments in
#'   the ellipsis that are valid for the target function, and respctive values
#'   provided by the user. THe resulting list can be used in the  form
#'   `do.call(fun, filter_args(fun, list(...)))`.
#'
#' @noRd
filter_args <- function(fun, dots) {
  # Helper function selecting from an ellipsis only the arguments matching those
  # of the interesting function
  valid <- names(formals(fun))
  dots[names(dots) %in% valid]
}

#' Merge a list of omics layers into a single data set
#'
#' @param layers A named R list of source matched omics layers.
#'
#' @returns `merge_layers()` returns a single data set merging all the omics
#'   layers.
#'
#' @noRd
merge_layers <- function(layers) {
  # Make a single data set out of layers list
  if (length(unique(sapply(layers, nrow))) == 1) {
    layers <- do.call(cbind, layers)
  } else {
    layers <- t(do.call(rbind, layers))
  }
  return(layers)
}

#' Print the omics types for which carmon is tailored on
#'
#' @returns `which_omics()` prints a message with the omics types that carmon is
#'   tailored on, together with a list of the synonyms of each type that carmon
#'   can comprehend, and the default tailored marginal distribution of each
#'   omics type.
#'
#' @export
#'
#' @examples
#' # See the compatible omics types, their synonyms, and their distributions
#' which_omics()
which_omics <- function() {
  cat('"rna-seq", also as "rnaseq", "gene counts", "transcriptomics"
   is modeled by default as count data with a negative-binomial marginal.\n
"proteomics", also as "protein fragments", "protein counts"
   is modeled by default as count data with a negative-binomial marginal.\n
"metabolomics", also as "lc-ms", "gc-ms", "ms"
   is modeled by default as positive continuous data with a log-normal marginal.\n
Anything else is modeled with the empirical marginal.')
}

#' Print the marginal distributions currently implemented inside carmon
#'
#' @returns `which_marginals()` prints a message with the marginal distributions
#'   implemented inside carmon in the current version of the package, together
#'   with the nomenclature and the abbreviations that carmon can comprehend when
#'   giving a custom input to the argument `marginals`.
#'
#' @export
#'
#' @examples
#' # See the marginal distributions you can select to model your omics data with carmon
#' which_marginals()
which_marginals <- function() {
  cat('"e" or "empirical" for using the empirical marginal distribution;
"n" or "normal" for using the normal marginal distribution;
"ln" or "lognormal" for using the log-normal marginal distribution;
"nb" or "negative binomial" for using the log-normal marginal distribution.')
}

#' Extract a `carmon` network
#'
#' `get_network()` extracts the reconstructed (and selected) network from a
#' `carmon` object.
#'
#' @encoding UTF-8
#' @param carmon_obj The object of `S3` class `carmon` or of `S3`
#'   class `carmon_rec`.
#' @param labels Optional. Used only when `carmon_obj` is of the `S3` class
#'   `carmon_rec`, a vector of strings containing all the variable names.
#'
#' @return `get_network()` returns the selected network, in the form of an
#' object of class `igraph`. If using an object of class `carmon_rec`, nodes
#' will have no labels, unless provided through `labels`.
#' @export
#'
#' @examples
#' c_obj <- carmon(multi_omics_small, net_method = "correlation",
#'                 cor_quant = 0.25, analysis = FALSE, plot = FALSE,
#'                 verbose = FALSE)
#' network <- get_network(c_obj)
#'
get_network <- function(carmon_obj, labels = NULL) {

  network <- igraph::graph_from_adjacency_matrix(carmon_obj$sel_adj, mode = "max")
  if (!inherits(carmon_obj, "carmon_rec")){
    igraph::V(network)$label <- unlist(lapply(carmon_obj$layers, colnames))
  } else if (!is.null(labels)) {
    igraph::V(network)$label <- labels
  }
  return(network)
}

#' Extract the weights of the edges of a `carmon` network
#'
#' `get_edge_weights()` extracts the weights of the edges from the reconstructed
#' network of a `carmon` object.
#'
#' It extracts by default the weights of the network reconstructed by carmon.
#' If the reconstruction method was "coglasso" or "glasso", it extracts the
#' partial correlation matrix of the network built through them that was
#' selected by the designed selection method. If instead the reconstruction
#' method was "correlation", it extracts the matrix of the correlation values
#' corresponding to the edges of the network. For the "mb" method, the extracted
#' weights are all equal to 1 the estimated edges, 0 otherwise.
#'
#' @param carmon_obj The object of `S3` class `carmon`
#'
#' @return `get_edge_weights()` returns the edge weights matrix of the carmon network.
#'
#' @noRd
#'
get_edge_weights <- function(carmon_obj){
  if(carmon_obj$net_method %in% c("coglasso", "glasso")){
    weights <- stats::cov2cor(carmon_obj$sel_icov)
    weights <- -weights
    diag(weights) <- 1
  } else if (carmon_obj$net_method == "correlation"){
    weights <- carmon_obj$cor
    weights[which(carmon_obj$sel_adj == 0, arr.ind = TRUE)] <- 0
    diag(weights) <- 1
  } else {
    #method is mb
    weights <- matrix(data = 1, nrow = ncol(carmon_obj$sel_adj), ncol = ncol(carmon_obj$sel_adj))
    weights[Matrix::which(carmon_obj$sel_adj == 0, arr.ind = TRUE)] <- 0
    diag(weights) <- 1
  }

  colnames(weights) <- rownames(weights) <- unlist(lapply(carmon_obj$layers, colnames))

  return(weights)
}

#' Assemble an object of class `carmon_cop` and one of class `carmon_rec` into a `carmon` object
#'
#' @param carmon_cop An object of `S3` class `carmon_cop`, resulting from `copulize()`.
#' @param carmon_rec An object of `S3` class `carmon_rec`, resulting from `reconstruct()`.
#'
#' @returns `assemble_carmon_obj()` returns an object of `S3` class `carmon`,
#'   assembled from the results of the `copulize()` and the `reconstruct()`
#'   functions. For the details of the composition of this object see
#'   [carmon()], keeping in mind that the assembled object will resemble the
#'   result of a `carmon()` run setting the argument `analysis = FALSE`.
#'
#' @export
#'
#' @examples
#' carmon_cop_obj <- copulize(multi_omics_micro, verbose = FALSE)
#' carmon_rec_obj <- reconstruct(carmon_cop_obj$layers,
#'                               net_method = "correlation",
#'                               cor_quant = 0.5, verbose = FALSE)
#' carmon_obj <- assemble_carmon_obj(carmon_cop_obj, carmon_rec_obj)
#'
assemble_carmon_obj <- function(carmon_cop, carmon_rec) {
  carmon_obj <- c(carmon_cop, carmon_rec)
  carmon_obj$network <- get_network(carmon_obj)
  class(carmon_obj) <- "carmon"
  return(carmon_obj)
}
