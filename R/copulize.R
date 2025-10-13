#' Use copulae to transfer every omics layer to the normal realm
#'
#' @param layers The omics layers to analyze. Preferably provided as a named R
#'   list of \strong{non-normalized} omics data sets. If possible the names of the list
#'   should correspond to the respective omics type. If not possible, for
#'   example because of two layers from the same technology, please provide the
#'   omics types with the parameter `omics`.  To see a list of available
#'   omics types run the function `which_omics()`.
#'   Each data set should be source-matched (same amount of matched samples or
#'   individuals across each data set). Placing the samples (or individuals)
#'   should also be consistent: either along the rows for \emph{all} the data
#'   sets, or along the columns for \emph{all} the data sets, nothing in
#'   between. All the samples (or individuals) should also have consistent
#'   naming across the data sets.\cr
#'   `layers` can also be a single unified data set, but then it is necessary to
#'   specify the argument `p`.
#' @param p Optional, to be specified only in case layers is a single data set.
#'   A vector with with the number of variables for each omic layer of the
#'   data set (e.g. the number of transcripts, metabolites etc.), in the same
#'   order the layers have in the data set. If given a single number, carmon
#'   assumes that the total of data sets is two, and that the number given is
#'   the dimension of the first one.
#' @param omics Highly recommended. A vector of as many elements as the number
#'   of layers, naming what omics each layer contains, in the same order as
#'   provided in the input `layers`
#'   (\emph{e.g.} `omics = c("RNA-seq", "proteomics", "metabolomics")`). To
#'   see a list of terms and omics technologies for which carmon is specifically
#'   tailored, use the function `which_omics()`.
#' @param marginals Optional, to be specified when the user prefers to use
#'   different marginal distributions than the default distribution carmon
#'   tailored for each omics layers. A vector of as many elements as the number
#'   of layers, specifying which marginal distribution should be used for each
#'   omics layer, in the same order as provided in the input `layers`. To see a
#'   list of available marginal distributions, use the function
#'   `which_marginals()`. For a more custom setting, place a `0` in the vector
#'   in the position corresponding to the omics layers for which the default
#'   distribution is desired. Otherwise, specify the desired marginal
#'   distribution.
#' @param noninv_method A placeholder for future functionalities of carmon, do
#'   not use.
#' @param copula A placeholder for future functionalities of carmon, do not use.
#' @param verbose The level of verbosity of the copulization process. `0`
#'   suppresses the information output, while `1` and `2` give progressively
#'   increasing amounts of information about the inner computations happening
#'   inside `copulize()`.
#'
#' @returns `copulize()` returns an object of `S3` class `carmon_cop`, having
#'   the following elements:
#' * `layers` is an R list, each element being a data set of the corresponding
#'   layer, already copulized and transferred to the normal realm.
#' * `omics` is a vector containing the omics type assigned to each omics layer.
#' * `marginals` is a vector containing the marginal distributions used to
#'   transfer each omics layer to the normal realm.
#' * `copulize_call` is the matched call.
#' @export
#'
#' @examples
#' # To apply the copula-based transition to the normal realm, it is sufficient
#' # to provide the input data as a named R list, with each element being the
#' # data set of an omics layer.
#' copulized <- copulize(multi_omics_micro, verbose = FALSE)
#'
copulize <- function(layers, p = NULL, omics = NULL, marginals = NULL,
                     noninv_method = NULL, copula = NULL, verbose = FALSE) {
  # Make layers either a data frame itself or a list of datasets. Make sure to
  # specify in help that if layers is a dataframe, samples MUST be on the rows.
  # Copulize will work a bit as a front door to carmon. Think with the logic of
  # having a GUI at a certain point, with the user able to separately add data
  # sets.
  #
  #
  # Add here a call to estimate the marginals parameters. If layers are a list,
  # marginals should be a vector of strings indicating the marginals to assume for the dataset. Add a list of supported marginals and "other" or "no assumption" in the help.
  #
  # Make a parameter for the omic technology given as input. Make it overwriteable by marginals when chosen as an option. Add a list of supported omics and "other" in the help.
  #
  # Regarding the GUI:
  # - Make it obligatory to add two data sets at least.
  # - Make sure to ask whether the row names are actually the sample\replicate names
  call <- match.call()
  if (verbose) {
    mes <- "Checking sample-matching, formatting the data, and other formalities...."
    cat(mes, "\r")
    flush.console()
  }
  if (inherits(layers, "list")) {
    # Check that layers match in sample and sample names
    layers <- check_layers_dims(layers)
    # Check consistency of layer number across arguments
    check_layers_num(layers = layers, omics = omics, marginals = marginals)
  } else if (!is.null(dim(layers))){ # layers is a matrix, table or data.frame
    if (is.null(p)){
      stop("When the multi-omics data set is given as a single matrix, please provide the number of features in each layer through the argument p.")
    } else {
      if (sum(p) != ncol(layers)) { # Assume that dimensions of the last layer are not specified
        p_last <- ncol(layers) - sum(p)
        p <- c(p, p_last)
      }

      # Check consistency of layer number across arguments
      check_layers_num(p = p, omics = omics, marginals = marginals)

      # Make the input multi-omics data set a list of layers
      layers <- split_layers(layers, p = p, omics = omics)
    }
  }

  # Define default parameters
  if (is.null(noninv_method)) {
    noninv_method <- "median"
  }
  if (is.null(copula)) {
    copula <- "gaussian"
  }
  if (is.null(marginals)) {
    marginals <- rep(0, length(layers))
  }
  if (is.null(omics)) {
    if(check_omics(names(layers), marginals)) {
      omics <- names(layers)
    } else {
      omics <- rep("other", length(layers))
    }
  }
  # Check and in case assign layer and variable names
  layers <- check_layers_names(layers, omics)
  layers <- check_gen_colnames(layers)

  # Manage default NULLs, NAs etc. for marginals
  marginals <- omics2marginals(omics, marginals)

  if (verbose) {
    mes <- "Checking sample-matching, formatting the data, and other formalities....done"
    cat(mes, "\r")
    cat("\n\n")
    cat("****************Beginning copulization****************")
    cat("\n")
    flush.console()
  }

  for (l in 1:length(layers)) {
    if (verbose) {
      mes <- paste(c("Copulizing layer ", l, " of ", length(layers), " (", omics[l], ")"), sep ="")
      cat(mes, "\r")
      cat("\n")
    }
    if (tolower(marginals[l]) %in% c("e", "empirical")) {
      if (verbose > 1) {
        mes <- "Using empirical marginals...."
        cat(mes, "\r")
      }
      marginals[l] <- "empirical"
      layers[[l]] <- e_m_copulizer(layers[[l]],
        noninv_method = noninv_method,
        copula = copula
      )
      if (verbose > 1) {
        mes <- "Using empirical marginals....done"
        cat(mes, "\r")
        cat("\n")
      }
    } else if (tolower(marginals[l]) %in% c("nb", "negative binomial")) {
      if (verbose > 1) {
        mes <- "Using negative binomial marginals...."
        cat(mes, "\r")
      }
      marginals[l] <- "negative binomial"
      # Add check for counts, NB distribution expects counts
      layers[[l]] <- nb_m_copulizer(layers[[l]],
        design_matrix = NULL,
        design_formula = NULL,
        noninv_method = noninv_method, copula = copula
      )
      if (verbose > 1) {
        mes <- "Using negative binomial marginals....done"
        cat(mes, "\r")
        cat("\n")
      }
#     } else if (tolower(marginals[l]) %in% c("bb", "beta binomial")) {
#       if (verbose > 1) {
#         mes <- "Using beta binomial marginals...."
#         cat(mes, "\r")
#       }
#       marginals[l] <- "beta binomial"
#       # The methylation layer expects a list of two matrices, the methylated
#       # counts and the site coverage
#       layers[[l]] <- bb_m_copulizer(layers[[l]],
#         design_matrix = NULL,
#         design_formula = NULL,
#         noninv_method = noninv_method, copula = copula
#       )
#       if (verbose > 1) {
#         mes <- "Using beta binomial marginals....done"
#         cat(mes, "\r")
#         cat("\n")
#       }
    } else if (tolower(marginals[l]) %in% c("ln", "lognormal")) {
      if (verbose > 1) {
        mes <- "Using log-normal marginals...."
        cat(mes, "\r")
      }
      marginals[l] <- "lognormal"
      layers[[l]] <- ln_m_copulizer(layers[[l]],
        design_matrix = NULL,
        design_formula = NULL, copula = copula
      )
      if (verbose > 1) {
        mes <- "Using log-normal marginals....done"
        cat(mes, "\r")
        cat("\n")
      }
    } else if (tolower(marginals[l]) %in% c("n", "normal")) {
      if (verbose > 1) {
        mes <- "Using normal marginals...."
        cat(mes, "\r")
      }
      marginals[l] <- "normal"
      layers[[l]] <- n_m_copulizer(layers[[l]],
        design_matrix = NULL,
        design_formula = NULL, copula = copula
      )
      if (verbose > 1) {
        mes <- "Using normal marginals....done"
        cat(mes, "\r")
        cat("\n")
      }
    }
  }

  copulized <- list()
  copulized$layers <- layers
  copulized$omics <- omics
  copulized$marginals <- marginals
  copulized$copulize_call <- call
  class(copulized) <- "carmon_cop"

  if (verbose) {
    cat("****************Copulization complete*****************")
    cat("\n\n")
  }
  return(copulized)
}

#' Obtain the default marginal distribution for each omics layer
#'
#' @inherit copulize
#'
#' @returns It returns a vector with the default or chosen marginal distrbutions
#'   for each omics layer.
#'
#' @noRd
omics2marginals <- function(omics, marginals) {
  for (l in 1:length(omics)) {
    if (marginals[l] == 0) {
      if (tolower(omics[l]) %in% c("rna-seq", "rnaseq", "gene counts", "transcriptomics", "proteomics", "protein fragments", "protein counts")) {
        marginals[l] <- "nb"
#       } else if (tolower(omics[l]) %in% c("bs-seq", "bsseq", "wgbs", "methylomics")) {
#         marginals[l] <- "bb"
      } else if (tolower(omics[l]) %in% c("metabolomics", "lc-ms", "gc-ms", "ms")) {
        marginals[l] <- "ln"
      } else {
        marginals[l] <- "e"
      }
    }
  }
  return(marginals)
}

#' Use the Negative-Binomial marginal distribution to convert a layer to normal
#'
#' @param layer The omics layer to transfer to the normal realm. Should be a
#'   data set with samples (or individuals) on the rows and omics features on
#'   the columns.
#' @param design_matrix A placeholder for future functionalities of carmon, do
#'   not use.
#' @param design_formula A placeholder for future functionalities of carmon, do
#'   not use.
#' @param noninv_method A placeholder for future functionalities of carmon, do
#'   not use.
#' @param copula A placeholder for future functionalities of carmon, do not use.
#'
#' @returns It returns the layer given as an input, but transferred to the
#'   normal realm through the negative binomial marginal distribution and the
#'   specified copula.
#'
#' @noRd
nb_m_copulizer <- function(layer, design_matrix = NULL,
                           design_formula = NULL, noninv_method = NULL,
                           copula = NULL) {
  if (is.null(design_matrix)) {
    dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(t(round(layer)),
                                          colData = data.frame(row.names(layer)),
                                          design = ~1))
  } else {
    # How to change to keep into account covariates/confounders?
    dds <- DESeq2::DESeqDataSetFromMatrix(t(round(layer)),
                                          colData = design_matrix,
                                          design = design_formula
    )
  }
  dds <- DESeq2::estimateSizeFactors(dds, quiet = TRUE)
  dds <- DESeq2::estimateDispersions(dds, quiet = TRUE)
  means <- t(dds@assays@data$mu)
  dispersions <- DESeq2::dispersions(dds)
  if (noninv_method == "median") {
    median_probabilities <- (stats::pnbinom(t(layer - 1), mu = t(means), size = 1/dispersions) +
                               stats::pnbinom(t(layer), mu = t(means), size = 1/dispersions)) / 2
    if (copula == "gaussian") {
      layer <- t(stats::qnorm(median_probabilities))
    }
  }
  layer <- hackInf(layer)
  return(layer)
}

#' Use the empirical marginal distribution to convert a layer to normal
#'
#' @inherit nb_m_copulizer
#'
#' @returns It returns the layer given as an input, but transferred to the
#'   normal realm through the empirical marginal distribution and the
#'   specified copula.
#'
#' @noRd
e_m_copulizer <- function(layer, noninv_method = NULL, copula = NULL) {
  for (j in 1:ncol(layer)) {
    empirical_cdf <- stats::ecdf(layer[, j])
    if (noninv_method == "median") {
      median_probabilities <- (empirical_cdf(layer[, j] -
        1) + empirical_cdf(layer[, j])) / 2
      if (copula == "gaussian") {
        layer[, j] <- stats::qnorm(median_probabilities)
      }
    }
  }
  layer <- hackInf(layer)
  return(layer)
}

#' Use the log-Normal marginal distribution to convert a layer to normal
#'
#' @inherit nb_m_copulizer
#'
#' @returns It returns the layer given as an input, but transferred to the
#'   normal realm through the log-normall marginal distribution and the
#'   specified copula.
#'
#' @noRd
ln_m_copulizer <- function(layer, design_matrix = NULL,
                           design_formula = NULL, copula = NULL) {
  layer <- log(layer)

  # Implement DESeq2-like normalization and variance shrinkage estimation, for
  # small number of replicates give biased MLE estimates of the variance

  infinites <- which(is.infinite(layer), arr.ind = T)
  layer[infinites] <- NA

  if (is.null(design_matrix)) {
    means <- apply(layer, 2, mean, na.rm = TRUE)
    sds <- apply(layer, 2, stats::sd, na.rm = TRUE)
  } else {
    # Needs to change to keep into account covariates/confounders
    means <- apply(layer, 2, mean, na.rm = TRUE)
    sds <- apply(layer, 2, stats::sd, na.rm = TRUE)
  }

  layer[infinites] <- -Inf

  probabilities <- stats::pnorm(t(layer), mean = means, sd = sds)

  if (copula == "gaussian") {
    layer <- t(stats::qnorm(probabilities))
  }

  layer <- hackInf(layer)
  return(layer)
}

#' Use the normal marginal distribution to convert a layer to normal
#'
#' @inherit nb_m_copulizer
#'
#' @returns It returns the layer given as an input, but transferred to the
#'   normal realm through the normal marginal distribution and the
#'   specified copula.
#'
#' @noRd
n_m_copulizer <- function(layer, design_matrix = NULL,
                          design_formula = NULL, copula = NULL) {
  # Implement DESeq2-like normalization and variance shrinkage estimation, for
  # small number of replicates give biased MLE estimates of the variance

  infinites <- which(is.infinite(layer), arr.ind = T)
  layer[infinites] <- NA

  if (is.null(design_matrix)) {
    means <- apply(layer, 2, mean, na.rm = TRUE)
    sds <- apply(layer, 2, stats::sd, na.rm = TRUE)
  } else {
    # Needs to change to keep into account covariates/confounders
    means <- apply(layer, 2, mean, na.rm = TRUE)
    sds <- apply(layer, 2, stats::sd, na.rm = TRUE)
  }

  layer[infinites] <- -Inf

  probabilities <- stats::pnorm(t(layer), mean = means, sd = sds)

  if (copula == "gaussian") {
    layer <- t(stats::qnorm(probabilities))
  }

  layer <- hackInf(layer)
  return(layer)
}

#' Replace infinities with values two standard deviations away from max and min values
#'
#' @param layer The omics layer whose infinities need to be converted to regular
#'   extreme values. Should be a data set with samples (or individuals) on the
#'   rows and omics features on the columns.
#'
#' @returns It returns the layer given as an input, but the +Inf values have
#'   been substituted by values two standard deviations higher than the maximum
#'   value of the feature they are measured in. The opposite holds true for the
#'   -Inf values, substitued by the minimum value in the column, minus two
#'   standard deviations.
#'
#' @noRd
hackInf <- function(layer) {
  # Replaces +Inf with col max + 2 and -Inf with col min -2
  # (2 standard deviations away from max or min)
  if (any(is.infinite(as.matrix(layer)))) {
    layer <- apply(layer, 2, function(x) {
      if (any(is.infinite(x))) {
        ind_pos_inf <- which(is.infinite(x) & x > 0)
        ind_neg_inf <- which(is.infinite(x) & x < 0)
        x[which(is.infinite(x))] <- NA
        x[ind_pos_inf] <- max(x, na.rm = TRUE) + 2
        x[ind_neg_inf] <- min(x, na.rm = TRUE) - 2
      }
      x
    })
  }
  as.matrix(layer)
}
