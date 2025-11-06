#' Use copulas to transfer all omics layers to the normal realm
#'
#' @param layers The omics layers to analyze. They can be provided in three
#' possible formats, and they should always contain source-matched,
#' \strong{non-normalized} data, to use \pkg{carmon} to its full potential.
#' \itemize{
#'   \item{} A named R list of omics data sets (\emph{recommended}). If possible
#'   the names of the list should correspond to the respective omics type. If
#'   not possible, for example because of two layers from the same technology,
#'   please provide the omics types with the parameter `omics`. To see a list of
#'   available omics types use the function `which_omics()`.\cr
#'   Each data set should be source-matched (same amount of matched samples or
#'   individuals across each data set). Placing of the samples (or individuals)
#'   should also be consistent: either along the rows for \emph{all} the data
#'   sets, or along the columns for \emph{all} the data sets, nothing in
#'   between. All the samples (or individuals) should also have consistent
#'   naming across the data sets.
#'   \item{} An object of `S4` class `MultiAssayExperiment` (more
#'   \link[MultiAssayExperiment:MultiAssayExperiment]{here}). In that case, we
#'   recommend using the `omics` parameter to specify which are the omics layers
#'   contained in the object, in the same order as presented by
#'   `MultiAssayExperiment::experiments(layers)`. This allows to use
#'   \pkg{carmon} to its full potential. To see a list of terms and omics
#'   technologies for which \pkg{carmon} is specifically tailored, use the
#'   function `which_omics()`. Please remember that `carmon()` expects data to
#'   be non-normalized, meaning that for RNA-seq data, for example, it will
#'   expect data to be in the form of counts.
#'   \item{} A single unified data set, but then it is
#'   necessary to specify the argument `p`. We also recommend using the
#'   `omics` parameter to specify which are the omics layers contained in the
#'   data set.}
#' @param p Optional, to be specified only in case layers is a single data set.
#'   A vector with with the number of variables for each omic layer of the
#'   data set (e.g. the number of transcripts, metabolites etc.), in the same
#'   order the layers have in the data set. If given a single number, carmon
#'   assumes that the total of data sets is two, and that the number given is
#'   the dimension of the first one.
#' @param omics Highly recommended. A vector of as many elements as the number
#'   of layers, naming what omics each layer contains, in the same order as
#'   provided in the input `layers`
#'   (\emph{e.g.} `omics = c('RNA-seq', 'proteomics', 'metabolomics')`). To
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
#' data(multi_omics_micro)
#' copulized <- copulize(multi_omics_micro, verbose = FALSE)
#'
copulize <- function(
    layers, p = NULL, omics = NULL, marginals = NULL, noninv_method = NULL,
    copula = NULL, verbose = FALSE) {
    checked <- check_formalities(layers, p, omics, marginals, noninv_method,
                                copula, verbose)
    layers <- checked$layers
    marginals <- checked$marginals
    omics <- checked$omics
    for (l in seq_len(length(layers))) {
        if (verbose) {
    message("Copulizing layer ", l, " of ", length(layers), " (", omics[l], ")")
        }
        if (tolower(marginals[l]) %in% c("e", "empirical")) {
            marginals[l] <- "empirical"
            copulization_deep_verbose(marginals[l], omics[l], FALSE, verbose)
            layers[[l]] <- e_m_copulizer(layers[[l]],
                noninv_method = checked$noninv_method, copula = checked$copula)
            copulization_deep_verbose(marginals[l], omics[l], TRUE, verbose)
        } else if (tolower(marginals[l]) %in% c("nb", "negative binomial")) {
            marginals[l] <- "negative binomial"
            copulization_deep_verbose(marginals[l], omics[l], FALSE, verbose)
            layers[[l]] <- nb_m_copulizer(layers[[l]], design_matrix = NULL,
                design_formula = NULL, noninv_method = checked$noninv_method,
                copula = checked$copula)
            copulization_deep_verbose(marginals[l], omics[l], TRUE, verbose)
        } else if (tolower(marginals[l]) %in% c("ln", "lognormal")) {
            marginals[l] <- "lognormal"
            copulization_deep_verbose(marginals[l], omics[l], FALSE, verbose)
            layers[[l]] <- ln_m_copulizer(layers[[l]], design_matrix = NULL,
                design_formula = NULL, copula = checked$copula)
            copulization_deep_verbose(marginals[l], omics[l], TRUE, verbose)
        } else if (tolower(marginals[l]) %in% c("n", "normal")) {
            marginals[l] <- "normal"
            copulization_deep_verbose(marginals[l], omics[l], FALSE, verbose)
            layers[[l]] <- n_m_copulizer(layers[[l]], design_matrix = NULL,
                design_formula = NULL, copula = checked$copula)
            copulization_deep_verbose(marginals[l], omics[l], TRUE, verbose)
        } # Add new `else if` statement here if expanding the marginals library
    }
    copulized <- list()
    copulized$layers <- layers
    copulized$omics <- omics
    copulized$marginals <- marginals
    copulized$copulize_call <- match.call()
    class(copulized) <- "carmon_cop"
    if (verbose) {
        message("****************Copulization complete*****************\n")
    }
    return(copulized)
}

#' Check formalities of input to carmon, set defaults
#'
#' @inherit copulize
#'
#' @returns This function returns a list containing the checked and re-formatted
#'     input arguments. Where the user did not provide custom values, the
#'     arguments are set to a default.
#'
#' @noRd
#'
check_formalities <- function(layers, p = NULL, omics, marginals, noninv_method,
                                copula, verbose) {
    if (verbose) {
        message("Checking sample-matching, formatting the data,
    and other formalities....")
        flush.console()
    }
    if (inherits(layers, "MultiAssayExperiment")) {
        layers <- layers_from_MAE(layers)
        check_layers_num(layers = layers, omics = omics, marginals = marginals)
    } else if (inherits(layers, "list")) {
        # Layers match in sample and sample names, and consistent layer number
        layers <- check_layers_dims(layers)
        check_layers_num(layers = layers, omics = omics, marginals = marginals)
    } else if (!is.null(dim(layers))) {
        # layers is a matrix, table or data.frame
        if (is.null(p)) {
            stop("When the multi-omics data set is given as a single matrix,
    please provide the number of features in each layer through the argument p")
        } else {
            if (sum(p) != ncol(layers)) {
                p_last <- ncol(layers) - sum(p)
                p <- c(p, p_last)
            }
            # Consistency of layer number across arguments, split layers in list
            check_layers_num(p = p, omics = omics, marginals = marginals)
            layers <- split_layers(layers, p = p, omics = omics)
        }
    }
    if (is.null(noninv_method)) noninv_method <- "median"
    if (is.null(copula))        copula <- "gaussian"
    if (is.null(marginals))     marginals <- rep(0, length(layers))
    if (is.null(omics)) {
        if (check_omics(names(layers), marginals)) {
            omics <- names(layers)
        } else {
            omics <- rep("other", length(layers))
        }
    }
    # Check and in case assign layer and variable names, and manage marginals
    layers <- check_layers_names(layers, omics)
    layers <- check_gen_colnames(layers)
    marginals <- omics2marginals(omics, marginals)
    if (verbose) {
    message("Done!\n\n****************Beginning copulization****************")
        flush.console()
    }
    return(list("layers" = layers, "omics" = omics, "marginals" = marginals,
                "noninv_method" = noninv_method, "copula" = copula))
}

#' Obtain the default marginal distribution for each omics layer
#'
#' @inherit copulize
#'
#' @returns It returns a vector with the default or chosen marginal
#'   distributions for each omics layer.
#'
#' @noRd
#'
omics2marginals <- function(omics, marginals) {
    for (l in seq_len(length(omics))) {
        if (marginals[l] == 0) {
            if (tolower(omics[l]) %in% c("rna-seq", "rnaseq", "rna",
                    "gene counts", "transcriptomics",
                    "mirna-seq", "mirnaseq", "microrna-seq", "micrornaseq",
                    "mirna", "microrna",
                    "proteomics", "protein fragments", "protein counts")) {
                # Negative Binomial marginal distribution
                marginals[l] <- "nb"
            # } else if (tolower(omics[l]) %in% c("bs-seq", "bsseq", "wgbs",
            #             "methylomics")) {
            #    # Beta Binomial marginal distribution
            #    marginals[l] <- 'bb'
            } else if (tolower(omics[l]) %in% c("metabolomics", "lc-ms",
                            "gc-ms", "ms")) {
                # log Normal marginal distribution
                marginals[l] <- "ln"
            } else if (tolower(omics[l]) %in% c("")) {
                # Normal marginal distribution
                # Nothing implemented yet as normally distributed by default
                marginals[l] <- "n"
            } else {
                # Empirical marginal distribution
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
nb_m_copulizer <- function(
    layer, design_matrix = NULL, design_formula = NULL,
    noninv_method = NULL, copula = NULL) {
    if (is.null(design_matrix)) {
        dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(t(round(layer)),
            colData = data.frame(row.names(layer)), design = ~1
        ))
    } else {
        # How to change to keep into account covariates/confounders?
        dds <- DESeq2::DESeqDataSetFromMatrix(t(round(layer)),
            colData = design_matrix, design = design_formula
        )
    }
    dds <- DESeq2::estimateSizeFactors(dds, quiet = TRUE)
    dds <- DESeq2::estimateDispersions(dds, quiet = TRUE)
    means <- t(dds@assays@data$mu)
    dispersions <- DESeq2::dispersions(dds)
    if (noninv_method == "median") {
        median_probabilities <- (stats::pnbinom(t(layer - 1),
            mu = t(means),
            size = 1 / dispersions
        ) + stats::pnbinom(t(layer),
            mu = t(means),
            size = 1 / dispersions
        )) / 2
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
#'
e_m_copulizer <- function(layer, noninv_method = NULL, copula = NULL) {
    for (j in seq_len(ncol(layer))) {
        empirical_cdf <- stats::ecdf(layer[, j])
        if (noninv_method == "median") {
            median_probabilities <- (empirical_cdf(layer[, j] - 1) +
                empirical_cdf(layer[, j])) / 2
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
#'
ln_m_copulizer <- function(
    layer, design_matrix = NULL, design_formula = NULL,
    copula = NULL) {
    layer <- log(layer)

    # Implement DESeq2-like normalization and variance shrinkage estimation,
    # for small number of replicates give biased MLE estimates of the variance

    infinites <- which(is.infinite(layer), arr.ind = TRUE)
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
#'
n_m_copulizer <- function(
    layer, design_matrix = NULL, design_formula = NULL,
    copula = NULL) {
    # Implement DESeq2-like normalization and variance shrinkage estimation,
    # for small number of replicates give biased MLE estimates of the variance

    infinites <- which(is.infinite(layer), arr.ind = TRUE)
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

#' Replace infinities with values two standard deviations away from max and min
#' values
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
#'
hackInf <- function(layer) {
    # Replaces +Inf with col max + 2 and -Inf with col min -2 (2 standard
    # deviations away from max or min)
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
