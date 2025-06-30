copulize <- function(data, p, omics = NULL, marginals = NULL,
                     noninv_method = NULL, copula = NULL) {
  # Make data either a data frame itself or a list of datasets. Copulize will
  # work a bit as a front door to rmonet. Think with the logic of having a GUI
  # at a certain point, with the user able to separately add data sets.
  #
  # Add here a call to estimate the marginals parameters. If data are a list,
  # marginals should be a vector of strings indicating the marginals to assume for the dataset. Add a list of supported marginals and "other" or "no assumption" in the help.
  #
  # Make a parameter for the omic technology given as input. Make it overwriteable by marginals when chosen as an option. Add a list of supported omics and "other" in the help.
  #
  # Regarding the GUI:
  # - Make it obligatory to add two data sets at least.
  # - Make sure to ask whether the row names are actually the sample\replicate names

  # Define default parameters
  if (is.null(noninv_method)) {
    noninv_method <- "median"
  }
  if (is.null(copula)) {
    copula <- "gaussian"
  }
  if (is.null(marginals)) {
    marginals <- rep(0, length(data))
  }
  if (is.null(omics)) {
    omics <- rep("other", length(data))
  }

  # Manage default NULLs, NAs etc. for marginals
  marginals <- omics2marginals(omics, marginals)

  for (d in 1:length(data)) {
    if (tolower(marginals[d]) %in% c("e", "empirical")) {
      data[[d]] <- e_m_copulizer(data[[d]],
        noninv_method = noninv_method,
        copula = copula
      )
    } else if (tolower(marginals[d]) %in% c("nb", "negative binomial")) {
      # Add check for counts, NB distribution expects counts
      data[[d]] <- nb_m_copulizer(data[[d]],
        design_matrix = NULL,
        design_formula = NULL,
        noninv_method = noninv_method, copula = copula
      )
    } else if (tolower(marginals[d]) %in% c("bb", "beta binomial")) {
      # The methylation layer expects a list of two matrices, the methylated
      # counts and the site coverage
      data[[d]] <- bb_m_copulizer(data[[d]],
        design_matrix = NULL,
        design_formula = NULL,
        noninv_method = noninv_method, copula = copula
      )
    } else if (tolower(marginals[d]) %in% c("ln", "lognormal")) {
      data[[d]] <- ln_m_copulizer(data[[d]],
        design_matrix = NULL,
        design_formula = NULL, copula = copula
      )
    } else if (tolower(marginals[d]) %in% c("n", "normal")) {
      data[[d]] <- n_m_copulizer(data[[d]],
        design_matrix = NULL,
        design_formula = NULL, copula = copula
      )
    }
  }

  return(data)
}

omics2marginals <- function(omics, marginals) {
  for (d in 1:length(omics)) {
    if (marginals[d] == 0) {
      if (tolower(omics[d]) %in% c("rna-seq", "rnaseq", "proteomics", "protein fragments")) {
        marginals[d] <- "nb"
      } else if (tolower(omics[d]) %in% c("bs-seq", "bsseq", "wgbs", "methylomics")) {
        marginals[d] <- "bb"
      } else if (tolower(omics[d]) %in% c("metabolomics", "lc-ms", "gc-ms", "ms")) {
        marginals[d] <- "ln"
      } else {
        marginals[d] <- "e"
      }
    }
  }
  return(marginals)
}

e_m_copulizer <- function(layer, noninv_method = NULL, copula = NULL) {
  for (j in 1:dim(layer)[2]) {
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

nb_m_copulizer <- function(layer, design_matrix = NULL,
                           design_formula = NULL, noninv_method = NULL,
                           copula = NULL) {
  if (is.null(design_matrix)) {
    dds <- DESeq2::DESeqDataSetFromMatrix(t(round(layer)),
      colData = data.frame(row.names(layer)),
      design = ~1
    )
  } else {
    # How to change to keep into account covariates/confounders?
    dds <- DESeq2::DESeqDataSetFromMatrix(t(round(layer)),
      colData = design_matrix,
      design = design_formula
    )
  }
  dds <- DESeq2::estimateSizeFactors(dds)
  dds <- DESeq2::estimateDispersions(dds)
  means <- t(dds@assays@data$mu)
  dispersions <- DESeq2::dispersions(dds)

  if (noninv_method == "median") {
    median_probabilities <- (stats::pnbinom(layer - 1, mu = means, size = dispersions) +
      stats::pnbinom(layer, mu = means, size = dispersions)) / 2
    if (copula == "gaussian") {
      layer <- stats::qnorm(median_probabilities)
    }
  }
  layer <- hackInf(layer)
  return(layer)
}

bb_m_copulizer <- function(layer, design_matrix = NULL,
                           design_formula = NULL, noninv_method = NULL,
                           copula = NULL) {
  # Implement also smoothed average
  # Implement also complex designs

  if (is.null(design_matrix)) {
    dss_parameters <- dss_parameters(t(layer[[1]]), t(layer[[2]]))
  } else {
    # lacking complex designs for the moment
    dss_parameters <- dss_parameters(t(layer[[1]]), t(layer[[2]]))
  }

  probs <- t(dss_parameters[["prob"]])
  dispersions <- dss_parameters[["phi"]]

  if (noninv_method == "median") {
    median_probabilities <- (pbetabinom(layer[[1]] - 1,
      n = layer[[2]],
      prob = probs, phi = dispersions
    ) +
      pbetabinom(layer[[1]],
        n = layer[[2]],
        prob = probs, phi = dispersions
      )) / 2
    if (copula == "gaussian") {
      layer <- stats::qnorm(median_probabilities)
    }
  }
  layer <- hackInf(layer)
  return(layer)
}

ln_m_copulizer <- function(layer, design_matrix = NULL,
                           design_formula = NULL, copula = NULL) {
  layer <- log(layer)

  # Implement DESeq2-like normalization and variance shrinkage estimation, for
  # small number of replicates give biased MLE estimates of the variance

  infinites <- which(is.infinite(layer), arr.ind = T)
  layer[infinites] <- NA

  if (is.null(design_matrix)) {
    means <- apply(layer, 2, mean, na.rm = TRUE)
    variances <- apply(layer, 2, stats::var, na.rm = TRUE)
  } else {
    # Needs to change to keep into account covariates/confounders
    means <- apply(layer, 2, mean, na.rm = TRUE)
    variances <- apply(layer, 2, stats::var, na.rm = TRUE)
  }

  layer[infinites] <- -Inf

  probabilities <- stats::pnorm(layer, mu = means, size = variances)

  if (copula == "gaussian") {
    layer <- stats::qnorm(probabilities)
  }

  layer <- hackInf(layer)
  return(layer)
}

n_m_copulizer <- function(layer, design_matrix = NULL,
                          design_formula = NULL, copula = NULL) {
  # Implement DESeq2-like normalization and variance shrinkage estimation, for
  # small number of replicates give biased MLE estimates of the variance

  infinites <- which(is.infinite(layer), arr.ind = T)
  layer[infinites] <- NA

  if (is.null(design_matrix)) {
    means <- apply(layer, 2, mean, na.rm = TRUE)
    variances <- apply(layer, 2, stats::var, na.rm = TRUE)
  } else {
    # Needs to change to keep into account covariates/confounders
    means <- apply(layer, 2, mean, na.rm = TRUE)
    variances <- apply(layer, 2, stats::var, na.rm = TRUE)
  }

  layer[infinites] <- -Inf

  probabilities <- stats::pnorm(layer, mu = means, size = variances)

  if (copula == "gaussian") {
    layer <- stats::qnorm(probabilities)
  }

  layer <- hackInf(layer)
  return(layer)
}

hackInf <- function(layer) {
  # Replaces +Inf with col max + 2 and -Inf with col min -2
  # (2 standard deviations away from max or min)
  if (any(is.infinite(as.matrix(layer)))) {
    data <- apply(layer, 2, function(x) {
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
