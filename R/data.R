#' Multi-omics dataset of sleep deprivation in mouse
#'
#' @description A dataset containing transcript and metabolite values analysed
#' in Albanese et al. 2024, subset of the multi-omics data set published in
#' Jan, M., Gobet, N., Diessler, S. et al. A multi-omics digital research object
#' for the genetics of sleep regulation. Sci Data 6, 258 (2019).
#' These data set and their generation follow closely what can be found in the
#' R package \pkg{coglasso}. Also see the file in the \pkg{carmon} R package,
#' under 'carmon/data-raw/multi-omics_data.R' to see how this data set is
#' generated.
#'
#' `multi_omics_small` is a smaller version, limited to the transcript Cirbp
#' and the transcripts and metabolites belonging to its neighborhood as
#' described in Albanese et al. 2024
#'
#' `multi_omics_micro` is a minimal version with Cirbp and a selection of its
#'  neighborhood.
#'
#' @format ## `multi_omics`
#' A named R list of two elements, each being a data set with 30 observations:
#' \describe{
#'   \item{rnaseq}{data set of the RNA-seq counts of 162 transcripts in mouse
#'     cortex under sleep deprivation (0--34991)}
#'   \item{metabolomics}{data set of the abundance values of 76 metabolites
#'     (0.02--1112.67)}
#' }
#' @source Jan, M., Gobet, N., Diessler, S. et al. A multi-omics digital
#' research object for the genetics of sleep regulation. Sci Data 6, 258
#' (2019) doi: \href{https://doi.org/10.1038/s41597-019-0171-x
#' }{10.1038/s41597-019-0171-x}
#' @source \href{https://figshare.com/articles/dataset/
#' Input_data_for_systems_genetics_of_sleep_regulation/7797434
#' }{Figshare folder of the original manuscript}
"multi_omics"

#' @rdname multi_omics
#' @format ## `multi_omics_small`
#' A named R list of two elements, each being a data set with 30 observations:
#' \describe{
#'   \item{rnaseq}{data set of the RNA-seq counts of 14 transcripts in mouse
#'     cortex under sleep deprivation (412--16235)}
#'   \item{metabolomics}{data set of the abundance values of 5 metabolites
#'   (0.17--145.33)}
#' }
"multi_omics_small"

#' @rdname multi_omics
#' @format ## `multi_omics_micro`
#' A named R list of two elements, each being a data set with 30 observations:
#' \describe{
#'   \item{rnaseq}{data set of the RNA-seq counts of 2 transcripts in mouse
#'     cortex under sleep deprivation (623--13304)}
#'   \item{metabolomics}{data set of the abundance values of 2 metabolites
#'   (58.80--145.33)}
#' }
"multi_omics_micro"
