test_that("copulize works with a dataset", {
  dataset <- merge_layers(multi_omics_micro)
  expect_error(copulize(dataset))
  expect_length(copulize(dataset, p = 2, omics = c("rnaseq", "metabolomics"))$layers, 2)
})

test_that("non-tailored marginals work (including verbosity)", {
  layers <- multi_omics_micro
  layers$metabolomics <- log(layers$metabolomics)
  expect_output(expect_length(copulize(layers, p = 2, marginals = c("e", "n"), verbose = 2)$layers, 2))
})

test_that("hackInf works", {
  copulized <- copulize(multi_omics_micro)
  max_trp <- max(copulized$layers$metabolomics[,1])
  min_trp <- min(copulized$layers$metabolomics[,1])
  copulized$layers$metabolomics[4,1] <- +Inf
  copulized$layers$metabolomics[3,1] <- -Inf
  copulized$layers$metabolomics <- hackInf(copulized$layers$metabolomics)
  expect_equal(max_trp, max(copulized$layers$metabolomics[,1])-2)
  expect_equal(min_trp, min(copulized$layers$metabolomics[,1])+2)
})
