test_that("reconstruction with glasso works, with verbosity 1", {
  copulized <- copulize(multi_omics_micro)
  expect_output(reconstruct(copulized$layers,
    net_method = "glasso", nlambda = 1,
    rep.num = 2, verbose = 1
  ))
  expect_no_error(reconstruct(copulized$layers,
    net_method = "glasso",
    sel_method = "ric", nlambda = 1, rep.num = 2, verbose = FALSE
  ))
})

test_that("sel_method errors, and 'criterion' instead of 'sel_method' work", {
  copulized <- copulize(multi_omics_micro)
  # by setting sel_method and criterion = foo all the three elements are triggered
  expect_error(suppressWarnings(reconstruct(copulized$layers,
    net_method = "glasso", sel_method = "foo",
    criterion = "foo", nlambda = 3, rep.num = 3
  )))
  expect_error(suppressWarnings(reconstruct(copulized$layers,
    net_method = "correlation", sel_method = "foo",
    criterion = "foo", nlambda = 3, rep.num = 3
  )))
  expect_error(suppressWarnings(reconstruct(copulized$layers,
    net_method = "mb", sel_method = "foo",
    criterion = "foo", nlambda = 3, rep.num = 3
  )))
})

# test_that("sel_method errors, warnings, and 'criterion' instead of 'sel_method' work", {
#  copulized <- copulize(multi_omics_micro)
#  # by setting sel_method and criterion = foo all the three elements are triggered
#  expect_error(reconstruct(copulized$layers, net_method = "glasso", sel_method = "foo",
#                           criterion = "foo", nlambda=3, rep.num = 3))
#  expect_error(reconstruct(copulized$layers, net_method = "correlation", sel_method = "foo",
#                           criterion = "foo", nlambda=3, rep.num = 3))
#  expect_error(reconstruct(copulized$layers, net_method = "mb", sel_method = "foo",
#                           criterion = "foo", nlambda=3, rep.num = 3))
# })

test_that("reconstruction with coglasso works with verbosity 1", {
  copulized <- copulize(multi_omics_micro)
  expect_output(reconstruct(copulized$layers,
    net_method = "coglasso",
    nlambda_w = 1, nlambda_b = 1, c = 1, rep.num = 2,
    verbose = 1
  ))
})

test_that("reconstruction with ct works, with verbosity 1", {
  copulized <- copulize(multi_omics_micro)
  expect_output(reconstruct(copulized$layers,
    net_method = "correlation",
    nlambda = 1, rep.num = 2, verbose = 1
  ))
})

test_that("reconstruction with mb works, with verbosity 1", {
  copulized <- copulize(multi_omics_micro)
  expect_output(reconstruct(copulized$layers,
    net_method = "mb", nlambda = 1,
    rep.num = 2, verbose = 1
  ))
})

test_that("minimal output works", {
  copulized <- copulize(multi_omics_micro)
  expect_length(reconstruct(copulized$layers,
    net_method = "coglasso",
    nlambda_w = 1, nlambda_b = 1, c = 1, rep.num = 2,
    minimal_output = TRUE
  ), 5)
})
