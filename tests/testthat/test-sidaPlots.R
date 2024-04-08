library(mvlearnR)

test_that("CorrelationPlots returns ggplot on two-view data", {
  expect_equal({
    load("../../chacon_scripts/sida_and_sidanet_vignette_models.Rdata")
    sida_filtered = filterOmics
    sida_fit = fit.cvsida
    x = CorrelationPlots(sida_filtered$X, 
                         Ytest=sida_filtered$Y, sida_fit$hatalpha)
    class(x)
  }, 
  c("gg", "ggplot"))
})

test_that("CorrelationPlots returns list of ggplots on three-view data", {
  expect_equal({
    load("../../chacon_scripts/three_view_data.Rdata")
    three_filtered = filterOmics
    three_fit = fit.cvsida
    x = CorrelationPlots(three_filtered$X, 
                         Ytest=three_filtered$Y, three_fit$hatalpha)
    c(class(x), class(x[[1]]))
  }, 
  c("list","gg", "ggplot"))
})

test_that("DiscriminantPlots returns two ggplots on two-view data", {
  expect_equal({
    load("../../chacon_scripts/sida_and_sidanet_vignette_models.Rdata")
    sida_filtered = filterOmics
    sida_fit = fit.cvsida
    x = DiscriminantPlots(sida_filtered$X, sida_filtered$Y, sida_fit$hatalpha)
    c(length(x), class(x), class(x[[1]]))
  }, 
  c(2, "list", "gg", "ggplot"))
})

test_that("DiscriminantPlots returns three ggplots on tthree-view data", {
  expect_equal({
    load("../../chacon_scripts/three_view_data.Rdata")
    three_filtered = filterOmics
    three_fit = fit.cvsida
    x = DiscriminantPlots(three_filtered$X, three_filtered$Y, three_fit$hatalpha)
    c(length(x), class(x), class(x[[1]]))
  }, 
  c(3, "list", "gg", "ggplot"))
})