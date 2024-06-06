
## Correlationplots

test_that("CorrelationPlots returns ggplot on two-view data", {
  expect_equal({
    data = system.file("extdata", "sida_and_sidanet_vignette_models.Rdata", package = "mvlearnR")
    load(data)
    sida_filtered = filterOmics
    sida_fit = fit.cvsida
    x = CorrelationPlots(sida_filtered$X, 
                         Ytest=sida_filtered$Y, sida_fit$hatalpha, plotIt = FALSE)
    class(x)
  }, 
  c("gg", "ggplot"))
})

test_that("CorrelationPlots returns list of ggplots on three-view data", {
  expect_equal({
    data = system.file("extdata", "three_view_data.Rdata", package = "mvlearnR")
    load(data)
    three_filtered = filterOmics
    three_fit = fit.cvsida
    x = CorrelationPlots(three_filtered$X, 
                         Ytest=three_filtered$Y, three_fit$hatalpha, plotIt = FALSE)
    c(class(x), class(x[[1]]))
  }, 
  c("list","gg", "ggplot"))
})

test_that("CorrelationPlots returns ggplot on two-view sidanet", {
  expect_equal({
    data = system.file("extdata", "sida_and_sidanet_vignette_models.Rdata", package = "mvlearnR")
    load(data)
    x = CorrelationPlots(sidanetData[[1]], 
                         Ytest=sidanetData[[2]], mycvsidanet$hatalpha, plotIt = FALSE)
    c(class(x))
  }, 
  c("gg", "ggplot"))
})

## Discriminant Plots

test_that("DiscriminantPlots returns two ggplots on two-view data", {
  expect_equal({
    data = system.file("extdata", "sida_and_sidanet_vignette_models.Rdata", package = "mvlearnR")
    load(data)
    sida_filtered = filterOmics
    sida_fit = fit.cvsida
    x = DiscriminantPlots(sida_filtered$X, sida_filtered$Y, sida_fit$hatalpha, plotIt = FALSE)
    c(length(x), class(x), class(x[[1]]))
  }, 
  c(2, "list", "gg", "ggplot"))
})

test_that("DiscriminantPlots returns three ggplots on tthree-view data", {
  expect_equal({
    data = system.file("extdata", "three_view_data.Rdata", package = "mvlearnR")
    load(data)
    three_filtered = filterOmics
    three_fit = fit.cvsida
    x = DiscriminantPlots(three_filtered$X, three_filtered$Y, three_fit$hatalpha, plotIt = FALSE)
    c(length(x), class(x), class(x[[1]]))
  }, 
  c(3, "list", "gg", "ggplot"))
})

test_that("DiscriminantPlots returns two ggplots on sidanet data with three classes", {
  expect_equal({
    data = system.file("extdata", "sida_and_sidanet_vignette_models.Rdata", package = "mvlearnR")
    load(data)
    x = DiscriminantPlots(sidanetData[[1]], sidanetData[[2]], mycvsidanet$hatalpha, plotIt = FALSE)
    c(length(x), class(x), class(x[[1]]))
  }, 
  c(2, "list", "gg", "ggplot"))
})

### LoadingsPlots

test_that("LoadingsPlots returns two ggplots on two-view sidanet data", {
  expect_equal({
    data = system.file("extdata", "sida_and_sidanet_vignette_models.Rdata", package = "mvlearnR")
    load(data)
    sidanet_fit = mycvsidanet
    x = LoadingsPlots(sidanet_fit, keep.loadings = c(7,7), plotIt = FALSE)
    c(length(x), class(x), class(x[[1]]))
  }, 
  c(2, "list", "gg", "ggplot"))
})

test_that("LoadingsPlots works without supplying keep.loadings", {
  expect_equal({
    data = system.file("extdata", "sida_and_sidanet_vignette_models.Rdata", package = "mvlearnR")
    load(data)
    sidanet_fit = mycvsidanet
    x = LoadingsPlots(sidanet_fit, plotIt = FALSE)
    c(length(x), class(x), class(x[[1]]))
  }, 
  c(2, "list", "gg", "ggplot"))
})

test_that("LoadingsPlots throws error on SIDA data with 1 discriminant vector", {
  expect_error({
    data = system.file("extdata", "sida_and_sidanet_vignette_models.Rdata", package = "mvlearnR")
    load(data)
    sidanet_fit = mycvsidanet
    x = LoadingsPlots(sida_fit, keep.loadings = c(7,7), plotIt = FALSE)
    })
})

test_that("Loadings returns a dataframe with four columns and 131 rows on sidanet e.g. data", {
  expect_equal({
    data = system.file("extdata", "sida_and_sidanet_vignette_models.Rdata", package = "mvlearnR")
    load(data)
    sidanet_fit = mycvsidanet
    x = Loadings(sidanet_fit)
    c(class(x), ncol(x), nrow(x))
  },
  c("data.frame", "4", "131"))
})


### WithinviewBiplots
test_that("WithinViewBiplot returns a list of two ggplots objects on sidanet eg", {
  expect_equal({
    data = system.file("extdata", "sida_and_sidanet_vignette_models.Rdata", package = "mvlearnR")
    load(data)
    sidanet_fit = mycvsidanet
    x = WithinViewBiplot(sidanet_fit,sidanetData[[2]],Xtest=NULL, 
                         keep.loadings = c(3,3), plotIt = FALSE)
    c(length(x), class(x), class(x[[1]]))
  }, 
  c(2, "list", "gg", "ggplot"))
})

test_that("WithinViewBiplot works without keep.loadings", {
  expect_equal({
    data = system.file("extdata", "sida_and_sidanet_vignette_models.Rdata", package = "mvlearnR")
    load(data)
    sidanet_fit = mycvsidanet
    x = WithinViewBiplot(sidanet_fit,sidanetData[[2]],Xtest=NULL, 
                         plotIt = FALSE)
    c(length(x), class(x), class(x[[1]]))
  }, 
  c(2, "list", "gg", "ggplot"))
})



### BetweenviewBiplots
test_that("BetweenViewBiplot returns a single ggplots objects on sidanet eg", {
  expect_equal({
    data = system.file("extdata", "sida_and_sidanet_vignette_models.Rdata", package = "mvlearnR")
    load(data)
    sidanet_fit = mycvsidanet
    x = BetweenViewBiplot(sidanet_fit,sidanetData[[2]],Xtest=NULL, 
                         keep.loadings = c(3,3), plotIt = FALSE)
    c(length(x), class(x), class(x[[1]]))
  }, 
  c(1, "list", "gg", "ggplot"))
})

test_that("BetweenViewBiplot works without keep.loadings", {
  expect_equal({
    data = system.file("extdata", "sida_and_sidanet_vignette_models.Rdata", package = "mvlearnR")
    load(data)
    sidanet_fit = mycvsidanet
    x = BetweenViewBiplot(sidanet_fit,sidanetData[[2]],Xtest=NULL, 
                         plotIt = FALSE)
    c(length(x), class(x), class(x[[1]]))
  }, 
  c(1, "list", "gg", "ggplot"))
})
