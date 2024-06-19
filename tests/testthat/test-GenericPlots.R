# VarImportancePlots
test_that("VarImportance returns dataframe on SIDA vignette", {
  expect_equal({
    data = system.file("extdata", "sida_and_sidanet_vignette_models.Rdata", package = "mvlearnR")
    load(data)
    sida_fit = fit.cvsida
    class(VarImportance(sida_fit))
  }, 
  c("data.frame"))
})


test_that("VarImportance returns dataframe on SIDANet vignette", {
  expect_equal({
    data = system.file("extdata", "sida_and_sidanet_vignette_models.Rdata", package = "mvlearnR")
    load(data)
    sidanet_fit = mycvsidanet
    class(VarImportance(sidanet_fit))
  }, 
  c("data.frame"))
})

test_that("VarImportance returns a list of two ggplots on SIDA vignette", {
  expect_equal({
    data = system.file("extdata", "sida_and_sidanet_vignette_models.Rdata", package = "mvlearnR")
    load(data)
    sida_fit = fit.cvsida
    x = VarImportancePlot(sida_fit, plotIt = FALSE)
    c(length(x), class(x), class(x[[1]]))
  }, 
  c(2, "list","gg", "ggplot"))
})

test_that("VarImportance returns nothing when plotIt == TRUE", {
  expect_equal({
    data = system.file("extdata", "sida_and_sidanet_vignette_models.Rdata", package = "mvlearnR")
    load(data)
    sida_fit = fit.cvsida
    x = VarImportancePlot(sida_fit)
    x
  }, 
  NULL)
})

test_that("VarImportance returns a list of three ggplots on threeview data", {
  expect_equal({
    data = system.file("extdata", "three_view_data.Rdata", package = "mvlearnR")
    load(data)    
    three_fit = fit.cvsida
    x = VarImportancePlot(three_fit, plotIt = FALSE)
    c(length(x), class(x), class(x[[1]]))
  }, 
  c(3, "list","gg", "ggplot"))
})