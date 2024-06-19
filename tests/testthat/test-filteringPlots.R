test_that("umapPlot returns a list of two ggplots on sida filtered ", {
  expect_equal({
    data = system.file("extdata", "sida_and_sidanet_vignette_models.Rdata", package = "mvlearnR")
    load(data)
    sida_filtered = filterOmics
    sida_fit = fit.cvsida
    x = umapPlot(sida_filtered, plotIt = FALSE)
    c(class(x), length(x), class(x[[1]]))
  }, 
  c("list", "2", "gg", "ggplot"))
})

test_that("volcanoPlot returns a list of two ggplots on sida filtered ", {
  expect_equal({
    data = system.file("extdata", "sida_and_sidanet_vignette_models.Rdata", package = "mvlearnR")
    load(data)
    sida_filtered = filterOmics
    sida_fit = fit.cvsida
    x = volcanoPlot(sida_filtered, plotIt = FALSE)
    c(class(x), length(x), class(x[[1]]))
  }, 
  c("list", "2", "gg", "ggplot"))
})
