
test_that("performanceMetrics returns a dataframe with 10 rows on binomial data", {
  expect_equal({
    data = system.file("extdata", "sida_and_sidanet_vignette_models.Rdata", package = "mvlearnR")
    load(data)
    sida_filtered = filterOmics
    sida_fit = fit.cvsida
    Y.pred=sida_fit$PredictedClass.train-1 #to get this in 0 and 1
    Y.train=sida_filtered$Y
    train.metrics=PerformanceMetrics(Y.pred,Y.train,family='binomial')
    c(nrow(train.metrics), class(train.metrics))
  }, 
  c("10", "data.frame"))
})


test_that("performanceMetricsPlot returns a ggplot", {
  expect_equal({
    data = system.file("extdata", "sida_and_sidanet_vignette_models.Rdata", package = "mvlearnR")
    load(data)
    sida_filtered = filterOmics
    sida_fit = fit.cvsida
    Y.pred=sida_fit$PredictedClass.train-1 #to get this in 0 and 1
    Y.train=sida_filtered$Y
    x=PerformanceMetricsPlot(Y.pred,Y.train,family='binomial')
    class(x)
  }, 
  c("gg", "ggplot"))
})