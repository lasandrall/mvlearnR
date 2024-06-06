test_that("DevianceTable returns data.frame on SIDA data", {
  expect_equal({
    data = system.file("extdata", "sida_and_sidanet_vignette_models.Rdata", package = "mvlearnR")
    load(data)
    data("sidaData")
    x = devianceTable(fit.cvsida, 
                      sidaData[[3]][[1]],
                      sidaData[[3]][[2]], 
                      sidaData[[4]])
    class(x)
  }, 
  c("data.frame"))
})