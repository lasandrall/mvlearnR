test_that("DevianceTable returns data.frame on selp-predict data", {
  expect_equal({
    load("../../data/sida_and_sidanet_vignette_models.Rdata")
    x = devianceTable(selp_predict, 
                      sp_testdata[["Xtestdata1"]],
                      sp_testdata[["Xtestdata2"]], 
                      sp_testdata[["Ytest"]])
    class(x)
  }, 
  c("data.frame"))
})