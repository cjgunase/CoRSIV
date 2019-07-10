data("kidney_probe_data")
test_that("use", {
  expect_output(getCorrelationMatrix(22,kidney_probe_data),"done")
})
