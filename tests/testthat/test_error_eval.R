test_that(".eval_error_str nicely formats multiple errors", {
  error_str <- c(
    "Error 1: This is the first error.\n",
    "Error 2: This is the second error.\n",
    "Error 3: This is the third error.\n"
  )

  expect_error(
    .eval_error_str(error_str),
    paste0(
      "Error 1: This is the first error.\n",
      "Additionally, the following errors occurred:\n",
      "Error 2: This is the second error.\n",
      "Error 3: This is the third error.\n"
    )
  )
})
