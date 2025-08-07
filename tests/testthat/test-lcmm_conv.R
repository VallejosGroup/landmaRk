test_that("lcmm convergence messages are raised", {
  example <- list(conv = 1)
  expect_message(
    .check_lcmm_convergence(example),
    "LCMM model converged successfully."
  )

  example <- list(conv = 2)
  expect_warning(
    .check_lcmm_convergence(example),
    "Maximum number of iterations reached without convergence."
  )

  example <- list(conv = 3)
  expect_message(
    .check_lcmm_convergence(example),
    "Convergence criteria satisfied with a partial Hessian matrix."
  )

  example <- list(conv = 4)
  expect_error(
    .check_lcmm_convergence(example),
    "Problem occurred during optimisation of the LCMM model."
  )
})
