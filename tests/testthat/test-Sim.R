test_that("Correct simulation", {
  skip_on_cran()
  set.seed(1)
  X <- readRDS(test_path("fixtures", "X.rds"))
  sim <- sim_DGP(N = 100, n_Periods = 200, drift = TRUE, drift_I1 = TRUE, r_I1 = 2, r_I0 = 2)
  expect_equal(sim, X)
  factors <- readRDS(test_path("fixtures", "factors.rds"))
  sim_factors <- sim_DGP(N = 100, n_Periods = 200, drift = TRUE, drift_I1 = TRUE, r_I1 = 2, r_I0 = 2, return_factor = TRUE)
  expect_equal(sim_factors, factors)
  X_0 <- readRDS(test_path("fixtures", "X_0.rds"))
  sim_0 <-  sim_DGP(N = 10, n_Periods = 20, drift = 0, drift_I1 = 0, r_I1 = 0, r_I0 = 0)
  expect_equal(sim_0, X_0)
  # Various factor combinations
  expect_no_error(sim_DGP(N = 10, n_Periods = 10, drift = TRUE, drift_I1 = FALSE, r_I1 = 0, r_I0 = 1))
  expect_no_error(sim_DGP(N = 10, n_Periods = 10, drift = FALSE, drift_I1 = 0, r_I1 = 0, r_I0 = 1))
  expect_no_error(sim_DGP(N = 10, n_Periods = 10, drift = FALSE, drift_I1 = TRUE, r_I1 = 3, r_I0 = 1))
  expect_no_error(sim_DGP(N = 10, n_Periods = 10, drift = FALSE, drift_I1 = TRUE, r_I1 = 3, r_I0 = 0))
})

test_that("Sim inputs", {
  skip_on_cran()
  expect_error(sim_DGP(N = 10, n_Periods = 0))
  expect_error(sim_DGP(N = 10, n_Periods = 10, r_I0 = -1))
  expect_error(sim_DGP(N = 10, n_Periods = 10, drift = TRUE, drift_I1 = TRUE, r_I1 = 0))
  expect_error(sim_DGP(N = 10, n_Periods = 10, drift = TRUE, drift_I1 = FALSE, r_I0 = 0))
})
