test_that("BT_test_result", {
  skip_on_cran()
  X <- readRDS(test_path("fixtures", "X.rds"))
  set.seed(1)
  # Output
  result <- result_num <- BTtest(X = X)
  names(result_num) <- NULL
  expect_identical(result_num, c(1, 1, 2))
  expect_identical(names(result), paste0("r_", 1:3, "_hat"))
  # Selecting a specific R
  expect_no_error(BTtest(X = X, R = 100))
  # BT2
  expect_no_error(BTtest(X = X, BT1 = FALSE))
  # Test very long panel
  X_long <- sim_DGP(N = 10, n_Periods = 200)
  expect_no_error(BTtest(X = X_long))
})


test_that("BT_test_inputs", {
  skip_on_cran()
  X <- readRDS(test_path("fixtures", "X.rds"))
  set.seed(1)
  # Incorrect inputs
  expect_error(BTtest(X = X, r_max = 0))
  expect_error(BTtest(X = X, alpha = 1))
  expect_error(BTtest(X = X, BT1 = 1))
  X_char <- matrix(as.character(X), ncol = ncol(X))
  expect_error(BTtest(X = X_char))
  # Correct inputs
  result <- BTtest(X = as.data.frame(X))
})
