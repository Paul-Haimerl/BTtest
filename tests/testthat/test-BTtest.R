test_that("BT_test_result", {
  skip_on_cran()
  X <- readRDS(test_path("fixtures", "X.rds"))
  set.seed(1)
  # Output
  result <- result_num <- BTtest(X = X, r_max = 10, alpha = 0.05, BT1 = TRUE)
  names(result_num) <- NULL
  expect_identical(result_num, c(1, 1, 2))
  expect_identical(names(result), paste0("r_", 1:3, "_hat"))
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
