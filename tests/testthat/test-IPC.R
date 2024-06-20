test_that("IPC_result", {
  skip_on_cran()
  X <- readRDS(test_path("fixtures", "X.rds"))
  set.seed(1)
  # Output
  result <- result_num <- BaiIPC(X = X, r_max = 10)
  names(result_num) <- NULL
  expect_identical(result_num, c(2, 2, 2))
  expect_identical(names(result), paste0("IPC_", 1:3))
})


test_that("IPC_inputs", {
  skip_on_cran()
  X <- readRDS(test_path("fixtures", "X.rds"))
  set.seed(1)
  # Incorrect inputs
  expect_error(BaiIPC(X = X, r_max = 0))
  X_char <- matrix(as.character(X), ncol = ncol(X))
  expect_error(BaiIPC(X = X_char))
  # Correct inputs
  result <- BaiIPC(X = as.data.frame(X))
})
