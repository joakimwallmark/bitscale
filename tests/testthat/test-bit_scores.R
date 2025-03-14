library(mirt)

# Create a small example mirt model for testing
set.seed(123)
a <- matrix(rlnorm(30, meanlog = 0, sdlog = 0.5))
d <- matrix(rnorm(30, 0, 1))
data <- simdata(a, d, 500, itemtype = "2PL")
mirt_model <- mirt(data, 1)
thetas <- fscores(mirt_model, full.scores.SE = FALSE)
thetas <- as.matrix(thetas)

test_that("bit_scores returns a numeric matrix with grid output", {
  result <- bit_scores(mirt_model, thetas, return_grid = TRUE)
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 2)
  expect_equal(colnames(result), c("theta", "bit_score"))
})

test_that("bit_scores returns a numeric matrix for input thetas", {
  result <- bit_scores(mirt_model, thetas, return_grid = FALSE)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), nrow(thetas))
  expect_equal(ncol(result), 1)
})

test_that("bit_scores returns a numeric matrix with SEs", {
  result <- bit_scores(mirt_model, thetas, return_grid = FALSE, compute_SEs = TRUE)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), nrow(thetas))
  expect_equal(ncol(result), 2)
})

test_that("bit_scores errors for non-mirt model", {
  expect_error(bit_scores("not a model", thetas))
})

test_that("bit_scores errors when items is not numeric", {
  expect_error(bit_scores(mirt_model, thetas, items = "1:5"))
})

test_that("bit_scores errors when grid_size is not a positive integer", {
  expect_error(bit_scores(mirt_model, thetas, grid_size = -100))
})

test_that("bit_scores errors when return_grid is not logical", {
  expect_error(bit_scores(mirt_model, thetas, return_grid = "yes"))
})

test_that("higher theta yields a higher computed bit score", {
  custom_thetas <- matrix(seq(-4, 4, length.out = 50), ncol = 1)
  computed_scores <- bit_scores(mirt_model, custom_thetas, return_grid = FALSE)
  computed_scores <- as.vector(computed_scores)
  # Check that each successive bit score is higher than the previous one.
  expect_true(all(diff(computed_scores) > 0))
})

