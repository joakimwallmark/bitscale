library(mirt)

# Create a small example mirt model for testing
set.seed(123)
a <- matrix(rlnorm(30, meanlog = 0, sdlog = 0.5))
d <- matrix(rnorm(30, 0, 1))
data <- simdata(a, d, 500, itemtype = "2PL")
mirt_model <- mirt(data, 1)
# Generate original thetas and their corresponding "true" bit scores
thetas_orig <- fscores(mirt_model, full.scores.SE = FALSE)
bit_scores_target <- bit_scores(mirt_model, thetas_orig) # Compute bit scores to use as target input

test_that("bit_scores_inverse returns numeric vector of correct length (linear)", {
  # Default is linear
  approx_thetas <- bit_scores_inverse(mirt_model, bit_scores_target)
  expect_true(is.numeric(approx_thetas))
  expect_equal(length(approx_thetas), nrow(bit_scores_target)) # nrow as bit_scores returns matrix
})

test_that("bit_scores_inverse approximates original thetas reasonably well (linear)", {
  # Use a sufficiently large grid for good approximation
  approx_thetas <- bit_scores_inverse(mirt_model, bit_scores_target, grid_size = 10000)
  # Check max absolute difference - should be small
  # Tolerance might need adjustment based on model/items/grid_size
  expect_lt(max(abs(approx_thetas - thetas_orig)), 0.01)
})

test_that("bit_scores_inverse works with spline interpolation", {
  approx_thetas_spline <- bit_scores_inverse(mirt_model, bit_scores_target,
                                             interpolation_method = "spline",
                                             grid_size = 10000)
  expect_true(is.numeric(approx_thetas_spline))
  expect_equal(length(approx_thetas_spline), nrow(bit_scores_target))
  # Spline might have slightly different accuracy, check it's also reasonable
  expect_lt(max(abs(approx_thetas_spline - thetas_orig)), 0.01) # Adjust tolerance if needed
})

test_that("bit_scores_inverse handles extrapolation with rule = 2 (default)", {
  # Generate a grid to find approximate score range
  # Using a smaller grid here to make extrapolation more likely/testable
  grid_data <- bit_scores(mirt_model, theta = 0, grid_size = 500, return_grid = TRUE)
  finite_bs <- grid_data[is.finite(grid_data[, "bit_score"]), "bit_score"]
  min_bs <- min(finite_bs, na.rm = TRUE)
  max_bs <- max(finite_bs, na.rm = TRUE)

  # Target scores slightly outside the calculated range
  outside_scores <- c(min_bs - 1, max_bs + 1)

  approx_thetas_extrap <- bit_scores_inverse(mirt_model, outside_scores, grid_size = 500)

  # Expect finite numeric values (thetas corresponding to min/max bit score)
  expect_true(is.numeric(approx_thetas_extrap))
  expect_true(all(is.finite(approx_thetas_extrap)))
  expect_equal(length(approx_thetas_extrap), 2)

  # Optionally check they are near the boundary thetas (can be approximate)
  theta_at_min_bs <- bit_scores_inverse(mirt_model, min_bs, grid_size = 500)
  theta_at_max_bs <- bit_scores_inverse(mirt_model, max_bs, grid_size = 500)
  expect_equal(approx_thetas_extrap[1], theta_at_min_bs, tolerance = 0.1)
  expect_equal(approx_thetas_extrap[2], theta_at_max_bs, tolerance = 0.1)
})

test_that("bit_scores_inverse handles extrapolation with rule = 1", {
  grid_data <- bit_scores(mirt_model, theta = 0, grid_size = 500, return_grid = TRUE)
  finite_bs <- grid_data[is.finite(grid_data[, "bit_score"]), "bit_score"]
  min_bs <- min(finite_bs, na.rm = TRUE)
  max_bs <- max(finite_bs, na.rm = TRUE)
  outside_scores <- c(min_bs - 1, max_bs + 1)

  approx_thetas_extrap_na <- bit_scores_inverse(mirt_model, outside_scores, grid_size = 500,
                                                interpolation_method = "linear",
                                                extrapolation_rule = 1)
  # Expect NA for scores outside the range when rule = 1
  expect_true(is.numeric(approx_thetas_extrap_na))
  expect_true(all(is.na(approx_thetas_extrap_na)))
  expect_equal(length(approx_thetas_extrap_na), 2)
})

test_that("bit_scores_inverse errors for non-mirt model (non-SingleGroupClass)", {
  expect_error(bit_scores_inverse("not a model", bit_scores_target),
               "The 'model' object must be a fitted mirt model object")
})

test_that("bit_scores_inverse errors for non-numeric bit_scores", {
  expect_error(bit_scores_inverse(mirt_model, c("score1", "score2")),
               "The 'bit_scores' argument must be a numeric vector.")
})

test_that("bit_scores_inverse errors for non-numeric items", {
  expect_error(bit_scores_inverse(mirt_model, bit_scores_target, items = "1,2,3"),
               "The 'items' argument must be a numeric vector.")
})

test_that("bit_scores_inverse errors for invalid grid_size", {
  expect_error(bit_scores_inverse(mirt_model, bit_scores_target, grid_size = -100),
               "The 'grid_size' argument must be a single positive integer.")
  expect_error(bit_scores_inverse(mirt_model, bit_scores_target, grid_size = 0),
               "The 'grid_size' argument must be a single positive integer.")
  expect_error(bit_scores_inverse(mirt_model, bit_scores_target, grid_size = 500.5),
               "The 'grid_size' argument must be a single positive integer.")
  expect_error(bit_scores_inverse(mirt_model, bit_scores_target, grid_size = c(100, 200)),
               "The 'grid_size' argument must be a single positive integer.")
})

test_that("bit_scores_inverse errors for invalid interpolation_method", {
  # match.arg should throw an error for non-matching arguments
  expect_error(bit_scores_inverse(mirt_model, bit_scores_target, interpolation_method = "nearest"))
})

test_that("bit_scores_inverse errors for invalid extrapolation_rule with linear", {
  expect_error(bit_scores_inverse(mirt_model, bit_scores_target,
                                  interpolation_method = "linear", extrapolation_rule = 0),
               "The 'extrapolation_rule' must be 1 or 2")
  expect_error(bit_scores_inverse(mirt_model, bit_scores_target,
                                  interpolation_method = "linear", extrapolation_rule = 3),
               "The 'extrapolation_rule' must be 1 or 2")
})

test_that("bit_scores_inverse handles empty input bit_scores", {
  approx_thetas_empty <- bit_scores_inverse(mirt_model, numeric(0))
  expect_equal(approx_thetas_empty, numeric(0))
})

test_that("bit_scores_inverse handles NA input in bit_scores", {
  scores_with_na <- c(bit_scores_target[1:2, 1], NA, bit_scores_target[4:5, 1])
  approx_thetas_na <- bit_scores_inverse(mirt_model, scores_with_na)

  expect_equal(length(approx_thetas_na), length(scores_with_na))
  # approxfun/splinefun propagate NAs
  expect_true(is.na(approx_thetas_na[3]))
  # Check that non-NA inputs produced finite outputs (assuming rule=2 or within range)
  expect_true(all(is.finite(approx_thetas_na[-3])))
})

test_that("bit_scores_inverse accuracy improves with larger grid_size", {
  # Test approximation with a small vs large grid
  approx_thetas_low_grid <- bit_scores_inverse(mirt_model, bit_scores_target, grid_size = 100)
  approx_thetas_high_grid <- bit_scores_inverse(mirt_model, bit_scores_target, grid_size = 10000)

  # Calculate max absolute errors
  max_diff_low <- max(abs(approx_thetas_low_grid - thetas_orig))
  max_diff_high <- max(abs(approx_thetas_high_grid - thetas_orig))

  # Expect the error to be significantly smaller with the higher grid size
  expect_lt(max_diff_high, max_diff_low)
  # Re-verify high grid accuracy meets expected tolerance
  expect_lt(max_diff_high, 0.01)
})

test_that("bit_scores_inverse handles subset of items", {
  items_subset <- 1:10
  bit_scores_subset <- bit_scores(mirt_model, thetas_orig, items = items_subset)
  approx_thetas_subset <- bit_scores_inverse(mirt_model, bit_scores_subset, items = items_subset, grid_size = 5000)

  expect_equal(length(approx_thetas_subset), nrow(bit_scores_subset))
  # Accuracy might be different with fewer items, adjust tolerance if needed
  expect_lt(max(abs(approx_thetas_subset - thetas_orig)), 0.05)
})
