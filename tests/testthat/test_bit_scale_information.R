library(mirt)

# Create a small example mirt model for testing
set.seed(123)
a <- matrix(rlnorm(30, meanlog = 0, sdlog = 0.5))
d <- matrix(rnorm(30, 0, 1))
data <- simdata(a, d, 500, itemtype = "2PL")
model <- mirt(data, 1, verbose = FALSE)
thetas <- as.matrix(matrix(seq(-3, 3, length.out = 7), ncol = 1))

test_that("returns fisher information using testinfo when no item specified", {
  info <- bit_score_information(model, thetas)

  expect_true(is.numeric(info))
  expect_equal(length(info), nrow(thetas))
})

test_that("returns fisher information for a specified item", {
  info_item <- bit_score_information(model, thetas, item = 1)

  expect_true(is.numeric(info_item))
  expect_equal(length(info_item), nrow(thetas))
})
