#' Convert theta scores to bit scores
#'
#' This function takes a matrix of theta scores from a fitted IRT model and converts them to bit scores.
#'
#' @importFrom mirt extract.mirt extract.item probtrace
#'
#' @param model An fitted mirt model object from the mirt package.
#' @param theta A vector or matrix with theta scores from which to compute bit scores. Typically the matrix returned from the mirt::fscores method. If a second column is present, it is assumed to be the standard errors of the theta scores and bit score standard errors are computed.
#' @param items A numeric vector indicating which items to use for computation. By default, all items are used.
#' @param grid_size An integer specifying the size of the theta grid used for bit score computation. A higher value leads to improved accuracy. Default is 10000.
#' @param return_grid Whether or not to return the bit score for each value in the grid used for computation or only the bit scores for the input thetas. Default is FALSE.
#'
#' @return A numeric matrix. If `return_grid` is TRUE, the matrix has two columns named "theta" and "bit_score", corresponding to the complete theta grid and the associated computed bit scores. If `return_grid` is FALSE, the matrix contains a single column of bit scores corresponding to the input theta values, preserving their original order.
#'
#' @examples
#' \dontrun{
#' library(mirt)
#' # As an example, we simulate a dataset with 30 items and 1000 respondents
#' # discrimination parameters from a log-normal distribution
#' a <- matrix(rlnorm(30, meanlog = 0, sdlog = 0.5))
#' # difficulty parameters from a standard normal distribution
#' d <- matrix(rnorm(30, 0, 1))
#' data <- simdata(a, d, 1000, itemtype = '2PL')
#' # Fit the model and compute theta scores
#' mirt_model <- mirt(data, 1)
#' # Compute the bit scores
#' thetas <- fscores(mirt_model, full.scores.SE = FALSE)
#' bit <- bit_scores(mirt_model, thetas)
#' # Compute the bit scores with standard errors
#' thetas <- fscores(mirt_model, method="ML", full.scores.SE = TRUE)
#' bit_with_se <- bit_scores(mirt_model, thetas)
#' }
#'
#' @export
bit_scores <- function(
    model,
    theta,
    items = 1:extract.mirt(model, "nitems"),
    grid_size = 10000,
    return_grid = FALSE
) {
  compute_se <- FALSE
  if (!"SingleGroupClass" %in% class(model)) {
    stop("The model object must be a fitted mirt model object.")
  }
  if (model@Model$model != 1) {
    stop("The model must be unidimensional (1 factor model).")
  }
  if (!is.numeric(theta) && !is.matrix(theta)) {
    stop("The 'thetas' argument must be a numeric vector or a one/two column matrix.")
  }
  if (!is.numeric(items)) {
    stop("The 'items' argument must be a numeric vector.")
  }
  if (!(is.numeric(grid_size) && grid_size > 0)) {
    stop("The 'grid_size' argument must be a positive integer.")
  }
  if (!is.logical(return_grid)) {
    stop("The 'return_grid' argument must be a logical.")
  }
  if (is.matrix(theta)) {
    if (ncol(theta) == 2) {
      compute_se <- TRUE
      theta_se <- theta[, 2]
    }
    theta <- theta[, 1]
  }
  else {
    theta <- c(theta)
  }

  theta_grid <- setdiff(seq(-10, 10, length.out = grid_size), theta) |> # remove potential duplicates from grid
    c(theta) |> # merge with sample thetas
    sort()
  theta_grid[1] <- -200
  theta_grid[length(theta_grid)] <- 200

  bit_scores <- vector("numeric", length = length(theta_grid))
  for (item in items) {
    mirt_item <- extract.item(model, item)
    probs <- probtrace(mirt_item, theta_grid)
    total_item_dist <- 0
    entropies <- rowSums(-probs * log(probs, base = 2))
    for (i in 2:length(entropies)) {
      # add up total item distance to theta corresponding to surp[i, ]
      total_item_dist <- total_item_dist + sum(abs(entropies[i] - entropies[i - 1]))
      # add this to the total for the same theta
      bit_scores[i] <- bit_scores[i] + total_item_dist
    }
  }
  if (return_grid) {
    mat <- as.matrix(cbind(theta = theta_grid, bit_score = bit_scores))
    colnames(mat) <- c("theta", "bit_score")
    return(mat)
  } else {
    # Extract the entropy scores from the grid
    matching_indices <- which(theta_grid %in% theta)
    ordered_bit_scores <- bit_scores[matching_indices]
    # Order the results based on the original order of the 'thetas' argument
    order_vec <- order(order(theta))
    return_matrix <- as.matrix(ordered_bit_scores[order_vec])

    if (compute_se) {
      bit_score_gradient <- bit_score_gradient(model, theta, items)
      return_matrix <- cbind(return_matrix, theta_se * bit_score_gradient)
    }
    return(return_matrix)
  }
}

