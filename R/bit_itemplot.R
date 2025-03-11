#' Plot item response curves on the bit score scale
#'
#' This function takes a fitted unidimensional mirt model object and plots the item
#' response curves (category probability curves) for a given item on the bit scale.
#'
#' @importFrom mirt extract.mirt extract.item probtrace
#' @importFrom ggplot2 ggplot aes geom_line labs theme_bw scale_x_continuous
#'
#' @param model A fitted unidimensional mirt model object from the mirt package.
#' @param item A single numeric value indicating which item to plot.
#' @param bit_grid_size An integer specifying the internal grid size for
#'   \code{\link{bit_scores}}. A larger value produces more accurate bit
#'   score estimates. Default is 10000.
#'
#' @return A ggplot2 object displaying the item response curves on the bit scale.
#'
#' @examples
#' \dontrun{
#' library(mirt)
#' # Example dataset
#' a <- matrix(rlnorm(30, meanlog = 0, sdlog = 0.5))
#' d <- matrix(rnorm(30, 0, 1))
#' data <- simdata(a, d, 1000, itemtype = '2PL')
#' # Fit the model
#' mirt_model <- mirt(data, 1)
#' # Plot the response curves for item 1 on bit scale
#' bit_itemplot(mirt_model, item = 1)
#' }
#'
#' @export
bit_itemplot <- function(
    model,
    item,
    bit_grid_size = 10000
) {
  # Basic checks
  if (!"SingleGroupClass" %in% class(model)) {
    stop("The model object must be a fitted mirt model object.")
  }
  if (model@Model$model != 1) {
    stop("The model must be unidimensional (1 factor model).")
  }
  if (!is.numeric(item) || length(item) != 1) {
    stop("The 'item' argument must be a single numeric value.")
  }
  if (!(is.numeric(bit_grid_size) && bit_grid_size > 0)) {
    stop("The 'bit_grid_size' argument must be a positive integer.")
  }

  theta_to_plot <- seq(-10, 10, length.out = 1000)
  theta_to_plot[1] <- -Inf
  theta_to_plot[length(theta_to_plot)] <- Inf
  mirt_item <- extract.item(model, item)
  probs <- probtrace(mirt_item, theta_to_plot)

  bit_vals <- bit_scores(
    model = model,
    thetas = matrix(theta_to_plot, ncol = 1),
    grid_size = bit_grid_size,
    return_grid = FALSE
  )

  bit_vals <- as.vector(bit_vals)
  prob_df <- data.frame(
    bit = rep(bit_vals, times = ncol(probs)),
    response = factor(rep(seq_len(ncol(probs)), each = length(theta_to_plot))),
    probability = as.vector(probs)
  )

  p <- ggplot(prob_df, aes(x = .data$bit, y = .data$probability, color = .data$response)) +
    geom_line(size = 1.1) +  # slightly thicker lines
    scale_x_continuous(
      limits = range(prob_df$bit),
      expand = c(0, 0)
    ) +
    theme_bw() +
    labs(
      x = "Bit Score",
      y = "Probability",
      color = "Response",
      title = paste("Item Response Curves for Item", item, "on Bit Scale")
    )

  return(p)
}
