#' Plot thetas against their bit transformed counterparts
#'
#' @importFrom mirt extract.mirt
#' @importFrom ggplot2 ggplot aes geom_line labs theme_bw scale_x_continuous scale_y_continuous
#'
#' @param model A fitted unidimensional mirt model object from the mirt package.
#' @param bit_grid_size An integer specifying the internal grid size for \code{\link{bit_scores}}. A larger value produces more accurate bit score estimates. Default is 10000.
#' @param bit_scale_items A numeric vector indicating which items to use for bit score computation. By default, all items are used.
#' @param thetas_to_plot The range of theta values to plot the test information curve. The bit scores and the information is computed for these values. Default is from -9 to 9 with 501 points.
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
#' bit_plot_transformation(mirt_model)
#' }
#'
#' @export
bit_plot_transformation <- function(
    model,
    bit_grid_size = 10000,
    bit_scale_items = 1:extract.mirt(model, "nitems"),
    thetas_to_plot = seq(-9, 9, length.out = 501)
) {
  # Basic checks
  if (!"SingleGroupClass" %in% class(model)) {
    stop("The model object must be a fitted mirt model object.")
  }
  if (model@Model$model != 1) {
    stop("The model must be unidimensional (1 factor model).")
  }
  if (!(is.numeric(bit_grid_size) && bit_grid_size > 0)) {
    stop("The 'bit_grid_size' argument must be a positive integer.")
  }

  bit_scores <- bit_scores(
    model = model,
    theta = matrix(thetas_to_plot, ncol = 1),
    items = bit_scale_items,
    grid_size = bit_grid_size,
    return_grid = FALSE
  )
  bit_scores <- as.vector(bit_scores)

  info_df <- data.frame(
    theta = thetas_to_plot,
    bit = bit_scores
  )

  p <- ggplot(info_df, aes(x = theta, y = bit)) +
    geom_line(linewidth = 1.1, color = "steelblue") +
    scale_x_continuous(
      limits = range(thetas_to_plot),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, max(bit_scores)),
      expand = c(0, 0)
    ) +
    theme_bw() +
    labs(
      x = expression(paste("Theta")),
      y = "Bit score",
      title = paste("Bit scale transformation")
    )

  return(p)
}
