#' Plot item response curves or item information on the bit score scale
#'
#' This function takes a fitted unidimensional mirt model object and plots the item
#' response curves (category probability curves) for a given item on the bit scale.
#'
#' @importFrom mirt extract.mirt extract.item probtrace
#' @importFrom ggplot2 ggplot aes geom_line labs theme_bw scale_x_continuous
#'
#' @param model A fitted unidimensional mirt model object from the mirt package.
#' @param item A single numeric value indicating which item to plot.
#' @param type A character string indicating the type of plot to produce. Can be item trace lines ('trace) or information ('info').
#' @param bit_grid_size An integer specifying the internal grid size for \code{\link{bit_scores}}. A larger value produces more accurate bit score estimates. Default is 10000.
#' @param bit_scale_items A numeric vector indicating which items to use for bit score computation. By default, all items are used.
#' @param thetas_to_plot The range of theta values to plot the test information curve. The bit scores and the item curves are computed at these values. When type = "trace", negative/positive infinity is always added to the beginning/end of the sequence. Default is from -9 to 9 with 501 points.
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
#' bit_plot_item(mirt_model, item = 1)
#' bit_plot_item(mirt_model, item = 1, type = "info")
#' }
#'
#' @export
bit_plot_item <- function(
    model,
    item,
    type = "trace",
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
  if (!is.numeric(item) || length(item) != 1) {
    stop("The 'item' argument must be a single numeric value.")
  }
  if (!(is.numeric(bit_grid_size) && bit_grid_size > 0)) {
    stop("The 'bit_grid_size' argument must be a positive integer.")
  }

  if (type != "info") {
    thetas_to_plot[1] <- -Inf
    thetas_to_plot[length(thetas_to_plot)] <- Inf
  }

  bit_vals <- bit_scores(
    model = model,
    thetas = matrix(thetas_to_plot, ncol = 1),
    items = bit_scale_items,
    grid_size = bit_grid_size,
    return_grid = FALSE
  )
  bit_vals <- as.vector(bit_vals)
  if (type == "trace") {
    mirt_item <- extract.item(model, item)
    probs <- probtrace(mirt_item, thetas_to_plot)

    prob_df <- data.frame(
      bit = rep(bit_vals, times = ncol(probs)),
      response = factor(rep(seq_len(ncol(probs)) - 1, each = length(thetas_to_plot))),
      probability = as.vector(probs)
    )

    p <- ggplot(prob_df, aes(x = bit, y = probability, color = response)) +
      geom_line(size = 1.1) +
      scale_x_continuous(
        limits = range(prob_df$bit),
        expand = c(0, 0)
      ) +
      theme_bw() +
      labs(
        x = "Bit Score",
        y = "Probability",
        color = "Response",
        title = paste("Item Response Curves for Item", item)
      )

  } else if (type == "info") {
    info <- bit_score_information(model, matrix(thetas_to_plot, ncol = 1), item = item, bit_scale_items = bit_scale_items)

    info_df <- data.frame(
      bit = bit_vals,
      information = info
    )

    p <- ggplot(info_df, aes(x = bit, y = information)) +
      geom_line(size = 1.1, color = "steelblue") +
      scale_x_continuous(
        limits = range(info_df$bit),
        expand = c(0, 0)
      ) +
      theme_bw() +
      labs(
        x = "Bit Score",
        y = "Information",
        title = paste("Information Curve for Item", item)
      )
  } else {
    stop("Invalid 'type' argument. Choose 'trace' or 'info'.")
  }

  return(p)
}
