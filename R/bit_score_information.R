#' Approximate the bit score Fisher information for a given mirt model at the bit scale locations corresponding to the supplied theta scores.
#'
#' @importFrom mirt extract.mirt probtrace extract.item testinfo iteminfo
#'
#' @param model An fitted mirt model object from the mirt package
#' @param thetas A one column matrix with theta scores. Typically returned from the mirt::fscores method.
#' @param item If provided, the information provided by the specified item is computed as opposed to all items.
#' @param bit_scale_items A numeric vector indicating which items to use for bit score computation. By default, all items are used.
#'
#' @return An information vector for each element in thetas.
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
#' thetas <- fscores(mirt_model, full.scores.SE = FALSE)
#' # Compute the bit score test information
#' test_info <- bit_score_information(mirt_model, thetas)
#' # Compute the bit score item information for item 3
#' item_info <- bit_score_information(mirt_model, thetas, item = 3)
#' }
#'
#' @export
bit_score_information <- function(
    model,
    thetas,
    item = NULL,
    bit_scale_items = 1:extract.mirt(model, "nitems")
) {
  if (!"SingleGroupClass" %in% class(model)) {
    stop("The model object must be a fitted mirt model object.")
  }
  if (model@Model$model != 1) {
    stop("The model must be unidimensional (1 factor model).")
  }
  if (!is.numeric(bit_scale_items)) {
    stop("The 'bit_scale_items' argument must be a numeric vector.")
  }

  squared_entropy_deriv_sum <- 0
  for (bit_scale_item in bit_scale_items) {
    item_deriv_list <- mirt:::DerivTheta(model@ParObjects$pars[[bit_scale_item]], thetas)$grad
    item_derivs <- do.call(cbind, item_deriv_list)
    item_probs <- probtrace(extract.item(model, bit_scale_item), Theta = thetas)
    item_entropy_derivs <- rowSums((log(item_probs) + 1) * item_derivs)
    squared_entropy_deriv <- item_entropy_derivs ** 2
    squared_entropy_deriv_sum <- squared_entropy_deriv_sum + squared_entropy_deriv
  }
  if (!is.null(item)) {
    fisher_info <- iteminfo(extract.item(model, item), Theta = thetas) / squared_entropy_deriv_sum
  } else {
    fisher_info <- testinfo(model, thetas) / squared_entropy_deriv_sum
  }
  return(fisher_info)
}
