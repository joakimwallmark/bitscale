#' Calculate derivative of the bit scores at the provided theta locations with respect to theta.
#'
#' This internal function computes the cumulative bit score derivative for a given
#' mirt model at specified theta values. It is used as a helper function in the
#' computation of bit score information.
#'
#' @param model A fitted mirt model object.
#' @param theta A matrix of theta scores
#' @param bit_scale_items A numeric vector indicating which items to use for the bit scale computation.
#'
#' @return A vector representing the cumulative bit score derivative.
#'
#' @keywords internal
#' @noRd
bit_score_gradient <- function(model, theta, bit_scale_items) {
  theta <- matrix(theta, ncol = 1)
  gradient <- 0
  for (bit_scale_item in bit_scale_items) {
    item_probs <- probtrace(extract.item(model, bit_scale_item), Theta = theta)
    item_deriv_list <- mirt:::DerivTheta(model@ParObjects$pars[[bit_scale_item]], theta)$grad
    item_derivs <- do.call(cbind, item_deriv_list)
    entropy_gradient_summands <- log2(item_probs) * item_derivs

    item_entropy_derivs <- -rowSums((log2(item_probs) + 1) * item_derivs)
    squared_entropy_deriv <- item_entropy_derivs ** 2
    gradient <- gradient + abs(rowSums(entropy_gradient_summands))
  }
  return(gradient)
}
