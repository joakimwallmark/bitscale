#' Approximate Theta Scores from Bit Scores
#'
#' This function approximates the original theta scores corresponding to a given
#' set of bit scores by inverting the relationship computed by the
#' `bit_scores` function using interpolation.
#'
#' @importFrom stats approxfun splinefun
#' @importFrom mirt extract.mirt extract.item probtrace
#'
#' @param model An fitted unidimensional mirt model object (class 'SingleGroupClass')
#'   from the mirt package.
#' @param bit_scores A numeric vector of bit scores for which the
#'   corresponding theta scores are to be approximated.
#' @param items A numeric vector indicating which items were used for the
#'   original bit score computation. Defaults to all items in the model.
#' @param grid_size An integer specifying the size of the theta grid used for
#'   generating the internal lookup table. A higher value leads to improved
#'   accuracy but increases computation time for the lookup table generation.
#'   Default is 10000.
#' @param interpolation_method The method used for interpolation. Options are
#'   "linear" (using `stats::approxfun`, default) or "spline" (using
#'   `stats::splinefun`).
#' @param extrapolation_rule Applicable only when `interpolation_method` is "linear".
#'   Specifies how extrapolation should be handled by `stats::approxfun`.
#'   Default is `2`, which returns the theta value at the nearest data extreme.
#'   Use `1` to return `NA` for points outside the grid range.
#'
#' @return A numeric vector of the same length as `bit_scores`, containing
#'   the approximated theta values. Returns `NA` for any `target_bit_scores` that
#'   are `NA` or fall outside the computed range if `extrapolation_rule = 1`.
#'
#' @examples
#' \dontrun{
#' library(mirt)
#' # Simulate data and fit model
#' set.seed(123)
#' a <- matrix(rlnorm(30, meanlog = 0, sdlog = 0.5))
#' d <- matrix(rnorm(30, 0, 1))
#' data <- simdata(a, d, 1000, itemtype = '2PL')
#' mirt_model <- mirt(data, 1)
#'
#' # Compute some bit scores
#' thetas_orig <- fscores(mirt_model, full.scores.SE = FALSE)
#' bit_scores_computed <- bit_scores(mirt_model, thetas_orig)
#' # Aproximate thetas back from bit scores
#' approximated_thetas <- bit_scores_inverse(
#'   model = mirt_model,
#'   bit_scores = bit_scores_computed,
#'   grid_size = 10000
#' )
#' # Should be only a small difference
#' (approximated_thetas-thetas_orig) |> abs() |> max()
#'
#' # Example with spline interpolation
#' approximated_thetas_spline <- bit_scores_inverse(
#'   model = mirt_model,
#'   bit_scores = bit_scores_computed,
#'   grid_size = 10000,
#'   interpolation_method = "spline"
#' )
#' (approximated_thetas_spline-thetas_orig) |> abs() |> max()
#' }
#'
#' @export
bit_scores_inverse <- function(
    model,
    bit_scores,
    items = 1:extract.mirt(model, "nitems"),
    grid_size = 10000,
    interpolation_method = c("linear", "spline"),
    extrapolation_rule = 2
) {
  if (!"SingleGroupClass" %in% class(model)) {
    stop("The 'model' object must be a fitted mirt model object (SingleGroupClass).")
  }
  if (!is.numeric(bit_scores)) {
    stop("The 'bit_scores' argument must be a numeric vector.")
  }
  if (!is.numeric(items)) {
    stop("The 'items' argument must be a numeric vector.")
  }
  if (!(is.numeric(grid_size) && length(grid_size) == 1 && grid_size > 0 && grid_size == floor(grid_size))) {
    stop("The 'grid_size' argument must be a single positive integer.")
  }
  interpolation_method <- match.arg(interpolation_method)
  if (interpolation_method == "linear" && !(extrapolation_rule %in% c(1, 2))) {
    stop("The 'extrapolation_rule' must be 1 or 2 when using linear interpolation.")
  }

  # --- Step 1: Generate the Lookup Table using the bit_scores function ---
  # Call the existing bit_scores function to get the theta <-> bit_score mapping
  # We pass a dummy theta value (like 0) just to satisfy the argument requirement,
  # as we are primarily interested in the grid itself.
  grid_data <- bit_scores(model = model,
                          theta = 0,
                          items = items,
                          grid_size = grid_size,
                          return_grid = TRUE)

  grid_thetas <- grid_data[, "theta"]
  grid_bit_scores <- grid_data[, "bit_score"]

  # --- Step 2: Prepare data for interpolation ---
  # Filter out non-finite values which might arise from -Inf/Inf thetas at edges
  # Interpolation functions generally require finite inputs.
  finite_indices <- is.finite(grid_thetas) & is.finite(grid_bit_scores)
  grid_thetas_finite <- grid_thetas[finite_indices]
  grid_bit_scores_finite <- grid_bit_scores[finite_indices]

  if (length(grid_thetas_finite) < 2) {
    stop("Insufficient finite points generated on the grid to perform interpolation. Try increasing grid_size or check the model/items.")
  }

  # --- Step 3: Create and Apply the Interpolation Function ---
  approx_theta <- NULL # Initialize

  if (interpolation_method == "linear") {
    # Create the inverse interpolation function: bit_score -> theta
    inverse_interpolator <- stats::approxfun(x = grid_bit_scores_finite,
                                             y = grid_thetas_finite,
                                             method = "linear",
                                             rule = extrapolation_rule)
    # Apply to the target bit scores
    approx_theta <- inverse_interpolator(bit_scores)

  } else if (interpolation_method == "spline") {
    # Create the inverse spline interpolation function: bit_score -> theta
    # Note: Spline interpolation might behave poorly if extrapolating far.
    # `splinefun` doesn't have a simple 'rule' argument like approxfun for extrapolation.
    # It extrapolates using the cubic polynomial from the edge interval.
    inverse_interpolator <- stats::splinefun(x = grid_bit_scores_finite,
                                             y = grid_thetas_finite,
                                             method = "fmm") # standard Catmull-Rom spline
    # Apply to the target bit scores
    approx_theta <- inverse_interpolator(bit_scores)
  }

  return(approx_theta)
}
