#' Group PALM/iPALM localizations based on their physical separation distance
#'
#' PALM/iPALM data results in a list of spatial coordinates for fluorophore
#' localizations. This function groups nearby localizations if they are within
#' the provided critical distance from each other.
#'
#' @title Return the group number for each localization.
#'
#' @param input A numeric matrix where each row is a localization and each
#' column is a spatial axis.
#' @param critDist The critical distance for which localizations nearer than
#' this distance are deemed part of the same group.
#'
#' @author Zach Colburn
#'
#' @examples
#' # Generate random data.
#' set.seed(10)
#' input <- as.matrix(data.frame(x=rnorm(10),y=rnorm(10)))
#'
#' # Perform linking.
#' euclidean_linker(input, 0.4)
#'
#' @export
#'
#' @importFrom assertthat assert_that is.number
euclidean_linker <- function(input, critDist) {
  assert_that(class(input) == "matrix")
  assert_that(length(input) > 1)
  assert_that(class(input[1]) %in% c("integer", "numeric"))
  assert_that(nrow(input) > 1)
  assert_that(ncol(input) > 0)
  assert_that(is.number(critDist))
  assert_that(critDist > 0)

  .euclidean_linker_cpp(input, critDist)
}

