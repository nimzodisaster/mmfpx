#' Simulated Longitudinal Neuroimaging Measurements
#'
#' A simulated dataset containing longitudinal measurements for multiple
#' subjects with several volumetric outcomes following non-linear growth
#' trajectories. The data are meant to illustrate the workflows in the
#' package, including fractional polynomial mixed-effects modelling and
#' normalization utilities.
#'
#' @format A tibble with 559 rows and 9 variables:
#' \describe{
#'   \item{subj\_id}{Integer subject identifier, ranging from 1 to 120.}
#'   \item{diagnosis}{Character indicator of group membership (e.g., patient vs control).}
#'   \item{visit}{Numeric visit number.}
#'   \item{age}{Positive numeric age measurement at each visit.}
#'   \item{volumeA}{Simulated volumetric outcome exhibiting moderate curvature.}
#'   \item{volumeB}{Outcome with stronger non-linear scaling.}
#'   \item{volumeC}{Outcome with rapidly changing magnitude across age.}
#'   \item{volumeD}{Outcome with mild non-linearity.}
#'   \item{volumeE}{Outcome with extreme scaling behaviour.}
#' }
#'
#' @source Simulated for package examples.
"simulated_longitudinal_data"
