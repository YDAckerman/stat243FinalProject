#' plot
#'
#' Create summary plots for the result of the multiple factor analysis
#' @param mfa - the mfa object
#' @param dims - the dimensions of the resulting matrices to plot
#' @return a three plots with the projection of each observation onto
#'         the first two extracted components, the partial factor scores,
#'         and the variable loadings
#' @export
plot.mfa <- function(mfa, dims = 1:2){
  check_mfa(mfa)
  check_plot_dims(length(attributes(mfa)$sets), dims)
  plot_compromise(mfa, dims)
  plot_partial_factor(mfa, dims)
  plot_variable_loadings(mfa, dims)
}

# plot(mfa1)
