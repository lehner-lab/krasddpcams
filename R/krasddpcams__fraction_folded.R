
#' krasddpcams__fraction_folded
#'
#' Calculate fraction folded in 2-state model.
#'
#' @param folding_energy folding energies for 2-state model (required)
#' @param RT constant (default:0.001987*(273+30))
#'
#' @return Fraction folded
#' @export
krasddpcams__fraction_folded <- function(
  folding_energy,
  RT = 0.001987*(273+30)
  ){
  return(1/(1+exp(folding_energy/RT)))
}