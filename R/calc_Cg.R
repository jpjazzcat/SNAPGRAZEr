#' @title Daily biomass consumption during the growing season
#'
#' @description Calculation of the parameter Cg, the daily biomass consumption rate during the growing season, from Ritchie (2020)
#' @param W Average animal body size (kg live weight)
#' @export

calc_Cg = function(W) {
  Cg = (5300+770*log(W))
  return(Cg)
}
