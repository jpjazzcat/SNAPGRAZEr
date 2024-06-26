#' @title Loss of biomass during the off season
#'
#' @description Estimating the amount of biomass lost during the dormant season (dry, non-growing season).
#' @param Co Daily biomass consumption rate during the offseason (g/animal/day).
#' @param Cg Daily biomass consumption rate during the growing season (g/animal/day).
#' @param Gdays Total number of days in the growing season.
#' @param DdaysO Number of days in offseason grazing episode. Default = 365 - Gdays.
#' @param d_off Stocking density during the offseason (head/ha)
#' @export

calc_Lo= function(Cg, Co = NA, Gdays, DdaysO = NA, d_off) {

  if(is.na(Co)) {
    Co = Cg/2
  }

  if(is.na(DdaysO)) {
    DdaysO = 365 - Gdays
  }

  Lo = (Co*DdaysO*d_off*10^(-4))

  return(Lo)

}
