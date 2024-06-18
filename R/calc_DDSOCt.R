#' @title Dung-Derived SOC in year t
#'
#' @description Dung-derived SOC inputs for year t (PDSOCt) (g/m2)
#' @param LIGCELL Lignin and cellulose fraction of plant biomass (0 - 1)
#' @param Ddays Number of days of grazing episode
#' @param Cg Daily consumption rate (g/animal/day). This is the average rate during the growing season - for the dormant season, the equations assume half this consumption rate.
#' @param n Number of "pastures" per total area, A.
#' @param d Stocking density (head/ha)
#' @param Lo Loss of biomass during the off season. Output of calc_Lo().
#' @export

calc_DDSOCt = function(LIGCELL, Ddays, Cg, n, d, Lo) {

  DDSOCt = LIGCELL*0.45*0.3*(Ddays*Cg*n*d*10^(-4)+Lo)
  return(DDSOCt)

}

