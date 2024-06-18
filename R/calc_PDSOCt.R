#' @title Plant-Derived SOC in year t
#'
#' @description Plant-derived SOC inputs for year t (PDSOCt) (g/m2)
#' @param BNPPt_est Output of calc_BNPPest()
#' @param Sf Biomass at the end of the growing season. Output of calc_Sf() (g/m2)
#' @param Lo Loss of biomass during the non-growing season. Output of calc_Lo() (g/m2)
#' @param LIGCELL Lignin and cellulose fraction of plant biomass (0 - 1)
#' @param FIRE Proportion of biomass consumed by fire during the non-growing season (0-1)
#' @export


calc_PDSOCt = function(BNPPt_est, Sf, Lo, LIGCELL, FIRE) {

  PDSOCt = 0.45*((LIGCELL*(Sf-Lo/2)*(1-FIRE)) + (LIGCELL+0.05)*BNPPt_est )
  return(PDSOCt)

}
