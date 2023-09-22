#' @title Delta SOC
#'
#' @description To calculate the change in SOC for year t (deltaSOCt), we need the combination of PDSOCt and DDSOCt, but also the maximum rate of microbial respiration for year t (MRESPt). A key input to MRESPt is WETDAYS, which is calculated as part of this function.
#' @param PDSOCt Output of calc_PDSOCt()
#' @param DDSOCt Output of calc_DDSOCt()
#' @param SAND Sand % in top 30 cm soil
#' @param RAIN MAP for year t (mm/year)
#' @param Gdays Total number of days in the growing season. Default = 153 (October to March-ish).
#' @param SOC Starting soil carbon stocks (g/m2).
#' @param lowSOC Default = FALSE. Different regression equation for respiration rate is applied for low and high SOC to avoid a negative respiration rate (which isn't physically possible). Threshold for what qualifies as "low SOC" is 4,600 gC/m^2 (i.e. 46 t/ha). Low SOC regression equation is applicable for higher SOC, but just with slightly lower R-squared.
#' @export deltaSOC Change in soil carbon stocks over year t (g/m2)

calc_deltaSOC = function(PDSOCt, DDSOCt, SAND, RAIN, Gdays, SOC, lowSOC = FALSE) {

  WETDAYS = (0.00044*RAIN-0.025)*Gdays

  MRESP = WETDAYS*(0.7+(0.3*SAND/100))*(0.00044*SOC-0.579)

  if(lowSOC) {

    MRESPt = (WETDAYS*(0.7+(0.3*SAND/100)))*(exp(-10.872)*SOC^1.296)

  } else {

    MRESPt = (WETDAYS*(0.7+(0.3*SAND/100)))*(0.00044*SOC-0.579)

  }

  deltaSOC = PDSOCt+DDSOCt-MRESPt

  return(deltaSOC)

}

