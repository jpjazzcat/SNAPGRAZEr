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
#' @param orig Default = FALSE. Use the original DMRESP equations from Ritchie 2020 or the updated ones from Ruan deWet
#' @export

calc_deltaSOC = function(PDSOCt, DDSOCt, SAND, RAIN, Gdays, SOC, lowSOC = FALSE, orig = FALSE) {

  WETDAYS = (0.00044*RAIN-0.025)*Gdays

  if(orig) {

    # Ritchie 2020 original equations

    if(lowSOC){

      ## SOC < 4600
      DMRESP = exp(-10.872+(1.296*log(SOC)))

    } else {

      ## SOC > 4600
      DMRESP = -0.579 + 0.0004 * SOC


    }

  } else {

    # Ruan deWet updated equations

    if(lowSOC){

      ## SOC < 4600
      DMRESP = exp(-10.872+(1.296*log(SOC)))

    } else {

      ## SOC > 4600
      DMRESP = -0.579 + 0.00036 * SOC


  }


  MRESPt = (WETDAYS*(0.7+(0.3*SAND/100)))*DMRESP
  deltaSOC = PDSOCt+DDSOCt-MRESPt

  return(deltaSOC)

  }
}

