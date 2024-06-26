#' @title Delta SOC
#'
#' @description Calculate the change in SOC for year t (deltaSOCt) in tSOC/ha using plant-derived soil organic carbon (PDSOCt), dung-derived soil organic carbon (DDSOCt), and microbial respiration (MRESPt). A key input to MRESPt is WETDAYS, which is calculated as part of this function.
#' @param PDSOCt Output of calc_PDSOCt()
#' @param DDSOCt Output of calc_DDSOCt()
#' @param SAND Sand fraction in the top 30 cm of soil (0 - 1)
#' @param RAIN MAP for year t (mm/year)
#' @param Gdays Total number of days in the growing season.
#' @param SOC Starting soil organic carbon stocks (tSOC/ha).
#' @param orig Default = FALSE. Use the original DMRESP equations from Ritchie 2020 or the updated ones from Ruan deWet.
#' @export

calc_deltaSOC = function(PDSOCt, DDSOCt, SAND, RAIN, Gdays, SOC, orig = FALSE) {

  WETDAYS = (0.00044*RAIN-0.025)*Gdays

  if(orig) {

    # Ritchie 2020 original equations

    if(SOC < 46){

      DMRESP = exp(-10.872+(1.296*log(SOC*100)))

      # Multiply SOC*100 for all DMRESP equations to convert units to g/m2

    } else {

      DMRESP = -0.579 + 0.0004 * (SOC*100)


    }

  } else {

    # Ruan deWet updated equations

    if(SOC < 46){

      DMRESP = exp(-10.872+(1.296*log(SOC*100)))

    } else {

      DMRESP = -0.579 + 0.00036 * (SOC*100)


  }


  MRESPt = (WETDAYS*(0.7+(0.3*SAND)))*DMRESP
  deltaSOC = (PDSOCt+DDSOCt-MRESPt)/100

  return(deltaSOC) # Units = tSOC/ha for single year

  }
}

