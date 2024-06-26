#' @title Equilibrium SOC
#'
#' @description Calculate the SOC stock in tSOC/ha at equilibrium under the specified conditions. This is the point at which SOC inputs from plant and dung equal SOC outputs from soil respiration (deltaSOC = 0).
#' @param PDSOCt Output of calc_PDSOCt() (g/m2)
#' @param DDSOCt Output of calc_DDSOCt() (g/m2)
#' @param DEPTH Soil depth at which equilibrium SOC stock is reported (cm)
#' @param SAND Sand fraction in the top 30 cm of soil (0 - 1)
#' @param RAIN Mean annual precipitation for year t (mm/year)
#' @param Gdays Total number of days in the growing season.
#' @param orig Default = FALSE. Use the original DMRESP equations from Ritchie 2020 or the updated ones from Ruan deWet
#' @export

calc_SOCeq = function(PDSOCt, DDSOCt, DEPTH, SAND, RAIN, Gdays, orig = FALSE) {

  WETDAYS = (0.00044*RAIN-0.025)*Gdays

  if(orig) {

    SOCeq = min(((PDSOCt+DDSOCt)/(WETDAYS*(0.7+0.3*SAND)*exp(-10.872)))^(1/1.296),
                ((((PDSOCt+DDSOCt)/(WETDAYS*(0.7+(0.3*SAND))))+0.579)/0.00044))

    SOCeq = SOCeq*(-0.35+0.37*log(DEPTH))/100

  } else {

    SOCeq = min(((PDSOCt+DDSOCt)/(WETDAYS*(0.7+0.3*SAND)*exp(-10.872)))^(1/1.296),
                ((((PDSOCt+DDSOCt)/(WETDAYS*(0.7+(0.3*SAND))))+0.579)/0.00036))

    SOCeq = SOCeq*(-0.35+0.37*log(DEPTH))/100


    # Minimum function reflects piecewise equation to calculate MRESP. See 'calc_deltaSOC' for more.
    # Depth function adjusts SOC stock to specified depth. Original SNAP model was calibrated to 40 cm. Not sure where coefficients come from, but they exist in older versions of code. Worth confirming.
    # Divide by 100 to convert units from g/m2 to t/ha

  }

  return(SOCeq)

}
