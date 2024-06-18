#' @title Maximum Above-ground Net Primary Production
#'
#' @description This function calculates the max ANPP  given mean annual precipitation and temperature.
#' @param RAIN Long-term Mean Annual Precipitation (mm/year).
#' @param MAT Long-term Mean Annual Precipitation (degrees Celsius).
#' @param SAND Sand fraction in the top 30 cm of soil (0 - 1)
#' @export


calc_ANPPmax = function(RAIN, MAT, SAND) {

  e_suffix = 12.039+0.718*log(RAIN)-(25.18/(0.00834*(273.15+MAT)))
  ANPPmax = exp(e_suffix)*(1.33-0.0075*(SAND*100))
  return(ANPPmax)

}
