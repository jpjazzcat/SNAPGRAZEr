#' @title Length of growing season
#'
#' @description Calculation of the length of the growing season in days from Ritchie (2020)
#' @param MAT Mean annual temperature (degrees Celsius)
#' @param RAIN Mean annual precipitation (mm/yr)
#' @export

calc_Gdays = function(RAIN, MAT) {
  Gdays = 22.99*MAT-0.94*MAT^2+0.073*RAIN
  return(Gdays)
}
