#' @title Biomass consumed during grazing episode
#'
#' @description Calculation of the parameter Lg, the biomass loss during grazing episode, from Ritchie (2020)
#' @param Ddays Number of days of grazing episode
#' @param d Stocking density (head/ha)
#' @param n Number of "pastures" per total area, A.
#' @param W Average animal body size (kg live weight)
#' @param Cg Daily consumption rate (g/animal/day)
#' @export

calc_Lg = function(Ddays, d, n, W, Cg = NA) {

  if(is.numeric(Cg)) {

    Cg = Cg

  } else if(is.numeric(W)) {

    Cg = calc_Cg(W)

  } else {print("You need to either provide Cg (daily consumption) or W (animal body size)")}

  Lg = Ddays*d*Cg*n*10^(-4)
  return(Lg)

}
