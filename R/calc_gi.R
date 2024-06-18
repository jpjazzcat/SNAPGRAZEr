#' @title Grazing intensity
#'
#' @description Calculation of the parameter gi, the grazing intensity of the specified grazing management practices
#' @param Sf Biomass at the end of the growing season. Output of calc_Sf.
#' @param Lo Biomass consumed during the offseason. Output of calc_Lo.
#' @param Sk The steady state of biomass in the absence of grazing for a given location. Either measured directly or calculated as output of calc_ANPPmax.
#' @param Lg Biomass consumed during the growing season. Output of calc_LG.
#' @param ANPPest Aboveground net primary productivity during the growing season. Output of calc_ANPPest.
#' @param method Method used to calculate grazing intensity. Options are 'Sk' or 'ANPPest'. 'Sk' calculates grazing intensity as the proportional difference in plant biomass between grazed and ungrazed treatments at the end of the offseason. 'ANPPest' calculates grazing intensity as the amount of biomass consumed during both the growing and offseason relative to the amount of biomass produced.
#' @export

calc_gi = function(Sf, Lo, Sk, Lg, ANPPest, method) {
  if(method == "Sk"){
    gi = 1-((Sf-Lo)/Sk) #stop("ERROR: 'method' must be either 'Sk' or 'ANPPest'.")
  }
  else if(method == "ANPPest"){
    gi = (Lg+Lo)/ANPPest
  }
  else{
    stop("ERROR: 'method' must be either 'Sk' or 'ANPPest'.")
  }
  return(gi)
}
