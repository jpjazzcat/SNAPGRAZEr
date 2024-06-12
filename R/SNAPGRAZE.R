#' @title SNAPGRAZE Wrapper Function
#'
#' @description This is a wrapper function around all the SNAPGRAZE equations to calculate SOCeq based on the below input parameters.
#' @param SAND Sand % in top 30 cm soil
#' @param RAIN Long-term Mean Annual Precipitation (mm/year).
#' @param MAT Long-term Mean Annual Precipitation (degrees Celsius).
#' @param FIRE Average number of fires per year (#/year)
#' @param LIGCELL Lignin and cellulose content of livestock feed for year t (%)
#' @param Sk The steady state of biomass in the absence of grazing for a given location. This should ideally be measured directly using grazing exclosures.
#' @param S0 Biomass condition prior to the growing season (at the end of the dry season) that is mostly comprised of carbon stores in rhizomes. Default is 0.1*SK.
#' @param Edays Number of days within the growing season prior to grazing episode
#' @param Ddays Number of days of grazing episode
#' @param Fdays Number of days left in the growing season after the grazing episode. Fdays = Gdays - Edays - Ddays
#' @param Gdays Total number of days in the growing season. Default = 153 (October to March-ish).
#' @param d Stocking density (head/ha)
#' @param n Number of "pastures" per total area, A.
#' @param W Average animal body size (kg live weight)
#' @param Cg Daily consumption rate (g/animal/day)
#' @param r Maximum relative growth rate of grass biomass. Can be calculated as the y-intercept of the relationship between the measured relative growth rate and biomass at a given time. Default for tropical grasslands is 0.05 and for temperate grasslands is 0.035.
#' @param APCcorrection Default = FALSE. A correction factor for the influence of annual versus perennial plant growth strategies on belowground production. If correction is applied, then APC = 0.291. This is where forage is dominated by annuals (or shrubs often associated with annuals in drier grasslands). Otherwise, APC = 1.
#' @param DEPTH Default = 30. Depth of soil sampling / estimation (cm). The original SNAP model was developed based on measurements to a depth of 40 cm.
#' @param orig Default = FALSE. Use the original DMRESP equations from Ritchie 2020 or the updated ones from Ruan deWet
#' @export

SNAPGRAZE = function(SAND, RAIN, MAT, FIRE, LIGCELL,
                     Sk = NA, S0 = 0.1*Sk, Edays, Ddays, Fdays, Gdays = NA, d_off,
                     d, n, W, Cg = NA, r = 0.05, APCcorrection = FALSE, DEPTH = 30, orig = FALSE) {

  if(SAND<0|SAND>100){
    stop("ERROR: SAND must be between 0 and 100")
  }

  if(RAIN<0|RAIN>10000){
    stop("ERROR: RAIN must be between 0 and 10000")
  }

  if(MAT < -20|MAT>40){
    stop("ERROR: MAT must be between -20 and 40")
  }

  if(FIRE<0|FIRE>1){
    stop("ERROR: FIRE must be between 0 and 1")
  }

  if(LIGCELL<0|LIGCELL>100){
    stop("ERROR: LIGCELL must be between 0 and 100")
  }

  if(Gdays!=sum(Edays, Ddays, Fdays)){
    stop("ERROR: Gdays must equal Edays + Ddays + Fdays")
  }

  if(d<0|d>100){
    stop("ERROR: d must be between 0 and 100")
  }

  if(W<0|W>5000){
    stop("ERROR: W must be between 0 and 5000")
  }

  if(r<=0|r>0.2){
    stop("ERROR: r must be between 0 and 0.2")
  }

  if(APCcorrection<0|APCcorrection>1){
    stop("ERROR: APCcorrection must be either TRUE, FALSE, or between 0 and 1")
  }

  if(DEPTH<0|DEPTH>200){
    stop("ERROR: DEPTH must be between 0 and 200")
  }

  if(is.na(Gdays)){
    Gdays = 22.99*MAT-0.94*MAT^2+0.073*RAIN
  }

  if(is.na(Sk)){
    Sk = calc_ANPPmax(RAIN, MAT,SAND)/0.9
  }

  if(is.na(Cg)) {
    Cg = (5300+770*log(W))
  }

  # Episodic Herbivory Model (EHM)
  Se = calc_SE(Sk, Edays, S0, r)
  Lg = calc_Lg(Ddays, d, n, W, Cg)
  Sg = calc_Sg(Sk, Se, Lg, Ddays, n, d, r, W, Cg)
  Sf = calc_Sf(Sk, Sg, r, Fdays)
  Pg = calc_Pg(Se, Sg, Sf, Sk, S0)
  Lo = calc_Lo(Cg, Gdays, d_off)

  # dmax = calc_dmax(Sf, Sk, Cg, Gdays)

  # if(dmax > d){
  #   print("Looks like stocking density is greater than what is theoretically sustainable.")
  #   } else {print("Stocking density A-O-K!")}

  # Productivity
  ANPPt_max = calc_ANPPmax(RAIN, MAT, SAND)
  ANPPt_est = calc_ANPPest(Se, Sg, Sf, Sk, S0)
  BNPPt_est = calc_BNPPest(RAIN, MAT, ANPPt_est, Sk, S0, APCcorrection, DEPTH)

  # SOC
  PDSOCt = calc_PDSOCt(BNPPt_est, Sf, Lo, LIGCELL, FIRE)
  DDSOCt = calc_DDSOCt(LIGCELL, Ddays, Cg, n, d, Lo)

  SOCeq = calc_SOCeq(PDSOCt, DDSOCt, DEPTH, SAND, RAIN, Gdays, orig)

  return(SOCeq)

}

