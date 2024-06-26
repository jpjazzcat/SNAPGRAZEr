#' @title SNAPGRAZE delta Wrapper Function
#'
#' @description This is a wrapper function around all the SNAPGRAZE equations to calculate SOC delta over a given number of years.
#' @param SAND Sand fraction in the top 30 cm of soil (0 - 1)
#' @param RAIN Long-term Mean Annual Precipitation (mm/year).
#' @param MAT Long-term Mean Annual Precipitation (degrees Celsius).
#' @param FIRE Average number of fires per year (#/year)
#' @param LIGCELL Lignin and cellulose fraction of plant biomass (0 - 1)
#' @param years Number of years over which to run the simulation.
#' @param SOC Starting soil carbon stocks (t/ha).
#' @param Sk The steady state of biomass in the absence of grazing for a given location. This should ideally be measured directly using grazing exclosures.
#' @param S0 Biomass condition prior to the growing season (at the end of the dry season) that is mostly comprised of carbon stores in rhizomes. Default is 0.1*SK.
#' @param Edays Number of days within the growing season prior to grazing episode
#' @param Ddays Number of days of grazing episode
#' @param Gdays Total number of days in the growing season. Default = 153 (October to March-ish).
#' @param DdaysO Number of days of offseason grazing episode. Default = 365 - Gdays.
#' @param d Stocking density during the growing season (head/ha)
#' @param d_off stocking density during the offseason (head/ha)
#' @param n Number of "pastures" per total area, A.
#' @param W Average animal body size (kg live weight)
#' @param Cg Daily biomass consumption rate during the growing season (g/animal/day). Default = calc_Cg(W)
#' @param Co Daily biomass consumption rate during the offseason (g/animal/day). Default = Cg/2
#' @param r Maximum relative growth rate of grass biomass. Can be calculated as the y-intercept of a relationship between the measured relative growth rate and biomass at a given time. Default for tropical grasslands is 0.05 and for temperate grasslands is 0.035.
#' @param APCcorrection Default = FALSE. A correction factor for the influence of annual versus perennial plant growth strategies on belowground production. If correction is applied, then APC = 0.291. This is where forage is dominated by annuals (or shrubs often associated with annuals in drier grasslands). Otherwise, APC = 1.
#' @param DEPTH Default = 30. Depth of soil sampling / estimation (cm). The original SNAP model was developed based on measurements to a depth of 40 cm.
#' @param orig Default = FALSE. Use the original DMRESP equations from Ritchie 2020 or the updated ones from Ruan deWet.
#' @export

SNAPGRAZE_delta = function(SAND, RAIN, MAT, FIRE, LIGCELL, years, SOC,
                     Sk = NA, S0 = 0.1*Sk, Edays, Ddays, Gdays = NA, DdaysO = NA, d, d_off,
                     n, W, Cg = NA, Co = NA, r, APCcorrection = FALSE, DEPTH = 30, orig = FALSE) {

  if(SAND<0|SAND>1){
    stop("ERROR: SAND must be between 0 and 1")
  } else if(RAIN<0|RAIN>10000){
    stop("ERROR: RAIN must be between 0 and 10000")
  } else if(MAT < -20|MAT>40){
    stop("ERROR: MAT must be between -20 and 40")
  } else if(FIRE<0|FIRE>1){
    stop("ERROR: FIRE must be between 0 and 1")
  } else if(LIGCELL<0|LIGCELL>1){
    stop("ERROR: LIGCELL must be between 0 and 1")
  } else if(Edays<0|Edays>365){
    stop("ERROR: Edays must be between 0 and 365")
  } else if(Ddays<0|Ddays>365){
    stop("ERROR: Ddays must be between 0 and 365")
  } else if(d<0|d>100){
    stop("ERROR: d must be between 0 and 100")
  } else if(d_off<0|d_off>100){
    stop("ERROR: d_off must be between 0 and 100")
  } else if(W<0|W>5000){
    stop("ERROR: W must be between 0 and 5000")
  } else if(r<=0|r>0.2){
    stop("ERROR: r must be between 0 and 0.2")
  } else if(APCcorrection<0|APCcorrection>1){
    stop("ERROR: APCcorrection must be either TRUE, FALSE, or between 0 and 1")
  } else if(DEPTH<0|DEPTH>200){
    stop("ERROR: DEPTH must be between 0 and 200")
  }

  if(is.na(Sk)){
    Sk = calc_ANPPmax(RAIN, MAT,SAND)/0.9
  }

  if(is.na(Gdays)){
    Gdays = calc_Gdays(RAIN, MAT)
  }

  if(Edays + Ddays > Gdays){
    stop("ERROR: Edays + Ddays must be <= Gdays")
  } else if (Gdays > 365){
    stop("ERROR: Gdays must be <= 365")
  }

  if(is.na(DdaysO)){
    DdaysO = (365 - Gdays)
  }

  if(DdaysO < 0 | DdaysO > 365){
    stop("ERROR: DdaysO must be between 0 and 365")
  } else if(DdaysO + Gdays > 365){
    stop("ERROR: DdaysO + Gdays must be <= 365")
  }

  if(is.na(Cg)) {
    Cg = calc_Cg(W)
  }

  if(is.na(Co)) {
    Co = Cg/2
  }


  # Episodic Herbivory Model (EHM)
  Se = calc_SE(Sk, Edays, S0, r)
  Lg = calc_Lg(Ddays, d, n, W, Cg)
  Sg = calc_Sg(Sk, Se, Lg, Ddays, n, d, r, W, Cg)
  Fdays = calc_Fdays(Gdays, Edays, Ddays)
  Sf = calc_Sf(Sk, Sg, r, Fdays)
  Pg = calc_Pg(Se, Sg, Sf, Sk, S0)
  Lo = calc_Lo(Cg, Co, Gdays, DdaysO, d_off)

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

  soc_list = vector("list", years)
  soc_list[[1]] <- (SOC) # Convert SOC from t/ha to g/m2

  for(i in 1:years){
    x <- soc_list[[i]]
    deltaSOC = calc_deltaSOC(PDSOCt, DDSOCt, SAND, RAIN, Gdays, SOC = x, orig)
    # SOC stock at the end of year i
    SOCi =  x+deltaSOC
    i <- i+1
    soc_list[[i]] <- SOCi

  }

  return(soc_list)

  }

