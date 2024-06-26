#' @title SNAPGRAZE delta Wrapper Function
#'
#' @description This is a wrapper function around all the SNAPGRAZE equations to calculate SOC delta over a given number of years. Management and weather inputs can be provided as either single values or as a vector/list of values representing the corresponding value for each year being modeled.
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
#' @param d_off Stocking density during the offseason (head/ha)
#' @param n Number of "pastures" per total area, A.
#' @param W Average animal body size (kg live weight)
#' @param Cg Daily biomass consumption rate during the growing season (g/animal/day). Default = calc_Cg(W)
#' @param Co Daily biomass consumption rate during the offseason (g/animal/day). Default = Cg/2
#' @param r Maximum relative growth rate of grass biomass. Can be calculated as the y-intercept of a relationship between the measured relative growth rate and biomass at a given time. Default for tropical grasslands is 0.05 and for temperate grasslands is 0.035.
#' @param APCcorrection Default = FALSE. A correction factor for the influence of annual versus perennial plant growth strategies on belowground production. If correction is applied, then APC = 0.291. This is where forage is dominated by annuals (or shrubs often associated with annuals in drier grasslands). Otherwise, APC = 1.
#' @param DEPTH Default = 30. Depth of soil sampling / estimation (cm). The original SNAP model was developed based on measurements to a depth of 40 cm.
#' @param orig Default = FALSE. Use the original DMRESP equations from Ritchie 2020 or the updated ones from Ruan deWet.
#' @export

SNAPGRAZE_delta_ann = function(SAND, RAIN, MAT, FIRE, LIGCELL, years, SOC,
                           Sk = NA, S0 = NA, Edays, Ddays, Gdays = NA, DdaysO = NA, d, d_off,
                           n, W, Cg = NA, Co = NA, r, APCcorrection = FALSE, DEPTH = 30, orig = FALSE) {


  soc_list = vector("list", (years+1))
  soc_list[[1]] <- (SOC)

  test_input <- function(x){
    if(length(unlist(x)) > 1){
      x <- unlist(x)[[i]]
    }else{
      x
    }
  }

  for(i in 1:years){

    RAIN_i <- test_input(RAIN)
    MAT_i <- test_input(MAT)
    FIRE_i <- test_input(FIRE)
    LIGCELL_i <- test_input(LIGCELL)
    Edays_i <- test_input(Edays)
    Ddays_i <- test_input(Ddays)
    DdaysO_i <- test_input(DdaysO)
    d_i <- test_input(d)
    d_off_i <- test_input(d_off)
    n_i <- test_input(n)
    APCcorrection_i <- test_input(APCcorrection)

    if(SAND<0|SAND>1){
      stop("ERROR: SAND must be between 0 and 1")
    } else if(RAIN_i<0|RAIN_i>10000){
      stop("ERROR: RAIN must be between 0 and 10000")
    } else if(MAT_i < -20|MAT_i>40){
      stop("ERROR: MAT must be between -20 and 40")
    } else if(FIRE_i<0|FIRE_i>1){
      stop("ERROR: FIRE must be between 0 and 1")
    } else if(LIGCELL_i<0|LIGCELL_i>1){
      stop("ERROR: LIGCELL must be between 0 and 1")
    } else if(Edays_i<0|Edays_i>365){
      stop("ERROR: Edays must be between 0 and 365")
    } else if(Ddays_i<0|Ddays_i>365){
      stop("ERROR: Ddays must be between 0 and 365")
    } else if(d_i<0|d_i>100){
      stop("ERROR: d must be between 0 and 100")
    } else if(d_off_i<0|d_off_i>100){
      stop("ERROR: d_off must be between 0 and 100")
    } else if(W<0|W>5000){
      stop("ERROR: W must be between 0 and 5000")
    } else if(r<=0|r>0.2){
      stop("ERROR: r must be between 0 and 0.2")
    } else if(APCcorrection_i<0|APCcorrection_i>1){
      stop("ERROR: APCcorrection must be either TRUE, FALSE, or between 0 and 1")
    } else if(DEPTH<0|DEPTH>200){
      stop("ERROR: DEPTH must be between 0 and 200")
    }


    if(is.na(Sk)){
      Sk_i = calc_ANPPmax(RAIN_i, MAT_i, SAND)/0.9
    }else{
      Sk_i = Sk
    }

    if(is.na(S0)){
      S0_i = Sk_i*0.1
    }else{
      S0_i = S0
    }

    if(is.na(Gdays)){
      Gdays_i = calc_Gdays(RAIN_i, MAT_i)
    }else{
      Gdays_i = Gdays
    }

    if(Edays_i + Ddays_i > Gdays_i){
      stop("ERROR: Edays + Ddays must be <= Gdays")
    } else if (Gdays_i > 365){
      stop("ERROR: Gdays must be <= 365")
    }

    if(is.na(DdaysO)){
      DdaysO_i = (365 - Gdays_i)
    }else{
      DdaysO_i = DdaysO
    }

    if(DdaysO_i < 0|DdaysO_i > 365){
      stop("ERROR: DdaysO must be between 0 and 365")
    } else if(DdaysO_i + Gdays_i > 365){
      stop("ERROR: DdaysO + Gdays must be <= 365")
    }

    if(is.na(Cg)) {
      Cg = (5300+770*log(W))
    }

    if(is.na(Co)) {
      Co = Cg/2
    }


    # Episodic Herbivory Model (EHM)
    Se = calc_SE(Sk = Sk_i, Edays = Edays_i, S0 = S0_i, r)
    Lg = calc_Lg(Ddays = Ddays_i, d = d_i, n = n_i, W, Cg)
    Sg = calc_Sg(Sk = Sk_i, Se, Lg, Ddays = Ddays_i, n = n_i, d = d_i, r, W, Cg)
    Fdays = calc_Fdays(Gdays = Gdays_i, Edays = Edays_i, Ddays = Ddays_i)
    Sf = calc_Sf(Sk = Sk_i, Sg, r, Fdays)
    Pg = calc_Pg(Se, Sg, Sf, Sk = Sk_i, S0 = S0_i)
    Lo = calc_Lo(Cg, Co, Gdays = Gdays_i, DdaysO = DdaysO_i, d_off = d_off_i)

    # dmax = calc_dmax(Sf, Sk, Cg, Gdays)

    # if(dmax > d){
    #   print("Looks like stocking density is greater than what is theoretically sustainable.")
    #   } else {print("Stocking density A-O-K!")}

    # Productivity
    ANPPt_max = calc_ANPPmax(RAIN = RAIN_i, MAT = MAT_i, SAND)
    ANPPt_est = calc_ANPPest(Se, Sg, Sf, Sk = Sk_i, S0 = S0_i)
    BNPPt_est = calc_BNPPest(RAIN = RAIN_i, MAT = MAT_i, ANPPt_est, Sk = Sk_i, S0 = S0_i, APCcorrection = APCcorrection_i, DEPTH)

    # SOC
    PDSOCt = calc_PDSOCt(BNPPt_est, Sf, Lo, LIGCELL = LIGCELL_i, FIRE = FIRE_i)
    DDSOCt = calc_DDSOCt(LIGCELL = LIGCELL_i, Ddays = Ddays_i, Cg, n = n_i, d = d_i, Lo)

    x <- soc_list[[i]]
    deltaSOC = calc_deltaSOC(PDSOCt, DDSOCt, SAND, RAIN = RAIN_i, Gdays = Gdays_i, SOC = x, orig=FALSE)
    # SOC stock at the end of year i
    SOCi_end =  x+deltaSOC
    soc_list[[i+1]] <- SOCi_end

  }

  return(soc_list)

}

