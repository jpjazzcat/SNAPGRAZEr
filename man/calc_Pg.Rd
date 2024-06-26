% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_Pg.R
\name{calc_Pg}
\alias{calc_Pg}
\title{Total seasonal grass production}
\usage{
calc_Pg(Se, Sg, Sf, Sk, S0 = 0.1 * Sk)
}
\arguments{
\item{Se}{Biomass at the start of the grazing episode. Output of calc_Se.}

\item{Sg}{Biomass at the end of grazing episode. Output of calc_Sg.}

\item{Sf}{Biomass at the end of the growing season. Output of calc_Sf.}

\item{Sk}{The steady state of biomass in the absence of grazing for a given location. This should ideally be measured directly using grazing exclosures.}

\item{S0}{Biomass condition prior to the growing season (at the end of the dry season) that is mostly comprised of carbon stores in rhizomes. Default is 0.1*SK.}
}
\description{
Calculation of the parameter Pg, the total seasonal grass production for the given grazing management system, from Ritchie (2020)
}
