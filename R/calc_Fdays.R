#' @title Number of growing season days after the grazing event.
#'
#' @description Number of growing season days after the grazing event.
#' @param Gdays Total number of days in the growing season.
#' @param Edays Number of growing season days before the grazing event.
#' @param Ddays Number of days of grazing event.
#' @export

calc_Fdays = function(Gdays, Edays, Ddays) {

  Fdays = (Gdays - Edays - Ddays)

  if(Fdays < 0){
    stop("ERROR: Fdays must be greater than zero")
  } else{
    return(Fdays)
  }
}
