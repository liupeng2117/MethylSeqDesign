#' Title Obtain predicted EDR
#'
#' @param result input result from the function Estimate.EDR.from.pilot 
#'
#' @return a list of results, including EDR, FDR etc
#' @export obtain.EDR
#'
#' @examples 
#' 
#' EDR<-obtain.EDR(result=power.result)
#' 
#' 
obtain.EDR<-function(result){
  y<-result[[1]] #return edr
  data.ref<-Reduce("+",y)/length(y)
  return(data.ref)
}
