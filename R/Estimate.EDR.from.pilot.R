#' Title Estimate EDR based on p values 
#'
#' @param p.values resulting p values from DMR analysis 
#' @param model results from DMR analysis
#' @param thresh.p threshold for CBUM mixture model fitting, default is set to 0.005
#' @param FDR False discovery rate, by default is set to 0.05
#' @param M Number of parametric bootstrap, default is 10
#' @param N0 pilot sample size per group
#' @param target.N target sample size per group
#' @param method by defualt "CBUM+CDD"
#' @param transform.null equals False, we do not transform test statistics for non-DMR genes
#'
#' @return a list of estimation result
#' @export Estimate.EDR.from.pilot
#'
#' @examples
Estimate.EDR.from.pilot <- function(p.values, model,thresh.p=0.005,
                                    FDR=0.05, M=10,
                                    N0, target.N,
                                    method="CBUM+CDD",transform.null = F){
  
  q.value = p.adjust(p.values, method = "BH") ##?
  
  Fitted.Model = switch(method, 
                        "CBUM+CDD" = MixtureModel.Fittting.pi02(p.values,thresh.p=thresh.p,restrict = F)
                        )
  
  p.values.mod <- p.values
  if (any(is.na(p.values.mod))) {
    p.values.mod <- p.values.mod[!is.na(p.values.mod)]
    model = model[!is.na(p.values.mod),]
  }
  if (min(p.values.mod) == 0) {
    min.nonzero <- min(p.values.mod[p.values.mod > 0])
    p.values.mod[p.values.mod == 0] <- min.nonzero/2
  }
  if (max(p.values.mod) == 1) {
    max.nonone <- max(p.values.mod[p.values.mod < 1])
    p.values.mod[p.values.mod == 1] <- 1-(1-max.nonone)/2
  }
  
  ngenes = length(p.values.mod)
  
  if(Fitted.Model[1]=="error"){
    return("error")
  }else{
    lambda = as.numeric(Fitted.Model[1])
    r=as.numeric(Fitted.Model[2])
    s=as.numeric(Fitted.Model[3])
    if(method == "TwoBeta"){
      r2=as.numeric(Fitted.Model[4])
      s2=as.numeric(Fitted.Model[5])
      posterior = (lambda*dbeta(p.values.mod,r2,s2))/(dbeta(p.values.mod,r,s)*(1-lambda)+lambda*dbeta(p.values.mod,r2,s2)) #prob be non-DE
    } else posterior = lambda/(dbeta(p.values.mod,r,s)*(1-lambda)+lambda) #prob be non-DE
    
    if(lambda>0.99|lambda<0.5){
      #return(Fitted.Model)
      print(paste("lambda >0.95 or <0.5,lambda=",Fitted.Model[1]))
    }
    else{
      parameter = list(n.old = N0, n.new= target.N)
      Resampling <- function(target.N){
        # DE_status_bootstrap = sample(c(FALSE,TRUE),ngenes,prob = c(lambda,1-lambda),replace=TRUE)
        DE_status_posterior = sapply(1:length(posterior),function(x) sample(c(FALSE,TRUE),1,prob = c(posterior[x],1-posterior[x]),replace=TRUE))
        
        transform.p.values.norm <- function(p.values.each, DE_status_each, parameter, #model, 
                                           transform.null = F){
          n.old = parameter[[1]]
          n.new = parameter[[2]]
          
          if(DE_status_each){
            #beta.estimate = model[1:2]
            #delta = model[3]
            statistic.old = qnorm(p.values.each/2,lower.tail=FALSE)
            old<-n.old
            new<-n.new
            
            statistic.new = statistic.old*sqrt(new/old)
            
            p.values.star.star = (1-pnorm(abs(statistic.new)))*2
            return(p.values.star.star)
          } else if(transform.null) return(rep(runif(1),length(n.new)))
          else return(rep(p.values.each, length(n.new)))
          # } else{return(rep(p.values.each,length(n.new)))}
        }
        
        

        # p.values.star.posterior.old = sapply(1:ngenes,function(x) transform.p.values.old(p.values.mod[x], DE_status_posterior[x], N0,
        #                                                                            sample.size, transform.null = transform.null))
        p.values.star.posterior = sapply(1:ngenes,function(x) 
                                                         transform.p.values.norm(p.values.mod[x], 
                                                         DE_status_posterior[x], 
                                                         parameter = parameter,
                                                         transform.null = transform.null))
        
        
        
        p.values.star.posterior = matrix(p.values.star.posterior, ncol = length(target.N), byrow = T)
        
        Estimate_Posterior <- function(p.values.star.star){
          if (min(p.values.star.star,na.rm=T) == 0) {
            min.nonzero <- min(p.values.star.star[p.values.star.star > 0])
            p.values.star.star[p.values.star.star == 0] <- min.nonzero/2
          }
          if (max(p.values.star.star,na.rm=T) == 1) {
            max.non1 <- max(p.values.star.star[p.values.star.star <1])
            p.values.star.star[p.values.star.star == 1] <- (max.non1+1)/2
          }
          ## empirical FDR control
          p.DE=p.values.star.star[DE_status_posterior]
          p.nonDE=p.values.star.star[!DE_status_posterior]
          p.unique=sort(unique(p.values.star.star))
          
          FDR.unique <- vector(length=length(p.unique))
          for(i in 1:length(p.unique)){
            FDR.unique[i]=sum(p.nonDE<=p.unique[i])/sum(p.values.star.star<=p.unique[i])
          }
          index = which(FDR.unique<=FDR)
          p.values.cut = if(length(index)){
            p.unique[max(index)]
          } else NA #everything is beyond FDR 0.05
          if(min(p.values.star.star)>p.values.cut | is.na(p.values.cut)){Declare_status=rep("nonDE",ngenes)}
          else{
            Declare_status = DE_status_posterior
            Declare_status[which(p.values.star.star<=p.values.cut)] = "DE"
            Declare_status[-which(p.values.star.star<=p.values.cut)] = "nonDE"
          }
          A = sum((Declare_status=="nonDE")*(!DE_status_posterior))
          B = sum((Declare_status=="nonDE")*(DE_status_posterior))
          C = sum((Declare_status=="DE")*(!DE_status_posterior))
          D = sum((Declare_status=="DE")*(DE_status_posterior))
          if((C+D)==0){
            TP_hat_post=0
            ## no declared genes
          }
          else{
            TP_hat_post = D/(C+D)
          }
          if((A+B)==0){
            TN_hat_post=0
            ## no declared genes
          }
          else{
            TN_hat_post = A/(A+B)
          }
          # if(mean(p.values)>0.5+1.96*sqrt(1/(12*length(p.values)))){EDR_post=D/length(p.values)}
          # else{EDR_post = D/(B+D)}
          EDR_post = D/(B+D)
          Declare_post  = sum(Declare_status=="DE")
          # return(c(TP_hat_post,TN_hat_post,EDR_post,Declare_post, B = B, D = D, p.values.cut = p.values.cut))
          return(c(TP_hat_post, TN_hat_post, EDR_post, Declare_post, A = A,B = B, C = C, D = D,
                   p.values.cut = p.values.cut, FDR = FDR.unique[max(index)]))
        }
        
        Estimate.Posterior.Result = matrix(apply(p.values.star.posterior,2,Estimate_Posterior), ncol = length(target.N)); round(Estimate.Posterior.Result, 3)
        
        row.names(Estimate.Posterior.Result) = c("TP", "TN", "EDR", "DeclareDE", "A", "B", "C", "D", "FDRcut", "FDR")
        colnames(Estimate.Posterior.Result) = paste(target.N[1], "x" ,"2")
        return(Estimate.Posterior.Result)
      }
      Result = lapply(1:M,function(x){
        #print(x)
        Resampling(target.N)
      })
      
      ave.result<-Reduce("+",Result)/length(Result)
      pred<-list(EDR=ave.result["EDR",],DeclareDE=ave.result["DeclareDE",],FDR=ave.result["FDR",])
      return(pred) 
    }
  }
}



# Use censored beta-uniform mixer model to estimate pi0(lambda)
MixtureModel.Fittting.pi02 <- function(p.values, restrict = TRUE, l.upper = 0.95,thresh.p=0.05){
  library(limma)
  if (all(is.na(p.values))) {
    stop("all p-values were NA; nothing to compute")
  }
  orig.pvals <- p.values
  if (any(is.na(p.values))) {
    p.values <- p.values[!is.na(pvals)]
  }
  if (min(p.values) == 0) {
    min.nonzero <- min(p.values[p.values > 0])
    p.values[p.values == 0] <- min.nonzero/2
  }
  
  fBUM <- function(z,d) {
    p.values = d
    r = z[1]
    s = z[2]
    Target.Func = -sum(log(lambda+(1-lambda)*p.values^(r-1)*(1-p.values)^(s-1)/beta(r,s)))
    return(Target.Func)
  }
  
  ## inital value(MME and non-parametric)
  mean.p = mean(p.values)
  var.p = var(p.values)
  init.r = ((1-mean.p)*mean.p^2-mean.p*var.p)/var.p
  init.s= ((1-mean.p)^2*mean.p-(1-mean.p)*var.p)/var.p
  init = c(max(0,min(init.r,.9)),max(1,init.s))
  
  if(restrict){
    lambda = max(min(CBUM(p.values,thresh.censor=thresh.p),l.upper), 0.8)
    r.upper = .9
    s.upper = Inf
    r.lower = 0
    s.lower = 1
  } else{
    lambda = min(CBUM(p.values,thresh.censor=thresh.p),0.99)
    r.upper = 1
    s.upper = Inf
    r.lower = 0
    s.lower = 1
  }
  if(unique(tryCatch(print(c(optim(init, fBUM,d=p.values,method= "L-BFGS-B",upper=c(r.upper,s.upper),lower=c(r.lower,s.lower))$par)),error = function(e) {print("error")})=="error")){return(init)}
  else{
    pars = optim(init, fBUM,d=p.values,method= "L-BFGS-B",upper=c(r.upper,s.upper),lower=c(r.lower,s.lower))$par
    LL=sum(log(lambda+(1-lambda)*p.values^(pars[1]-1)*(1-p.values)^(pars[2]-1)/beta(pars[1],pars[2])))
    return(c(lambda,pars,LL))
  }
}



# Get dmr when controling true FDR at 0.05 level
get.dm.regions<-function(p,level=0.05,dmr){
  p.de<-p[dmr]
  p.ordered<-p[order(p)]
  tdr<-rep(0,length(p))
  for(i in 1:length(p)){
    tdr[i]=sum(p.de<=p.ordered[i])/sum(p.ordered<=p.ordered[i])
  }
  id<-which(tdr>=1-level)
  if(length(id)==0){
    return(print("0 DMR found"))
  }else num.de<-max(id)
  claim.dmr<-order(p)[1:num.de]
  return(claim.dmr)
}



