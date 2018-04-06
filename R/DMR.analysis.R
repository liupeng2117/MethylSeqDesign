#' Differential methylation analysis on pilot data
#'
#' @param x a vector containing study design information, the number of each element is 
#' the sample size for each group, the length of x is the number of groups. E.g., if 
#' there are two groups with equal sample size 5, x=c(5,5)
#' @param cov.matrix a matrix of coverage, the rows are genes and columns are samples
#' @param methyl.matrix a matrix of methylated read count, the rows are genes and columns are samples. A one to one correspondence between the matrix of methylated read count and the matrix of coverage is required.
#' @return a list of two, one is the p values of each gene, the second element is a vector of 
#' model parameters 
#' @export DMR.analysis
#'
#' @examples
DMR.analysis<-function(x,cov.matrix,methyl.matrix){
  fai.est<-estimate.fai2(x=x, cov.matrix, methyl.matrix)
  test.result<-Z.wald(x, cov.matrix, methyl.matrix, fai=fai.est,a=1)
  p.values<-test.result[[1]]
  model<-cbind(beta0=test.result[[1]],beta1=test.result[[2]],fai=fai.est)
  return(list(p.values,model))
}


# Z-value + wald-test
Z.wald<-function(x,cov.matrix,sim.meth,fai,a=1){
  G<-dim(cov.matrix)[1]
  D<-sum(x)
  n.groups<-length(x)
  # Assign each group size
  
  for(i in 1:n.groups){
    assign(paste0("D",i),x[i])
  }
  
  ### calculate beta1-hat
  chi2<-rep(0,G)
  
  for(i in 1:G){
    tryCatch({
      X<-matrix(c(rep(1,D),rep(1,D1),rep(0,D-D1)),ncol=2)
      
      if(n.groups>2){
        c<-D1
        for(j in 2:(n.groups-1)){
          col<-c(rep(0,c),rep(1,get(paste0("D",j))),rep(0,D-c-get(paste0("D",j))))
          X<-cbind(X,col)
          c=c+get(paste0("D",j))
        }        
      }
      
      V.inverse<-diag(cov.matrix[i,]/(1+(cov.matrix[i,]-1)*fai))
      Z<-asin(2*(sim.meth[i,]+a)/(cov.matrix[i,]+2*a)-1)
      sigma<-solve(t(X)%*%V.inverse%*%X)
      beta.hat<-sigma %*% t(X) %*% V.inverse %*% Z
      
      C<-diag(rep(1,n.groups))[,-1]
      var<-t(C) %*% sigma %*% C
      chi2[i]<-t(t(C) %*% beta.hat) %*% solve(var) %*% (t(C) %*% beta.hat)
    }, error = function(err) {
      # error handler picks up where error was generated
      print(paste("MY_ERROR:  ",err))
    })
  }
  # Calculate test statistic, p and q values
  p.Zw<-pchisq(chi2,df=n.groups-1,lower.tail=F)
  return(list(p.Zw=p.Zw,z.statistic=chi2))
}



## Estimate fai using park's method
estimate.fai2<-function(cov.matrix,methyl.matrix,a=1,x){
  G<-dim(cov.matrix)[1]
  D<-dim(cov.matrix)[2]
  n.groups<-length(x)
  # Assign each group size
  for(i in 1:n.groups){
    assign(paste0("D",i),x[i])
  }
  
  if(D==2) fai.est=0.001 else{
    fai<-rep(0,G)
    for (i in 1:G){
      tryCatch({
        V0.inverse<-diag(cov.matrix[i,])
        # X matrix
        X<-matrix(c(rep(1,D),rep(1,D1),rep(0,D-D1)),ncol=2)
        if(n.groups>2){
          c<-D1
          for(j in 2:(n.groups-1)){
            col<-c(rep(0,c),rep(1,get(paste0("D",j))),rep(0,D-c-get(paste0("D",j))))
            X<-cbind(X,col)
            c=c+get(paste0("D",j))
          }
        }
        Z<-asin(2*(methyl.matrix[i,]+a)/(cov.matrix[i,]+2*a)-1)
        sigma<-solve(t(X)%*%V0.inverse%*%X)
        beta.hat0<-sigma%*%t(X)%*%V0.inverse%*%Z
        chi2<-sum(cov.matrix[i,]*(Z-X %*% beta.hat0)^2)
        sigma2.hat<-chi2/(D-n.groups)
        fai[i]<-(D*(sigma2.hat-1))/sum(cov.matrix[i,]-1)
      }, error = function(err) {
        # error handler picks up where error was generated
        print(paste("MY_ERROR:  ",err))
      })
    }
    fai<-fai[complete.cases(fai)]
    fai<-ifelse(fai<0,10^(-6),fai)
    fai<-ifelse(fai>1,1-10^(-6),fai)
    fai.est<-mean(fai)
  }
  return(fai.est)
}