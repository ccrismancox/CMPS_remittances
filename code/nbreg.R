
NBgrad <- function(b, X, Y){
     ## Gradient function for negative binomial 
     
     eta <- b[length(b)]
     alpha<-exp(eta)
     a2 <- exp(2*eta)
     beta<-b[-length(b)]
     XB <- as.numeric(X %*% beta)
     e.XB<-exp(XB)
     e.2XB<-exp(2*XB)
     eXBA <- exp(XB+eta)
     
     e.XB[e.XB<.Machine$double.eps] <- .Machine$double.eps
     e.2XB[e.2XB<.Machine$double.eps] <- .Machine$double.eps
     eXBA[eXBA<.Machine$double.eps] <- .Machine$double.eps
     
     
     dB <- X* (Y-e.XB)/(eXBA+1)
     dA <-( exp(-eta)/(eXBA+1)) * (Y * alpha + 
                                        eXBA * log(eXBA+1) - 
                                        eXBA * digamma(Y + 1/alpha) +
                                        eXBA * digamma(1/alpha) -
                                        eXBA+
                                        log(eXBA + 1)-
                                        digamma(Y + 1/alpha) + 
                                        digamma(1/alpha))
     return(colSums(cbind(dB,dA)))
     
}


NBreg<-function(b, X, Y, single=TRUE){
     ## loglikelihood function for negative binomial 
  
     alpha<-exp(b[length(b)])
     beta<-b[-length(b)]
     e.XB<-exp(X %*% beta)
     e.XB[e.XB<.Machine$double.eps] <- .Machine$double.eps
     LL<-Y*log((alpha*e.XB)/(1+alpha*e.XB))-(1/alpha)*log(1+alpha*e.XB) + 
          lgamma(Y + 1/alpha)-lgamma(Y+1)-lgamma(1/alpha)
     if(single){
          return(sum(LL))
     }else{
          return(LL)
     }
}


