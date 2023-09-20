library(data.table)
library(sandwich)
library(margins)
library(MASS)
library(car)
library(maxLik)
rm(list=ls())

source("nbreg.R")
load("../data/mainRemittanceData.rdata")


m0.dummies <- glm.nb(domAttacks~ lremitPC*ldem + lremitPC*lano
                     + factor(ccode)+factor(year)-1,
                     data=remitNEW, subset=sam==1, maxit=25, x=TRUE)
while(m0.dummies$converged == FALSE || !is.null( m0.dummies$th.warn) ){
  m0.dummiesa <- maxLik(NBreg, NBgrad,
                        start=c(m0.dummies$coef, log(1/m0.dummies$theta)),
                        X=m0.dummies$x, Y=m0.dummies$y,
                        method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  m0.dummies <- glm.nb(domAttacks~ lremitPC*ldem + lremitPC*lano
                       + factor(ccode)+factor(year)-1,
                       data=remitNEW, subset=sam==1, maxit=25, x=TRUE,
                       start=m0.dummiesa$est[-length(m0.dummiesa$est)], init.theta=1/exp(m0.dummiesa$est[length(m0.dummiesa$est)]))
  
}
df.margins <- as.data.frame(remitNEW[sam==1])

v0.dummies <- vcovCL(m0.dummies, remitNEW[sam==1]$ccode)
mar0 <- summary(margins(m0.dummies,
                        "lremitPC",
                        vcov=v0.dummies,
                        data = df.margins,
                        at=data.frame(ldem=c(1,0), lano=c(0,1))))




m2.dummies <- glm.nb(domAttacks~ lremitPC*ldem + lremitPC*lano + llmil + llpop + lgdpgrowth + llgdppc + lpress
                     +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                     + factor(ccode)+factor(year)-1,
                     data=remitNEW, subset=sam==1, maxit=25, x=TRUE)
while(m2.dummies$converged == FALSE){
  m2.dummiesa <- maxLik(NBreg, NBgrad,
                        start=c(m2.dummies$coef, log(1/m2.dummies$theta)),
                        X=m2.dummies$x, Y=m2.dummies$y,
                        method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  m2.dummies <- glm.nb(domAttacks~  lremitPC*ldem + lremitPC*lano+ llmil + llpop + lgdpgrowth + llgdppc + lpress
                       +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                       + factor(ccode)+factor(year)-1,
                       data=remitNEW, subset=sam==1, maxit=25, x=TRUE,
                       start=m2.dummiesa$est[-length(m2.dummiesa$est)],
                       init.theta=1/exp(m2.dummiesa$est[length(m2.dummiesa$est)]))
  
}
v2.dummies <- vcovCL(m2.dummies, remitNEW[sam==1]$ccode)
mar2 <- summary(margins(m2.dummies,
                        "lremitPC",
                        vcov=v2.dummies,
                        data = df.margins,
                        at=data.frame(ldem=c(1,0,0), lano=c(0,0,1))))







m3.dummies <- glm.nb(domAttacks~ lremitPC*ldem + lremitPC*lano + llmil + llpop + lgdpgrowth + llgdppc + lpress
                     +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                     + lag_domAttacks
                     + factor(ccode) + factor(year)-1,
                     data=remitNEW, subset=sam==1, maxit=25, x=TRUE)
while(m3.dummies$converged == FALSE){
  m3.dummiesa <- maxLik(NBreg, NBgrad,
                        start=c(m3.dummies$coef, log(1/m3.dummies$theta)),
                        X=m3.dummies$x, Y=m3.dummies$y,
                        method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  m3.dummies <- glm.nb(domAttacks~  lremitPC*ldem + lremitPC*lano+ llmil + llpop + lgdpgrowth + llgdppc + lpress
                       +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                       + lag_domAttacks
                       + factor(ccode) + factor(year)-1,
                       data=remitNEW, subset=sam==1, maxit=25, x=TRUE,
                       start=m3.dummiesa$est[-length(m3.dummiesa$est)],
                       init.theta=1/exp(m3.dummiesa$est[length(m3.dummiesa$est)]))
  
}
v3.dummies <- vcovCL(m3.dummies, remitNEW[sam==1]$ccode)
mar3 <- summary(margins(m3.dummies,
                        "lremitPC",
                        vcov=v3.dummies,
                        data = df.margins,
                        at=data.frame(ldem=c(1,0,0), lano=c(0,0,1))))



m4.dummies <- glm.nb(domAttacks~ lremitPC*ldem*llgdppc + lremitPC*lano*llgdppc + llmil + llpop + lgdpgrowth  + lpress
                     +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                     + factor(ccode) + factor(year)-1,
                     data=remitNEW, subset=sam==1,  maxit=25, x=TRUE)
while(m4.dummies$converged == FALSE |!is.null(m4.dummies$th.warn) ){
  m4.dummiesa <- maxLik(NBreg, NBgrad,
                        start=c(m4.dummies$coef, log(1/m4.dummies$theta)),
                        X=m4.dummies$x, Y=m4.dummies$y,
                        method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  m4.dummies <- glm.nb(domAttacks~ lremitPC*ldem*llgdppc + lremitPC*lano*llgdppc+ llmil + llpop + lgdpgrowth + lpress
                       +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                       + factor(ccode) + factor(year)-1,
                       data=remitNEW[sam==1], maxit=25, x=TRUE,
                       start=m4.dummiesa$est[-length(m4.dummiesa$est)],
                       init.theta=1/exp(m4.dummiesa$est[length(m4.dummiesa$est)]))
  
}
v4.dummies <- vcovCL(m4.dummies, remitNEW[sam==1]$ccode)

# Quantiles by regime
newDat <- remitNEW[sam==1,quantile(llgdppc, c(0.25,.75))]
newDat <- expand.grid(ldem=c(0,1),llgdppc=newDat)
coefs4 <- do.call(rbind,
                  apply(newDat, 1, \(x){deltaMethod(m4.dummies,
                                                    paste("lremitPC+",
                                                          x[1],"*`lremitPC:ldem`+",
                                                          x[2],"*`lremitPC:llgdppc`+",
                                                          x[1]*x[2],"*`lremitPC:ldem:llgdppc`"),
                                                    vcov=v4.dummies)})
)
rownames(coefs4) <- c("auto-lo", "dem-lo", "auto-hi", "dem-hi")
mar4 <- summary(margins(m4.dummies,
                        "lremitPC",
                        vcov=v4.dummies,
                        data = df.margins,
                        at=newDat))
cat("GDP pc quantiles\n")
print(exp(newDat))



save.image("../output/mainResultsRemittances.rdata")
