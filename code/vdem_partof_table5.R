library(data.table)
library(lmtest)
library(car)
library(sandwich)
library(MASS)
library(xtable)
library(stringr)
library(stargazer)
library(maxLik)
library(nonnest2)
rm(list=ls())
load("../data/mainRemittanceData.rdata")


setkey(remitNEW,"ccode", "year")

remitNEW[,lv2x_partipdem := shift(v2x_partipdem,1), by=ccode]
remitNEW[,lv2xlg_legcon:= shift(v2xlg_legcon,1), by=ccode]
remitNEW[,lv2x_cspart := shift(v2x_cspart,1), by=ccode]
remitNEW[,lv2psoppaut := shift(v2psoppaut,1), by=ccode]

remitNEW.1a <- copy(remitNEW)
remitNEW.1a$sam <- 1-apply(remitNEW.1a[,list(domAttacks,lremitPC, 
                                             lv2x_partipdem,lv2xlg_legcon,lv2x_cspart,
                                             llmil, llpop, lgdpgrowth, llgdppc, lpress, year)], 1, anyNA)
# Then we flag all 0 countries
remitNEW.1a[,sumAttacksDom:= sum(domAttacks*sam, na.rm=T) , by=ccode ]
remitNEW.1a[,sam:=sam*(sumAttacksDom>0)] #remove countries with no attacks in sample (FE -> -infty)




source("nbreg.R")
#### participatory democracy #####
r1.dummies <- glm(domAttacks~ lremitPC*lv2x_partipdem 
                  + llmil + llpop + lgdpgrowth + llgdppc + lpress 
                  +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                  +factor(ccode)+factor(year)-1, family="poisson",
                  data=remitNEW.1a, subset=sam==1, maxit=25, x=TRUE)
while(r1.dummies$converged == FALSE|  r1.dummies$family$family=="poisson"){
  r1.dummiesa <- maxLik(NBreg, NBgrad,
                        start=c(r1.dummies$coef, 0), 
                        X=r1.dummies$x, Y=r1.dummies$y,
                        method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  r1.dummies <- glm.nb(domAttacks~  lremitPC*lv2x_partipdem
                       + llmil + llpop + lgdpgrowth + llgdppc + lpress 
                       +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                       +factor(ccode)+factor(year)-1-1, 
                       data=remitNEW.1a, subset=sam==1, maxit=25, x=TRUE,
                       start=r1.dummiesa$est[-length(r1.dummiesa$est)], init.theta=1/exp(r1.dummiesa$est[length(r1.dummiesa$est)]))
  
}
v1.dummies <- vcovCL(r1.dummies, remitNEW.1a[sam==1]$ccode)
df.margins <- as.data.frame(remitNEW.1a[as.numeric(rownames(r1.dummies$model)),])
P <- matrix(quantile(remitNEW.1a[sam==1]$lv2x_partipdem,c(0.1,0.9)), ncol=1)
tests <- 
  cbind(P,
        r1.dummies$coef["lremitPC"] 
        +  P[,1] *r1.dummies$coef["lremitPC:lv2x_partipdem"],
        t(sapply(paste("lremitPC+",
                       P[,1], "*lremitPC:lv2x_partipdem"), 
                 function(x){
                   H <- linearHypothesis(r1.dummies, x, vcov=v1.dummies)
                   return(c(H$Chisq[2], H$Pr[2]))})))
tests <- as.data.frame(tests)
colnames(tests) <- c("lv2x_partipdem", "coef","chi", "p")
tests <- tests[order(tests$lv2x_partipdem),]
rownames(tests) <- NULL
tests$chi <- formatC(tests$chi,digits=2 ,format='f')
tests$coef <- formatC(tests$coef,digits=2 ,format='f')
tests$coef <- ifelse(tests$p < 0.05,
                     paste(tests$coef, "^{**}", sep=""),
                     ifelse(tests$p <0.1, 
                            paste(tests$coef, "^*", sep=""),
                            paste(tests$coef)))

tests$p <- NULL

#### participation and constraints ####
r2.dummies <- glm(domAttacks~ lremitPC*lv2x_partipdem +lremitPC*lv2xlg_legcon
                  + llmil + llpop + lgdpgrowth + llgdppc + lpress 
                  +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                  +factor(ccode)+factor(year)-1, family="poisson",
                  data=remitNEW.1a, subset=sam==1, maxit=25, x=TRUE)
while(r2.dummies$converged == FALSE|  r2.dummies$family$family=="poisson"){
  r2.dummiesa <- maxLik(NBreg, NBgrad,
                        start=c(r2.dummies$coef, 0), 
                        X=r2.dummies$x, Y=r2.dummies$y,
                        method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  r2.dummies <- glm.nb(domAttacks~ lremitPC*lv2x_partipdem + lremitPC*lv2xlg_legcon
                       + llmil + llpop + lgdpgrowth + llgdppc + lpress 
                       +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                       +factor(ccode)+factor(year)-1, 
                       data=remitNEW.1a, subset=sam==1, maxit=25, x=TRUE,
                       start=r2.dummiesa$est[-length(r2.dummiesa$est)], init.theta=1/exp(r2.dummiesa$est[length(r2.dummiesa$est)]))
  
}
v2.dummies <- vcovCL(r2.dummies, remitNEW.1a[sam==1]$ccode)


P2 <- expand.grid(lv2x_partipdem=quantile(remitNEW.1a[sam==1]$lv2x_partipdem,c(0.1, .9)),
                  lv2xlg_legcon=median(remitNEW.1a[sam==1]$lv2xlg_legcon))
tests2 <-
  cbind(P2,
        r2.dummies$coef["lremitPC"]
        +  P2[,1] *r2.dummies$coef["lremitPC:lv2x_partipdem"]
        +  P2[,2] *r2.dummies$coef["lremitPC:lv2xlg_legcon"],
        t(sapply(paste("lremitPC+",
                       P2[,1], "*lremitPC:lv2x_partipdem+",
                       P2[,2], "*lremitPC:lv2xlg_legcon"),
                 function(x){
                   H <- linearHypothesis(r2.dummies, x, vcov=v2.dummies)
                   return(c(H$Chisq[2], H$Pr[2]))})))
tests2 <- as.data.frame(tests2)
colnames(tests2) <- c("lv2x_partipdem", "lv2xlg_legcon", "coef","chi", "p")
tests2 <- tests2[order(tests2$lv2x_partipdem, tests2$lv2xlg_legcon),]
rownames(tests2) <- NULL
tests2$chi <- formatC(tests2$chi,digits=2 ,format='f')
tests2$coef <- formatC(tests2$coef,digits=2 ,format='f')
tests2$chi <- ifelse(tests2$p < 0.05,
                     paste(tests2$chi, "^{**}", sep=""),
                     ifelse(tests2$p <0.1,
                            paste(tests2$chi, "^*", sep=""),
                            paste(tests2$chi)))
tests2$p <- NULL




######## r3 through r5 mentioned, but not presented in text #############
#### Constraints only ####
r3.dummies <- glm(domAttacks~ lremitPC*lv2xlg_legcon
                  + llmil + llpop + lgdpgrowth + llgdppc + lpress 
                  +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                  +factor(ccode)+factor(year)-1, family="poisson",
                  data=remitNEW.1a, subset=sam==1, maxit=25, x=TRUE)

while(r3.dummies$converged == FALSE|r3.dummies$family$family=="poisson"){
  r3.dummiesa <- maxLik(NBreg, NBgrad,
                        start=c(r3.dummies$coef, 0), 
                        X=r3.dummies$x, Y=r3.dummies$y,
                        method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  r3.dummies <- glm.nb(domAttacks~  lremitPC*lv2xlg_legcon
                       + llmil + llpop + lgdpgrowth + llgdppc + lpress 
                       +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                       +factor(ccode)+factor(year)-1, 
                       data=remitNEW.1a, subset=sam==1, maxit=25, x=TRUE,
                       start=r3.dummiesa$est[-length(r3.dummiesa$est)], init.theta=1/exp(r3.dummiesa$est[length(r3.dummiesa$est)]))
  
}
v3.dummies <- vcovCL(r3.dummies, remitNEW.1a[sam==1]$ccode)



#### participation other measure: civil participation##### 
r4.dummies <- glm(domAttacks~ lremitPC*lv2x_cspart
                  + llmil + llpop + lgdpgrowth + llgdppc + lpress 
                  +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                  +factor(ccode)+factor(year)-1, family="poisson",
                  data=remitNEW.1a, subset=sam==1, maxit=25, x=TRUE)
# r3.dummies$converged <-FALSE
# r3.dummies$theta <- .5
while(r4.dummies$converged == FALSE|r4.dummies$family$family=="poisson"){
  r4.dummiesa <- maxLik(NBreg, NBgrad,
                        start=c(r4.dummies$coef, 0), 
                        X=r4.dummies$x, Y=r4.dummies$y,
                        method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  r4.dummies <- glm.nb(domAttacks~  lremitPC*lv2x_cspart
                       + llmil + llpop + lgdpgrowth + llgdppc + lpress 
                       +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                       +factor(ccode)+factor(year)-1, 
                       data=remitNEW.1a, subset=sam==1, maxit=25, x=TRUE,
                       start=r4.dummiesa$est[-length(r4.dummiesa$est)], init.theta=1/exp(r4.dummiesa$est[length(r4.dummiesa$est)]))
  
}
v4.dummies <- vcovCL(r4.dummies, remitNEW.1a[sam==1]$ccode)
## Combined coefs
deltaMethod(r4.dummies, "lremitPC+0.25*`lremitPC:lv2x_cspart`", vcov=v4.dummies)
deltaMethod(r4.dummies, "lremitPC+0.93*`lremitPC:lv2x_cspart`", vcov=v4.dummies)


#### participation other measure: opposition autonomy #### 
remitNEW.1a[,sam2 := as.numeric(!is.na(lv2psoppaut))]
r5.dummies <- glm(domAttacks~ lremitPC*lv2psoppaut
                  + llmil + llpop + lgdpgrowth + llgdppc + lpress 
                  +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                  +factor(ccode)+factor(year)-1, family="poisson",
                  data=remitNEW.1a, subset=sam==1& sam2==1, maxit=25, x=TRUE)
while(r5.dummies$converged == FALSE|r5.dummies$family$family=="poisson"){
  r5.dummiesa <- maxLik(NBreg, NBgrad,
                        start=c(r5.dummies$coef, 0), 
                        X=r5.dummies$x, Y=r5.dummies$y,
                        method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  r5.dummies <- glm.nb(domAttacks~  lremitPC*lv2psoppaut
                       + llmil + llpop + lgdpgrowth + llgdppc + lpress 
                       +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                       +factor(ccode)+factor(year)-1, 
                       data=remitNEW.1a, subset=sam==1 & sam2==1, maxit=25, x=TRUE,
                       start=r5.dummiesa$est[-length(r5.dummiesa$est)], init.theta=1/exp(r5.dummiesa$est[length(r5.dummiesa$est)]))
  
}
v5.dummies <- vcovCL(r5.dummies, remitNEW.1a[sam==1& sam2==1]$ccode)
## Combined coefs
deltaMethod(r5.dummies, "lremitPC+-1.436*`lremitPC:lv2psoppaut`", vcov=v5.dummies)
deltaMethod(r5.dummies, "lremitPC+2.761*`lremitPC:lv2psoppaut`", vcov=v5.dummies)


cat("Vuong test for the V-Dem measures: competition v. constraints\n")
vuongtest(r1.dummies, r3.dummies) 
cat("Inconclusive, but suggestive in favor of competition\n")

# model fit stats favor competition
cat("AIC selection for the V-Dem measures: competition v. constraints\n")
which.min(c(AIC(r1.dummies),AIC(r3.dummies)))
cat("Selects compeition\n")


vdem1.dummies <- r1.dummies
vdem2.dummies <- r2.dummies 

vdemV1.dummies <- v1.dummies
vdemV2.dummies <- v2.dummies

# combined coef test
test.vdem1 <- tests
test.vdem2 <- tests2

VdemDataset <- copy(remitNEW.1a)
save(list=c("vdem1.dummies", "vdem2.dummies",
            "vdemV1.dummies", "vdemV2.dummies", "VdemDataset",
            "test.vdem1", "test.vdem2"), file="../output/vdem_competition.rdata")

