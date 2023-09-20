library(data.table)
library(car)
library(sandwich)
library(maxLik)
library(stargazer)
library(stringr)
library(MASS)
rm(list=ls())
source("nbreg.R")

load("../data/mainRemittanceData.rdata")


#### State side ####
remitNEW[,lledspending := shift(log(edspending+1)), by=ccode]
remitNEW[,llmilspending := shift(log(milspending+1)), by=ccode]
remitNEW[,llinfant := shift(log(inf.mortality+1)), by=ccode]
remitNEW[,llmilex := shift(log(milex/(pop)+1)), by=ccode]

m1.state <- glm.nb(domAttacks~ lremitPC*ldem + lremitPC*lano + llmil + llpop + lgdpgrowth + llgdppc + lpress  + llmilspending
                   +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                   +factor(ccode)+factor(year)-1,
                   data=remitNEW, subset=sam==1, maxit=25, x=TRUE)

while(m1.state$converged == FALSE | !is.null(m1.state$th.warn)){
  m1.statea <- maxLik(NBreg, NBgrad,
                      start=c(m1.state$coef, log(1/m1.state$theta)), 
                      X=m1.state$x, Y=m1.state$y,
                      method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  m1.state <- glm.nb(domAttacks~ lremitPC*ldem + lremitPC*lano + llmil + llpop + lgdpgrowth + llgdppc + lpress +llmilspending
                     +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                     +factor(ccode)+factor(year)-1,
                     data=remitNEW, subset=sam==1, maxit=25, x=TRUE,
                     start=m1.statea$est[-length(m1.statea$est)], init.theta=1/exp(m1.statea$est[length(m1.statea$est)]))
  
}
v1.state <- vcovCL(m1.state, remitNEW[sam==1]$ccode)


#### State side COW####
m1.stateB <- glm.nb(domAttacks~ lremitPC*ldem + lremitPC*lano + llmil + llpop + lgdpgrowth + llgdppc + lpress  + llmilex
                    +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                    +factor(ccode)+factor(year)-1,
                    data=remitNEW, subset=sam==1, maxit=25, x=TRUE)

while(m1.stateB$converged == FALSE | !is.null(m1.stateB$th.warn)){
  m1.stateBa <- maxLik(NBreg, NBgrad,
                       start=c(m1.stateB$coef, log(1/m1.stateB$theta)), 
                       X=m1.stateB$x, Y=m1.stateB$y,
                       method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  m1.stateB <- glm.nb(domAttacks~ lremitPC*ldem + lremitPC*lano + llmil + llpop + lgdpgrowth + llgdppc + lpress +llmilex
                      +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                      +factor(ccode)+factor(year)-1, 
                      data=remitNEW, subset=sam==1, maxit=25, x=TRUE,
                      start=m1.stateBa$est[-length(m1.stateBa$est)], init.theta=1/exp(m1.stateBa$est[length(m1.stateBa$est)]))
  
}
v1.stateB <- vcovCL(m1.stateB, remitNEW[sam==1]$ccode)



#### State side only infant####
m1.stateA <- glm.nb(domAttacks~ lremitPC*ldem + lremitPC*lano + llmil + llpop + lgdpgrowth + llgdppc + lpress  + llinfant
                    +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                    +factor(ccode)+factor(year)-1, 
                    data=remitNEW, subset=sam==1, maxit=25, x=TRUE)

while(m1.stateA$converged == FALSE | !is.null(m1.stateA$th.warn)){
  m1.stateAa <- maxLik(NBreg, NBgrad,
                       start=c(m1.stateA$coef, log(1/m1.stateA$theta)), 
                       X=m1.stateA$x, Y=m1.stateA$y,
                       method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  m1.stateA <- glm.nb(domAttacks~ lremitPC*ldem + lremitPC*lano + llmil + llpop + lgdpgrowth + llgdppc + lpress+llinfant
                      +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                      +factor(ccode)+factor(year)-1,
                      data=remitNEW, subset=sam==1, maxit=25, x=TRUE,
                      start=m1.stateAa$est[-length(m1.stateAa$est)], init.theta=1/exp(m1.stateAa$est[length(m1.stateAa$est)]))
  
}
v1.stateA <- vcovCL(m1.stateA, remitNEW[sam==1]$ccode)


#### State side infant and military####
m1.stateC <- glm.nb(domAttacks~ lremitPC*ldem + lremitPC*lano + llmil + llpop + lgdpgrowth + llgdppc + lpress  + llmilex+ llinfant
                    +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                    +factor(ccode)+factor(year)-1, 
                    data=remitNEW, subset=sam==1, maxit=25, x=TRUE)

while(m1.stateC$converged == FALSE | !is.null(m1.stateC$th.warn)){
  m1.stateCa <- maxLik(NBreg, NBgrad,
                       start=c(m1.stateC$coef, log(1/m1.stateC$theta)), 
                       X=m1.stateC$x, Y=m1.stateC$y,
                       method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  m1.stateC <- glm.nb(domAttacks~ lremitPC*ldem + lremitPC*lano + llmil + llpop + lgdpgrowth + llgdppc + lpress +llmilex+llinfant
                      +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                      +factor(ccode)+factor(year)-1,
                      data=remitNEW, subset=sam==1, maxit=25, x=TRUE,
                      start=m1.stateCa$est[-length(m1.stateCa$est)], init.theta=1/exp(m1.stateCa$est[length(m1.stateCa$est)]))
  
}
v1.stateC <- vcovCL(m1.stateC, remitNEW[sam==1]$ccode)




#### State side D####
m1.stateD <- glm(domAttacks~ lremitPC*ldem + lremitPC*lano + llmil + llpop + lgdpgrowth + llgdppc + lpress
                    + llmilex
                    +llinfant
                    +lledspending
                 +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                 +factor(ccode)+factor(year)-1,family="poisson",
                    data=remitNEW, subset=sam==1, maxit=25, x=TRUE)

while(m1.stateD$converged == FALSE | !is.null(m1.stateD$th.warn)|m1.stateD$family$family=="poisson"){
  m1.stateDa <- maxLik(NBreg, NBgrad,
                       start=c(m1.stateD$coef, 0), 
                       X=m1.stateD$x, Y=m1.stateD$y,
                       method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  m1.stateD <- glm.nb(domAttacks~ lremitPC*ldem + lremitPC*lano + llmil + llpop + lgdpgrowth + llgdppc + lpress    + llmilex
                      +llinfant
                      +lledspending
                      +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                      +factor(ccode)+factor(year)-1,
                      data=remitNEW, subset=sam==1, maxit=25, x=TRUE,
                      start=m1.stateDa$est[-length(m1.stateDa$est)], init.theta=1/exp(m1.stateDa$est[length(m1.stateDa$est)]))
  
}
v1.stateD <- vcovCL(m1.stateD, remitNEW[sam==1]$ccode)



cat("\nTable D.5\n")
stargazer(m1.stateA, m1.stateB, m1.stateC,
          title = "Effect of remittances after controlling for spending",
          label="tab:alternative",
          dep.var.labels = "Domestic terrorist attacks",
          omit=c("factor*", "*.bar",
                 "lremitPC:lano", "lano", "^ldem$",
                 "^llmil$", "llpop", "lgdpgrowth", "llgdppc", "lpress", 
                 "lgini", "lmaxlowx", "lwexcluded", "lnum_ongoing_ucdp"),
          order=c("lremitPC"),
          se=list(sqrt(diag(v1.stateA)), sqrt(diag(v1.stateB)), 
                  sqrt(diag(v1.stateC))),
          align=TRUE,       
          no.space=TRUE,
          digits=2,
          covariate.labels = c("Remittances per capita",
                               "Remittances per capita $\\times$ Dem.",
                               "Infant mortality",
                               "Military expenditures"),
          omit.stat = "AIC",
          star.cutoffs = c(0.1, .05, NA),
          notes=c("$p<[0.*], p<[0.**]$. Coefficients from negative binomial models. Clustered standard errors in parentheses."),
          notes.append=FALSE
)

est <- c(m1.stateA$coefficients["lremitPC"]+m1.stateA$coefficients["lremitPC:ldem"],
         m1.stateB$coefficients["lremitPC"]+m1.stateB$coefficients["lremitPC:ldem"],
         m1.stateC$coefficients["lremitPC"]+m1.stateC$coefficients["lremitPC:ldem"])
se <- sqrt(c(v1.stateA[1, 1] + v1.stateA["lremitPC:ldem","lremitPC:ldem"]+2* v1.stateA[1,"lremitPC:ldem"],
             v1.stateB[1, 1] + v1.stateB["lremitPC:ldem","lremitPC:ldem"]+2* v1.stateB[1,"lremitPC:ldem"],
             v1.stateC[1, 1] + v1.stateC["lremitPC:ldem","lremitPC:ldem"]+2* v1.stateC[1,"lremitPC:ldem"]))
p.vals <- pnorm(abs(est/se), lower=F)*2
cat("Combined coefficients for table D5\n")
print(round(rbind(est,se,p.vals), 2))



