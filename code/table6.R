library(data.table)
library(sandwich)
library(MASS)
library(modelsummary)
library(car)
library(maxLik)
rm(list=ls())


#### Set up ####
load("../data/mainRemittanceData.rdata")
source("nbreg.R")
setkey(remitNEW, "ccode", "year")



remitNEW[,lPR:=shift(pr,1), by=ccode]
remitNEW.dem <- copy(remitNEW)
remitNEW.dem <- subset(remitNEW.dem, polity2>=5)




remitNEW.dem[,lparties := (shift(ENP_nat)), by=ccode]
remitNEW.dem[,llparties := log(shift(ENP_nat)), by=ccode]
remitNEW.dem[,lwexcluded := shift(weightedExclude), by=ccode]
remitNEW.dem[,lnexcluded := shift(nExclude), by=ccode]


##### PR ####
remitNEW.dem$sam <- 1-apply(remitNEW.dem[,list(domAttacks,lremitPC, lPR, llmil, llpop,
                                               lgdpgrowth, llgdppc, lpress,
                                               lgini, lnum_ongoing_ucdp,
                                               lfrac, lmaxlowx, ethfrac,
                                               year)], 1, anyNA)
## Then we flag all 0 countries
remitNEW.dem[,sumAttacksDom:= sum(domAttacks*sam, na.rm=T) , by=ccode ]
remitNEW.dem[,sam:=sam*(sumAttacksDom>0)] #remove countries with no attacks in sample (FE -> -infty)

r3.dummies <- glm(domAttacks~ lremitPC*lPR 
                  + llmil + llpop + lgdpgrowth + llgdppc + lpress
                  +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                  +factor(ccode)+factor(year)-1, family="poisson", 
                  data=remitNEW.dem, subset=sam==1, maxit=25, x=TRUE)
r3.dummies$converged <-FALSE
r3.dummies$theta <- .5
while(r3.dummies$converged == FALSE){
  r3.dummiesa <- maxLik(NBreg, NBgrad,
                        start=c(r3.dummies$coef, log(1/r3.dummies$theta)), 
                        X=r3.dummies$x, Y=r3.dummies$y,
                        method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  r3.dummies <- glm.nb(domAttacks~ lremitPC*lPR 
                       + llmil + llpop + lgdpgrowth + llgdppc + lpress
                       +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                       +factor(ccode)+factor(year)-1, 
                       data=remitNEW.dem, subset=sam==1, maxit=25, x=TRUE,
                       start=r3.dummiesa$est[-length(r3.dummiesa$est)],
                       init.theta=1/exp(r3.dummiesa$est[length(r3.dummiesa$est)]))
  
}
v3.dummies <- vcovCL(r3.dummies, remitNEW.dem[sam==1]$ccode)
df.margins <- as.data.frame(remitNEW[as.numeric(rownames(r3.dummies$model)),])
Est1 <- r3.dummies$coef['lremitPC']+r3.dummies$coef['lremitPC:lPR']
SE1 <- sqrt(v3.dummies['lremitPC', 'lremitPC']+
              v3.dummies['lremitPC:lPR', 'lremitPC:lPR'] +
              2*v3.dummies['lremitPC:lPR', 'lremitPC'])


#### remittances and frac/parties #### 

r3A.dummies <- glm(domAttacks~ lremitPC*lfrac
                   + llmil + llpop + lgdpgrowth + llgdppc + lpress 
                   +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                   +factor(ccode)+factor(year)-1, family="poisson", 
                   data=remitNEW.dem, subset=sam==1, maxit=25, x=TRUE)
r3A.dummies$converged <-FALSE
r3A.dummies$theta <- .5
while(r3A.dummies$converged == FALSE){
  r3A.dummiesa <- maxLik(NBreg, NBgrad,
                         start=c(r3A.dummies$coef, log(1/r3A.dummies$theta)), 
                         X=r3A.dummies$x, Y=r3A.dummies$y,
                         method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  r3A.dummies <- glm.nb(domAttacks~ lremitPC*lfrac
                        + llmil + llpop + lgdpgrowth + llgdppc + lpress 
                        +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                        +factor(ccode)+factor(year)-1, 
                        data=remitNEW.dem, subset=sam==1 & polity2 >=5, maxit=25, x=TRUE,
                        start=r3A.dummiesa$est[-length(r3A.dummiesa$est)],
                        init.theta=1/exp(r3A.dummiesa$est[length(r3A.dummiesa$est)]))
  
}
v3A.dummies <- vcovCL(r3A.dummies, remitNEW.dem[sam==1&polity2>=5]$ccode)
# cat("lfrac quantiles \n ")
# print(round(quantile(remitNEW.dem[sam==1]$lfrac),2))

lincom3A <- deltaMethod(r3A.dummies, 
                        "lremitPC+0.57*`lremitPC:lfrac`", 
                        vcov=v3A.dummies)
est3A <- lincom3A$Estimate
SE3A <- lincom3A$SE
lincom3A <- deltaMethod(r3A.dummies, 
                        "lremitPC+0.76*`lremitPC:lfrac`", 
                        vcov=v3A.dummies)
est3A2 <- lincom3A$Estimate
SE3A2 <- lincom3A$SE
equalH <- linearHypothesis(r3A.dummies, 
                 "lremitPC+0.57*lremitPC:lfrac=lremitPC+0.76*lremitPC:lfrac", 
                 vcov=v3A.dummies)


##### leg and eth frac groups #####
r3D.dummies <- glm(domAttacks~ lremitPC*ethfrac*lfrac-ethfrac
                   + llmil + llpop + lgdpgrowth + llgdppc + lpress 
                   +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                   +factor(ccode)+factor(year)-1, family="poisson", 
                   data=remitNEW.dem, subset=sam==1, maxit=25, x=TRUE)
r3D.dummies$converged <-FALSE
r3D.dummies$theta <- .5
while(r3D.dummies$converged == FALSE){
  r3D.dummiesa <- maxLik(NBreg, NBgrad,
                         start=c(r3D.dummies$coef, log(1/r3D.dummies$theta)), 
                         X=r3D.dummies$x, Y=r3D.dummies$y,
                         method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  r3D.dummies <- glm.nb(domAttacks~  lremitPC*ethfrac*lfrac-ethfrac
                        + llmil + llpop + lgdpgrowth + llgdppc + lpress 
                        +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                        +factor(ccode)+factor(year)-1, 
                        data=remitNEW.dem, subset=sam==1 , maxit=25, x=TRUE,
                        start=r3D.dummiesa$est[-length(r3D.dummiesa$est)],
                        init.theta=1/exp(r3D.dummiesa$est[length(r3D.dummiesa$est)]))
  
}
v3D.dummies <- vcovCL(r3D.dummies, remitNEW.dem[sam==1]$ccode)
# cat("Quantiles for eth and l frac\n")
# print(round(quantile(remitNEW.dem[sam==1]$ethfrac, na.rm=T), 2))
# print(round(quantile(remitNEW.dem[sam==1]$lfrac, na.rm=T), 2))

lincom3D <- deltaMethod(r3D.dummies, 
                        "lremitPC+0.08*`lremitPC:ethfrac` + 
                        .57*`lremitPC:lfrac` + 
                        0.08*0.57*`lremitPC:ethfrac:lfrac`", 
                        vcov=v3D.dummies)
est3D <- lincom3D$Estimate
SE3D <- lincom3D$SE
lincom3D <- deltaMethod(r3D.dummies, 
                        "lremitPC+0.08*`lremitPC:ethfrac` + 
                        .76*`lremitPC:lfrac` + 
                        0.08*0.76*`lremitPC:ethfrac:lfrac`", 
                        vcov=v3D.dummies)
est3D2 <- lincom3D$Estimate
SE3D2 <- lincom3D$SE



lincom3D <- deltaMethod(r3D.dummies, 
                        "lremitPC+.53*`lremitPC:ethfrac` + .57*`lremitPC:lfrac` +
                        .53*.57*`lremitPC:ethfrac:lfrac`", 
                        vcov=v3D.dummies)
est3D3 <- lincom3D$Estimate
SE3D3 <- lincom3D$SE
lincom3D <- deltaMethod(r3D.dummies, 
                        "lremitPC+.53*`lremitPC:ethfrac` + 
                        .76*`lremitPC:lfrac` + .53*.76*`lremitPC:ethfrac:lfrac`", 
                        vcov=v3D.dummies)
est3D4 <- lincom3D$Estimate
SE3D4 <- lincom3D$SE

eqH1 <- linearHypothesis(r3D.dummies, 
                        "lremitPC+0.08*lremitPC:ethfrac + 
                        .57*lremitPC:lfrac + 
                        0.0456*lremitPC:ethfrac:lfrac=
lremitPC+0.08*lremitPC:ethfrac + 
                        .76*lremitPC:lfrac + 
                        0.0608*lremitPC:ethfrac:lfrac",
                 vcov=v3D.dummies)
eqH2 <- linearHypothesis(r3D.dummies, 
                        "lremitPC+0.53*lremitPC:ethfrac + 
                        .57*lremitPC:lfrac + 
                        .3021*lremitPC:ethfrac:lfrac=
lremitPC+0.53*lremitPC:ethfrac + 
                        .76*lremitPC:lfrac + 
                        0.4028*lremitPC:ethfrac:lfrac",
                 vcov=v3D.dummies)


###### frac and PR ####
r3C.dummies <- glm(domAttacks~ lremitPC*ethfrac*lPR-ethfrac
                   + llmil + llpop + lgdpgrowth + llgdppc + lpress 
                   +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                   +factor(ccode)+factor(year)-1, family="poisson", 
                   data=remitNEW.dem, subset=sam==1, maxit=25, x=TRUE)
r3C.dummies$converged <-FALSE
r3C.dummies$theta <- .5
while(r3C.dummies$converged == FALSE){
  r3C.dummiesa <- maxLik(NBreg, NBgrad,
                         start=c(r3C.dummies$coef, log(1/r3C.dummies$theta)), 
                         X=r3C.dummies$x, Y=r3C.dummies$y,
                         method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  r3C.dummies <- glm.nb(domAttacks~ lremitPC*ethfrac*lPR-ethfrac
                        + llmil + llpop + lgdpgrowth + llgdppc + lpress 
                        +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                        +factor(ccode)+factor(year)-1, 
                        data=remitNEW.dem, subset=sam==1 , maxit=25, x=TRUE,
                        start=r3C.dummiesa$est[-length(r3C.dummiesa$est)],
                        init.theta=1/exp(r3C.dummiesa$est[length(r3C.dummiesa$est)]))
  
}
v3C.dummies <- vcovCL(r3C.dummies, remitNEW.dem[sam==1]$ccode)
# cat("quantiles for ethfrac\n")
# print(round(quantile(remitNEW.dem[sam==1]$ethfrac, na.rm=T),2))



lincom3C <- deltaMethod(r3C.dummies, 
                        "lremitPC+0.08*`lremitPC:ethfrac`", 
                        vcov=v3C.dummies)
est3C <- lincom3C$Estimate
SE3C <- lincom3C$SE

lincom3C <- deltaMethod(r3C.dummies, 
                        "lremitPC+0.08*`lremitPC:ethfrac` + 
                        `lremitPC:lPR` + 0.08*`lremitPC:ethfrac:lPR`", 
                        vcov=v3C.dummies)
est3C2 <- lincom3C$Estimate
SE3C2 <- lincom3C$SE


lincom3C <- deltaMethod(r3C.dummies, 
                        "lremitPC+.53*`lremitPC:ethfrac` ", 
                        vcov=v3C.dummies)
est3C3 <- lincom3C$Estimate
SE3C3 <- lincom3C$SE

lincom3C <- deltaMethod(r3C.dummies, 
                        "lremitPC+.53*`lremitPC:ethfrac` +
                        `lremitPC:lPR` + .53*`lremitPC:ethfrac:lPR`", 
                        vcov=v3C.dummies)
est3C4 <- lincom3C$Estimate
SE3C4 <- lincom3C$SE



##### REMITTANCES AND COMPETITION #####

remitNEW.leg <- subset(remitNEW.dem, legelec.lag==1)
r4.dummies <- lm(oppvote~ lremitPC+lPR+ lfrac
                 +llpop  + llgdppc + lgdpgrowth + lpress
                 +lgini +  lmaxlowx + lwexcluded 
                 +factor(ccode)+factor(year)-1,
                 data=remitNEW.leg, #subset=sam==1 ,
                 x=TRUE)

v4.dummies <- vcovCL(r4.dummies, remitNEW.leg$ccode)



cat("Table 6\n")
glance_custom.negbin <- function(x, ...) {
  theta <- x$theta
  out <- data.frame("theta" = theta)
  return(out)
}
rows <- rbind.data.frame(c("Country fixed effects", rep("\\multicolumn{1}{c}{Yes}", 5)),
                         c("Year fixed effects", rep("\\multicolumn{1}{c}{Yes}", 5)),
                         c("Controls", rep("\\multicolumn{1}{c}{Yes}", 5)))
attr(rows, 'position') <- c(15:17)

cat(modelsummary(list(r3.dummies,
                  r3C.dummies,
                  r3A.dummies,
                  r3D.dummies, r4.dummies),
             output="latex",
             fmt=2,
             vcov=list(v3.dummies,
                       v3C.dummies,
                        v3A.dummies,
                       v3D.dummies,
                       v4.dummies),
             stars=c("*"=.1, "**"=0.05),
             coef_map=c("lremitPC"="Remittances",
                        "lremitPC:lPR" ="Remittances $\\times$ PR",
                        "lremitPC:lfrac" ="Remittances $\\times$ legislative frac.",
                        "lremitPC:ethfrac"  = "Remittances $\\times$ ELF",
                        "lremitPC:ethfrac:lPR" = "Remittances $\\times$ ELF  $\\times$ PR",
                        "lremitPC:ethfrac:lfrac" = "Remittances $\\times$ ELF $\\times$ leg. frac."),
             escape=FALSE,
             gof_map=data.frame(raw=c("nobs", "logLik", "theta", "r.squared"),
                                clean=c("Observations", "Log Likelihood", "$\\theta$", "$R^2$"),
                                fmt=c(0,2,2,2)),
             add_rows = rows,
             align="lddddd",
             notes="\\\\footnotesize $^*p<0.1$, $^{**}p<0.05$. Regression coefficients. Standard errors in parentheses, clustered on country.",
             title = "Political competition within democracies\\label{tab:demMech}"))


cat("Table 6 combined coefficients (est se p-val) \n")
fulComb <- cbind(c(r3.dummies$coef[1], sqrt(v3.dummies[1,1]),
                   pnorm(abs(r3.dummies$coef[1]/sqrt(v3.dummies[1,1])), lower=F)*2, 
                   rep(NA, 6),
                   Est1, SE1, pnorm(abs(Est1/SE1), lower=FALSE)*2,
                   rep(NA, 6)),
                 c(NA, NA, NA,
                   est3C, SE3C, pnorm(abs(est3C/SE3C), lower=FALSE)*2,
                   est3C3, SE3C3, pnorm(abs(est3C3/SE3C3), lower=FALSE)*2,
                   NA, NA, NA,
                   est3C2, SE3C2, pnorm(abs(est3C2/SE3C2), lower=FALSE)*2,
                   est3C4, SE3C4, pnorm(abs(est3C4/SE3C4), lower=FALSE)*2),
                 c(est3A, SE3A, pnorm(abs(est3A/SE3A), lower=FALSE)*2,      
                   rep(NA, 6),
                   est3A2, SE3A2, pnorm(abs(est3A2/SE3A2), lower=FALSE)*2,
                   rep(NA, 6)),
                 c(rep(NA, 3),
                   est3D, SE3D, pnorm(abs(est3D/SE3D), lower=FALSE)*2,
                   est3D3, SE3D3, pnorm(abs(est3D3/SE3D3), lower=FALSE)*2,
                   rep(NA, 3),
                   est3D2, SE3D2, pnorm(abs(est3D2/SE3D2), lower=FALSE)*2,
                   est3D4, SE3D4, pnorm(abs(est3D4/SE3D4), lower=FALSE)*2),
                 NA)
print(round(fulComb,2))
