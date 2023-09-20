library(data.table)
library(car)
library(sandwich)
library(MASS)
library(modelsummary)
library(nonnest2)
library(maxLik)
rm(list=ls())
load("../data/mainRemittanceData.rdata")


## set up ##
source("nbreg.R")
remitNEW.1a <- copy(remitNEW)
remitNEW.1a$sam <- 1-apply(remitNEW.1a[,list(domAttacks,lremitPC, lxrcomp,lxconst, 
                                             llmil, llpop, lgdpgrowth, llgdppc, lpress,
                                             lgini,lnum_ongoing_ucdp,lmaxlowx,
                                             year)], 1, anyNA)
# Then we flag all 0 countries
remitNEW.1a[,sumAttacksDom:= sum(domAttacks*sam, na.rm=T) , by=ccode ]
remitNEW.1a[,sam:=sam*(sumAttacksDom>0)] #remove countries with no attacks in sample (FE -> -infty)




### competition ## 
r1.dummies <- glm.nb(domAttacks~ lremitPC*lxrcomp 
                     + llmil + llpop + lgdpgrowth + llgdppc + lpress 
                     +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                     +factor(ccode)+factor(year)-1, 
                     data=remitNEW.1a, subset=sam==1, maxit=25, x=TRUE)
while(r1.dummies$converged == FALSE){
  r1.dummiesa <- maxLik(NBreg, NBgrad,
                        start=c(r1.dummies$coef, log(1/r1.dummies$theta)), 
                        X=r1.dummies$x, Y=r1.dummies$y,
                        method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  r1.dummies <- glm.nb(domAttacks~  lremitPC*lxrcomp
                       + llmil + llpop + lgdpgrowth + llgdppc + lpress 
                       +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                       +factor(ccode)+factor(year)-1-1, 
                       data=remitNEW.1a, subset=sam==1, maxit=25, x=TRUE,
                       start=r1.dummiesa$est[-length(r1.dummiesa$est)], init.theta=1/exp(r1.dummiesa$est[length(r1.dummiesa$est)]))
  
}
v1.dummies <- vcovCL(r1.dummies, remitNEW.1a[sam==1]$ccode)
df.margins <- as.data.frame(remitNEW.1a[as.numeric(rownames(r1.dummies$model)),])

P <- as.matrix(na.omit(unique(remitNEW.1a[sam==1,list(lxrcomp)])))
tests <- 
  cbind(P,
        r1.dummies$coef["lremitPC"] 
        +  P[,1] *r1.dummies$coef["lremitPC:lxrcomp"],
        t(sapply(paste("lremitPC+",
                       P[,1], "*lremitPC:lxrcomp"), 
                 function(x){
                   H <- linearHypothesis(r1.dummies, x, vcov=v1.dummies)
                   return(c(H$Chisq[2], H$Pr[2]))})))
tests <- as.data.frame(tests)
colnames(tests) <- c("lxrcomp", "coef","chi", "p")
tests <- tests[order(tests$lxrcomp),]
rownames(tests) <- NULL
tests$lxrcomp <- c("Closed selection", "No rules state","Transition", "Competitive")
tests$chi <- formatC(tests$chi,digits=2 ,format='f')
tests$coef <- formatC(tests$coef,digits=2 ,format='f')
tests$chi <- ifelse(tests$p < 0.05,
                    paste(tests$chi, "^{**}", sep=""),
                    ifelse(tests$p <0.1, 
                           paste(tests$chi, "^*", sep=""),
                           paste(tests$chi)))

tests$p <- NULL




r2.dummies <- glm(domAttacks~ lremitPC*lxrcomp +lremitPC*lxconst
                     + llmil + llpop + lgdpgrowth + llgdppc + lpress 
                  +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                  +factor(ccode)+factor(year)-1, family="poisson",
                     data=remitNEW.1a, subset=sam==1, maxit=25, x=TRUE)

while(r2.dummies$converged == FALSE| r2.dummies$family$family=="poisson"){
  r2.dummiesa <- maxLik(NBreg, NBgrad,
                        start=c(r2.dummies$coef, 0), 
                        X=r2.dummies$x, Y=r2.dummies$y,
                        method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  r2.dummies <- glm.nb(domAttacks~ lremitPC*lxrcomp + lremitPC*lxconst
                       + llmil + llpop + lgdpgrowth + llgdppc + lpress 
                       +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                       +factor(ccode)+factor(year)-1, 
                       data=remitNEW.1a, subset=sam==1, maxit=25, x=TRUE,
                       start=r2.dummiesa$est[-length(r2.dummiesa$est)], init.theta=1/exp(r2.dummiesa$est[length(r2.dummiesa$est)]))
  
}
v2.dummies <- vcovCL(r2.dummies, remitNEW.1a[sam==1]$ccode)


P2 <- expand.grid(lxrcomp=c(min(remitNEW.1a[sam==1]$lxrcomp),max(remitNEW.1a[sam==1]$lxrcomp)),
                  lxconst=c(median(remitNEW.1a[sam==1]$lxconst)))
tests2 <-
  cbind(P2,
        r2.dummies$coef["lremitPC"]
        +  P2[,1] *r2.dummies$coef["lremitPC:lxrcomp"]
        +  P2[,2] *r2.dummies$coef["lremitPC:lxconst"],
        t(sapply(paste("lremitPC+",
                       P2[,1], "*lremitPC:lxrcomp+",
                       P2[,2], "*lremitPC:lxconst"),
                 function(x){
                   H <- linearHypothesis(r2.dummies, x, vcov=v2.dummies)
                   return(c(H$Chisq[2], H$Pr[2]))})))
tests2 <- as.data.frame(tests2)
colnames(tests2) <- c("lxrcomp", "lxconst", "coef","chi", "p")
tests2 <- tests2[order(tests2$lxrcomp, tests2$lxconst),]
rownames(tests2) <- NULL
tests2$chi <- formatC(tests2$chi,digits=2 ,format='f')
tests2$coef <- formatC(tests2$coef,digits=2 ,format='f')
tests2$chi <- ifelse(tests2$p < 0.05,
                     paste(tests2$chi, "^{**}", sep=""),
                     ifelse(tests2$p <0.1,
                            paste(tests2$chi, "^*", sep=""),
                            paste(tests2$chi)))
tests2$p <- NULL





load("../output/vdem_competition.rdata")
glance_custom.negbin <- function(x, ...) {
  theta <- x$theta
  out <- data.frame("theta" = theta,
                    "logLik"=logLik(x))
  return(out)
}

cat(modelsummary(list("(5)"=r1.dummies,
                  "(6)"=r2.dummies, 
                  "(7)"=vdem1.dummies, 
                  "(8)"=vdem2.dummies),
             title = "Exploring the political competition mechanism \\label{tab:mechTable}",
             output="latex",
             fmt=2,
             vcov=list(v1.dummies,
                       v2.dummies,
                       vdemV1.dummies,
                       vdemV2.dummies),
             stars=c("*"=.1, "**"=0.05),
             coef_map=c("lremitPC"="Remittances",
                        "lremitPC:lxrcomp"="Remittances $\\times$ Competition",
                        "lremitPC:lv2x_partipdem"="Remittances $\\times$ Competition",
                        "lremitPC:lxconst"="Remittances $\\times$ Constraints",
                        "lremitPC:lv2xlg_legcon"="Remittances $\\times$ Constraints"),
             escape=FALSE,
             gof_map=data.frame(raw=c("nobs", "logLik", "theta"),
                                clean=c("Observations", "Log Likelihood", "$\\theta$"),
                                fmt=c(0,2,2)),
             align="ldddd",
             notes="\\\\footnotesize $^*p<0.1$, $^{**}p<0.05$. Coefficients from negative binomial models. Standard errors in parentheses, clustered on country. For the combined coefficients, the high (low) values are the 90th (10th) percentile of the competition measure. For models with constraints, the constraints variable is fixed to its median value when combined.  Coefficients for the control variables are suppressed for space."))
cat("\n")
## test on whether constraints explain anything
cat("Wald tests for model 6 v model 5\n")
print(linearHypothesis(r2.dummies, c("lremitPC:lxconst", "lxconst"), vcov=v2.dummies )) 
cat("Wald tests for model 8 v model 7\n")
print(linearHypothesis(vdem2.dummies, c("lremitPC:lv2xlg_legcon", "lv2xlg_legcon"), vcov=vdemV2.dummies )) 
cat("With polity, we fail to reject the null that constrains add nothing to the model\n")
cat("With V-Dem, we  reject the null that constrains add nothing to the model\n")


XR <- cbind(1,matrix(quantile(remitNEW.1a[sam==1]$lxrcomp,c(0.1,0.9)), ncol=1))
bR <- r1.dummies$coef[c("lremitPC", "lremitPC:lxrcomp")]
vR <- v1.dummies[c("lremitPC", "lremitPC:lxrcomp"), c("lremitPC", "lremitPC:lxrcomp")]
est1 <- XR %*% bR
var1 <- XR %*% vR %*% t(XR)
se1 <- sqrt(diag(var1))

XR <- cbind(1,
            matrix(quantile(remitNEW.1a[sam==1]$lxrcomp,c(0.1,0.9)), ncol=1),
            rep(quantile(remitNEW.1a[sam==1]$lxconst,c(.5)),2))
bR <- r2.dummies$coef[c("lremitPC", "lremitPC:lxrcomp",  "lremitPC:lxconst")]
vR <- v2.dummies[c("lremitPC", "lremitPC:lxrcomp","lremitPC:lxconst"),
                 c("lremitPC", "lremitPC:lxrcomp", "lremitPC:lxconst")]
est2 <- XR %*% bR
var2 <- XR %*% vR %*% t(XR)
se2 <- sqrt(diag(var2))


XR <- cbind(1,matrix(quantile(VdemDataset[sam==1]$lv2x_partipdem,c(0.1,0.9)), ncol=1))
bR <- vdem1.dummies$coef[c("lremitPC", "lremitPC:lv2x_partipdem")]
vR <- vdemV1.dummies[c("lremitPC", "lremitPC:lv2x_partipdem"), c("lremitPC", "lremitPC:lv2x_partipdem")]
estV1 <- XR %*% bR
varV1 <- XR %*% vR %*% t(XR)
seV1 <- sqrt(diag(varV1))


XR <- cbind(1,
            matrix(quantile(VdemDataset[sam==1]$lv2x_partipdem,c(0.1,0.9)), ncol=1),
            rep(quantile(VdemDataset[sam==1]$lv2xlg_legcon,c(.5)),2))
bR <- vdem2.dummies$coef[c("lremitPC", "lremitPC:lv2x_partipdem", "lremitPC:lv2xlg_legcon")]
vR <- vdemV2.dummies[c("lremitPC", "lremitPC:lv2x_partipdem","lremitPC:lv2xlg_legcon"),
                     c("lremitPC", "lremitPC:lv2x_partipdem", "lremitPC:lv2xlg_legcon")]
estV2 <- XR %*% bR
varV2 <- XR %*% vR %*% t(XR)
seV2 <- sqrt(diag(varV2))
cat("Combined coefficients for Table 5\n")
combined.betas <- 
  rbind(
    c(est1, est2, estV1, estV2),
    c(se1, se2, seV1, seV2),
    c(pnorm(abs(est1/se1), lower=FALSE)*2,
      pnorm(abs(est2/se2), lower=FALSE)*2,
      pnorm(abs(estV1/seV1), lower=FALSE)*2,
      pnorm(abs(estV2/seV2), lower=FALSE)*2))
combined.betas <- matrix(combined.betas,ncol=4)
rownames(combined.betas) <- rep(c("est", "s.e.", "p-value"),2)
print(round(combined.betas,2))





r3.dummies <- glm(domAttacks~ lremitPC*lxconst
                  + llmil + llpop + lgdpgrowth + llgdppc + lpress 
                  +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                  +factor(ccode)+factor(year)-1, family="poisson",
                  data=remitNEW.1a, subset=sam==1, maxit=25, x=TRUE)

while(r3.dummies$converged == FALSE| r3.dummies$family$family=="poisson"){
  r3.dummiesa <- maxLik(NBreg, NBgrad,
                        start=c(r3.dummies$coef, 0), 
                        X=r3.dummies$x, Y=r3.dummies$y,
                        method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  r3.dummies <- glm.nb(domAttacks~  lremitPC*lxconst
                       + llmil + llpop + lgdpgrowth + llgdppc + lpress 
                       +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                       +factor(ccode)+factor(year)-1, 
                       data=remitNEW.1a, subset=sam==1, maxit=25, x=TRUE,
                       start=r3.dummiesa$est[-length(r3.dummiesa$est)], init.theta=1/exp(r3.dummiesa$est[length(r3.dummiesa$est)]))
  
}
v3.dummies <- vcovCL(r3.dummies, remitNEW.1a[sam==1]$ccode)

# Vuong test (POLITY)
# Fails the first test but some evidence in 
# favor of r1 (competition)
cat("Vuong test for the polity measures: competition v. constraints\n")
vuongtest(r1.dummies, r3.dummies) 
cat("Inconclusive, but suggestive in favor of competition\n")

# model fit stats favor competition
cat("AIC selection for the polity measures: competition v. constraints\n")
which.min(c(AIC(r1.dummies),AIC(r3.dummies)))
cat("Selects compeition\n")



























