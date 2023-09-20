library(data.table)
library(car)
library(sandwich)
library(MASS)
library(maxLik)
library(splines)
library(modelsummary)
library(ggplot2)
rm(list=ls())
source("nbreg.R")


## For older drafts revert to pre 9/21/2021
load("../output/mainResultsRemittances.rdata")
load("../data/mainRemittanceData.rdata")

remitNEW[,`:=`(lvd1 = as.numeric(lvd==0),
               lvd2 = as.numeric(lvd==1),
               lvd3 = as.numeric(lvd==2),
               lvd4 = as.numeric(lvd==3))]



remitNEW[, remitPCT := (remittance_received_2010usd/gdp_2010usd)*100]
remitNEW[, lremitPCT := shift(remitPCT), by=ccode]
remitNEW[, remitPC_PCT := (remitPC/gdp_2010usd)*100]
remitNEW[, lremitPC_PCT := shift(remitPC_PCT), by=ccode]
##### Differerent measurements: Table B1 #####

## log ##

m1.log <- glm.nb(domAttacks~ log(lremitPC)*ldem + log(lremitPC)*lano 
                 + llmil + llpop + lgdpgrowth + llgdppc + lpress
                 +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                 + factor(ccode)+factor(year)-1,
                 data=remitNEW, subset=sam==1, maxit=25, x=TRUE)
while(m1.log$converged == FALSE){
  m1.loga <- maxLik(NBreg, NBgrad,
                    start=c(m1.log$coef, log(1/m1.log$theta)),
                    X=m1.log$x, Y=m1.log$y,
                    method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  m1.log <- glm.nb(domAttacks~ log(lremitPC)*ldem + log(lremitPC)*lano 
                   + llmil + llpop + lgdpgrowth + llgdppc + lpress
                   +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                   + factor(ccode)+factor(year)-1,
                   data=remitNEW, subset=sam==1, maxit=25, x=TRUE,
                   start=m1.loga$est[-length(m1.loga$est)], init.theta=1/exp(m1.loga$est[length(m1.loga$est)]))
  
}
v1.log <- vcovCL(m1.log, remitNEW[sam==1]$ccode)


## sqrt ##

m1.sqrt <- glm.nb(domAttacks~ sqrt(lremitPC)*ldem + sqrt(lremitPC)*lano 
                    + llmil + llpop + lgdpgrowth + llgdppc + lpress
                  +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                  + factor(ccode)+factor(year)-1,
                  data=remitNEW, subset=sam==1, maxit=25, x=TRUE)
while(m1.sqrt$converged == FALSE){
  m1.sqrta <- maxLik(NBreg, NBgrad,
                     start=c(m1.sqrt$coef, sqrt(1/m1.sqrt$theta)),
                     X=m1.sqrt$x, Y=m1.sqrt$y,
                     method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  m1.sqrt <- glm.nb(domAttacks~ sqrt(lremitPC)*ldem + sqrt(lremitPC)*lano
                      + llmil + llpop + lgdpgrowth + llgdppc + lpress
                    +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                    + factor(ccode)+factor(year)-1,
                    data=remitNEW, subset=sam==1, maxit=25, x=TRUE,
                    start=m1.sqrta$est[-length(m1.sqrta$est)], init.theta=1/exp(m1.sqrta$est[length(m1.sqrta$est)]))
  
}
v1.sqrt <- vcovCL(m1.sqrt, remitNEW[sam==1]$ccode)



## detrend ##
detrendModels <- list()
for(d in 1:3){
  detrended <- lm(remitPC~bs(year,degree=d)*factor(ccode)-1,  data=remitNEW)$resid
  remitNEW[as.numeric(names(detrended)),detrend:=detrended]
  
  remitNEW[,detrend:=shift(detrend), by=ccode]
  
  
  ## Randomly check some countries to see if it worked ###
  # plot.df <- melt(remitNEW[country %in% sample(remitNEW[sam==1]$country, size=4),
  #                          .(country,lremitPC, detrend, year)],
  #                 id.var=c("country", "year"))
  # print( ggplot(plot.df)+
  #          geom_line(aes(x=year, y=value))+
  #          stat_smooth(aes(x=year, y=value), method='lm', formula=y~x)+
  #          facet_grid(country~variable, scales="free"))
  
  m1.detrend <- glm(domAttacks~ detrend*ldem + detrend*lano 
                    + llmil + llpop + lgdpgrowth + llgdppc + lpress
                    +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                    + factor(ccode)+factor(year)-1, family="poisson",
                    data=remitNEW, subset=sam==1, maxit=25, x=TRUE)
  
  while(m1.detrend$converged == FALSE|| m1.detrend$family$family=="poisson"){
    m1.detrenda <- maxLik(NBreg, NBgrad,
                          start=c(m1.detrend$coef, 0),
                          X=m1.detrend$x, Y=m1.detrend$y,
                          method="NR", control=list(tol=1e-6, gradtol=1e-6))
    
    m1.detrend <- glm.nb(domAttacks~ detrend*ldem + detrend*lano 
                         + llmil + llpop + lgdpgrowth + llgdppc + lpress
                         +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                         + factor(ccode)+factor(year)-1,
                         data=remitNEW, subset=sam==1, maxit=25, x=TRUE,
                         start=m1.detrenda$est[-length(m1.detrenda$est)], init.theta=1/exp(m1.detrenda$est[length(m1.detrenda$est)]))
    
  }
  detrendModels[[d]] <- m1.detrend
  
}
k <- which.min(sapply(detrendModels, AIC))
cat("detrend model picks ", k, "\n")
### Degree 2 chosen based on BIC comparison for 1:5

m1.detrend <- detrendModels[[k]]
v1.detrend <- vcovCL(m1.detrend, remitNEW[sam==1]$ccode)

### remit as pct of gdp ## 

m1.gdp.dummies <- glm.nb(domAttacks~ log(lremitPCT)*ldem +
                           log(lremitPCT)*lano
                         + llmil + llpop + lgdpgrowth + llgdppc + lpress
                         +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                         +factor(ccode)+factor(year)-1,
                         data=remitNEW, subset=sam==1, maxit=25, x=TRUE)
while(m1.gdp.dummies$converged == FALSE){
  m1.gdpA.dummies <- maxLik(NBreg, NBgrad,
                            start=c(m1.gdp.dummies$coef, log(1/m1.gdp.dummies$theta)),
                            X=m1.gdp.dummies$x, Y=m1.gdp.dummies$y,
                            method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  m1.gdp.dummies <- glm.nb(domAttacks~log(lremitPCT)*ldem +
                             log(lremitPCT)*lano  
                             + llmil + llpop + lgdpgrowth + llgdppc + lpress
                           +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                           +factor(ccode)+factor(year)-1,
                           data=remitNEW, subset=sam==1, maxit=25, x=TRUE,
                           start=m1.gdpA.dummies$est[-length(m1.gdpA.dummies$est)], init.theta=1/exp(m1.gdpA.dummies$est[length(m1.gdpA.dummies$est)]))
  
}
v1.gdp.dummies <- vcovCL(m1.gdp.dummies, remitNEW[sam==1]$ccode)


### remit/capita as pct of gdp ## 
m1.gdp.dummies2 <- glm.nb(domAttacks~ log(lremitPC_PCT)*ldem +
                           log(lremitPC_PCT)*lano
                         + llmil + llpop + lgdpgrowth + llgdppc + lpress
                         +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                         +factor(ccode)+factor(year)-1,
                         data=remitNEW, subset=sam==1, maxit=25, x=TRUE)
while(m1.gdp.dummies2$converged == FALSE){
  m1.gdpA.dummies2 <- maxLik(NBreg, NBgrad,
                            start=c(m1.gdp.dummies2$coef, log(1/m1.gdp.dummies2$theta)),
                            X=m1.gdp.dummies2$x, Y=m1.gdp.dummies2$y,
                            method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  m1.gdp.dummies2 <- glm.nb(domAttacks~log(lremitPC_PCT)*ldem +
                             log(lremitPC_PCT)*lano  
                           + llmil + llpop + lgdpgrowth + llgdppc + lpress
                           +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                           +factor(ccode)+factor(year)-1,
                           data=remitNEW, subset=sam==1, maxit=25, x=TRUE,
                           start=m1.gdpA.dummies2$est[-length(m1.gdpA.dummies2$est)], 
                           init.theta=1/exp(m1.gdpA.dummies2$est[length(m1.gdpA.dummies$est2)]))
  
}
v1.gdp.dummies2 <- vcovCL(m1.gdp.dummies2, remitNEW[sam==1]$ccode)



## IMF ##
m1.imf <- glm(domAttacks~ limf.remitPC*ldem + limf.remitPC*lano 
              + llmil + llpop + lgdpgrowth + llgdppc + lpress
              +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
              + factor(ccode)+factor(year)-1,  family="poisson",
              data=remitNEW, subset=sam.imf==1, maxit=25, x=TRUE)
while(m1.imf$converged == FALSE| m1.imf$family$family=="poisson"){
  m1.imfa <- maxLik(NBreg, NBgrad,
                    start=c(m1.imf$coef, 0),
                    X=m1.imf$x, Y=m1.imf$y,
                    method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  m1.imf <- glm.nb(domAttacks~ limf.remitPC*ldem + limf.remitPC*lano 
                   + llmil + llpop + lgdpgrowth + llgdppc + lpress
                   +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                   + factor(ccode)+factor(year)-1,
                   data=remitNEW, subset=sam.imf==1, maxit=25, x=TRUE,
                   start=m1.imfa$est[-length(m1.imfa$est)], init.theta=1/exp(m1.imfa$est[length(m1.imfa$est)]))
  
}
v1.imf <- vcovCL(m1.imf, remitNEW[sam.imf==1]$ccode)



## Comparing main to other specifications
which.min(c(AIC(m2.dummies),
            AIC(m1.log),
            AIC(m1.sqrt),
            AIC(m1.detrend),
            AIC(m1.gdp.dummies),
            AIC(m1.gdp.dummies2)))
## OG is the winner here 


## Check the main against IMF with the same samples
m2.refit.imf <- update(m2.dummies, start=NULL, subset=sam.imf==1 & sam==1)
m.refit.imf <- update(m1.imf, start=NULL, subset=sam.imf==1 & sam==1)
all(names(m2.refit.imf$resid) == names(m.refit.imf$resid) )
which.min(c(AIC(m2.refit.imf),AIC(m.refit.imf)))
## OG specification not preferred, but lose a lot of data (more than 20%)


glance_custom.negbin <- function(x, ...) {
  theta <- x$theta
  out <- data.frame("theta" = theta)
  return(out)
}







est <- c(m2.dummies$coefficients["lremitPC"]+m2.dummies$coefficients["lremitPC:ldem"],
         m1.log$coefficients["log(lremitPC)"]+m1.log$coefficients["log(lremitPC):ldem"],
         m1.sqrt$coefficients["sqrt(lremitPC)"]+m1.sqrt$coefficients["sqrt(lremitPC):ldem"],
         m1.detrend$coefficients["detrend"]+m1.detrend$coefficients["detrend:ldem"],
         m1.gdp.dummies$coefficients["log(lremitPCT)"]+m1.gdp.dummies$coefficients["log(lremitPCT):ldem"],
         m1.gdp.dummies2$coefficients["log(lremitPC_PCT)"]+m1.gdp.dummies2$coefficients["log(lremitPC_PCT):ldem"],
         m1.imf$coefficients["limf.remitPC"]+m1.imf$coefficients["limf.remitPC:ldem"])
se <- sqrt(c(v2.dummies["lremitPC","lremitPC"] + v2.dummies["lremitPC:ldem","lremitPC:ldem"]+2* v2.dummies["lremitPC","lremitPC:ldem"],
             v1.log["log(lremitPC)","log(lremitPC)"] + v1.log["log(lremitPC):ldem","log(lremitPC):ldem"]+2* v1.log["log(lremitPC)","log(lremitPC):ldem"],
             v1.sqrt["sqrt(lremitPC)","sqrt(lremitPC)"] + v1.sqrt["sqrt(lremitPC):ldem","sqrt(lremitPC):ldem"]+2* v1.sqrt["sqrt(lremitPC)","sqrt(lremitPC):ldem"],
             v1.detrend["detrend", "detrend"] + v1.detrend["detrend:ldem","detrend:ldem"]+2* v1.detrend["detrend","detrend:ldem"],
             v1.gdp.dummies["log(lremitPCT)","log(lremitPCT)"] + v1.gdp.dummies["log(lremitPCT):ldem","log(lremitPCT):ldem"]+2* v1.gdp.dummies["log(lremitPCT)","log(lremitPCT):ldem"],
             v1.gdp.dummies2["log(lremitPC_PCT)","log(lremitPC_PCT)"] + v1.gdp.dummies2["log(lremitPC_PCT):ldem","log(lremitPC_PCT):ldem"]+2* v1.gdp.dummies2["log(lremitPC_PCT)","log(lremitPC_PCT):ldem"],
             v1.imf["limf.remitPC", "limf.remitPC"] + v1.imf["limf.remitPC:ldem","limf.remitPC:ldem"]+2* v1.imf["limf.remitPC","limf.remitPC:ldem"]))
p.vals <- pnorm(abs(est/se), lower=F)*2

betas <- rbind(c("$\\hat{\\beta}_{\\text{Remittances}}  + \\hat{\\beta}_{\\text{Remittances $\\times$ Democracy}}$",
                 paste0(formatC(est,2, format="f"), ifelse(p.vals < .05, "^{**}", ifelse(p.vals < .1, "^*", "")))),
               c("", paste0("(", formatC(se,2, format="f"), ")")))


rows <- rbind.data.frame(betas, 
                         c("\\midrule Country fixed effects", rep("\\multicolumn{1}{c}{Yes}", 7)),
                         c("Year fixed effects", rep("\\multicolumn{1}{c}{Yes}", 7)),
                         c("Controls", rep("\\multicolumn{1}{c}{Yes}", 7)))
attr(rows, 'position') <- c(11:15)

cat("Table D1\n")
cat(modelsummary(list("Remittances/capita (Model 2)"=m2.dummies,
                  "Log remittances/capita"=m1.log,
                  "Sqrt remittances/capita" = m1.sqrt,
                  "Detrended remittances/capita"= m1.detrend,
                  "Log Remittances/GDP"=m1.gdp.dummies,
                  "Log (Remittances/capita)/GDP"=m1.gdp.dummies2,
                  "IMF Remittances/capita"=m1.imf),
             title = "Negative binomials with different measures and transformations of remittances \\label{tab:diff.remit}",
             output="latex",
             fmt=2,
             vcov=list(v2.dummies,
                       v1.log,
                       v1.sqrt,
                       v1.detrend,
                       v1.gdp.dummies,
                       v1.gdp.dummies2,
                       v1.imf),
             stars=c("*"=.1, "**"=0.05),
             coef_map=c("lremitPC"="Remittances",
                        "log(lremitPC)"="Remittances",
                        "sqrt(lremitPC)"="Remittances",
                        "detrend"="Remittances",
                        "log(lremitPCT)"="Remittances",
                        "log(lremitPC_PCT)"="Remittances",
                        "limf.remitPC"="Remittances",
                        "lremitPC:ldem" ="Remittances $\\times$ Dem.",
                        "log(lremitPC):ldem" ="Remittances $\\times$ Dem.",
                        "sqrt(lremitPC):ldem" ="Remittances $\\times$ Dem.",
                        "detrend:ldem" ="Remittances $\\times$ Dem.",
                        "log(lremitPCT):ldem" ="Remittances $\\times$ Dem.",
                        "log(lremitPC_PCT):ldem" ="Remittances $\\times$ Dem.",
                        "limf.remitPC:ldem" ="Remittances $\\times$ Dem."),
             escape=FALSE,
             gof_map=data.frame(raw=c("nobs", "logLik", "theta"),
                                clean=c("Observations", "Log Likelihood", "$\\theta$"),
                                fmt=c(0,2,2)),
             add_rows = rows,
             align="lddddddd",
             notes="\\\\footnotesize $^*p<0.1$, $^{**}p<0.05$. Regression coefficients. Standard errors in parentheses, clustered on country."))

             

cat("Combined coefs for Table D1\n")
round(rbind(est,se,p.vals), 2)



