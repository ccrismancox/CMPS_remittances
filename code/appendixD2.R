library(data.table)
library(car)
library(sandwich)
library(MASS)
library(maxLik)
library(modelsummary)
rm(list=ls())
source("nbreg.R")


load("../output/mainResultsRemittances.rdata")
load("../data/mainRemittanceData.rdata")

remitNEW[,`:=`(lvd1 = as.numeric(lvd==0),
               lvd2 = as.numeric(lvd==1),
               lvd3 = as.numeric(lvd==2),
               lvd4 = as.numeric(lvd==3))]

###### Linear polity #####
mpolity.dummies <- glm.nb(domAttacks~ lremitPC*lpolity2 + llmil + llpop + lgdpgrowth + llgdppc + lpress
                          +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                          + factor(ccode) + factor(year)-1,
                     data=remitNEW, subset=sam==1,  maxit=25, x=TRUE)
while(mpolity.dummies$converged == FALSE |!is.null(mpolity.dummies$th.warn) ){
  mpolity.dummiesa <- maxLik(NBreg, NBgrad,
                        start=c(mpolity.dummies$coef, log(1/mpolity.dummies$theta)),
                        X=mpolity.dummies$x, Y=mpolity.dummies$y,
                        method="NR", control=list(tol=1e-6, gradtol=1e-6))

  mpolity.dummies <- glm.nb(domAttacks~  lremitPC*lpolity2+ llmil + llpop + lgdpgrowth + llgdppc + lpress
                            +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                            + factor(ccode) + factor(year)-1,
                       data=remitNEW[sam==1],maxit=25, x=TRUE,
                       start=mpolity.dummiesa$est[-length(mpolity.dummiesa$est)],
                       init.theta=1/exp(mpolity.dummiesa$est[length(mpolity.dummiesa$est)]))

}
vpolity.dummies <- vcovCL(mpolity.dummies, remitNEW[sam==1]$ccode)



#### Polity quadratic #####

#### polity quad####
remitNEW[,lpolity2.sq := lpolity2^2]
m1.polity.sq <- glm(domAttacks~ lremitPC*lpolity2 + lremitPC*I(lpolity2^2)
                    + llmil + llpop + lgdpgrowth + llgdppc + lpress
                    +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                    + factor(ccode)+factor(year)  -1,
                    data=remitNEW, subset=sam==1, maxit=25,x=TRUE, family="poisson")
while(m1.polity.sq$converged == FALSE| m1.polity.sq$family$family=="poisson"){
  m1.polity.sqA <- maxLik(NBreg, NBgrad,
                          start=c(m1.polity.sq$coef, 0),
                          X=m1.polity.sq$x, Y=m1.polity.sq$y,
                          method="NR", control=list(tol=1e-6, gradtol=1e-6))

  m1.polity.sq <- glm.nb(domAttacks~ lremitPC*lpolity2 + lremitPC*I(lpolity2^2)
                         + llmil + llpop + lgdpgrowth + llgdppc + lpress
                         +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                         + factor(ccode)+factor(year)  -1,
                         data=remitNEW, subset=sam==1, maxit=25, x=TRUE,
                         start=m1.polity.sqA$est[-length(m1.polity.sqA$est)], init.theta=1/exp(m1.polity.sqA$est[length(m1.polity.sqA$est)]))

}
v1.polity.sq <- vcovCL(m1.polity.sq, remitNEW[sam==1]$ccode)




#### polity alternative cut points #####

#### using 5 as the cutpoint
remitNEW[,ldem5 := ifelse(lpolity2 >= 5, 1,0)]
remitNEW[,lano5 := ifelse(lpolity2 > -5 & lpolity2 < 5, 1,0)]
m1.dummies5 <- glm.nb(domAttacks~ lremitPC*ldem5 + lremitPC*lano5 + llmil + llpop +
                          lgdpgrowth + llgdppc + lpress
                      +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                      + factor(ccode)+factor(year)-1,
                      data=remitNEW, subset=sam==1, maxit=25, x=TRUE)
while(m1.dummies5$converged == FALSE){
    m1.dummies5a <- maxLik(NBreg, NBgrad,
                           start=c(m1.dummies5$coef, log(1/m1.dummies5$theta)),
                           X=m1.dummies5$x, Y=m1.dummies5$y,
                           method="NR", control=list(tol=1e-6, gradtol=1e-6))

    m1.dummies5 <- glm.nb(domAttacks~ lremitPC*ldem5 + lremitPC*lano5
                          + llmil + llpop + lgdpgrowth + llgdppc + lpress
                          +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                          + factor(ccode)+factor(year)-1,
                          data=remitNEW, subset=sam==1, maxit=25, x=TRUE,
                          start=m1.dummies5a$est[-length(m1.dummies5a$est)], init.theta=1/exp(m1.dummies5a$est[length(m1.dummies5a$est)]))

}
v1.dummies5 <- vcovCL(m1.dummies5, remitNEW[sam==1]$ccode)

#### using 6 as the cutpoint
remitNEW[,ldem6 := ifelse(lpolity2 >= 6, 1,0)]
remitNEW[,lano6 := ifelse(lpolity2 > -6 & lpolity2 < 6, 1,0)]
m1.dummies6 <- glm.nb(domAttacks~ lremitPC*ldem6 + lremitPC*lano6 + llmil + llpop +
                        lgdpgrowth + llgdppc + lpress
                      +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                      + factor(ccode)+factor(year)-1,
                      data=remitNEW, subset=sam==1, maxit=25, x=TRUE)
while(m1.dummies6$converged == FALSE){
  m1.dummies6a <- maxLik(NBreg, NBgrad,
                         start=c(m1.dummies6$coef, log(1/m1.dummies6$theta)),
                         X=m1.dummies6$x, Y=m1.dummies6$y,
                         method="NR", control=list(tol=1e-6, gradtol=1e-6))

  m1.dummies6 <- glm.nb(domAttacks~ lremitPC*ldem6 + lremitPC*lano6
                        + llmil + llpop + lgdpgrowth + llgdppc + lpress
                        +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                        + factor(ccode)+factor(year)-1,
                        data=remitNEW, subset=sam==1, maxit=25, x=TRUE,
                        start=m1.dummies6a$est[-length(m1.dummies6a$est)], init.theta=1/exp(m1.dummies6a$est[length(m1.dummies6a$est)]))

}
v1.dummies6 <- vcovCL(m1.dummies6, remitNEW[sam==1]$ccode)




m1.vdem <- glm(domAttacks~ lremitPC*lvd2 + lremitPC*I(lvd3 +lvd4)
               + llmil + llpop + lgdpgrowth + llgdppc + lpress
               +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
               + factor(ccode)+factor(year)-1,
               data=remitNEW, subset=sam==1, maxit=25,x=TRUE, family=poisson)
while(m1.vdem$converged == FALSE |  m1.vdem$family$family=="poisson"){
  m1.vdemA <- maxLik(NBreg, NBgrad,
                     start=c(m1.vdem$coef,0),
                     X=m1.vdem$x, Y=m1.vdem$y,
                     method="NR", control=list(tol=1e-6, gradtol=1e-6))

  m1.vdem <- glm.nb(domAttacks~ lremitPC*lvd2 + lremitPC*I(lvd3 +lvd4)
                    + llmil + llpop + lgdpgrowth + llgdppc + lpress
                    +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                    + factor(ccode)+factor(year)-1,
                    data=remitNEW, subset=sam==1, maxit=25, x=TRUE,
                    start=m1.vdemA$est[-length(m1.vdemA$est)],
                    init.theta=1/exp(m1.vdemA$est[length(m1.vdemA$est)]))

}
v1.vdem <- vcovCL(m1.vdem, remitNEW[sam==1]$ccode)



#### V-Dem's Electoral democracy index ####
m1.vdem2 <- glm(domAttacks~ lremitPC*lpolyarch +
                    + llmil + llpop + lgdpgrowth + llgdppc + lpress
                +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                +factor(ccode)+factor(year)-1,
                data=remitNEW, subset=sam==1, maxit=25,x=TRUE, family=poisson)
while(m1.vdem2$converged == FALSE |  m1.vdem2$family$family=="poisson"){
  m1.vdem2A <- maxLik(NBreg, NBgrad,
                      start=c(m1.vdem2$coef,0),
                      X=m1.vdem2$x, Y=m1.vdem2$y,
                      method="NR", control=list(tol=1e-6, gradtol=1e-6))

  m1.vdem2 <- glm.nb(domAttacks~ lremitPC*lpolyarch +
                         + llmil + llpop + lgdpgrowth + llgdppc + lpress
                     +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                     +factor(ccode)+factor(year)-1,
                     data=remitNEW, subset=sam==1, maxit=25, x=TRUE,
                     start=m1.vdem2A$est[-length(m1.vdem2A$est)],
                     init.theta=1/exp(m1.vdem2A$est[length(m1.vdem2A$est)]))

}
v1.vdem2 <- vcovCL(m1.vdem2, remitNEW[sam==1]$ccode)







RHO <- cor(remitNEW[sam==1,
             .(2*ldem+lano, 
               2*ldem6 + lano6,
               2*ldem5 + lano5,
               lpolity2,
               lvd, 
               lpolyarch)], method = "spearman")
cat("Correlation among democracy measures\n")
print(max(RHO[lower.tri(RHO)]))
print(min(RHO[lower.tri(RHO)]))
glance_custom.negbin <- function(x, ...) {
  theta <- x$theta
  out <- data.frame("theta" = theta)
  return(out)
}



Combined.auto <- rbind(deltaMethod(m2.dummies,"lremitPC", v2.dummies),
                  deltaMethod(m1.dummies6,"lremitPC", v1.dummies6),
                  deltaMethod(m1.dummies5,"lremitPC", v1.dummies5),
                  deltaMethod(mpolity.dummies,"lremitPC-9*`lremitPC:lpolity2`", vpolity.dummies),
                  deltaMethod(m1.polity.sq,"lremitPC-9*`lremitPC:lpolity2`+81*`lremitPC:I(lpolity2^2)`", v1.polity.sq),
                  deltaMethod(m1.vdem,"lremitPC", v1.vdem),
                  deltaMethod(m1.vdem2,"lremitPC+.14*`lremitPC:lpolyarch`", v1.vdem2)
              )
est <- Combined.auto[,1]
se <- Combined.auto[,2]
p.vals <- pnorm(abs(est/se), lower=F)*2




betasA <- rbind(c("Combined coefficient on remittances",
                 paste0(formatC(est,2, format="f"), ifelse(p.vals < .05, "^{**}", ifelse(p.vals < .1, "^*", "")))),
               c("(Autocracy)", paste0("(", formatC(se,2, format="f"), ")")))
Combined <- rbind(deltaMethod(m2.dummies,"lremitPC+`lremitPC:ldem`", v2.dummies),
                  deltaMethod(m1.dummies6,"lremitPC+`lremitPC:ldem6`", v1.dummies6),
                  deltaMethod(m1.dummies5,"lremitPC+`lremitPC:ldem5`", v1.dummies5),
                  deltaMethod(mpolity.dummies,"lremitPC+9*`lremitPC:lpolity2`", vpolity.dummies),
                  deltaMethod(m1.polity.sq,"lremitPC+9*`lremitPC:lpolity2`+81*`lremitPC:I(lpolity2^2)`", v1.polity.sq),
                  deltaMethod(m1.vdem,"lremitPC+`lremitPC:I(lvd3 + lvd4)`", v1.vdem),
                  deltaMethod(m1.vdem2,"lremitPC+.88*`lremitPC:lpolyarch`", v1.vdem2)
)
est <- Combined[,1]
se <- Combined[,2]
p.vals <- pnorm(abs(est/se), lower=F)*2

betasD <- rbind(c("Combined coefficient on remittances",
                  paste0(formatC(est,2, format="f"), ifelse(p.vals < .05, "^{**}", ifelse(p.vals < .1, "^*", "")))),
                c("(Democracy)", paste0("(", formatC(se,2, format="f"), ")")))

rows <- rbind.data.frame(betasA, betasD ,
                         c("\\midrule Country fixed effects", rep("\\multicolumn{1}{c}{Yes}", 7)),
                         c("Year fixed effects", rep("\\multicolumn{1}{c}{Yes}", 7)),
                         c("Controls", rep("\\multicolumn{1}{c}{Yes}", 7)))
attr(rows, 'position') <- c(9:16)

cat("Table D2\n")
cat(modelsummary(list("\\begin{tabular}{c}Polity \\\\ cutpoints at $\\pm 7$\\\\ (Model 2)\\end{tabular}"=m2.dummies,
                "\\begin{tabular}{c}Polity \\\\ cutpoints at $\\pm 6$\\end{tabular}"=m1.dummies6,
                  "\\begin{tabular}{c}Polity \\\\ cutpoints at $\\pm 5$\\end{tabular}" = m1.dummies5,
                  "\\begin{tabular}{c}Polity \\\\ linear \\end{tabular}"= mpolity.dummies,
                  "\\begin{tabular}{c}Polity \\\\ quadratic \\end{tabular}"=m1.polity.sq,
                  "\\begin{tabular}{c}V-DEM \\\\ categorical \\end{tabular}"=m1.vdem,
                "\\begin{tabular}{c}V-DEM \\\\ Polyarchy index \\end{tabular}"=m1.vdem2
                ),
             title = "Negative binomials with different measures of democracy \\label{tab:diff.dem}",
             output="latex",
             fmt=2,
             vcov=list(v2.dummies,
                       v1.dummies6,
                       v1.dummies5,
                       vpolity.dummies,
                       v1.polity.sq,
                       v1.vdem,
                       v1.vdem2),
             stars=c("*"=.1, "**"=0.05),
             coef_map=c("lremitPC"="Remittances",
                        "lremitPC:ldem" ="Remittances $\\times$ Dem.",
                        "lremitPC:ldem6" ="Remittances $\\times$ Dem.",
                        "lremitPC:ldem5" ="Remittances $\\times$ Dem.",
                        "lremitPC:lpolity2" ="Remittances $\\times$ Dem.",
                        "lremitPC:I(lpolity2^2)" ="Remittances $\\times$ Dem. sq.",
                        "lremitPC:I(lvd3 + lvd4)" ="Remittances $\\times$ Dem.",
                        "lremitPC:lpolyarch" ="Remittances $\\times$ Dem.",
                        "lremitPC:lvd2" ="Remittances $\\times$ Electoral Autocracy"),
             escape=FALSE,
             gof_map=data.frame(raw=c("nobs", "logLik", "theta"),
                                clean=c("Observations", "Log Likelihood", "$\\theta$"),
                                fmt=c(0,2,2)),
             add_rows = rows,
             align="lddddddd",
             notes="\\\\footnotesize $^*p<0.1$, $^{**}p<0.05$. Regression coefficients. Standard errors in parentheses, clustered on country."))

cat("\n")

