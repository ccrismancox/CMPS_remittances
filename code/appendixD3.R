library(data.table)
library(car)
library(sandwich)
library(MASS)
library(pscl)
library(maxLik)
library(lme4)
library(stargazer)

rm(list=ls())
source("stargazerNoteCorrection.r")
source("nbreg.R")


load("../output/mainResultsRemittances.rdata")
load("../data/mainRemittanceData.rdata")

remitNEW[,`:=`(lvd1 = as.numeric(lvd==0),
               lvd2 = as.numeric(lvd==1),
               lvd3 = as.numeric(lvd==2),
               lvd4 = as.numeric(lvd==3))]


############ ALT Samples/DV ##############
### NOT OECD ##
oecd <- fread("../data/oecdmembers.csv")
oecd[,ccode := countrycode::countrycode(country, "country.name", "cown")]
oecd[,date:=as.Date(date, "%d %B %Y")]
oecd[, member71 := ifelse(date <= 1971, 1,0)]
oecd[, `:=`(country=NULL, date =NULL)]
remitNEW <- merge(remitNEW, oecd, by="ccode", all.x=TRUE)
remitNEW[, member71:= ifelse(is.na(member71), 0, member71)]
# with(remitNEW[sam==1], table(country, member71)) ## check that it worked

m1.notOECD <- glm.nb(domAttacks~ lremitPC*ldem + lremitPC*lano + llmil + llpop + lgdpgrowth + llgdppc + lpress
                     +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                     + factor(ccode)+factor(year)-1,
                     data=remitNEW, subset=sam==1 & member71 != 1, maxit=25, x=TRUE)
while(m1.notOECD$converged == FALSE){
  m1.notOECDa <- maxLik(NBreg, NBgrad,
                        start=c(m1.notOECD$coef, log(1/m1.notOECD$theta)),
                        X=m1.notOECD$x, Y=m1.notOECD$y,
                        method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  m1.notOECD <- glm.nb(domAttacks~ lremitPC*ldem + lremitPC*lano + llmil + llpop + lgdpgrowth + llgdppc + lpress
                       +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                       + factor(ccode)+factor(year)-1,
                       data=remitNEW,subset=sam==1 & member71 != 1, maxit=25, x=TRUE,
                       start=m1.notOECDa$est[-length(m1.notOECDa$est)], init.theta=1/exp(m1.notOECDa$est[length(m1.notOECDa$est)]))
  
}
v1.notOECD <- vcovCL(m1.notOECD, remitNEW[sam==1&member71 != 1]$ccode)



#### different DV ####


m1.esg.dummies <- glm.nb(esgAttacks~ lremitPC*ldem + lremitPC*lano + llmil + llpop + lgdpgrowth + llgdppc + lpress
                         +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                         + factor(ccode)+factor(year)-1,
                         data=remitNEW[sam.esg==1], subset=, maxit=25, x=TRUE)
while(m1.esg.dummies$converged == FALSE){
  m1.esgA.dummies <- maxLik(NBreg, NBgrad,
                            start=c(m1.esg.dummies$coef, log(1/m1.esg.dummies$theta)),
                            X=m1.esg.dummies$x, Y=m1.esg.dummies$y,
                            method="NR", control=list(tol=1e-6, gradtol=1e-6))
  
  m1.esg.dummies <- glm.nb(esgAttacks~ lremitPC*ldem + lremitPC*lano + llmil + llpop + lgdpgrowth + llgdppc + lpress
                           +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                           + factor(ccode)+factor(year)-1,
                           data=remitNEW[sam.esg==1], maxit=25, x=TRUE,
                           start=m1.esgA.dummies$est[-length(m1.esgA.dummies$est)], init.theta=1/exp(m1.esgA.dummies$est[length(m1.esgA.dummies$est)]))
  
}
v1.esg.dummies <- vcovCL(m1.esg.dummies, remitNEW[sam.esg==1]$ccode)



cat("Table D3\n")

stargazer(m1.notOECD,m1.esg.dummies,
          title = "Alternative dependent variables and samples",
          label="tab:alt_sam",
          model.names = TRUE, model.numbers = TRUE,
          omit=c("factor*", "*.bar",
                 "lremitPC:lano", "lano", "^ldem$",
                 "llmil", "llpop", "lgdpgrowth", "llgdppc", "lpress", 
                 "lgini", "lmaxlowx", "lwexcluded", "lnum_ongoing_ucdp"),
          dep.var.labels = c("Domestic Attacks (main)", "Domestic Attacks (ESG)"),
          order=c("lremitPC","lremitPC:ldem"),
          se=list(sqrt(diag(v1.notOECD)), sqrt(diag(v1.esg.dummies))),
          align=TRUE,
          no.space=TRUE,
          digits=2,
          omit.stat = c("AIC", "BIC", "f","rsq", "ser"),
          star.cutoffs = c(0.1, .05, NA),
          covariate.labels = c("Remittances per capita",
                               "Remittances per capita $\\times$ dem."),
          notes=c("$p<[0.*], p<[0.**]$. Coefficients from negative binomial models. Clustered standard errors in parentheses."),
          notes.append=FALSE
)

est <- c(m1.notOECD$coefficients['lremitPC']+m1.notOECD$coefficients['lremitPC:ldem'],
         m1.esg.dummies$coefficients['lremitPC']+m1.esg.dummies$coefficients['lremitPC:ldem'])
se <- sqrt(c(v1.notOECD['lremitPC','lremitPC']+
               v1.notOECD['lremitPC:ldem','lremitPC:ldem']+
               2*v1.notOECD['lremitPC','lremitPC:ldem'],
             v1.esg.dummies['lremitPC', 'lremitPC'] + v1.esg.dummies['lremitPC:ldem','lremitPC:ldem']+
               2* v1.esg.dummies['lremitPC', 'lremitPC:ldem'] ))
p.vals <- pnorm(abs(est/se), lower=F)*2
cat("Combined coefs, se's, and p-vals for Table D3\n")
print(round(rbind(est,se,p.vals), 2))


############## Structure #####################
##### Zero inflation ####


##Try it the Drakos and Gofas way##
remit.zero <- copy(remitNEW)
remit.zero <- subset(remit.zero, sam==1)

logit.cols <- c("ldem","lano","lpress", "lag_domAttacks", "lpolity2")
remit.zero[, (paste0(logit.cols,".bar")) := lapply(.SD, mean, na.rm=T), by=.(ccode),
           .SDcols=logit.cols]
m1.zero.dummies2 <- zeroinfl(domAttacks~ lremitPC*ldem + lremitPC*lano
                             + llmil + llpop + lgdpgrowth + llgdppc + lpress
                             +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                             + factor(ccode)+factor(year)-1|
                               ldem+lano+ lpress+ lag_domAttacks+
                               ldem.bar+ lano.bar+lpress.bar +lag_domAttacks.bar
                             ,
                             data=remit.zero, dist="negbin", x=TRUE)
df.margins <- as.data.frame(remitNEW[as.numeric(rownames(m1.zero.dummies2$model)),])
v1.zero <- vcovCL(m1.zero.dummies2, remit.zero$ccode)






#### not fixed effects ####


m1.dummies.pooled <- glm.nb(domAttacks~ lremitPC*ldem + lremitPC*lano + llmil +
                              llpop + lgdpgrowth + llgdppc + lpress
                            +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                            ,
                            data=remitNEW, maxit=25, x=TRUE)
v1.dummies.pooled <- vcovCL(m1.dummies.pooled, remitNEW$ccode)


m1.dummies.re <- glmer.nb(domAttacks~ lremitPC*ldem + lremitPC*lano +
                            llmil + llpop + lgdpgrowth + llgdppc + lpress
                          +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                          +(1|ccode),
                          data=remitNEW)




remitNEW.cre <- copy(remitNEW)

remitNEW.cre[,lremitPC.ldem := lremitPC*ldem]
remitNEW.cre[,lremitPC.lano := lremitPC*lano]
remitNEW.cre[,lremitPC.bar := mean(lremitPC,na.rm=T), by=ccode]
remitNEW.cre[,ldem.bar := mean(ldem,na.rm=T), by=ccode]
remitNEW.cre[,lano.bar := mean(lano,na.rm=T), by=ccode]
remitNEW.cre[,lremitPC.ldem.bar := mean(lremitPC.ldem,na.rm=T), by=ccode]
remitNEW.cre[,lremitPC.lano.bar := mean(lremitPC.lano,na.rm=T), by=ccode]
remitNEW.cre[,llmil.bar := mean(llmil,na.rm=T), by=ccode]
remitNEW.cre[,llpop.bar := mean(llpop,na.rm=T), by=ccode]
remitNEW.cre[,lgdpgrowth.bar := mean(lgdpgrowth,na.rm=T), by=ccode]
remitNEW.cre[,lpress.bar := mean(lpress,na.rm=T), by=ccode]
remitNEW.cre[,llgdppc.bar := mean(llgdppc,na.rm=T), by=ccode]
remitNEW.cre[,lmaxlowx.bar := mean(lmaxlowx,na.rm=T), by=ccode]
remitNEW.cre[,lgini.bar := mean(lgini,na.rm=T), by=ccode]
remitNEW.cre[,lnum_ongoing_ucdp.bar := mean(lnum_ongoing_ucdp,na.rm=T), by=ccode]
remitNEW.cre[,lwexcluded.bar := mean(lwexcluded,na.rm=T), by=ccode]





m1.dummies.cre <- glm(domAttacks~ lremitPC*ldem + lremitPC*lano
                      + llmil + llpop     + lgdpgrowth + llgdppc + lpress+
                        +lgini +  lmaxlowx + lwexcluded+ lnum_ongoing_ucdp 
                      +factor(year)-1
                      +lremitPC.bar+ldem.bar+lano.bar+
                        +lremitPC.ldem.bar+lremitPC.lano.bar+
                        llmil.bar+ llpop.bar+lgdpgrowth.bar+lpress.bar+llgdppc.bar
                      +lgini.bar+  lmaxlowx.bar + lwexcluded.bar+ lnum_ongoing_ucdp.bar, 
                      data=remitNEW.cre, family="poisson", x=T)
v1.dummies.cre <- vcovCL(m1.dummies.cre, cluster = remitNEW.cre$ccode)
df.margins <- m1.dummies.cre$model
colnames(df.margins)[grep(colnames(df.margins), pattern="factor")] <- "year"

mar.cre <- summary(margins(m1.dummies.cre,
                        "lremitPC",
                        vcov=v1.dummies.cre,
                        data = df.margins,
                        at=data.frame(ldem=c(1,0), lano=c(0,0))))



m1.zero.dummies2$vcov <- v1.zero
cat("Table D4\n")
stargazer(m1.zero.dummies2, m1.dummies.pooled, m1.dummies.re, m1.dummies.cre, #m1.poisson,
          title = "Alternative model structure",
          label="tab:diff.models", #float.env = "sidewaystable",
          model.names = TRUE, model.numbers = TRUE,
          omit=c("factor*", "*.bar",
                 "lremitPC:lano", "lano", "^ldem$",
                 "llmil", "llpop", "lgdpgrowth", "llgdppc", "lpress", 
                 "lgini", "lmaxlowx", "lwexcluded", "lnum_ongoing_ucdp"),
          dep.var.labels = c("Domestic Attacks"),
          order=c("lremitPC","ldem", "lano"),
          se=list(NULL,
                  sqrt(diag(v1.dummies.pooled)),
                  NULL, sqrt(diag(v1.dummies.cre))),
          align=TRUE,
          no.space=TRUE,
          digits=2,
          omit.stat = c("AIC", "BIC"),
          star.cutoffs = c(0.1, .05, NA),
          covariate.labels = c("Remittances per capita",
                               "Remittances per capita $\\times$ dem."),
          notes=c("$p<[0.*], p<[0.**]$. Coefficients from negative binomial models. Clustered standard errors in parentheses (except for RE model)."),
          notes.append=FALSE
)
cat("Theta values for zero-inflated and RE negative binomials for Table D4\n")
print(round(m1.zero.dummies2$theta,2))
print(round(getME(m1.dummies.re, "glmer.nb.theta"),2))

cat("Combined coefs, se's, and p-vals for Table D4\n")
  est <- c(m1.zero.dummies2$coefficients$count['lremitPC']+m1.zero.dummies2$coefficients$count['lremitPC:ldem'],
         m1.dummies.pooled$coefficients['lremitPC']+m1.dummies.pooled$coefficients['lremitPC:ldem'],
         m1.dummies.re@beta[2]+m1.dummies.re@beta[14],
         m1.dummies.cre$coefficients['lremitPC']+m1.dummies.cre$coefficients['lremitPC:ldem'])
se <- sqrt(c(v1.zero['count_lremitPC','count_lremitPC']+
               v1.zero['count_lremitPC:ldem','count_lremitPC:ldem']+
               2*v1.zero['count_lremitPC','count_lremitPC:ldem'],
             v1.dummies.pooled['lremitPC', 'lremitPC'] +
               v1.dummies.pooled['lremitPC:ldem','lremitPC:ldem']+
               2* v1.dummies.pooled['lremitPC', 'lremitPC:ldem'] ,
             vcov(m1.dummies.re)[2,2] + vcov(m1.dummies.re)[14,14]
             +2*vcov(m1.dummies.re)[2,14],
             v1.dummies.cre['lremitPC', 'lremitPC'] +
               v1.dummies.cre['lremitPC:ldem','lremitPC:ldem']+
               2* v1.dummies.cre['lremitPC', 'lremitPC:ldem'] ))
p.vals <- pnorm(abs(est/se), lower=F)*2
print(round(rbind(est,se,p.vals), 2))
