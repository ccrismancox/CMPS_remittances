library(data.table)
library(car)
library(margins)
library(MASS)
library(stargazer)
library(xtable)
library(ggplot2)
library(matrixStats)
rm(list=ls())
source('stargazerNoteCorrection.r')
load("../output/mainResultsRemittances.rdata")

Est0 <- m0.dummies$coef['lremitPC']+m0.dummies$coef['lremitPC:ldem']
SE0 <- sqrt(v0.dummies['lremitPC', 'lremitPC']+v0.dummies['lremitPC:ldem', 'lremitPC:ldem'] + 2*v0.dummies['lremitPC:ldem', 'lremitPC'])

Est2 <- m2.dummies$coef['lremitPC']+m2.dummies$coef['lremitPC:ldem']
SE2 <- sqrt(v2.dummies['lremitPC', 'lremitPC']+v2.dummies['lremitPC:ldem', 'lremitPC:ldem'] + 2*v2.dummies['lremitPC:ldem', 'lremitPC'])
round(c(Est2, SE2, pnorm(abs(Est2/SE2), lower=FALSE)*2),3)


Est3 <- m3.dummies$coef['lremitPC']+m3.dummies$coef['lremitPC:ldem']
SE3 <- sqrt(v3.dummies['lremitPC', 'lremitPC']+v3.dummies['lremitPC:ldem', 'lremitPC:ldem'] + 2*v3.dummies['lremitPC:ldem', 'lremitPC'])
round(c(Est3, SE3, pnorm(abs(Est3/SE3), lower=FALSE)*2),3)
m3.boot <- mvrnorm(5000, m3.dummies$coef, v3.dummies)

cat("Information in Tables 3 and B1\n")
stargazer( m0.dummies,
           m2.dummies, m3.dummies, m4.dummies,
           title = "Negative Binomial Results", label="tab:regTable",
           dep.var.labels = "Domestic terrorist attacks",
           no.space=TRUE, digits=2,
           omit.stat="AIC",
           order =c("lremitPC"),
           covariate.labels = c("Remittances",
                                "Remittances $\\times$ Democracy",
                                "Remittances $\\times$ GDP pc",
                                "Remittances $\\times$ Anocracy",
                                "Remittances $\\times$ Democracy $\\times$ GDP pc",
                                "Remittances $\\times$ Anocracy $\\times$ GPP pc",
                                "Democracy",
                                "Anocracy",
                                "Military personnel",
                                "Population",
                                "GDP growth",
                                "GDP per capita",
                                "Free Press",
                                "GINI",
                                "Horizontal inequality",
                                "Excluded population",
                                "\\# of ongoing civil conflicts",
                                "Lag attacks",
                                "Democracy $\\times$ GDP pc",
                                "Anocracy $\\times$ GPP pc"),
           omit=c("factor*"),
           se=list(sqrt(diag(v0.dummies)),
                   # sqrt(diag(v1.dummies)), 
                   sqrt(diag(v2.dummies)), 
                   sqrt(diag(v3.dummies)), sqrt(diag(v4.dummies))),
           align=TRUE,      
           star.cutoffs = c(0.1, .05, NA),
           notes=c("$p<[0.*]$", "$p<[0.**]$", "Coefficients from negative binomial models. Clustered stanard errors in parenthesis."),
           notes.append=FALSE
)
combined.betas <- 
  rbind(
    c(Est0, Est2, Est3),
    c(SE0,SE2, SE3),
    c(pnorm(abs(Est0/SE0), lower=FALSE)*2,
      pnorm(abs(Est2/SE2), lower=FALSE)*2,
      pnorm(abs(Est3/SE3), lower=FALSE)*2))
rownames(combined.betas) <- c("est", "se", "p value")
cat("Combined coefficients Table 3\n")
print(round(combined.betas,2))
print(round(cbind(coefs4[,1:2],
                  pnorm(abs(coefs4[,1]/coefs4[,2]), lower=F)*2)[c(1,3,2,4),],
                  2))


marg.out <- data.frame(est = mar2$AME[c(1,3)],
                       lo = mar2$lower[c(1,3)],
                       hi = mar2$upper[c(1,3)])


marg.tab <- rbind(paste("$", sprintf("%.2f", marg.out$est), "$"),
                  paste("$(", sprintf("%.2f", marg.out$lo),
                        ", ",sprintf("%.2f", marg.out$hi), ")$", sep=""))
colnames(marg.tab) <- c("Autocracy", "Democracy")
rownames(marg.tab) <- c("Change in attacks", "")
cat("Table 4\n")
print(xtable(marg.tab, align="rcc",
             caption="Average marginal effect  of remittances on domestic terrorist attacks (Model 2)", 
             label="tab.mar"), 
      caption.placement="top",
      booktabs=T,
      sanitize.text.function=function(x){x})





m2.boot <- mvrnorm(5000, mu=m2.dummies$coefficients, Sigma=v2.dummies)
exp.remit.auto1 <- exp.remit.dem1 <- exp.remit.ano1 <- exp.remit.diff <- matrix(0, nrow=50, ncol=3)
RSeq <- seq(0, 3, length=50)

for(i in 1:50){
  aX <- dX <-  anX <- m2.dummies$x
  aX[,"ldem"] <- 0
  aX[,"lano"] <- 0
  aX[,"lremitPC:ldem"] <- 0
  aX[,"lremitPC:lano"] <- 0
  aX[,"lremitPC"] <- RSeq[i]
  exp.remit.auto1[i,] <- quantile(rowMeans(t(exp(aX %*% t(m2.boot)))), probs=c(0.025,.5, .975))
  dX[,"ldem"] <- 1
  dX[,"lano"] <- 0
  dX[,"lremitPC:lano"] <- 0
  dX[,"lremitPC"] <- dX[,"lremitPC:ldem"] <-  RSeq[i]
  anX[,"ldem"] <- 0
  anX[,"lano"] <- 1
  anX[,"lremitPC:lano"] <-   anX[,"lremitPC"]  <- RSeq[i]
  anX[,"lremitPC:ldem"] <- 0
  exp.remit.dem1[i,] <- quantile(rowMeans(t(exp(dX %*% t(m2.boot)))), probs=c(0.025,.5, .975))
  exp.remit.ano1[i,] <- quantile(rowMeans(t(exp(anX %*% t(m2.boot)))), probs=c(0.025,.5, .975))
  exp.remit.diff[i,] <- quantile(rowMeans(t(exp(aX %*% t(m2.boot)))- t(exp(dX %*% t(m2.boot)))), probs=c(0.025,.5, .975))
  
}
plot.df <- rbind.data.frame(exp.remit.auto1, exp.remit.dem1)
colnames(plot.df) <- c("lo", "Est", "hi")
plot.df$Regime <- c(rep(c("Autocracy", "Democracy"), each=50))
plot.df$Remittances <- RSeq

ExpAttacks <- ggplot(plot.df)+
  geom_ribbon(aes(x=Remittances, ymin=lo, ymax=hi), alpha=.2)+
  geom_line(aes(x=Remittances, y=Est))+
  geom_rug(aes(lremitPC), data=remitNEW[lremitPC <= max(RSeq)& sam==1])+
  facet_grid(~Regime)+
  theme_bw(14)+
  xlab("Remittances (Hundreds of USD per person)")+
  ylab("Expected number of\ndomestic terrorist attacks")

print(ExpAttacks)

ggsave(ExpAttacks, file="../output/Figure1.pdf", width=8,height=4)


plot.df <- rbind.data.frame(exp.remit.diff)
colnames(plot.df) <- c("lo", "Est", "hi")
plot.df$Remittances <- RSeq

DiffAttacks <- ggplot(plot.df)+
  geom_ribbon(aes(x=Remittances, ymin=lo, ymax=hi), alpha=.2)+
  geom_line(aes(x=Remittances, y=Est))+
  geom_rug(aes(lremitPC), data=remitNEW[lremitPC <= max(RSeq) & sam==1])+
  theme_bw(14)+
  xlab("Remittances (Hundreds of USD per person)")+
  ylab("Expected difference in\ndomestic terrorist attacks")

print(DiffAttacks)

ggsave(DiffAttacks, file="../output/FigureC1.pdf", width=8,height=4)


