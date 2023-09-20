#######Libraries, Functions, and Data#####
library(xtable)
library(data.table)
rm(list=ls())

load("../output/mainResultsRemittances.rdata")
load("../data/mainRemittanceData.rdata")

num2str <- function(x){
  formatC(x, digits=2, format='f')
}

MODEL <- m2.dummies$model[,grep("factor",colnames(m2.dummies$model), invert = TRUE)]

results <- matrix("0", ncol=6, nrow=ncol(MODEL))



results[,2]<-num2str(apply(MODEL, 2, min, na.rm=TRUE))
results[,3]<-num2str(apply(MODEL, 2, mean, na.rm=TRUE))
results[,4]<-num2str(apply(MODEL, 2, sd, na.rm=TRUE))
results[,5]<-num2str(apply(MODEL, 2, max, na.rm=TRUE))
length(unique(remitNEW[sam==1]$ccode))
range(unique(remitNEW[sam==1]$year))

results[,1]<-c("Domestic attacks",
               "Remittances per capita",
               "Democracy",
               "Anocracy",
               "Mil. per. pc (logged)",
               "Population (logged)",
               "Economic growth",
               "GDP per capita (logged)",
               "Free press",
               "GINI",
               "NHI",
               "Political inequality",
               "\\# of ongoing civil conflicts" )

results[, 6] <- c("GTD",
                  "World Bank",
                  "Polity IV",
                  "Polity IV",
                  "COW-NMC",
                  "World Bank",
                  "World Bank",
                  "World Bank",
                  "Li (2005)/Freedom House",
                  "World Bank",
                  "\\citet{BCG2014}",
                  "EPR",
                  "UCDP")


colnames(results) <- paste("\\multicolumn{1}{c}",
                           c("{Variable}",
                             "{Min}",
                             "{Mean}",
                             "{St.~Dev.}",
                             "{Max}",
                             "{Source/Measurement}"),
                           sep="")

cat("Table 1\n")
print(xtable(results,
             caption="Summary statistics for main variables",
             label="tab:sumStat",
             align="rrddddl"),
      include.rownames=FALSE,
      sanitize.text.function=function(x){x},
      booktabs = TRUE,
      caption.placement="top",
      table.placement="ht!")



length(unique( remitNEW[sam==1]$ccode))
summary(unique( remitNEW[sam==1]$year))






attackCols <- c("assault",
                "assassination", 
                "bombing", 
                "infrastructure",
                "hijack",
                "hostage",
                "unarmed",
                "other",
                "nkill")


build.tests <- function(x){
  with(x,  c(estimate, statistic, p.value))
}
tab <- 
  do.call(rbind,
          lapply(
            list(t.test(I(assault/domAttacks)~demo, data=remitNEW[sam==1& anoc==0 & year != 1993]),
                 t.test(I(assassination/domAttacks)~demo, data=remitNEW[sam==1& anoc==0 & year != 1993]),
                 t.test(I(hijack/domAttacks)~demo, data=remitNEW[sam==1& anoc==0 & year != 1993]),
                 t.test(I(unarmed/domAttacks)~demo, data=remitNEW[sam==1& anoc==0 & year != 1993]),
                 t.test(I(bombing/domAttacks)~demo, data=remitNEW[sam==1& anoc==0 & year != 1993]),
                 t.test(I(infrastructure/domAttacks)~demo, data=remitNEW[sam==1& anoc==0 & year != 1993]),
                 t.test(I(hostage/domAttacks)~demo, data=remitNEW[sam==1& anoc==0 & year != 1993]),
                 t.test(I(other/domAttacks)~demo, data=remitNEW[sam==1& anoc==0 & year != 1993]),
                 t.test(domAttacks~demo, data=remitNEW[sam==1& anoc==0 & year != 1993]),
                 t.test(I(nkill/domAttacks)~demo, data=remitNEW[sam==1& anoc==0 & year != 1993])),
            build.tests)
  )
rownames(tab) <- c("Assault", 
                   "Assassination",
                   "Hijacking",
                   "Unarmed",
                   "Bombing",
                   "Infrastructure",
                   "Hostage",
                   "Other",
                   "Total attacks",
                   "Fatalities/attack")
colnames(tab) <- c("Autocracy mean", "Democracy mean", "$t$ stat.", "$p$-value")


cat("Table 2\n")

print(xtable(tab,
             caption = "Domestic terrorist attack descriptions by regime type",
             label = "tab:terror_diffs",
             align="lrrrr"),
      include.rownames=TRUE,
      sanitize.text.function=function(x){x},
      caption.placement="top",
      hline.after = c(-1,0,8,10),
      booktabs = TRUE,
      table.placement="ht!")


cat("Context for remittances \n")
#context
print(lm(remittance_received_2010usd~remitPC, data=remitNEW)$coef["remitPC"]/1000000)
print(lm(I(remittance_received_2010usd*weightedExclude)~remitPC, data=remitNEW)$coef["remitPC"]/1000000)




### Appendix A
countryList <-  remitNEW[sam==1][as.numeric(rownames(m4.dummies$model)),][,list(ccode, stateabb, country,year, domAttacks)]
countryList[, `:=`(nyears=length(year),
                   nAttacks=as.integer(sum(domAttacks))),
            by=ccode]
countryList[,year:=NULL]
countryList[,domAttacks:=NULL]
cat("Table A1\n")
print(xtable(unique(countryList)[,list(country, nyears, nAttacks)][order(country)]),
      include.rownames=FALSE, booktabs = T)
