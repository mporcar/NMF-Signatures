library("data.table")
library("ggplot2")
setwd("~/Desktop/DistributionTime")
source("Munge/Transform.R")


ggplot(dataset, aes(delay_part, fill = interval, colour = interval)) +
  geom_density(alpha = 0.1) +
  xlim(0, 150)


ggplot(dataset[!is.na(p_edad)], aes(delay_part, fill = p_edad   , colour = p_edad   )) +
  geom_density(alpha = 0.1) +
  xlim(0, 150)


ggplot(dataset[!is.na(CampaignType)], aes(delay_part, fill = CampaignType    , colour = CampaignType           )) +
  geom_density(alpha = 0.1) +
  xlim(0, 150)

ggplot(dataset, aes(delay_part, fill = Nsurvey    , colour = Nsurvey           )) +
  geom_density(alpha = 0.1) +
  xlim(0, 150)

ggplot(dataset, aes(delay_part, fill = claseSocial    , colour = claseSocial           )) +
  geom_density(alpha = 0.1) +
  xlim(0, 150)

dataset[,CS_EDAD:=paste0(interval," // ",p_edad)]

ggplot(dataset, aes(delay_part, fill = CS_EDAD    , colour = CS_EDAD)) +
  geom_density(alpha = 0.1) +
  xlim(0, 150)

emp2535=ecdf(dataset[p_edad == "Entre 25 y 35"]$delay_part)
del2535<-as.numeric(dataset[p_edad == "Entre 25 y 35"]$delay_part)
del5665<-as.numeric(dataset[p_edad == "Entre 56 y 65"]$delay_part)
quantile(del2535,0.20)
quantile(del5665,0.20)
quantile(del2535,0.80)
quantile(del5665,0.80)
?quantile
emp=ecdf(delay3)
save("emp.f", file="empFunction.Rdata")
save("emp", file="empFunction.Rdata")
load("empFunction.Rdata")
emp
emp(1)
?ecdf
Sys.time(emp(130))
system.time(emp(130))
length(unique(delay3))
quantile(emp.f,0.6)
library(fitdistrplus)
descdist(b, discrete = FALSE)
dlnorm
library(logspline)
?
# fit.weibull = fitdist(dt.par$diffpar, "weibull")
fit.gamma <- fitdist(b, "gamma")
fit.dlnorm <- fitdist(as.numeric(dataset$delay_part), "lnorm")
fit.dlnorm
?fitdist
plot(density(rlnorm(40000,meanlog=2.262182,sdlog=2.079152)),xlim=c(-1,200))
qlnorm(0.85,meanlog=2.262182,sdlog=2.079152)
plot(fit.dlnorm,
     dat <- data.frame(cond = factor(rep(c("A"), each=40000)), 
                       rating = c(rlnorm(40000,meanlog=2.262182,sdlog=2.079152)))
     
     plnorm(96,meanlog=2.262182,sdlog=2.079152)
     ggplot(dat, aes(x=rating)) +
       geom_density(alpha=.2, fill="#FF6666")+  # Overlay with transparent density plot+
       xlim(-1,100)
     ?f
     plot(fit.lnorm)
     ?fitdist
     # fit.weibull <- fitdist(dt.par$diffpar, distr = "weibull", method = "mle", lower = c(0, 0))
     # fit.gamma <- fitdist(dt.par$diffpar, distr = "gamma", method = "mle", lower = c(0, 0), start = list(scale = 1, shape = 1))
     # ?fitdist
     sum(dt.par[diffpar>10,])