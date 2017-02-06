install.packages("OIsurv")
install.packages("survsim")
install.packages("broom")
library(OIsurv)
library(dplyr)
library(broom)
library(ggplot2)
library(survsim)
library(KMsurv)

data(tongue); head(tongue) # real dataset (from KMSurv package)
summary(tongue)
help(tongue)

# Analyzing just one type of tumor (Aneuploid Tumor)
tongue2 = tongue %>% 
  filter(type == 1) 

# Converting into a Surv object 
tongue2_surv <- Surv(tongue2$time, tongue2$delta)

# Getting KM estimator
tongue2_survfit = survfit(tongue2_surv ~ 1)
plot(tongue2_survfit)

glance(tongue2_survfit) # median shows 50% of subjects died within 93 weeks of diagnosis

# Graphically Comparing KM estimator for 2 tumors
tongue_survfit = survfit(Surv(time = time, event = delta) ~ type, data = tongue)
plot(tongue_survfit, lty = 2:3, xlab = "weeks", ylab = "Proporation Survival")
legend(100, .8, c("Aneuploid", "Diploid"), lty = 2:3) 

# ggplot2 version of the plot
tongue_tidy = tidy(tongue_survfit)
mx = max(tongue_tidy$n.censor)
ggplot(tongue_tidy, aes(time, estimate, fill = strata)) + 
  geom_line()+
  geom_point(aes(shape = as.factor(n.censor)), size = 3)+
  scale_shape_manual(values=c(NA, 1:mx))+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.25)+
  xlab("weeks")+
  ylab("Proportion Survival") # confidence intervals overlaps for both type of tumors

# Statistical tests for comparison, whether there is significant difference between types of tumors
survdiff(Surv(time = time, event = delta) ~ type, data = tongue, rho = 0) # log-rank, default

# positive rho implies higher weight to initial part
survdiff(Surv(time = time, event = delta) ~ type, data = tongue, rho = 1)

# Simulating dataset
set.seed(2365)
d1 = simple.surv.sim(n = 500, foltime = 1000, dist.ev = 'weibull', anc.ev = 1.5, beta0.ev = 2, # event dist - weibull with p = 1.5
                     dist.cens = 'weibull', anc.cens = 10, beta0.cens = 2.01, # censoring dist - weibull with p = 10
                     #z = list(c("weibull", 1)), #assuming independent observations i.e. no within subject correlation 
                     beta = list(c(.6)), x = list(c("bern", .7)) # 1 binary covariate with p=0.7 and beta=0.6
)

# View(d1)
summary(d1)

summary(as.data.frame(d1))

d1_survfit = survfit(Surv(time = start, time2 = stop, event = status) ~ x, data = d1)
d1_survfit = survfit(Surv(time = stop, event = status) ~ x, data = d1) # same as above

plot(d1_survfit, lty = 2:3, xlab = "time", ylab = "Proporation Survival")
legend(100, .8, c("x=0", "x=1"), lty = 2:3) 

# ggplot2 version of the plot
d1_tidy = tidy(d1_survfit)
mx = max(d1_tidy$n.censor)
ggplot(d1_tidy, aes(time, estimate, fill = strata)) + 
  geom_line() +
  geom_point(aes(shape = as.factor(n.censor)), size = 2) + 
  scale_shape_manual(values=c(NA, 1:mx))+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.25) + 
  xlab("time") + 
  ylab("Proportion Survival")

# Statistical tests for comparison
survdiff(Surv(time = stop, event = status) ~ x, data = d1, rho = 0)
survdiff(Surv(time = stop, event = status) ~ x, data = d1, rho = 1)
survdiff(Surv(time = stop, event = status) ~ x, data = d1, rho = -1)

# TBC: https://www.r-bloggers.com/survival-analysis-2/
