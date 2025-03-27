# Regression discontinuity model code
# Timothy Kuiper timothykuiper@gmail.com
# March 2025
# For Science article: "Dehorning reduces rhino poaching" (Kuiper et al. 2025)

#Explanation

#Code used to run the Regression Discontinuity in Time model 
#(see main article Methods).

# Data set used and sensitivity ---------------------------------------------

#This code uses the data set "reg.disc.data" which consists of the following columns
#(the unit of analysis is the reserve-quarter)

#reserve ID
#year
#year quarter 
#rhino population size
#number of rhinos poached - n.poached
#number of quarters relative to dehorning event (negative or positive)
#the above is the "running variable" for regression discontinuity
#before or after dehorning dummy variable

#As explained in the main article, this data are sensitive. Queries related to 
#data access can be directed to the Greater Kruger Environmental Protection Foundation 
#info@gkepf.org”


# load packages -----------------------------------------------------------

library(brms)
library(tidyverse)

# Bayesian reg discontinuity with kink ------------------------------------

# Tests for change in slope AND treatment effect

bayes.rdd2<- brm(poached ~ quarters.rel.deh + bef_aft+
                   quarters.rel.deh*bef_aft+
                   (1|reserve)+(1|quarter)+
                   offset(log(pop)), 
                 family = negbinomial(), 
                 data = reg.disc.data,
                 chains = 4, 
                 iter = 3000, 
                 warmup = 1500, 
                 thin = 1,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.99),
                 cores = parallel::detectCores())
summary(bayes.rdd2)

# Plot the effect ---------------------------------------------------------

#Figure 4 in the main article  

ce <- conditional_effects(bayes.rdd2, conditions = data.frame(pop = 100),
                          probs = c(0.05, 0.95))

dp<-tibble(ce[["quarters.since.deh:bef_aft"]])
names(dp)[c(10,12,13)]<-c("est","low","high")

##curve before
dp.b<-dp |> filter(bef_aft=="Before\ndehorning",
                   quarters.since.deh<0)
ggplot(dp.b, aes(x=quarters.since.deh,y=est))+
  geom_line() +
  geom_ribbon(aes(ymin = low, ymax = high,alpha=0.8) )

##curve after
dp.a<-dp |> filter(bef_aft=="After\ndehorning",
                   quarters.since.deh>0)
ggplot(dp.a, aes(x=quarters.since.deh,y=est))+
  geom_line() +
  geom_ribbon(aes(ymin = low, ymax = high,alpha=0.8) )

#combine
dfull<-rbind(dp.b,dp.a)

ggplot(dfull, aes(x=quarters.since.deh,y=est,
                  col=bef_aft,fill=bef_aft))+
  geom_line() +
  geom_ribbon(aes(ymin = low, ymax = high),alpha=0.4)+
  scale_colour_manual(values=c("#E69F00","#56B4E9"),
                      name="")+
  scale_fill_manual(values=c("#E69F00","#56B4E9"),
                    name="")+
  geom_vline(xintercept = 0, linetype="dotted",
             linewidth=0.8)+
  geom_point(data=allpd.q.deh,aes(x=quarters.since.deh,y=p.rate*100,
                                  col=bef_aft,
                                  shape=at.least.1.horned),
             position = position_jitter(height = 0.15,width=0.15),
             alpha=0.4,size=2)+
  scale_shape_manual(values = c(16, 8))+
  xlab("Quarter relative to dehorning")+
  ylab("Mean quarterly poaching rate (%)")+
  theme_classic()+
  theme(legend.position="top")

