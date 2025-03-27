# Baseline all intervention model code
# Timothy Kuiper timothykuiper@gmail.com
# March 2025
# For Science article: "Dehorning reduces rhino poaching" (Kuiper et al. 2025)

#Explanation

#Code used to run the various versions of the baseline all intervention model 
#(see main article Methods). Models the effect of dehorning as well all the other 
#anti-poaching interventions

# Data set used and sensitivity ---------------------------------------------

#This code uses the data set "glmm.data" which consists of the following columns
#(the unit of analysis is the reserve-quarter)

#reserve ID
#year
#quarter
#rhino population size
#res size in km2
#number of poacher incursions - n.inc
#number of rhinos poached - n.poached
#proportion of the rhino population dehorned - p.dehorn (see Methods)
#a column for each intervention index (see paper Methods)
#exposure index (how exposed reserve is - see Methods)

#As explained in the main article, this data are sensitive. Queries related to 
#data access can be directed to the Greater Kruger Environmental Protection Foundation 
#info@gkepf.org”

# load packages -----------------------------------------------------------

library(brms)
library(tidyverse)

# Baseline all intervention model -----------------------------------------

mod.nb.main<- brm(n.poached ~ p.dehorn+ #proportion dehorned in each reserve-quarter
                    
                    #Other intervention indices
                    rangers.100km2+cr.rrt.index+
                    sec.test.roll+fail.prop.roll3+
                    exposure+
                    det.zones.std+fence.ind.std+
                    track.dogs.100km2+det.dogs.gate.std+
                    access.index.std + cam.tech.ZAR.km2
                    rhino.mon.index + fa.index+
                    
                    #Offset for population size
                    offset(log(tot.pop))+
                    
                    #Random effects
                    (1|reserve)+(1|quarter),
                    
                  #Prior - horsheoe for regularisation
                  set_prior("horseshoe(1)"),
                  
                  #Model parameters
                  family = negbinomial(), 
                  data = glmm.data,
                  chains = 4, # nb of chains
                  iter = 2000, # nb of iterations, including burnin
                  warmup = 1000, # burnin
                  thin = 1,
                  save_pars = save_pars(all = TRUE),
                  control = list(adapt_delta = 0.99))
summary(mod.nb.main)

#Posterior checks
pp_check(mod.nb.main)+xlim(0,40) #distribution of simulated versus observed data
pp_check(mod.nb.main, type='error_scatter_avg')
pp_check(mod.nb.main, type='stat', stat='mean')

# Model with no regularisation - flat prior -------------------------------
#Default prior for brm used

mod.no.hs<- brm(n.poached ~ p.dehorn+ #proportion dehorned in each reserve-quarter
                    
                    #Other intervention indices
                    rangers.100km2+cr.rrt.index+
                    sec.test.roll+fail.prop.roll3+
                    exposure+
                    det.zones.std+fence.ind.std+
                    track.dogs.100km2+det.dogs.gate.std+
                    access.index.std + cam.tech.ZAR.km2
                  rhino.mon.index + fa.index+
                    
                    #Offset for population size
                    offset(log(tot.pop))+
                    
                    #Random effects
                    (1|reserve)+(1|quarter),
                  
                  #Model parameters
                  family = negbinomial(), 
                  data = glmm.data,
                  chains = 4, # nb of chains
                  iter = 2000, # nb of iterations, including burnin
                  warmup = 1000, # burnin
                  thin = 1,
                  save_pars = save_pars(all = TRUE),
                  control = list(adapt_delta = 0.99))

# Model with LASSO regularisation -------------------------------

mod.lasso<- brm(n.poached ~ p.dehorn+ #proportion dehorned in each reserve-quarter
                  
                  #Other intervention indices
                  rangers.100km2+cr.rrt.index+
                  sec.test.roll+fail.prop.roll3+
                  exposure+
                  det.zones.std+fence.ind.std+
                  track.dogs.100km2+det.dogs.gate.std+
                  access.index.std + cam.tech.ZAR.km2
                rhino.mon.index + fa.index+
                  
                  #Offset for population size
                  offset(log(tot.pop))+
                  
                  #Random effects
                  (1|reserve)+(1|quarter),
                
                #LASSO
                prior(lasso(), class = "b"),
                
                #Model parameters
                family = negbinomial(), 
                data = glmm.data,
                chains = 4, # nb of chains
                iter = 2000, # nb of iterations, including burnin
                warmup = 1000, # burnin
                thin = 1,
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99))

# Original Poisson model to compare fit -----------------------------------

mod.poisson<- brm(n.poached ~ p.dehorn+ #proportion dehorned in each reserve-quarter
               
               #Other intervention indices
               rangers.100km2+cr.rrt.index+
               sec.test.roll+fail.prop.roll3+
               exposure+
               det.zones.std+fence.ind.std+
               track.dogs.100km2+det.dogs.gate.std+
               access.index.std + cam.tech.ZAR.km2
             rhino.mon.index + fa.index+
               
               #Offset for population size
               offset(log(tot.pop))+
               
               #Random effects
               (1|reserve)+(1|quarter),
             
             #Prior - horsheoe for regularisation
             set_prior("horseshoe(1)"),
             
             #Model parameters
             family = poisson(), 
             data = glmm.data,
             chains = 4, # nb of chains
             iter = 2000, # nb of iterations, including burnin
             warmup = 1000, # burnin
             thin = 1,
             save_pars = save_pars(all = TRUE),
             control = list(adapt_delta = 0.99))

#Negative binomial fits better
#Based on leave one out cross validation

LOO(mod.poiss, mod.nb.main, moment_match = TRUE,reloo = TRUE)
(mw <- model_weights(mod.poiss, mod.nb.main))

mod.poiss<- add_criterion(mod.poiss, "waic")
mod.nb.main<- add_criterion(mod.nb.main, "waic")

loo_compare(mod.poiss, mod.nb.main, criterion = "waic")

##Calculate jkapp value of NB - overdispersion
# Extract posterior samples of phi (shape parameter)
post <- posterior_samples(mod.nb.main, pars = "shape")

# Compute kappa (1/phi)
post$kappa <- 1 / post$shape

# Summary statistics of kappa
summary(post$kappa)
summary(post$shape)

kappa_estimate <- posterior_summary(mod.nb.main, pars = "shape")
kappa_value <- 1 / kappa_estimate[, 1]  # Compute kappa from shape
kappa_value

# Model with poacher incursions as response -------------------------------

mod.inc<- brm(n.inc ~ p.dehorn+ #proportion dehorned in each reserve-quarter
                    
                    #Other intervention indices
                    rangers.100km2+cr.rrt.index+
                    sec.test.roll+fail.prop.roll3+
                    exposure+
                    det.zones.std+fence.ind.std+
                    track.dogs.100km2+det.dogs.gate.std+
                    access.index.std + cam.tech.ZAR.km2
                  rhino.mon.index + fa.index+
                    
                    #Offset for reserve size
                    offset(log(res.km2))+
                    
                    #Random effects
                    (1|reserve)+(1|quarter),
                  
                  #Prior - horsheoe for regularisation
                  set_prior("horseshoe(1)"),
                  
                  #Model parameters
                  family = negbinomial(), 
                  data = glmm.data,
                  chains = 4, # nb of chains
                  iter = 2000, # nb of iterations, including burnin
                  warmup = 1000, # burnin
                  thin = 1,
                  save_pars = save_pars(all = TRUE),
                  control = list(adapt_delta = 0.99))


# Model with lagged effects -----------------------------------------------

##Lag by two quarters
##specify columns to lag and not

sel<-c(7,9,10,16,18,20,22,26)
t<-1:ncol(glmm.data); 
cols.not<-t[-which(t %in% sel)]

#Create lagged dataset
r.lst<-glmm.data %>% group_split(reserve)
df<-r.lst[[1]];lag=2
modify.df<-function(df,lag){
  r<-nrow(df)
  d1<-df[-(1:lag),cols.not];d1
  d2<-df[-((r-lag+1):r),sel];d2
  d3<-cbind(d1,d2)
  return(d3)
}
modify.df(r.lst[[1]],2)
r.lst.lagged<-lapply(r.lst,modify.df,lag=1)

glmm.data.lagged<-do.call("rbind",r.lst.lagged)

#Run lagged model

mod.lagged<- brm(n.poached ~ p.dehorn+ #proportion dehorned in each reserve-quarter
                
                #Other intervention indices
                rangers.100km2+cr.rrt.index+
                sec.test.roll+fail.prop.roll3+
                exposure+
                det.zones.std+fence.ind.std+
                track.dogs.100km2+det.dogs.gate.std+
                access.index.std + cam.tech.ZAR.km2
              rhino.mon.index + fa.index+
                
                #Offset for population size
                offset(log(tot.pop))+
                
                #Random effects
                (1|reserve)+(1|quarter),
              
              #Prior - horsheoe for regularisation
              set_prior("horseshoe(1)"),
              
              #Model parameters
              family = negbinomial(), 
              data = glmm.data,
              chains = 4, # nb of chains
              iter = 2000, # nb of iterations, including burnin
              warmup = 1000, # burnin
              thin = 1,
              save_pars = save_pars(all = TRUE),
              control = list(adapt_delta = 0.99))



# Models with key interactions --------------------------------------------

#See Methods for explanation

mod.interactions<- brm(n.poached ~ p.dehorn+ #proportion dehorned in each reserve-quarter
                   
                   #Main interactions
                     
                     det.zones.std*track.dogs.100km2+
                     cam.tech.ZAR.km2*track.dogs.100km2+
                     det.zones.std*rangers.100km2+
                     cam.tech.ZAR.km2*rangers.100km2
                   
                   #Offset for population size
                   offset(log(tot.pop))+
                   
                   #Random effects
                   (1|reserve)+(1|quarter),
                 
                 #Prior - horsheoe for regularisation
                 set_prior("horseshoe(1)"),
                 
                 #Model parameters
                 family = negbinomial(), 
                 data = glmm.data,
                 chains = 4, # nb of chains
                 iter = 2000, # nb of iterations, including burnin
                 warmup = 1000, # burnin
                 thin = 1,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.99))

# Notes on supplementary models due to correlated predictors --------------------------------

#Two models were run to account for correlated predictors

#Note tracking tech was correlated 0.76 with rangers, 
#0.84 with cam tech, 0.74 with ranger equipment
#ranger equipment correlated with rangers at 0.87
#Exposure and rhino monitoring have 0.70 correlations

#Therefore ran two supplementary models
#(1) remove rangers and rhino monitoring, track tech
#(2) remove rte, rangers, cam tech, and exposure
