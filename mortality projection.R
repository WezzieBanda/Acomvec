library(rstan)
options(mc.cores=4)
rstan_options(auto_write = TRUE)
library(dplyr)
library(tidyr)
library(readxl)
library(adegenet)

#setwd("Wezzie-Tom/EHT DATA/Benin Trial")

df<- read_excel("Benin_NNP_hut_trial_raw data (2).xlsx", 
                sheet = "combined2")

df2 <- df %>%
  filter(tot_dead !="NA" & total!="0" & total!="NA" & Treatment!="UT" ) #& Hut!="Ifakara")

#data with age 0
#dfa <- df2 %>%
 # filter(age == 0, Treatment %in% c("UT", "IG1", "IG2", "P3", "RG"))

dfa <- df2 %>%
  filter(age == 0, Treatment %in% c("IG1"))


dfa20 <- df2 %>%
  filter(age == 3, Treatment %in% c( "IG1"))

dfa_match = merge(dfa,dfa20,by="Treatment")

dfa_match = dfa_match %>% drop_na(total.x)
dfa_match = dfa_match %>% drop_na(tot_dead.x)
dfa_match = dfa_match %>% drop_na(total.y)
dfa_match = dfa_match %>% drop_na(tot_dead.y)

dfa_match = subset(dfa_match,dfa_match$total.x != 0)


## Simple model first
setup_inputs = function(df2){
  
  # Prep data
  S <- as.numeric(length(df2$total.x))
  prop_dead <- c(df2$tot_dead.x)/df2$total.x
  X_halflife <- c(df2$total.y - df2$tot_dead.y)
  X_caught_halflife <- df2$total.y
  
  
  # Need to pass it the data:
  data_stan <- list(S=S, 
                    prop_dead=prop_dead,
                    X_halflife=X_halflife,
                    N_caught_halflife=X_caught_halflife)#,
  # nsite = nsite,
  # site = site)
  
  
  return(data_stan)
  
}







# Compile modelC:\Users\user\Documents\Wezzie-Tom\EHT DATA\Benin Trial

full_model<- stan_model("C:/Users/user/Documents/Wezzie-Tom/EHT DATA/Benin Trial/gamma.stan") # flex params


stan_base <- rstan::stan(file="C:/Users/user/Documents/Wezzie-Tom/EHT DATA/Benin Trial/gamma.stan", 
                         data=setup_inputs(dfa_match), 
                         warmup=1000,
                         control = list(adapt_delta = 0.8,
                                        max_treedepth = 20),
                         iter=2000, chains=1)
base <- rstan::extract(stan_base)
median(base$a)
median(base$b)

run_f = function(df2,dataset){
  data_stan = setup_inputs(df2)
  
  # Run mode
  fit_full <- sampling(full_model, data=data_stan, iter=2000, chains=1)
  # params_a = extract(fit_full)
  
  # save fit
  saveRDS(fit_full, paste0("C:/Users/user/Documents/Wezzie-Tom/EHT DATA/Benin Trial/stan model outputs/half_life",dataset,".rds"))
}
#C:\Users\user\Documents\Wezzie-Tom\EHT DATA\Benin Trial\stan model outputs

run_f(dfa_match,"df2")#all data

# load fit
fit_a <- readRDS("C:/Users/user/Documents/Wezzie-Tom/EHT DATA/Benin Trial/stan model outputs/half_lifedf2.rds")

fit_b <- readRDS("C:/Users/user/Documents/Wezzie-Tom/EHT DATA/Benin Trial/ento_fit1_half_life_pyrethroid_nets (4).rds")

############################
##############################




prop_killed_0_wash = setup_inputs(dfa_match)$prop_dead
prop_killed20_wash = setup_inputs(dfa_match)$X_halflife/setup_inputs(dfa_match)$N_caught_halflife

test = base
names(test)
sim_X = seq(0,1,0.01)
P_hl1 <- 1/(1 + exp(-(median(test$a) + median(test$b) * sim_X)))
P_hl1_low <- 1/(1 + exp(-(quantile(test$a,0.025) + quantile(test$b,0.025) * sim_X)))
P_hl1_upp <- 1/(1 + exp(-(quantile(test$a,0.975) + quantile(test$b,0.975) * sim_X)))


par(mfrow=c(1,1))

lp_tau = seq(0,1,0.01)-0.5
mort1 = seq(0,1,length=length(lp_tau))
mu1 = -2.36
rho1 = -3.05
original_eLife = 1/(1 + exp(-(-mu1 - rho1 * lp_tau)))
median(original_eLife)
hw=log(2)/0.9137258

plot(original_eLife ~ c(1-mort1),ylim=c(0,1),xlim=c(0,1),
     ylab = "EHT Mortality after 20 washes (%)",
     xlab = "EHT Mortality after 0 washes (%)",
     xaxt="n",yaxt="n")
axis(1,at=seq(0,1,0.2),labels=seq(0,100,20))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))

lines(P_hl1 ~ sim_X,col="red")
polygon(c(sim_X,rev(sim_X)),
        c(P_hl1_low,rev(P_hl1_upp)),border=NA,col=adegenet::transp("red",0.4))
points(prop_killed_0_wash, prop_killed20_wash,
       ylim=c(0,1),xlim=c(0,1),
       col=adegenet::transp("darkred",0.6),pch=19)

Hw = log(2)/(1/(1+exp(-(median(test$a) + median(test$b) * sim_X))))



