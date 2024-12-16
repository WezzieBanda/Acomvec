library(lme4)
library(devtools)
library(MASS)
library(RColorBrewer)
library(cowplot)
library(ggplot2)
library(dplyr)
library(reshape2)
library(fitdistrplus)
library(lattice)
library(performance)
library(gtsummary)
library(readxl)
library(patchwork)
library(gridExtra)
library(writexl)


########;lets start with unwashed_funestus

Malawi_EHT_data <- read_excel("Malawi_EHT_data.xlsx")


##########lets model total bf#####################



total_mosquitoes_by_treatment <- Malawi_EHT_data %>%
  group_by(treatment, net_status, hut,species, code_captureur, compartment) %>%
  summarise(
    total_mosquitoes = sum(Total, na.rm = TRUE),
    total_unf_live = sum(unf_live, na.rm = TRUE),
    total_unf_dead = sum(unf_dead, na.rm = TRUE),
    total_bf_live=sum(bf_live,   na.rm=TRUE),
    total_bf_dead=sum(bf_dead,   na.rm=TRUE),
    total_gravid_live= sum(gravid_live, na.rm=TRUE),
    total_gravid_dead= sum(gravid_dead, na.rm=TRUE),
    total_dead= sum(tot_dead, na.rm=TRUE),
  )

total_mosquitoes_by_treatment$treatment <- as.factor(total_mosquitoes_by_treatment$treatment)

total_mosquitoes_by_treatment$treatment <- relevel(total_mosquitoes_by_treatment$treatment, ref = "UT")


#######now lets model the total nmber of mosquitoes caught 


fit1 <- glmer.nb(
  total_mosquitoes ~ treatment + net_status + species + species*net_status + 
    (1 | hut),
  data = total_mosquitoes_by_treatment
)
summary(fit1)

tbl_regression(fit1, exponentiate = TRUE)


#unwashed funestus
Total_mosqutioes <- c(
  exp(coef(summary(fit1))["(Intercept)", "Estimate"]) *108 ,
  exp(coef(summary(fit1))["(Intercept)", "Estimate"] + coef(summary(fit1))["treatmentIG1", "Estimate"]) *108 ,
  exp(coef(summary(fit1))["(Intercept)", "Estimate"] + coef(summary(fit1))["treatmentIG2", "Estimate"]) *108 ,
  exp(coef(summary(fit1))["(Intercept)", "Estimate"] + coef(summary(fit1))["treatmentOP", "Estimate"]) *108 ,
  exp(coef(summary(fit1))["(Intercept)", "Estimate"] + coef(summary(fit1))["treatmentP2", "Estimate"]) *108 ,
  exp(coef(summary(fit1))["(Intercept)", "Estimate"] + coef(summary(fit1))["treatmentRG", "Estimate"])*108 ,
#Washed funestus
  exp(coef(summary(fit1))["(Intercept)", "Estimate"] + coef(summary(fit1))["net_statusWashed", "Estimate"]) *108 ,
  exp(coef(summary(fit1))["(Intercept)", "Estimate"] + coef(summary(fit1))["treatmentIG1", "Estimate"]+coef(summary(fit1))["net_statusWashed", "Estimate"]) *108 ,
  exp(coef(summary(fit1))["(Intercept)", "Estimate"] + coef(summary(fit1))["treatmentIG2", "Estimate"]+coef(summary(fit1))["net_statusWashed", "Estimate"]) *108 ,
  exp(coef(summary(fit1))["(Intercept)", "Estimate"] + coef(summary(fit1))["treatmentOP", "Estimate"]+coef(summary(fit1))["net_statusWashed", "Estimate"]) *108 ,
  exp(coef(summary(fit1))["(Intercept)", "Estimate"] + coef(summary(fit1))["treatmentP2", "Estimate"]+coef(summary(fit1))["net_statusWashed", "Estimate"]) *108 ,
  exp(coef(summary(fit1))["(Intercept)", "Estimate"] + coef(summary(fit1))["treatmentRG", "Estimate"]+coef(summary(fit1))["net_statusWashed", "Estimate"]) *108 ,
#unwashed gambiae
  exp(coef(summary(fit1))["(Intercept)", "Estimate"] + coef(summary(fit1))["speciesan_gambiae_sl", "Estimate"])  *108 ,
  exp(coef(summary(fit1))["(Intercept)", "Estimate"] + coef(summary(fit1))["treatmentIG1", "Estimate"]+ coef(summary(fit1))["speciesan_gambiae_sl", "Estimate"]) *108 ,
  exp(coef(summary(fit1))["(Intercept)", "Estimate"] + coef(summary(fit1))["treatmentIG2", "Estimate"]+ coef(summary(fit1))["speciesan_gambiae_sl", "Estimate"]) *108 ,
  exp(coef(summary(fit1))["(Intercept)", "Estimate"] + coef(summary(fit1))["treatmentOP", "Estimate"] + coef(summary(fit1))["speciesan_gambiae_sl", "Estimate"]) *108 ,
  exp(coef(summary(fit1))["(Intercept)", "Estimate"] + coef(summary(fit1))["treatmentP2", "Estimate"] + coef(summary(fit1))["speciesan_gambiae_sl", "Estimate"]) *108 ,
  exp(coef(summary(fit1))["(Intercept)", "Estimate"] + coef(summary(fit1))["treatmentRG", "Estimate"] + coef(summary(fit1))["speciesan_gambiae_sl", "Estimate"]) *108 ,
#washed gambiae
  exp(coef(summary(fit1))["(Intercept)", "Estimate"] + coef(summary(fit1))["speciesan_gambiae_sl", "Estimate"] + coef(summary(fit1))["net_statusWashed", "Estimate"]+ coef(summary(fit1))["net_statusWashed:speciesan_gambiae_sl", "Estimate"]) *108,
  exp(coef(summary(fit1))["(Intercept)", "Estimate"] + coef(summary(fit1))["speciesan_gambiae_sl", "Estimate"] + coef(summary(fit1))["treatmentIG1", "Estimate"]     +coef(summary(fit1))["net_statusWashed", "Estimate"]+ coef(summary(fit1))["net_statusWashed:speciesan_gambiae_sl", "Estimate"]) *108 ,
  exp(coef(summary(fit1))["(Intercept)", "Estimate"] + coef(summary(fit1))["speciesan_gambiae_sl", "Estimate"] + coef(summary(fit1))["treatmentIG2", "Estimate"]     +coef(summary(fit1))["net_statusWashed", "Estimate"]+ coef(summary(fit1))["net_statusWashed:speciesan_gambiae_sl", "Estimate"]) *108,
  exp(coef(summary(fit1))["(Intercept)", "Estimate"] + coef(summary(fit1))["speciesan_gambiae_sl", "Estimate"] + coef(summary(fit1))["treatmentOP", "Estimate"]      +coef(summary(fit1))["net_statusWashed", "Estimate"]+ coef(summary(fit1))["net_statusWashed:speciesan_gambiae_sl", "Estimate"]) *108,
  exp(coef(summary(fit1))["(Intercept)", "Estimate"] + coef(summary(fit1))["speciesan_gambiae_sl", "Estimate"] + coef(summary(fit1))["treatmentP2", "Estimate"]      +coef(summary(fit1))["net_statusWashed", "Estimate"]+ coef(summary(fit1))["net_statusWashed:speciesan_gambiae_sl", "Estimate"]) *108,
  exp(coef(summary(fit1))["(Intercept)", "Estimate"] + coef(summary(fit1))["speciesan_gambiae_sl", "Estimate"] + coef(summary(fit1))["treatmentRG", "Estimate"]      +coef(summary(fit1))["net_statusWashed", "Estimate"]+ coef(summary(fit1))["net_statusWashed:speciesan_gambiae_sl", "Estimate"]) *108

)


Total_mosqutioes


# Create a data frame to save the results
results <- data.frame(
  Treatment = c("UTU", "IG1U", "IG2U", "OPU", "P2U", "RGU", "UTW", "IG1W", "IG2W", "OPW", "P2W", "RGW"),
  Total_mosqutioes= Total_mosqutioes
)


##########################now lets model total bf live mosquitoes using a logistic reg



fit2 <-
  glmer(
    cbind( total_bf_live, total_mosquitoes- total_bf_live ) ~
      treatment + net_status + species + species* net_status +(1|hut),
    family = binomial, data = total_mosquitoes_by_treatment) 
summary(fit2)


tbl_regression(fit2, exponentiate = TRUE)

InvLogit <- function(x) {
  exp(x) / (1 + exp(x))
}


Total_bflive <- c(
  InvLogit(coef(summary(fit2))["(Intercept)", "Estimate"]) ,
  InvLogit(coef(summary(fit2))["(Intercept)", "Estimate"] + coef(summary(fit2))["treatmentIG1", "Estimate"]) ,
  InvLogit(coef(summary(fit2))["(Intercept)", "Estimate"] + coef(summary(fit2))["treatmentIG2", "Estimate"])  ,
  InvLogit(coef(summary(fit2))["(Intercept)", "Estimate"] + coef(summary(fit2))["treatmentOP", "Estimate"]) ,
  InvLogit(coef(summary(fit2))["(Intercept)", "Estimate"] + coef(summary(fit2))["treatmentP2", "Estimate"])  ,
  InvLogit(coef(summary(fit2))["(Intercept)", "Estimate"] + coef(summary(fit2))["treatmentRG", "Estimate"]),
  #Washed funestus
  InvLogit(coef(summary(fit2))["(Intercept)", "Estimate"] + coef(summary(fit2))["net_statusWashed", "Estimate"]),
  InvLogit(coef(summary(fit2))["(Intercept)", "Estimate"] + coef(summary(fit2))["treatmentIG1", "Estimate"]+coef(summary(fit2))["net_statusWashed", "Estimate"]) ,
  InvLogit(coef(summary(fit2))["(Intercept)", "Estimate"] + coef(summary(fit2))["treatmentIG2", "Estimate"]+coef(summary(fit2))["net_statusWashed", "Estimate"]) ,
  InvLogit(coef(summary(fit2))["(Intercept)", "Estimate"] + coef(summary(fit2))["treatmentOP", "Estimate"]+coef(summary(fit2))["net_statusWashed", "Estimate"]),
  InvLogit(coef(summary(fit2))["(Intercept)", "Estimate"] + coef(summary(fit2))["treatmentP2", "Estimate"]+coef(summary(fit2))["net_statusWashed", "Estimate"]) ,
  InvLogit(coef(summary(fit2))["(Intercept)", "Estimate"] + coef(summary(fit2))["treatmentRG", "Estimate"]+coef(summary(fit2))["net_statusWashed", "Estimate"]) ,
  #unwashed gambiae
  InvLogit(coef(summary(fit2))["(Intercept)", "Estimate"] + coef(summary(fit2))["speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit2))["(Intercept)", "Estimate"] + coef(summary(fit2))["treatmentIG1", "Estimate"]+ coef(summary(fit2))["speciesan_gambiae_sl", "Estimate"]) ,
  InvLogit(coef(summary(fit2))["(Intercept)", "Estimate"] + coef(summary(fit2))["treatmentIG2", "Estimate"]+ coef(summary(fit2))["speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit2))["(Intercept)", "Estimate"] + coef(summary(fit2))["treatmentOP", "Estimate"] + coef(summary(fit2))["speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit2))["(Intercept)", "Estimate"] + coef(summary(fit2))["treatmentP2", "Estimate"] + coef(summary(fit2))["speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit2))["(Intercept)", "Estimate"] + coef(summary(fit2))["treatmentRG", "Estimate"] + coef(summary(fit2))["speciesan_gambiae_sl", "Estimate"]) ,
  #washed gambiae
  InvLogit(coef(summary(fit2))["(Intercept)", "Estimate"] + coef(summary(fit2))["speciesan_gambiae_sl", "Estimate"] + coef(summary(fit2))["net_statusWashed", "Estimate"]+ coef(summary(fit2))["net_statusWashed:speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit2))["(Intercept)", "Estimate"] + coef(summary(fit2))["speciesan_gambiae_sl", "Estimate"] + coef(summary(fit2))["treatmentIG1", "Estimate"]     +coef(summary(fit2))["net_statusWashed", "Estimate"]+ coef(summary(fit2))["net_statusWashed:speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit2))["(Intercept)", "Estimate"] + coef(summary(fit2))["speciesan_gambiae_sl", "Estimate"] + coef(summary(fit2))["treatmentIG2", "Estimate"]     +coef(summary(fit2))["net_statusWashed", "Estimate"]+ coef(summary(fit2))["net_statusWashed:speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit2))["(Intercept)", "Estimate"] + coef(summary(fit2))["speciesan_gambiae_sl", "Estimate"] + coef(summary(fit2))["treatmentOP", "Estimate"]      +coef(summary(fit2))["net_statusWashed", "Estimate"]+ coef(summary(fit2))["net_statusWashed:speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit2))["(Intercept)", "Estimate"] + coef(summary(fit2))["speciesan_gambiae_sl", "Estimate"] + coef(summary(fit2))["treatmentP2", "Estimate"]      +coef(summary(fit2))["net_statusWashed", "Estimate"]+ coef(summary(fit2))["net_statusWashed:speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit2))["(Intercept)", "Estimate"] + coef(summary(fit2))["speciesan_gambiae_sl", "Estimate"] + coef(summary(fit2))["treatmentRG", "Estimate"]      +coef(summary(fit2))["net_statusWashed", "Estimate"]+ coef(summary(fit2))["net_statusWashed:speciesan_gambiae_sl", "Estimate"])
  
)

Total_bflive

# Create a data frame to save the results
results2 <- data.frame(
  Treatment = c("UTU", "IG1U", "IG2U", "OPU", "P2U", "RGU", "UTW", "IG1W", "IG2W", "OPW", "P2W", "RGW"),
  Total_bflive= Total_bflive
)
results2




########################################################now lets do bffed dead############################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################




fit3<-
  glmer(
    cbind(total_unf_live, total_mosquitoes- total_unf_live ) ~
      treatment + net_status + species + species* net_status +(1|hut),
    family = binomial, data = total_mosquitoes_by_treatment) 
summary(fit3)


tbl_regression(fit3, exponentiate = TRUE)

InvLogit <- function(x) {
  exp(x) / (1 + exp(x))
}



Total_unfedlive <- c(
  InvLogit(coef(summary(fit3))["(Intercept)", "Estimate"]) ,
  InvLogit(coef(summary(fit3))["(Intercept)", "Estimate"] + coef(summary(fit3))["treatmentIG1", "Estimate"]) ,
  InvLogit(coef(summary(fit3))["(Intercept)", "Estimate"] + coef(summary(fit3))["treatmentIG2", "Estimate"])  ,
  InvLogit(coef(summary(fit3))["(Intercept)", "Estimate"] + coef(summary(fit3))["treatmentOP", "Estimate"]) ,
  InvLogit(coef(summary(fit3))["(Intercept)", "Estimate"] + coef(summary(fit3))["treatmentP2", "Estimate"])  ,
  InvLogit(coef(summary(fit3))["(Intercept)", "Estimate"] + coef(summary(fit3))["treatmentRG", "Estimate"]),
  #Washed funestus
  InvLogit(coef(summary(fit3))["(Intercept)", "Estimate"] + coef(summary(fit3))["net_statusWashed", "Estimate"]),
  InvLogit(coef(summary(fit3))["(Intercept)", "Estimate"] + coef(summary(fit3))["treatmentIG1", "Estimate"]+coef(summary(fit3))["net_statusWashed", "Estimate"]) ,
  InvLogit(coef(summary(fit3))["(Intercept)", "Estimate"] + coef(summary(fit3))["treatmentIG2", "Estimate"]+coef(summary(fit3))["net_statusWashed", "Estimate"]) ,
  InvLogit(coef(summary(fit3))["(Intercept)", "Estimate"] + coef(summary(fit3))["treatmentOP", "Estimate"]+coef(summary(fit3))["net_statusWashed", "Estimate"]),
  InvLogit(coef(summary(fit3))["(Intercept)", "Estimate"] + coef(summary(fit3))["treatmentP2", "Estimate"]+coef(summary(fit3))["net_statusWashed", "Estimate"]) ,
  InvLogit(coef(summary(fit3))["(Intercept)", "Estimate"] + coef(summary(fit3))["treatmentRG", "Estimate"]+coef(summary(fit3))["net_statusWashed", "Estimate"]) ,
  #unwashed gambiae
  InvLogit(coef(summary(fit3))["(Intercept)", "Estimate"] + coef(summary(fit3))["speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit3))["(Intercept)", "Estimate"] + coef(summary(fit3))["treatmentIG1", "Estimate"]+ coef(summary(fit3))["speciesan_gambiae_sl", "Estimate"]) ,
  InvLogit(coef(summary(fit3))["(Intercept)", "Estimate"] + coef(summary(fit3))["treatmentIG2", "Estimate"]+ coef(summary(fit3))["speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit3))["(Intercept)", "Estimate"] + coef(summary(fit3))["treatmentOP", "Estimate"] + coef(summary(fit3))["speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit3))["(Intercept)", "Estimate"] + coef(summary(fit3))["treatmentP2", "Estimate"] + coef(summary(fit3))["speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit3))["(Intercept)", "Estimate"] + coef(summary(fit3))["treatmentRG", "Estimate"] + coef(summary(fit3))["speciesan_gambiae_sl", "Estimate"]) ,
  #washed gambiae
  InvLogit(coef(summary(fit3))["(Intercept)", "Estimate"] + coef(summary(fit3))["speciesan_gambiae_sl", "Estimate"] + coef(summary(fit3))["net_statusWashed", "Estimate"]+ coef(summary(fit3))["net_statusWashed:speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit3))["(Intercept)", "Estimate"] + coef(summary(fit3))["speciesan_gambiae_sl", "Estimate"] + coef(summary(fit3))["treatmentIG1", "Estimate"]     +coef(summary(fit3))["net_statusWashed", "Estimate"]+ coef(summary(fit3))["net_statusWashed:speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit3))["(Intercept)", "Estimate"] + coef(summary(fit3))["speciesan_gambiae_sl", "Estimate"] + coef(summary(fit3))["treatmentIG2", "Estimate"]     +coef(summary(fit3))["net_statusWashed", "Estimate"]+ coef(summary(fit3))["net_statusWashed:speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit3))["(Intercept)", "Estimate"] + coef(summary(fit3))["speciesan_gambiae_sl", "Estimate"] + coef(summary(fit3))["treatmentOP", "Estimate"]      +coef(summary(fit3))["net_statusWashed", "Estimate"]+ coef(summary(fit3))["net_statusWashed:speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit3))["(Intercept)", "Estimate"] + coef(summary(fit3))["speciesan_gambiae_sl", "Estimate"] + coef(summary(fit3))["treatmentP2", "Estimate"]      +coef(summary(fit3))["net_statusWashed", "Estimate"]+ coef(summary(fit3))["net_statusWashed:speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit3))["(Intercept)", "Estimate"] + coef(summary(fit3))["speciesan_gambiae_sl", "Estimate"] + coef(summary(fit3))["treatmentRG", "Estimate"]      +coef(summary(fit3))["net_statusWashed", "Estimate"]+ coef(summary(fit3))["net_statusWashed:speciesan_gambiae_sl", "Estimate"])
  
)
Total_unfedlive

# Create a data frame to save the results
results3 <- data.frame(
  Treatment = c("UTU", "IG1U", "IG2U", "OPU", "P2U", "RGU", "UTW", "IG1W", "IG2W", "OPW", "P2W", "RGW"),
  Total_unflive= Total_unfedlive
)
results3




#######################3lastly let do mortality ###################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################





fit4<-
  glmer(
    cbind(total_dead, total_mosquitoes- total_dead ) ~
      treatment + net_status + species + species* net_status +(1|hut),
    family = binomial, data = total_mosquitoes_by_treatment) 
summary(fit4)


tbl_regression(fit4, exponentiate = TRUE)

InvLogit <- function(x) {
  exp(x) / (1 + exp(x))
}



Total_dead <- c(
  InvLogit(coef(summary(fit4))["(Intercept)", "Estimate"]) ,
  InvLogit(coef(summary(fit4))["(Intercept)", "Estimate"] + coef(summary(fit4))["treatmentIG1", "Estimate"]) ,
  InvLogit(coef(summary(fit4))["(Intercept)", "Estimate"] + coef(summary(fit4))["treatmentIG2", "Estimate"])  ,
  InvLogit(coef(summary(fit4))["(Intercept)", "Estimate"] + coef(summary(fit4))["treatmentOP", "Estimate"]) ,
  InvLogit(coef(summary(fit4))["(Intercept)", "Estimate"] + coef(summary(fit4))["treatmentP2", "Estimate"])  ,
  InvLogit(coef(summary(fit4))["(Intercept)", "Estimate"] + coef(summary(fit4))["treatmentRG", "Estimate"]),
  #Washed funestus
  InvLogit(coef(summary(fit4))["(Intercept)", "Estimate"] + coef(summary(fit4))["net_statusWashed", "Estimate"]),
  InvLogit(coef(summary(fit4))["(Intercept)", "Estimate"] + coef(summary(fit4))["treatmentIG1", "Estimate"]+coef(summary(fit4))["net_statusWashed", "Estimate"]) ,
  InvLogit(coef(summary(fit4))["(Intercept)", "Estimate"] + coef(summary(fit4))["treatmentIG2", "Estimate"]+coef(summary(fit4))["net_statusWashed", "Estimate"]) ,
  InvLogit(coef(summary(fit4))["(Intercept)", "Estimate"] + coef(summary(fit4))["treatmentOP", "Estimate"]+coef(summary(fit4))["net_statusWashed", "Estimate"]),
  InvLogit(coef(summary(fit4))["(Intercept)", "Estimate"] + coef(summary(fit4))["treatmentP2", "Estimate"]+coef(summary(fit4))["net_statusWashed", "Estimate"]) ,
  InvLogit(coef(summary(fit4))["(Intercept)", "Estimate"] + coef(summary(fit4))["treatmentRG", "Estimate"]+coef(summary(fit4))["net_statusWashed", "Estimate"]) ,
  #unwashed gambiae
  InvLogit(coef(summary(fit4))["(Intercept)", "Estimate"] + coef(summary(fit4))["speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit4))["(Intercept)", "Estimate"] + coef(summary(fit4))["treatmentIG1", "Estimate"]+ coef(summary(fit4))["speciesan_gambiae_sl", "Estimate"]) ,
  InvLogit(coef(summary(fit4))["(Intercept)", "Estimate"] + coef(summary(fit4))["treatmentIG2", "Estimate"]+ coef(summary(fit4))["speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit4))["(Intercept)", "Estimate"] + coef(summary(fit4))["treatmentOP", "Estimate"] + coef(summary(fit4))["speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit4))["(Intercept)", "Estimate"] + coef(summary(fit4))["treatmentP2", "Estimate"] + coef(summary(fit4))["speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit4))["(Intercept)", "Estimate"] + coef(summary(fit4))["treatmentRG", "Estimate"] + coef(summary(fit4))["speciesan_gambiae_sl", "Estimate"]) ,
  #washed gambiae
  InvLogit(coef(summary(fit4))["(Intercept)", "Estimate"] + coef(summary(fit4))["speciesan_gambiae_sl", "Estimate"] + coef(summary(fit4))["net_statusWashed", "Estimate"]+ coef(summary(fit4))["net_statusWashed:speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit4))["(Intercept)", "Estimate"] + coef(summary(fit4))["speciesan_gambiae_sl", "Estimate"] + coef(summary(fit4))["treatmentIG1", "Estimate"]     +coef(summary(fit4))["net_statusWashed", "Estimate"]+ coef(summary(fit4))["net_statusWashed:speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit4))["(Intercept)", "Estimate"] + coef(summary(fit4))["speciesan_gambiae_sl", "Estimate"] + coef(summary(fit4))["treatmentIG2", "Estimate"]     +coef(summary(fit4))["net_statusWashed", "Estimate"]+ coef(summary(fit4))["net_statusWashed:speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit4))["(Intercept)", "Estimate"] + coef(summary(fit4))["speciesan_gambiae_sl", "Estimate"] + coef(summary(fit4))["treatmentOP", "Estimate"]      +coef(summary(fit4))["net_statusWashed", "Estimate"]+ coef(summary(fit4))["net_statusWashed:speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit4))["(Intercept)", "Estimate"] + coef(summary(fit4))["speciesan_gambiae_sl", "Estimate"] + coef(summary(fit4))["treatmentP2", "Estimate"]      +coef(summary(fit4))["net_statusWashed", "Estimate"]+ coef(summary(fit4))["net_statusWashed:speciesan_gambiae_sl", "Estimate"]),
  InvLogit(coef(summary(fit4))["(Intercept)", "Estimate"] + coef(summary(fit4))["speciesan_gambiae_sl", "Estimate"] + coef(summary(fit4))["treatmentRG", "Estimate"]      +coef(summary(fit4))["net_statusWashed", "Estimate"]+ coef(summary(fit4))["net_statusWashed:speciesan_gambiae_sl", "Estimate"])
  
)
Total_dead

# Create a data frame to save the results
results4 <- data.frame(
  Treatment = c("UTU", "IG1U", "IG2U", "OPU", "P2U", "RGU", "UTW", "IG1W", "IG2W", "OPW", "P2W", "RGW"),
  Total_dead = Total_dead
)
results4



######################333combine all the 4 results

# Combine the data frames into a list
results_list <- list(
  "Results1" = results,
  "Results2" = results2,
  "Results3" = results3,
  "Results4" = results4
)

# Write to an Excel file with each data frame in a separate sheet
write_xlsx(results_list, "Combined_Results.xlsx")
