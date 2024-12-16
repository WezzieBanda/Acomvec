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
library(nnet)
library(VGAM)

########;lets start with unwashed_funestus

Malawi_EHT_data <- read_excel("Malawi_EHT_data.xlsx")
Malawi_EHT_data <- read_excel("Wezzie-Tom/EHT DATA/Malawi_EHT_data.xlsx")


##########lets model total bf#####################



total_mosquitoes_by_treatment <- Malawi_EHT_data %>%
  group_by(treatment, net_status, hut,species, code_captureur, compartment) %>%
  summarise(
    total_mosquitoes = sum(Total, na.rm = TRUE),
    total_unf_live = sum(UNF_live, na.rm = TRUE),
    total_unf_dead = sum(unf_dead, na.rm = TRUE),
    total_bf_live=sum(Bflive,   na.rm=TRUE),
    total_bf_dead=sum(bf_dead,   na.rm=TRUE),
    total_gravid_live= sum(gravid_live, na.rm=TRUE),
    total_gravid_dead= sum(gravid_dead, na.rm=TRUE),
    total_dead= sum(Total_dead, na.rm=TRUE),
  )

write_xlsx(total_mosquitoes_by_treatment, "total_mosquitoes_by_treatment.xlsx")

total_mosquitoes_by_treatment$total_bf_live[total_mosquitoes_by_treatment$total_bf_live == 2] <- 1
total_mosquitoes_by_treatment$total_unf_live[total_mosquitoes_by_treatment$total_unf_live == 2] <- 1

total_mosquitoes_by_treatment$treatment <- as.factor(total_mosquitoes_by_treatment$treatment)
total_mosquitoes_by_treatment$treatment <- relevel(total_mosquitoes_by_treatment$treatment, ref = "UT")


unique(total_mosquitoes_by_treatment$total_unf_live)
unique(total_mosquitoes_by_treatment$total_bf_live)

#total_mosquitoes_by_treatment <- total_mosquitoes_by_treatment[total_mosquitoes_by_treatment$total_bf_live != 2, ]
#total_mosquitoes_by_treatment <- total_mosquitoes_by_treatment[total_mosquitoes_by_treatment$total_unf_live != 2, ]


#######now lets model the total nmber of mosquitoes caught 


fit1 <- glmer.nb(
  total_mosquitoes ~ treatment + net_status + species + species*net_status + 
    (1 | hut),
  data = total_mosquitoes_by_treatment
)
summary(fit1)

tbl<- tbl_regression(fit1, exponentiate = TRUE)
tbl_df <- as_tibble(tbl)
write_xlsx(tbl_df, "regression_table.xlsx")
View(tbl_df) 
getwd()


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

results



##########################now lets model total bf live mosquitoes using a logistic reg



fit2 <-
  glmer( total_bf_live ~
      treatment + net_status + species + species* net_status +(1|hut),
    family = binomial, data = total_mosquitoes_by_treatment) 
summary(fit2)


tbl_regression(fit2, exponentiate = TRUE)

tbl2<- tbl_regression(fit2, exponentiate = TRUE)
tbl2_df <- as_tibble(tbl2)
write_xlsx(tbl2_df, "regression_table2.xlsx")
View(tbl2_df) 
getwd()

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




########################################################now lets do uffed live############################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################




fit3<-
  glmer(total_unf_live ~
      treatment + net_status + species + species* net_status +(1|hut),
    family = binomial, data = total_mosquitoes_by_treatment) 
summary(fit3)


tbl_regression(fit3, exponentiate = TRUE)


tbl3<- tbl_regression(fit3, exponentiate = TRUE)
tbl3_df <- as_tibble(tbl3)
write_xlsx(tbl3_df, "regression_table3.xlsx")
View(tbl3_df) 
getwd()

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
  glmer(total_dead ~
      treatment + net_status + species + species* net_status +(1|hut),
    family = binomial, data = total_mosquitoes_by_treatment) 
summary(fit4)


tbl_regression(fit4, exponentiate = TRUE)

tbl4<- tbl_regression(fit4, exponentiate = TRUE)
tbl4_df <- as_tibble(tbl4)
write_xlsx(tbl4_df, "regression_table4.xlsx")
View(tbl4_df) 
getwd()

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













##############lets compute confidence intervals for all the fits

# Load necessary library for exporting to Excel
library(writexl)

# Function to compute confidence intervals directly
compute_ci <- function(terms, multiplier = 108) {
  coef_estimates <- coef(summary(fit1))[, "Estimate"]
  vcov_matrix <- vcov(fit1)
  
  # Get indices for the terms
  coef_indices <- match(terms, names(coef_estimates))
  
  # If any terms are missing, return NULL
  if (any(is.na(coef_indices))) {
    return(NULL)
  }
  
  coefs <- coef_estimates[coef_indices]
  vcov_subset <- vcov_matrix[coef_indices, coef_indices]
  
  # Point estimate
  estimate <- exp(sum(coefs)) * multiplier
  
  # Standard error for combined coefficients
  se <- sqrt(sum(vcov_subset))
  
  # Confidence intervals
  lower_ci <- exp(sum(coefs) - 1.96 * se) * multiplier
  upper_ci <- exp(sum(coefs) + 1.96 * se) * multiplier
  
  return(c(Estimate = estimate, CI_Lower = lower_ci, CI_Upper = upper_ci))
}

# Define scenarios for all treatments and conditions
scenarios <- list(
  c("(Intercept)"),
  c("(Intercept)", "treatmentIG1"),
  c("(Intercept)", "treatmentIG2"),
  c("(Intercept)", "treatmentOP"),
  c("(Intercept)", "treatmentP2"),
  c("(Intercept)", "treatmentRG"),
  c("(Intercept)", "net_statusWashed"),
  c("(Intercept)", "treatmentIG1", "net_statusWashed"),
  c("(Intercept)", "treatmentIG2", "net_statusWashed"),
  c("(Intercept)", "treatmentOP", "net_statusWashed"),
  c("(Intercept)", "treatmentP2", "net_statusWashed"),
  c("(Intercept)", "treatmentRG", "net_statusWashed"),
  c("(Intercept)", "speciesan_gambiae_sl"),
  c("(Intercept)", "treatmentIG1", "speciesan_gambiae_sl"),
  c("(Intercept)", "treatmentIG2", "speciesan_gambiae_sl"),
  c("(Intercept)", "treatmentOP", "speciesan_gambiae_sl"),
  c("(Intercept)", "treatmentP2", "speciesan_gambiae_sl"),
  c("(Intercept)", "treatmentRG", "speciesan_gambiae_sl"),
  c("(Intercept)", "speciesan_gambiae_sl", "net_statusWashed", "net_statusWashed:speciesan_gambiae_sl"),
  c("(Intercept)", "speciesan_gambiae_sl", "treatmentIG1", "net_statusWashed", "net_statusWashed:speciesan_gambiae_sl"),
  c("(Intercept)", "speciesan_gambiae_sl", "treatmentIG2", "net_statusWashed", "net_statusWashed:speciesan_gambiae_sl"),
  c("(Intercept)", "speciesan_gambiae_sl", "treatmentOP", "net_statusWashed", "net_statusWashed:speciesan_gambiae_sl"),
  c("(Intercept)", "speciesan_gambiae_sl", "treatmentP2", "net_statusWashed", "net_statusWashed:speciesan_gambiae_sl"),
  c("(Intercept)", "speciesan_gambiae_sl", "treatmentRG", "net_statusWashed", "net_statusWashed:speciesan_gambiae_sl")
)

# Compute confidence intervals for all scenarios
results21 <- lapply(scenarios, compute_ci)

# Remove NULL values from the results21 list (in case some scenarios didn't return valid results)
results21 <- results21[!sapply(results21, is.null)]

# Combine into a data frame
results_df21 <- do.call(rbind, results21)

# Ensure the results are in data frame format
results_df21 <- as.data.frame(results_df21)

# Set column names for the results data frame
colnames(results_df21) <- c("Estimate", "CI_Lower", "CI_Upper")

# Add Scenario names to the data frame
results_df21$Scenario <- names(results21)

# Print results to the console
print(results_df21)

# Save the results to an Excel file
write_xlsx(results_df21, "confidence_intervals_results.xlsx")

# Confirm that the file is saved
cat("Results have been saved to 'confidence_intervals_results.xlsx'.\n")




##########################################good now computing intervals for fit2 bf_live######################################################



# Function to compute confidence intervals for each coefficient
compute_ci_glmer <- function(model, terms) {
  coef_estimates <- coef(summary(model))[, "Estimate"]
  vcov_matrix <- vcov(model)
  
  # Get indices for the terms
  coef_indices <- match(terms, names(coef_estimates))
  
  # If any terms are missing, return NULL
  if (any(is.na(coef_indices))) {
    return(NULL)
  }
  
  coefs <- coef_estimates[coef_indices]
  vcov_subset <- vcov_matrix[coef_indices, coef_indices]
  
  # Point estimate
  estimate <- InvLogit(sum(coefs))
  
  # Standard error for combined coefficients
  se <- sqrt(sum(vcov_subset))
  
  # Confidence intervals
  lower_ci <- InvLogit(sum(coefs) - 1.96 * se)
  upper_ci <- InvLogit(sum(coefs) + 1.96 * se)
  
  return(c(Estimate = estimate, CI_Lower = lower_ci, CI_Upper = upper_ci))
}

# Define scenarios for all treatments and conditions
scenarios <- list(
  c("(Intercept)"),
  c("(Intercept)", "treatmentIG1"),
  c("(Intercept)", "treatmentIG2"),
  c("(Intercept)", "treatmentOP"),
  c("(Intercept)", "treatmentP2"),
  c("(Intercept)", "treatmentRG"),
  c("(Intercept)", "net_statusWashed"),
  c("(Intercept)", "treatmentIG1", "net_statusWashed"),
  c("(Intercept)", "treatmentIG2", "net_statusWashed"),
  c("(Intercept)", "treatmentOP", "net_statusWashed"),
  c("(Intercept)", "treatmentP2", "net_statusWashed"),
  c("(Intercept)", "treatmentRG", "net_statusWashed"),
  c("(Intercept)", "speciesan_gambiae_sl"),
  c("(Intercept)", "treatmentIG1", "speciesan_gambiae_sl"),
  c("(Intercept)", "treatmentIG2", "speciesan_gambiae_sl"),
  c("(Intercept)", "treatmentOP", "speciesan_gambiae_sl"),
  c("(Intercept)", "treatmentP2", "speciesan_gambiae_sl"),
  c("(Intercept)", "treatmentRG", "speciesan_gambiae_sl"),
  c("(Intercept)", "speciesan_gambiae_sl", "net_statusWashed", "net_statusWashed:speciesan_gambiae_sl"),
  c("(Intercept)", "speciesan_gambiae_sl", "treatmentIG1", "net_statusWashed", "net_statusWashed:speciesan_gambiae_sl"),
  c("(Intercept)", "speciesan_gambiae_sl", "treatmentIG2", "net_statusWashed", "net_statusWashed:speciesan_gambiae_sl"),
  c("(Intercept)", "speciesan_gambiae_sl", "treatmentOP", "net_statusWashed", "net_statusWashed:speciesan_gambiae_sl"),
  c("(Intercept)", "speciesan_gambiae_sl", "treatmentP2", "net_statusWashed", "net_statusWashed:speciesan_gambiae_sl"),
  c("(Intercept)", "speciesan_gambiae_sl", "treatmentRG", "net_statusWashed", "net_statusWashed:speciesan_gambiae_sl")
)

# Compute confidence intervals for all scenarios
results_glmer <- lapply(scenarios, compute_ci_glmer, model = fit2)

# Remove NULL values from the results_glmer list (in case some scenarios didn't return valid results)
results_glmer <- results_glmer[!sapply(results_glmer, is.null)]

# Combine into a data frame
results_df_glmer <- do.call(rbind, results_glmer)

# Ensure the results are in data frame format
results_df_glmer <- as.data.frame(results_df_glmer)

# Set column names for the results data frame
colnames(results_df_glmer) <- c("Estimate", "CI_Lower", "CI_Upper")

# Add Scenario names to the data frame
results_df_glmer$Scenario <- names(results_glmer)

# Print results to the console
print(results_df_glmer)

# Save the results to an Excel file
write_xlsx(results_df_glmer, "confidence_intervals_results_k.xlsx")

# Confirm that the file is saved
cat("Results have been saved to 'glmer_confidence_intervals_results.xlsx'.\n")


##############################################################################################################################################################
#######################################good now lets compute CI for unfedlive#############################################################




# Function to compute confidence intervals for each coefficient in fit3 model
compute_ci_fit3 <- function(model, terms) {
  coef_estimates <- coef(summary(model))[, "Estimate"]
  vcov_matrix <- vcov(model)
  
  # Get indices for the terms
  coef_indices <- match(terms, names(coef_estimates))
  
  # If any terms are missing, return NULL
  if (any(is.na(coef_indices))) {
    return(NULL)
  }
  
  coefs <- coef_estimates[coef_indices]
  vcov_subset <- vcov_matrix[coef_indices, coef_indices]
  
  # Point estimate using InvLogit
  estimate <- InvLogit(sum(coefs))
  
  # Standard error for combined coefficients
  se <- sqrt(sum(vcov_subset))
  
  # Confidence intervals
  lower_ci <- InvLogit(sum(coefs) - 1.96 * se)
  upper_ci <- InvLogit(sum(coefs) + 1.96 * se)
  
  return(c(Estimate = estimate, CI_Lower = lower_ci, CI_Upper = upper_ci))
}

# Define scenarios for all treatments and conditions in fit3 model
scenarios <- list(
  c("(Intercept)"),
  c("(Intercept)", "treatmentIG1"),
  c("(Intercept)", "treatmentIG2"),
  c("(Intercept)", "treatmentOP"),
  c("(Intercept)", "treatmentP2"),
  c("(Intercept)", "treatmentRG"),
  c("(Intercept)", "net_statusWashed"),
  c("(Intercept)", "treatmentIG1", "net_statusWashed"),
  c("(Intercept)", "treatmentIG2", "net_statusWashed"),
  c("(Intercept)", "treatmentOP", "net_statusWashed"),
  c("(Intercept)", "treatmentP2", "net_statusWashed"),
  c("(Intercept)", "treatmentRG", "net_statusWashed"),
  c("(Intercept)", "speciesan_gambiae_sl"),
  c("(Intercept)", "treatmentIG1", "speciesan_gambiae_sl"),
  c("(Intercept)", "treatmentIG2", "speciesan_gambiae_sl"),
  c("(Intercept)", "treatmentOP", "speciesan_gambiae_sl"),
  c("(Intercept)", "treatmentP2", "speciesan_gambiae_sl"),
  c("(Intercept)", "treatmentRG", "speciesan_gambiae_sl"),
  c("(Intercept)", "speciesan_gambiae_sl", "net_statusWashed", "net_statusWashed:speciesan_gambiae_sl"),
  c("(Intercept)", "speciesan_gambiae_sl", "treatmentIG1", "net_statusWashed", "net_statusWashed:speciesan_gambiae_sl"),
  c("(Intercept)", "speciesan_gambiae_sl", "treatmentIG2", "net_statusWashed", "net_statusWashed:speciesan_gambiae_sl"),
  c("(Intercept)", "speciesan_gambiae_sl", "treatmentOP", "net_statusWashed", "net_statusWashed:speciesan_gambiae_sl"),
  c("(Intercept)", "speciesan_gambiae_sl", "treatmentP2", "net_statusWashed", "net_statusWashed:speciesan_gambiae_sl"),
  c("(Intercept)", "speciesan_gambiae_sl", "treatmentRG", "net_statusWashed", "net_statusWashed:speciesan_gambiae_sl")
)

# Compute confidence intervals for all scenarios using the fit3 model
results_fit3 <- lapply(scenarios, compute_ci_fit3, model = fit3)

# Remove NULL values from the results_fit3 list (in case some scenarios didn't return valid results)
results_fit3 <- results_fit3[!sapply(results_fit3, is.null)]

# Combine into a data frame
results_df_fit3 <- do.call(rbind, results_fit3)

# Ensure the results are in data frame format
results_df_fit3 <- as.data.frame(results_df_fit3)

# Set column names for the results data frame
colnames(results_df_fit3) <- c("Estimate", "CI_Lower", "CI_Upper")

# Add Scenario names to the data frame
results_df_fit3$Scenario <- names(results_fit3)

# Print results to the console
print(results_df_fit3)

# Save the results to an Excel file
write_xlsx(results_df_fit3, "confidence_intervals_results_fit3.xlsx")

# Confirm that the file is saved
cat("Results have been saved to 'confidence_intervals_results_j.xlsx'.\n")


########################################################################################################################################################
#########################################intervals for tot dead ##################################################################################################

# Function to compute confidence intervals for each coefficient in fit4 model
compute_ci_fit4 <- function(model, terms) {
  coef_estimates <- coef(summary(model))[, "Estimate"]
  vcov_matrix <- vcov(model)
  
  # Get indices for the terms
  coef_indices <- match(terms, names(coef_estimates))
  
  # If any terms are missing, return NULL
  if (any(is.na(coef_indices))) {
    return(NULL)
  }
  
  coefs <- coef_estimates[coef_indices]
  vcov_subset <- vcov_matrix[coef_indices, coef_indices]
  
  # Point estimate using InvLogit
  estimate <- InvLogit(sum(coefs))
  
  # Standard error for combined coefficients
  se <- sqrt(sum(vcov_subset))
  
  # Confidence intervals
  lower_ci <- InvLogit(sum(coefs) - 1.96 * se)
  upper_ci <- InvLogit(sum(coefs) + 1.96 * se)
  
  return(c(Estimate = estimate, CI_Lower = lower_ci, CI_Upper = upper_ci))
}

# Define scenarios for all treatments and conditions in fit4 model
scenarios <- list(
  c("(Intercept)"),
  c("(Intercept)", "treatmentIG1"),
  c("(Intercept)", "treatmentIG2"),
  c("(Intercept)", "treatmentOP"),
  c("(Intercept)", "treatmentP2"),
  c("(Intercept)", "treatmentRG"),
  c("(Intercept)", "net_statusWashed"),
  c("(Intercept)", "treatmentIG1", "net_statusWashed"),
  c("(Intercept)", "treatmentIG2", "net_statusWashed"),
  c("(Intercept)", "treatmentOP", "net_statusWashed"),
  c("(Intercept)", "treatmentP2", "net_statusWashed"),
  c("(Intercept)", "treatmentRG", "net_statusWashed"),
  c("(Intercept)", "speciesan_gambiae_sl"),
  c("(Intercept)", "treatmentIG1", "speciesan_gambiae_sl"),
  c("(Intercept)", "treatmentIG2", "speciesan_gambiae_sl"),
  c("(Intercept)", "treatmentOP", "speciesan_gambiae_sl"),
  c("(Intercept)", "treatmentP2", "speciesan_gambiae_sl"),
  c("(Intercept)", "treatmentRG", "speciesan_gambiae_sl"),
  c("(Intercept)", "speciesan_gambiae_sl", "net_statusWashed", "net_statusWashed:speciesan_gambiae_sl"),
  c("(Intercept)", "speciesan_gambiae_sl", "treatmentIG1", "net_statusWashed", "net_statusWashed:speciesan_gambiae_sl"),
  c("(Intercept)", "speciesan_gambiae_sl", "treatmentIG2", "net_statusWashed", "net_statusWashed:speciesan_gambiae_sl"),
  c("(Intercept)", "speciesan_gambiae_sl", "treatmentOP", "net_statusWashed", "net_statusWashed:speciesan_gambiae_sl"),
  c("(Intercept)", "speciesan_gambiae_sl", "treatmentP2", "net_statusWashed", "net_statusWashed:speciesan_gambiae_sl"),
  c("(Intercept)", "speciesan_gambiae_sl", "treatmentRG", "net_statusWashed", "net_statusWashed:speciesan_gambiae_sl")
)

# Compute confidence intervals for all scenarios using the fit4 model
results_fit4 <- lapply(scenarios, compute_ci_fit4, model = fit4)

# Remove NULL values from the results_fit4 list (in case some scenarios didn't return valid results)
results_fit4 <- results_fit4[!sapply(results_fit4, is.null)]

# Combine into a data frame
results_df_fit4 <- do.call(rbind, results_fit4)

# Ensure the results are in data frame format
results_df_fit4 <- as.data.frame(results_df_fit4)

# Set column names for the results data frame
colnames(results_df_fit4) <- c("Estimate", "CI_Lower", "CI_Upper")

# Add Scenario names to the data frame
results_df_fit4$Scenario <- names(results_fit4)

# Print results to the console (check if any estimates seem off)
print(results_df_fit4)

# Save the results to an Excel file (adjust path if necessary)
write_xlsx(results_df_fit4, "C:/Users/user/Documents/Wezzie-Tom/EHT DATA/confidence_intervals_results_l.xlsx")

# Confirm that the file is saved
cat("Results have been saved to 'confidence_intervals_results_l.xlsx'.\n")













######################################################now lets plot the graphs for these estimates including the confidence intervals #######################################################################################
######################################################now lets plot the graphs for these estimates including the confidence intervals #######################################################################################
######################################################now lets plot the graphs for these estimates including the confidence intervals #######################################################################################
######################################################now lets plot the graphs for these estimates including the confidence intervals #######################################################################################

library(dplyr)
library(tidyr)


modelling_estimates <- read_excel("modelling estimates.xlsx")

# Read the data (assuming you have it loaded already)
data_funestus <- read_excel("modelling estimates.xlsx", sheet = "funestus")



# Define the specific order for treatments
treatment_order <- c("UTU", "UTW", "IG1U", "IG1W", "IG2U", "IG2W", "OPU", "OPW", "P2U", "P2W", "RGU", "RGW")

# Convert the Treatment column to a factor with the specified order
data_funestus$Treatment <- factor(data_funestus$Treatment, levels = treatment_order)



# Convert the Treatment column to a factor with the specified order
data_funestus$Treatment <- factor(data_funestus$Treatment, levels = treatment_order)

# Manually prepare data for plotting
plot_data <- data.frame(
  Treatment = rep(data_funestus$Treatment, 3),  # Replicate treatments for each category
  Category = rep(c("bflive", "unflive", "dead"), each = nrow(data_funestus)),  # Define category labels
  Estimate = c(data_funestus$Total_bflive, data_funestus$Total_unflive, data_funestus$Total_dead),  # Combine estimates
  CI_Lower = c(data_funestus$CI_Lower_bflive, data_funestus$CI_Lower_unflive, data_funestus$CI_Lower_totdead),  # Combine lower CI
  CI_Upper = c(data_funestus$CI_Upper_bflive, data_funestus$CI_Upper_unflive, data_funestus$CI_Upper_totdead)  # Combine upper CI
)





ggplot(plot_data, aes(x = Treatment, y = Estimate, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +  # Bar plot with dodged bars for each category
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, position = position_dodge(0.6)) +  # Add error bars for CI
  scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits from 0 to 1
  labs(title = "Funestus Mosquito Estimates with Confidence Intervals",
       y = "Estimate", x = "Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels





plot1<-ggplot(plot_data, aes(x = Treatment, y = Estimate, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +  # Bar plot with dodged bars for each category
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, position = position_dodge(0.6)) +  # Add error bars for CI
  scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits from 0 to 1
  labs(title = "Funestus",
       y = "Probability Estimate", x = "Treatment") +
  theme_minimal() +  # Use minimal theme for a clean background
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.line = element_line(color = "black", linewidth = 1),  # Add border line around the plot
    axis.ticks = element_line(color = "black", size = 1),  # Add ticks on axis
    panel.border = element_blank(),  # Remove the default border
    panel.background = element_rect(fill = "white", color = "black", size = 1),  # Background with border (whole panel)
    panel.grid = element_blank(),  # Remove the gridlines to make the borders stand out
    legend.background = element_rect(fill = "white", color = "black"),  # Box around the legend
    legend.box.background = element_rect(fill = "lightgray", color = "black"),  # Background of the legend box
    legend.title = element_text(size = 12, face = "bold"),  # Bold legend title
    legend.text = element_text(size = 10),  # Size of the legend text
    #legend.position = c(0.9, 0.8)  # Position the legend inside the plot (right-top corner)
    legend.position = "none"  # Remove the legend
  )



##########################now lets do for gambiae 








library(dplyr)
library(tidyr)


modelling_estimates <- read_excel("modelling estimates.xlsx")

# Read the data (assuming you have it loaded already)
data_funestus <- read_excel("modelling estimates.xlsx", sheet = "gambiae")



# Define the specific order for treatments
treatment_order <- c("UTU", "UTW", "IG1U", "IG1W", "IG2U", "IG2W", "OPU", "OPW", "P2U", "P2W", "RGU", "RGW")

# Convert the Treatment column to a factor with the specified order
data_funestus$Treatment <- factor(data_funestus$Treatment, levels = treatment_order)



# Convert the Treatment column to a factor with the specified order
data_funestus$Treatment <- factor(data_funestus$Treatment, levels = treatment_order)

# Manually prepare data for plotting
plot_data <- data.frame(
  Treatment = rep(data_funestus$Treatment, 3),  # Replicate treatments for each category
  Category = rep(c("bflive", "unflive", "dead"), each = nrow(data_funestus)),  # Define category labels
  Estimate = c(data_funestus$Total_bflive, data_funestus$Total_unflive, data_funestus$Total_dead),  # Combine estimates
  CI_Lower = c(data_funestus$CI_Lower_bflive, data_funestus$CI_Lower_unflive, data_funestus$CI_Lower_totdead),  # Combine lower CI
  CI_Upper = c(data_funestus$CI_Upper_bflive, data_funestus$CI_Upper_unflive, data_funestus$CI_Upper_totdead)  # Combine upper CI
)





ggplot(plot_data, aes(x = Treatment, y = Estimate, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +  # Bar plot with dodged bars for each category
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, position = position_dodge(0.6)) +  # Add error bars for CI
  scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits from 0 to 1
  labs(title = "Funestus Mosquito Estimates with Confidence Intervals",
       y = "Estimate", x = "Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels





plot2<- ggplot(plot_data, aes(x = Treatment, y = Estimate, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +  # Bar plot with dodged bars for each category
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, position = position_dodge(0.6)) +  # Add error bars for CI
  scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits from 0 to 1
  labs(title = "gambiae",
       y = "Probability Estimate", x = "Treatment") +
  theme_minimal() +  # Use minimal theme for a clean background
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.line = element_line(color = "black", linewidth = 1),  # Add border line around the plot
    axis.ticks = element_line(color = "black", size = 1),  # Add ticks on axis
    panel.border = element_blank(),  # Remove the default border
    panel.background = element_rect(fill = "white", color = "black", size = 1),  # Background with border (whole panel)
    panel.grid = element_blank(),  # Remove the gridlines to make the borders stand out
    legend.background = element_rect(fill = "white", color = "black"),  # Box around the legend
    legend.box.background = element_rect(fill = "lightgray", color = "black"),  # Background of the legend box
    legend.title = element_text(size = 10), #, face = "bold"),  # Bold legend title
    legend.text = element_text(size = 10),  # Size of the legend text
    legend.position = c(0.9, 0.8)  # Position the legend inside the plot (right-top corner)
  )





plot2 <- ggplot(plot_data, aes(x = Treatment, y = Estimate, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +  # Bar plot with dodged bars for each category
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, position = position_dodge(0.6)) +  # Add error bars for CI
  scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits from 0 to 1
  labs(title = "gambiae",
       y = "Probability Estimate", x = "Treatment") +
  theme_minimal() +  # Use minimal theme for a clean background
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.line = element_line(color = "black", linewidth = 1),  # Add border line around the plot
    axis.ticks = element_line(color = "black", size = 1),  # Add ticks on axis
    panel.border = element_blank(),  # Remove the default border
    panel.background = element_rect(fill = "white", color = "black", size = 1),  # Background with border (whole panel)
    panel.grid = element_blank(),  # Remove the gridlines to make the borders stand out
    legend.background = element_rect(fill = "white", color = "black"),  # Box around the legend
    legend.box.background = element_rect(fill = "lightgray", color = "black"),  # Background of the legend box
    legend.title = element_text(size = 12, face = "bold"),  # Bold legend title
    legend.text = element_text(size = 10),  # Size of the legend text
    legend.position = "bottom"  # Move the legend to the bottom
  )






plot2 <- ggplot(plot_data, aes(x = Treatment, y = Estimate, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +  # Bar plot with dodged bars for each category
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, position = position_dodge(0.6)) +  # Add error bars for CI
  scale_y_continuous(limits = c(0, 1)) +  # Set y-axis limits from 0 to 1
  labs(title = "gambiae",
       x = "Treatment") +  # Remove the y-axis label here
  theme_minimal() +  # Use minimal theme for a clean background
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.line = element_line(color = "black", linewidth = 1),  # Add border line around the plot
    axis.ticks = element_line(color = "black", size = 1),  # Add ticks on axis
    panel.border = element_blank(),  # Remove the default border
    panel.background = element_rect(fill = "white", color = "black", size = 1),  # Background with border (whole panel)
    panel.grid = element_blank(),  # Remove the gridlines to make the borders stand out
    legend.background = element_rect(fill = "white", color = "black"),  # Box around the legend
    legend.box.background = element_rect(fill = "lightgray", color = "black"),  # Background of the legend box
    legend.title = element_text(size = 12, face = "bold"),  # Bold legend title
    legend.text = element_text(size = 10),  # Size of the legend text
    legend.position = "bottom",  # Move the legend to the bottom
    axis.title.y = element_blank()  # Remove the y-axis label entirely
  )



# Combine the plots side by side
plot1 + plot2



###########################now we can also do the same for total mosqutioes caught

# Read the data from the Excel file (make sure you have the correct path)
modelling_estimates <- read_excel("modelling estimates.xlsx")

# Read data for Funestus and Gambiae
data_funestus <- read_excel("modelling estimates.xlsx", sheet = "funestus")
data_gambiae <- read_excel("modelling estimates.xlsx", sheet = "gambiae")





# Define the specific order for treatments
treatment_order <- c("UTU", "UTW", "IG1U", "IG1W", "IG2U", "IG2W", "OPU", "OPW", "P2U", "P2W", "RGU", "RGW")

# Convert the Treatment column to a factor with the specified order
data_funestus$Treatment <- factor(data_funestus$Treatment, levels = treatment_order)
data_gambiae$Treatment <- factor(data_gambiae$Treatment, levels = treatment_order)


# Manually prepare data for Total_mosquitoes plotting
plot_data_funestus <- data.frame(
  Treatment = data_funestus$Treatment,  # Use the treatments as they are
  Estimate = data_funestus$Total_mosquitoes,  # Use the Total_mosquitoes estimate
  CI_Lower = data_funestus$CI_Lower_total, # Use the lower confidence interval for Total_mosquitoes
  CI_Upper = data_funestus$CI_Upper_total  # Use the upper confidence interval for Total_mosquitoes
)

plot_data_gambiae <- data.frame(
  Treatment = data_gambiae$Treatment,  # Use the treatments as they are
  Estimate = data_gambiae$Total_mosquitoes,  # Use the Total_mosquitoes estimate
  CI_Lower = data_gambiae$CI_Lower_total, # Use the lower confidence interval for Total_mosquitoes
  CI_Upper = data_gambiae$CI_Upper_total  # Use the upper confidence interval for Total_mosquitoes
)


print(head(plot_data_funestus))  # Check the data for Funestus plot

# Plot for Funestus (with borders)
plot_funestus <- ggplot(plot_data_funestus, aes(x = Treatment, y = Estimate)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6, fill = "gray") +  # Bar plot for Total_mosquitoes
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, position = position_dodge(0.6)) +  # Add error bars for CI
  scale_y_continuous(limits = c(0,400)) +  # Set y-axis limits from 0 to 1
  labs(title = "Funestus",
       y = "Total Mosquito Counts", x = "Treatment") +
  theme_minimal() +  # Use minimal theme for a clean background
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.line = element_line(color = "black", linewidth = 1),  # Add border line around the plot
    axis.ticks = element_line(color = "black", size = 1),  # Add ticks on axis
    #panel.border = element_rect(color = "black", size = 1),  # Explicit border for the panel
    panel.background = element_rect(fill = "white", color = "black", size = 1),  # Set background color with border
    panel.grid = element_blank(),  # Remove gridlines for cleaner borders
    legend.position = "none"  # Remove the legend
  )
plot_funestus



# Plot for Gambiae (with borders and no y-axis label)
plot_gambiae <- ggplot(plot_data_gambiae, aes(x = Treatment, y = Estimate)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6, fill = "gray") +  # Bar plot for Total_mosquitoes
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, position = position_dodge(0.6)) +  # Add error bars for CI
  scale_y_continuous(limits = c(0, 400)) +  # Set y-axis limits from 0 to 400
  labs(title = "Gambiae", x = "Treatment") +  # No y-axis label specified here
  theme_minimal() +  # Use minimal theme for a clean background
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.line = element_line(color = "black", linewidth = 1),  # Add border line around the plot
    axis.ticks = element_line(color = "black", size = 1),  # Add ticks on axis
    panel.background = element_rect(fill = "white", color = "black", size = 1),  # Set background color with border
    panel.grid = element_blank(),  # Remove gridlines for cleaner borders
    legend.position = "none",  # Remove the legend
    axis.title.y = element_blank()  # Remove the y-axis label completely
  )




plot_gambiae

plot_funestus+plot_gambiae








# Aggregate 'unf_live' (unfed live mosquitoes) for each 'treatment' and 'net_status' group
funestus_summary_unfed <- Malawi_EHT_data_funestus%>%
  filter(treatment == "P2") %>%  # Filter to include only 'IG1' treatment
  group_by(treatment, net_status) %>%
  summarise(Unfed_Live_sum = sum(Total, na.rm = TRUE))  # Sum 'unf_live' per treatment and net status

# View the summary
print(funestus_summary_unfed)

Malawi_EHT_data_funestus <- read_excel("Wezzie-Tom/EHT DATA/Malawi_EHT_data_gambiae.xlsx")



total_mosquitoes_by_treatment <- Malawi_EHT_data











###################################3now lets model the compartment##############################

library(nnet)
library(VGAM)



Malawi_EHT_data$compartment <- as.factor(Malawi_EHT_data$compartment)
Malawi_EHT_data$Total <- as.numeric(Malawi_EHT_data$Total)

str(Malawi_EHT_data)
Malawi_EHT_data$compartment <- factor(Malawi_EHT_data$compartment)




fit.ms <- vglm(compartment ~ Total  , multinomial(refLevel = 1),
               data = Malawi_EHT_data)



summary(fit.ms)



# Get the fitted values from the model
fitted_values <- fitted(fit.ms)

# Order the data based on the 'Total' variable
ordered_data <- order(Malawi_EHT_data$Total)

# Plot the fitted probabilities
matplot(Malawi_EHT_data$Total[ordered_data], fitted_values[ordered_data, ], 
        type = "l", las = 1, lwd = 2, ylim = c(0, 1),
        ylab = "Fitted probabilities", 
        xlab = "Total",
        col = c("blue", "red", "green"),
        main = "Fitted Probabilities by Compartment")  # Adding a title)  # Change colors as needed


legend("bottomright", 
       col = c("blue", "red", "green"),  # Match colors from the plot
       lty = 1:3,                       # Line types (one for each compartment)
       legend = c("Net", "Room", "Veranda"),  # Labels for the compartments
       lwd = 2,
       xpd = TRUE,                      # Allow legend to be drawn outside the plot area
       inset = c(0, 0.4))

# Add a legend
#legend(x = 52.5, y = 0.62, 
   #    col = c("blue", "red", "green"),  # Match colors from the plot
    #   lty = 1:3,                      # Line types
     #  legend = colnames(fitted_values),  # Labels for the categories
      # lwd = 2)

# Add dashed gridlines
abline(v = seq(min(Malawi_EHT_data$Total), max(Malawi_EHT_data$Total), by = 5), 
       h = seq(0, 1, by = 0.1), 
       col = "gray", lty = "dashed")






#########################
#now iwant to see the difference in gambiae and funestus both for washed andunwashed
###############################

funestus_unwashed <- read_excel("Wezzie-Tom/EHT DATA/Malawi_EHT_data_funestus.xlsx", 
                                       sheet = "unwashed")




funestus_unwashed$compartment <- as.factor(funestus_unwashed$compartment)
funestus_unwashed$Total <- as.numeric(funestus_unwashed$Total)

str(funestus_unwashed)
funestus_unwashed$compartment <- factor(funestus_unwashed$compartment)




fit3.ms <- vglm(compartment ~ Total  , multinomial(refLevel = 1),
               data = funestus_unwashed)



summary(fit3.ms)





# Get the fitted values from the model
fitted_values <- fitted(fit3.ms)

# Order the data based on the 'Total' variable
ordered_data <- order(funestus_unwashed$Total)

# Plot the fitted probabilities
matplot(funestus_unwashed$Total[ordered_data], fitted_values[ordered_data, ], 
        type = "l", las = 1, lwd = 2, ylim = c(0, 1),
        ylab = "Fitted probabilities", 
        xlab = "Total",
        col = c("blue", "red", "green"),
        main = "Fitted Probabilities by Compartment")  # Adding a title)  # Change colors as needed


legend("bottomright", 
       col = c("blue", "red", "green"),  # Match colors from the plot
       lty = 1:3,                       # Line types (one for each compartment)
       legend = c("Net", "Room", "Veranda"),  # Labels for the compartments
       lwd = 2,
       xpd = TRUE,                      # Allow legend to be drawn outside the plot area
       inset = c(0, 0.6))

# Add a legend
#legend(x = 52.5, y = 0.62, 
#    col = c("blue", "red", "green"),  # Match colors from the plot
#   lty = 1:3,                      # Line types
#  legend = colnames(fitted_values),  # Labels for the categories
# lwd = 2)

# Add dashed gridlines
abline(v = seq(min(funestus_unwashed$Total), max(funestus_unwashed$Total), by = 5), 
       h = seq(0, 1, by = 0.1), 
       col = "gray", lty = "dashed")



###############funestus unwashed########
##########################################

gambiae_unwashed <- read_excel("Wezzie-Tom/EHT DATA/Malawi_EHT_data_gambiae.xlsx", 
                              sheet = "unwashed")






gambiae_unwashed$compartment <- as.factor(gambiae_unwashed$compartment)
gambiae_unwashed$Total <- as.numeric(gambiae_unwashed$Total)

str(gambiae_unwashed)
gambiae_unwashed$compartment <- factor(gambiae_unwashed$compartment)




fit4.ms <- vglm(compartment ~ Total  , multinomial(refLevel = 1),
                data = gambiae_unwashed)



summary(fit4.ms)





# Get the fitted values from the model
fitted_values <- fitted(fit4.ms)

# Order the data based on the 'Total' variable
ordered_data <- order(gambiae_unwashed$Total)

# Plot the fitted probabilities
matplot(gambiae_unwashed$Total[ordered_data], fitted_values[ordered_data, ], 
        type = "l", las = 1, lwd = 2, ylim = c(0, 1),
        ylab = "Fitted probabilities", 
        xlab = "Total",
        col = c("blue", "red", "green"),
        main = "Fitted Probabilities by Compartment")  # Adding a title)  # Change colors as needed


legend("bottomright", 
       col = c("blue", "red", "green"),  # Match colors from the plot
       lty = 1:3,                       # Line types (one for each compartment)
       legend = c("Net", "Room", "Veranda"),  # Labels for the compartments
       lwd = 2,
       xpd = TRUE,                      # Allow legend to be drawn outside the plot area
       inset = c(0, 0.4))

# Add a legend
#legend(x = 52.5, y = 0.62, 
#    col = c("blue", "red", "green"),  # Match colors from the plot
#   lty = 1:3,                      # Line types
#  legend = colnames(fitted_values),  # Labels for the categories
# lwd = 2)

# Add dashed gridlines
abline(v = seq(min(gambiae_unwashed$Total), max(gambiae_unwashed$Total), by = 5), 
       h = seq(0, 1, by = 0.1), 
       col = "gray", lty = "dashed")














############combined plots for unwashed#############



# Set up the plotting area into 1 row and 2 columns
par(mfrow = c(1, 2)) 

# First plot (Funestus Unwashed)
funestus_unwashed <- read_excel("Wezzie-Tom/EHT DATA/Malawi_EHT_data_funestus.xlsx", sheet = "unwashed")
funestus_unwashed$compartment <- as.factor(funestus_unwashed$compartment)
funestus_unwashed$Total <- as.numeric(funestus_unwashed$Total)
fit3.ms <- vglm(compartment ~ Total, multinomial(refLevel = 1), data = funestus_unwashed)
fitted_values <- fitted(fit3.ms)
ordered_data <- order(funestus_unwashed$Total)

matplot(funestus_unwashed$Total[ordered_data], fitted_values[ordered_data, ], 
        type = "l", las = 1, lwd = 2, ylim = c(0, 1),
        ylab = "Fitted probabilities", 
        xlab = "Total",
        col = c("blue", "red", "green"),
        main = "Funestus Unwashed")  # Title for the first plot

legend("bottomright", 
       col = c("blue", "red", "green"), 
       lty = 1:3,                       
       legend = c("Net", "Room", "Veranda"),  
       lwd = 2, 
       xpd = TRUE, 
       inset = c(0, 0.6))

# Second plot (Gambiae Unwashed)
gambiae_unwashed <- read_excel("Wezzie-Tom/EHT DATA/Malawi_EHT_data_gambiae.xlsx", sheet = "unwashed")
gambiae_unwashed$compartment <- as.factor(gambiae_unwashed$compartment)
gambiae_unwashed$Total <- as.numeric(gambiae_unwashed$Total)
fit4.ms <- vglm(compartment ~ Total, multinomial(refLevel = 1), data = gambiae_unwashed)
fitted_values <- fitted(fit4.ms)
ordered_data <- order(gambiae_unwashed$Total)

matplot(gambiae_unwashed$Total[ordered_data], fitted_values[ordered_data, ], 
        type = "l", las = 1, lwd = 2, ylim = c(0, 1),
        ylab = "Fitted probabilities", 
        xlab = "Total",
        col = c("blue", "red", "green"),
        main = "Gambiae Unwashed")  # Title for the second plot

legend("bottomright", 
       col = c("blue", "red", "green"), 
       lty = 1:3,                       
       legend = c("Net", "Room", "Veranda"),  
       lwd = 2, 
       xpd = TRUE, 
       inset = c(0, 0.4))

# Reset the plot layout to the default (1 plot at a time)
par(mfrow = c(1, 1))




##############################now for washed





# Set up the plotting area into 1 row and 2 columns
par(mfrow = c(1, 2)) 

# First plot (Funestus Unwashed)
funestus_unwashed <- read_excel("Wezzie-Tom/EHT DATA/Malawi_EHT_data_funestus.xlsx", sheet = "washed")
funestus_unwashed$compartment <- as.factor(funestus_unwashed$compartment)
funestus_unwashed$Total <- as.numeric(funestus_unwashed$Total)
fit3.ms <- vglm(compartment ~ Total, multinomial(refLevel = 1), data = funestus_unwashed)
fitted_values <- fitted(fit3.ms)
ordered_data <- order(funestus_unwashed$Total)

matplot(funestus_unwashed$Total[ordered_data], fitted_values[ordered_data, ], 
        type = "l", las = 1, lwd = 2, ylim = c(0, 1),
        ylab = "Fitted probabilities", 
        xlab = "Total",
        col = c("blue", "red", "green"),
        main = "Funestus washed")  # Title for the first plot

legend("bottomright", 
       col = c("blue", "red", "green"), 
       lty = 1:3,                       
       legend = c("Net", "Room", "Veranda"),  
       lwd = 2, 
       xpd = TRUE, 
       inset = c(0, 0.3))

# Second plot (Gambiae Unwashed)
gambiae_unwashed <- read_excel("Wezzie-Tom/EHT DATA/Malawi_EHT_data_gambiae.xlsx", sheet = "washed")
gambiae_unwashed$compartment <- as.factor(gambiae_unwashed$compartment)
gambiae_unwashed$Total <- as.numeric(gambiae_unwashed$Total)
fit4.ms <- vglm(compartment ~ Total, multinomial(refLevel = 1), data = gambiae_unwashed)
fitted_values <- fitted(fit4.ms)
ordered_data <- order(gambiae_unwashed$Total)

matplot(gambiae_unwashed$Total[ordered_data], fitted_values[ordered_data, ], 
        type = "l", las = 1, lwd = 2, ylim = c(0, 1),
        ylab = "Fitted probabilities", 
        xlab = "Total",
        col = c("blue", "red", "green"),
        main = "Gambiae washed")  # Title for the second plot

legend("bottomright", 
       col = c("blue", "red", "green"), 
       lty = 1:3,                       
       legend = c("Net", "Room", "Veranda"),  
       lwd = 2, 
       xpd = TRUE, 
       inset = c(0, 0.3))

# Reset the plot layout to the default (1 plot at a time)
par(mfrow = c(1, 1))













