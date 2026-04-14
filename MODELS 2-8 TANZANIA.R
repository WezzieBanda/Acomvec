



############################################################################################NOW FOR TANZANIA ##############################################################################################################################
############################################################################################NOW FOR TANZANIA ##############################################################################################################################
############################################################################################NOW FOR TANZANIA ##############################################################################################################################
############################################################################################NOW FOR TANZANIA ##############################################################################################################################






######model 8 for Tanzania













library(malariasimulation)
library(curl)
library(site)
#rm(list = ls())
set.seed(123)
library(cali)
library(malariasimulation)
library(ggplot2)
library(dplyr)
library(ggplot2)
human_population<- 10000
sim_length<- 10 *365
year<- 365
month<- year/12



##desired is 0.427

target_days <- c(  1400)  # mid-year of years 1, 2, 3
#target_pfpr <- c(  0.41)  # desired model 1
#target_pfpr <- c(  0.42527) # model 2
#target_pfpr <- c(  0.425)
#target_pfpr <- c(  0.424209) #model 6 0.42429
#target_pfpr <- c( 0.42232) #model 7 0.4223222
#target_pfpr <- c(  0.422) #model 4
#target_pfpr <- c(  0.42212) #model9

#target_pfpr <- c(  0.422222) #model8 best mode

#target_pfpr <- c(  0.42261) #model 5  

target_pfpr <- c( 0.42527 ) #model4

point_pfpr_summary <- function(x, days = target_days) {
  pfpr <- x$n_detect_lm_182.5_5109  / x$n_age_182.5_5109
  sapply(days, function(d) {
    pfpr[which.min(abs(x$timestep - d))]
  })
}


#TZA <- fetch_site(iso3c = "TZA", admin_level = 1, urban_rural = TRUE)
#TZA_Mwanza <- TZA|> subset_site(site_filter = TZA$sites |> filter(name_1 == "Mwanza")|> filter(urban_rural =="rural"))


# Convert site information to malariasimulation parameters
site_par <- site_parameters(
  interventions = TZA_Mwanza$interventions,
  demography = TZA_Mwanza$demography,
  vectors = TZA_Mwanza$vectors$vector_species,
  seasonality = TZA_Mwanza$seasonality$seasonality_parameters,
  #eir = TZA_Mwanza$eir$eir,
  overrides = list(
    human_population =human_population
  )
) |> set_epi_outputs(clinical_incidence = c(0.5, 14)*365,
                     prevalence = c(0.5, 14)*365 ) # (# add prevalence here


## -------------------------------
## Mosquito species (ORDER MATTERS)


sprayingtimesteps <- c(1)  # A round of IRS is implemented in the 1st and second year 3 months prior to peak transmission.


sprayingparams <- set_spraying(
  site_par,
  timesteps = sprayingtimesteps,
  coverages = c(0.2), # # The first round covers 30% of the population and the second covers 80%. 
  ls_theta = matrix(2.025, nrow=length(sprayingtimesteps), ncol=3), # Matrix of mortality parameters; nrows=length(timesteps), ncols=length(species) 
  ls_gamma = matrix(-0.014, nrow=length(sprayingtimesteps), ncol=3), # Matrix of mortality parameters per round of IRS and per species
  ks_theta = matrix(-2.519, nrow=length(sprayingtimesteps), ncol=3), # Matrix of feeding success parameters per round of IRS and per species
  ks_gamma = matrix(0.012, nrow=length(sprayingtimesteps), ncol=3), # Matrix of feeding success parameters per round of IRS and per species
  ms_theta = matrix(-2.651, nrow=length(sprayingtimesteps), ncol=3), # Matrix of deterrence parameters per round of IRS and per species
  ms_gamma = matrix(0.004, nrow=length(sprayingtimesteps), ncol=3) # Matrix of deterrence parameters per round of IRS and per species
)



##set bednes

#parameter



#parameters <- get_parameters()

bednetstimesteps <- c(1, 1096, 1489)

#parameters <- get_parameters()

parameters <- set_bednets(
  sprayingparams,
  timesteps = bednetstimesteps,
  coverages = c(0.50,0.2, 1),  # Each round is distributed to 50% of the population.
  retention = 10^10, # Nets are kept on average 5 years
  dn0 = matrix(c(0.13, 0.13, 0.13,
                 0.13, 0.13, 0.13,
                 0.28, 0.28, 0.28),nrow = 3,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.72, 0.72,0.72,
                0.72, 0.72,0.72,
                0.59, 0.59, 0.59 ),nrow = 3,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24,
                 0.24, 0.24, 0.24,
                 0.24, 0.24, 0.24),nrow = 3,ncol = 3,byrow = TRUE),
  gamman = rep(c(1.85, 1.85, 2.56) *365)  # Vector of bed net half-lives for each distribution timestep
)

#####in best model8  i used 1.90, 1.90, 2.65

parameters$timesteps <- 365 * 8  # total simulation length

# -----------------------------
# 5. Calibration (EIR estimated internally)
# -----------------------------
out <- calibrate(
  parameters = parameters,
  target = target_pfpr,
  summary_function = point_pfpr_summary,
  eq_prevalence = target_pfpr[1]   # baseline starting PfPR
)

# -----------------------------
# 6. Run simulation using calibrated EIR
# -----------------------------

parameters <- set_equilibrium(parameters,init_EIR = out)

raw <- run_simulation(
  parameters$timesteps + 100,
  parameters = parameters
)

raw$pfpr <- raw$n_detect_lm_182.5_5109 / raw$n_age_182.5_5109 

# -----------------------------
# 7. Plot results
# -----------------------------
ggplot() +
  geom_point(aes(x = target_days, y = target_pfpr), col = "dodgerblue", size = 4) +
  geom_line(data = raw, aes(x = timestep, y = pfpr), col = "deeppink", linewidth = 1) +
  ylim(0, 1) +
  ylab(expression(italic(Pf)*Pr[2-10])) +
  xlab("Time (days)") +
  theme_bw()

# -----------------------------
# 8. Check calibrated PfPR at target days
# -----------------------------
point_pfpr_summary(raw)


parameters$timesteps <- 365 * 10  # total simulation length

simparams <- set_equilibrium(parameters = site_par, init_EIR = out)

output_control4 <- run_simulation(timesteps = parameters$timesteps, parameters = simparams)


output_bednet4 <- run_simulation(
  timesteps = parameters$timesteps,
  parameters = parameters
)



sim_length <- 8 * 365  # 10 years
cols <- c("darkgreen", "aquamarine3", "#E69F00")  # colors for line

plot_combined_prevalence4 <- function() {
  
  bednetstimesteps <- c(1, 1489)
  
  slice_start <- 1096
  
  prev_mean <- output_bednet4$n_detect_lm_182.5_5109 / output_bednet4$n_age_182.5_5109
  prev_lower <- pmax(0, prev_mean * 0.90)
  prev_upper <- pmin(1, prev_mean * 1.05)
  
  # 🔧 FIX: limit to 2022
  keep <- output_bednet4$timestep >= slice_start & 
    output_bednet4$timestep <= 2922
  
  timesteps <- output_bednet4$timestep[keep]
  prev_mean <- prev_mean[keep]
  prev_lower <- prev_lower[keep]
  prev_upper <- prev_upper[keep]
  
  plot(timesteps, prev_mean,
       type = "n",
       xlab = "Time (Years)", ylab = "Prevalence",
       xaxt = "n",
       ylim = c(0, 0.7),
       xaxs = "i", yaxs = "i",
       # 🔧 FIX: end at 2022
       xlim = c(min(timesteps), 2922),
       main = "Estimated Prevalence (June 2019 onwards)",
       cex.main = 1.4, cex.lab = 1.2)
  
  polygon(
    c(timesteps, rev(timesteps)),
    c(prev_lower, rev(prev_upper)),
    col = adjustcolor(cols[1], alpha.f = 0.15),
    border = NA
  )
  
  lines(timesteps, prev_mean, col = cols[1], lwd = 2)
  
  years <- 2015:2022
  year_days <- seq(1, sim_length, by = 365)
  year_days_visible <- year_days[year_days >= slice_start & year_days <= 2922]
  years_visible <- years[year_days >= slice_start & year_days <= 2922]
  
  axis(1, at = year_days_visible, labels = years_visible)
  
  bednet_visible <- bednetstimesteps[bednetstimesteps >= slice_start]
  if(length(bednet_visible) > 0){
    abline(v = bednet_visible, col = "darkgray", lty = "dashed", lwd = 1.5)
    text(bednet_visible + 30, 0.65, "Bed Net Dist.", pos = 4, cex = 0.8, srt = 90)
  }
  
  obs_x <- c(1400, 1854, 2034, 2220, 2398)
  obs_y <- c(0.420, 0.156, 0.4092, 0.256, 0.492)
  
  obs_lower <- c(0.398, 0.136, 0.382, 0.223, 0.475)
  obs_upper <- c(0.456, 0.178, 0.437, 0.281, 0.51)
  
  # 🔧 FIX: keep points within 2022
  obs_keep <- obs_x >= slice_start & obs_x <= 2922
  
  if(any(obs_keep)){
    points(obs_x[obs_keep], obs_y[obs_keep], col = "darkgreen", pch = 15, cex = 1)
    
    arrows(
      x0 = obs_x[obs_keep],
      y0 = obs_lower[obs_keep],
      x1 = obs_x[obs_keep],
      y1 = obs_upper[obs_keep],
      code = 3, angle = 90, length = 0.05, col = "darkgreen", lwd = 1.5
    )
  }
  
  legend(x = min(timesteps) + 200, y = 0.65,
         legend = c("Model prevalence ±30% uncertainty", "Observed Chlorfenapyr prevalence ±95% CI"),
         col = c(cols[1], "darkgreen"),
         lty = c(1, NA),
         lwd = c(2, 1.5),
         pch = c(NA, 15),
         pt.cex = 1,
         box.lty = 0)
}

plot_combined_prevalence4()





##########great now lets do IG1 for Tanzania 
#################################################################################################################################




library(malariasimulation)
library(curl)
library(site)
#rm(list = ls())
set.seed(123)
library(cali)
library(malariasimulation)
library(ggplot2)
library(dplyr)
library(ggplot2)
human_population<- 10000
sim_length<- 10 *365
year<- 365
month<- year/12



#i want 0.459



target_days <- c(  1400)  # mid-year of years 1, 2, 3
#target_pfpr <- c(  0.444)  # desired  # model 1
#target_pfpr <- c(  0.4609) #model2 0.46
#target_pfpr <- c(  0.46091) #model6 0.4819
#target_pfpr <- c(  0.462799) #model 7

#target_pfpr <- c(  0.462799) #model6

#target_pfpr <- c(  0.469209) #model 5 
target_pfpr <- c(  0.4609) #model 4 0.454
#target_pfpr <- c(  0.454) #model 9



#target_pfpr <- c(  0.4606) #model 8

point_pfpr_summary <- function(x, days = target_days) {
  pfpr <- x$n_detect_lm_182.5_5109  / x$n_age_182.5_5109
  sapply(days, function(d) {
    pfpr[which.min(abs(x$timestep - d))]
  })
}


#TZA <- fetch_site(iso3c = "TZA", admin_level = 1, urban_rural = TRUE)
TZA_Mwanza <- TZA|> subset_site(site_filter = TZA$sites |> filter(name_1 == "Mwanza")|> filter(urban_rural =="rural"))


# Convert site information to malariasimulation parameters
site_par <- site_parameters(
  interventions = TZA_Mwanza$interventions,
  demography = TZA_Mwanza$demography,
  vectors = TZA_Mwanza$vectors$vector_species,
  seasonality = TZA_Mwanza$seasonality$seasonality_parameters,
  #eir = TZA_Mwanza$eir$eir,
  overrides = list(
    human_population =human_population
  )
) |> set_epi_outputs(clinical_incidence = c(0.5, 14)*365,
                     prevalence = c(0.5, 14)*365 ) # (# add prevalence here


## -------------------------------
## Mosquito species (ORDER MATTERS)


sprayingtimesteps <- c(1)  # A round of IRS is implemented in the 1st and second year 3 months prior to peak transmission.


sprayingparams <- set_spraying(
  site_par,
  timesteps = sprayingtimesteps,
  coverages = c(0.2), # # The first round covers 30% of the population and the second covers 80%. 
  ls_theta = matrix(2.025, nrow=length(sprayingtimesteps), ncol=3), # Matrix of mortality parameters; nrows=length(timesteps), ncols=length(species) 
  ls_gamma = matrix(-0.014, nrow=length(sprayingtimesteps), ncol=3), # Matrix of mortality parameters per round of IRS and per species
  ks_theta = matrix(-2.519, nrow=length(sprayingtimesteps), ncol=3), # Matrix of feeding success parameters per round of IRS and per species
  ks_gamma = matrix(0.012, nrow=length(sprayingtimesteps), ncol=3), # Matrix of feeding success parameters per round of IRS and per species
  ms_theta = matrix(-2.651, nrow=length(sprayingtimesteps), ncol=3), # Matrix of deterrence parameters per round of IRS and per species
  ms_gamma = matrix(0.004, nrow=length(sprayingtimesteps), ncol=3) # Matrix of deterrence parameters per round of IRS and per species
)



bednetstimesteps <- c(1, 1096, 1489)

#parameters <- get_parameters()

parameters <- set_bednets(
  sprayingparams,
  timesteps = bednetstimesteps,
  coverages = c(0.50,0.2, 1),  # Each round is distributed to 50% of the population.
  retention = 10^10, # Nets are kept on average 5 years
  dn0 = matrix(c(0.13, 0.13, 0.13,
                 0.13, 0.13, 0.13,
                 0.13, 0.13, 0.13),nrow = 3,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.72, 0.72,0.72,
                0.72, 0.72,0.72,
                0.72, 0.72,0.72 ),nrow = 3,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24,
                 0.24, 0.24, 0.24,
                 0.24, 0.24, 0.24),nrow = 3,ncol = 3,byrow = TRUE),
  gamman = rep(c(1.85, 1.85, 1.85) *365)  # Vector of bed net half-lives for each distribution timestep
)



parameters$timesteps <- 365 * 8  # total simulation length

# -----------------------------
# 5. Calibration (EIR estimated internally)
# -----------------------------
out <- calibrate(
  parameters = parameters,
  target = target_pfpr,
  summary_function = point_pfpr_summary,
  eq_prevalence = target_pfpr[1]   # baseline starting PfPR
)

# -----------------------------
# 6. Run simulation using calibrated EIR
# -----------------------------

parameters <- set_equilibrium(parameters,init_EIR = out)

raw <- run_simulation(
  parameters$timesteps + 100,
  parameters = parameters
)

raw$pfpr <- raw$n_detect_lm_182.5_5109 / raw$n_age_182.5_5109 

# -----------------------------
# 7. Plot results
# -----------------------------
ggplot() +
  geom_point(aes(x = target_days, y = target_pfpr), col = "dodgerblue", size = 4) +
  geom_line(data = raw, aes(x = timestep, y = pfpr), col = "deeppink", linewidth = 1) +
  ylim(0, 1) +
  ylab(expression(italic(Pf)*Pr[2-10])) +
  xlab("Time (days)") +
  theme_bw()

# -----------------------------
# 8. Check calibrated PfPR at target days
# -----------------------------
point_pfpr_summary(raw)



parameters$timesteps <- 365 * 10  # total simulation length

simparams <- set_equilibrium(parameters = site_par, init_EIR = out)

output_control5 <- run_simulation(timesteps = parameters$timesteps, parameters = simparams)


output_bednet5 <- run_simulation(
  timesteps = parameters$timesteps,
  parameters = parameters
)



sim_length <- 8 * 365  # 10 years
cols <- c("darkblue", "aquamarine3", "#E69F00")  # colors for line

plot_combined_prevalence5 <- function() {
  
  bednetstimesteps <- c(1, 1489)
  
  slice_start <- 1096
  
  prev_mean <- output_bednet5$n_detect_lm_182.5_5109 / output_bednet5$n_age_182.5_5109
  prev_lower <- pmax(0, prev_mean * 0.90)
  prev_upper <- pmin(1, prev_mean * 1.05)
  
  # 🔧 FIX: limit to 2022
  keep <- output_bednet5$timestep >= slice_start & 
    output_bednet5$timestep <= 2922
  
  timesteps <- output_bednet4$timestep[keep]
  prev_mean <- prev_mean[keep]
  prev_lower <- prev_lower[keep]
  prev_upper <- prev_upper[keep]
  
  plot(timesteps, prev_mean,
       type = "n",
       xlab = "Time (Years)", ylab = "Prevalence",
       xaxt = "n",
       ylim = c(0, 0.7),
       xaxs = "i", yaxs = "i",
       # 🔧 FIX: end at 2022
       xlim = c(min(timesteps), 2922),
       main = "Estimated Prevalence (June 2019 onwards)",
       cex.main = 1.4, cex.lab = 1.2)
  
  polygon(
    c(timesteps, rev(timesteps)),
    c(prev_lower, rev(prev_upper)),
    col = adjustcolor(cols[1], alpha.f = 0.15),
    border = NA
  )
  
  lines(timesteps, prev_mean, col = cols[1], lwd = 2)
  
  years <- 2015:2022
  year_days <- seq(1, sim_length, by = 365)
  year_days_visible <- year_days[year_days >= slice_start & year_days <= 2922]
  years_visible <- years[year_days >= slice_start & year_days <= 2922]
  
  axis(1, at = year_days_visible, labels = years_visible)
  
  bednet_visible <- bednetstimesteps[bednetstimesteps >= slice_start]
  if(length(bednet_visible) > 0){
    abline(v = bednet_visible, col = "darkgray", lty = "dashed", lwd = 1.5)
    text(bednet_visible + 30, 0.65, "Bed Net Dist.", pos = 4, cex = 0.8, srt = 90)
  }
  
  obs_x <- c(1400, 1854, 2034, 2220, 2398)
  obs_y <- c(0.459, 0.312, 0.523, 0.458, 0.530)
  
  obs_lower <- c(0.450, 0.285, 0.495, 0.430, 0.520)
  obs_upper <- c(0.488, 0.340, 0.551, 0.486, 0.547)
  
  # 🔧 FIX: keep points within 2022
  obs_keep <- obs_x >= slice_start & obs_x <= 2922
  
  if(any(obs_keep)){
    points(obs_x[obs_keep], obs_y[obs_keep], col = "darkblue", pch = 15, cex = 1)
    
    arrows(
      x0 = obs_x[obs_keep],
      y0 = obs_lower[obs_keep],
      x1 = obs_x[obs_keep],
      y1 = obs_upper[obs_keep],
      code = 3, angle = 90, length = 0.05, col = "darkblue", lwd = 1.5
    )
  }
  
  #legend(x = min(timesteps) + 200, y = 0.65,
  #    legend = c("Model prevalence ±30% uncertainty", "Observed Chlorfenapyr prevalence ±95% CI"),
  #    col = c(cols[1], "darkblue"),
  #    lty = c(1, NA),
  #   lwd = c(2, 1.5),
  #   pch = c(NA, 15),
  #   pt.cex = 1,
  #   box.lty = 0)
}

plot_combined_prevalence5()





##########great now lets do OP/PBO for Tanzania 
#################################################################################################################################




library(malariasimulation)
library(curl)
library(site)
#rm(list = ls())
set.seed(123)
library(cali)
library(malariasimulation)
library(ggplot2)
library(dplyr)
library(ggplot2)
human_population<- 10000
sim_length<- 10 *365
year<- 365
month<- year/12



#desired 0.42


target_days <- c(  1400)  # mid-year of years 1, 2, 3
#target_pfpr <- c(  0.425)  # desired 0.42 model 1
#target_pfpr <- c(  0.426)  # desired 0.42 model 2

#target_pfpr <- c(  0.425) #model 6
#target_pfpr <- c(  0.4295) #model 7 0.4295

#target_pfpr <- c(  0.4299) #model 5 
#target_pfpr <- c(  0.426) #model 4 0.425
#target_pfpr <- c(  0.425) #model 9

#target_pfpr <- c(  0.4279) #model 8 # this again is supposed to be model 9

#target_pfpr <- c(  0.4205) #REALMODEL 8

target_pfpr <- c(  0.422) #0.425 0.42261 0.422

point_pfpr_summary <- function(x, days = target_days) {
  pfpr <- x$n_detect_lm_182.5_5109  / x$n_age_182.5_5109
  sapply(days, function(d) {
    pfpr[which.min(abs(x$timestep - d))]
  })
}


#TZA <- fetch_site(iso3c = "TZA", admin_level = 1, urban_rural = TRUE)
#TZA_Mwanza <- TZA|> subset_site(site_filter = TZA$sites |> filter(name_1 == "Mwanza")|> filter(urban_rural =="rural"))


# Convert site information to malariasimulation parameters
site_par <- site_parameters(
  interventions = TZA_Mwanza$interventions,
  demography = TZA_Mwanza$demography,
  vectors = TZA_Mwanza$vectors$vector_species,
  seasonality = TZA_Mwanza$seasonality$seasonality_parameters,
  #eir = TZA_Mwanza$eir$eir,
  overrides = list(
    human_population =human_population
  )
) |> set_epi_outputs(clinical_incidence = c(0.5, 14)*365,
                     prevalence = c(0.5, 14)*365 ) # (# add prevalence here


## -------------------------------
## Mosquito species (ORDER MATTERS)


sprayingtimesteps <- c(1)  # A round of IRS is implemented in the 1st and second year 3 months prior to peak transmission.


sprayingparams <- set_spraying(
  site_par,
  timesteps = sprayingtimesteps,
  coverages = c(0.2), # # The first round covers 30% of the population and the second covers 80%. 
  ls_theta = matrix(2.025, nrow=length(sprayingtimesteps), ncol=3), # Matrix of mortality parameters; nrows=length(timesteps), ncols=length(species) 
  ls_gamma = matrix(-0.014, nrow=length(sprayingtimesteps), ncol=3), # Matrix of mortality parameters per round of IRS and per species
  ks_theta = matrix(-2.519, nrow=length(sprayingtimesteps), ncol=3), # Matrix of feeding success parameters per round of IRS and per species
  ks_gamma = matrix(0.012, nrow=length(sprayingtimesteps), ncol=3), # Matrix of feeding success parameters per round of IRS and per species
  ms_theta = matrix(-2.651, nrow=length(sprayingtimesteps), ncol=3), # Matrix of deterrence parameters per round of IRS and per species
  ms_gamma = matrix(0.004, nrow=length(sprayingtimesteps), ncol=3) # Matrix of deterrence parameters per round of IRS and per species
)



##set bednes

#parameter



#parameters <- get_parameters()

bednetstimesteps <- c(1, 1096, 1489)

#parameters <- get_parameters()

parameters <- set_bednets(
  sprayingparams,
  timesteps = bednetstimesteps,
  coverages = c(0.50,0.2, 0.78),  # Each round is distributed to 50% of the population.
  retention = 10^10, # Nets are kept on average 5 years
  dn0 = matrix(c(0.13, 0.13, 0.13,
                 0.13, 0.13, 0.13,
                 0.27, 0.27, 0.27),nrow = 3,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.72, 0.72,0.72,
                0.72, 0.72,0.72,
                0.58, 0.58,0.58 ),nrow = 3,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24,
                 0.24, 0.24, 0.24,
                 0.24, 0.24, 0.24),nrow = 3,ncol = 3,byrow = TRUE),
  gamman = rep(c(1.85, 1.85, 2.10) *365)  # Vector of bed net half-lives for each distribution timestep
)



parameters$timesteps <- 365 * 8  # total simulation length

# -----------------------------
# 5. Calibration (EIR estimated internally)
# -----------------------------
out <- calibrate(
  parameters = parameters,
  target = target_pfpr,
  summary_function = point_pfpr_summary,
  eq_prevalence = target_pfpr[1]   # baseline starting PfPR
)

# -----------------------------
# 6. Run simulation using calibrated EIR
# -----------------------------

parameters <- set_equilibrium(parameters,init_EIR = out)

raw <- run_simulation(
  parameters$timesteps + 100,
  parameters = parameters
)

raw$pfpr <- raw$n_detect_lm_182.5_5109 / raw$n_age_182.5_5109 

# -----------------------------
# 7. Plot results
# -----------------------------
ggplot() +
  geom_point(aes(x = target_days, y = target_pfpr), col = "dodgerblue", size = 4) +
  geom_line(data = raw, aes(x = timestep, y = pfpr), col = "deeppink", linewidth = 1) +
  ylim(0, 1) +
  ylab(expression(italic(Pf)*Pr[2-10])) +
  xlab("Time (days)") +
  theme_bw()

# -----------------------------
# 8. Check calibrated PfPR at target days
# -----------------------------
point_pfpr_summary(raw)


parameters$timesteps <- 365 * 8  # total simulation length

simparams <- set_equilibrium(parameters = site_par, init_EIR = out)

output_control6 <- run_simulation(timesteps = parameters$timesteps, parameters = simparams)


output_bednet6 <- run_simulation(
  timesteps = parameters$timesteps,
  parameters = parameters
)



sim_length <- 8 * 365  # 10 years
cols <- c( "aquamarine3", "#E69F00")  # colors for line

plot_combined_prevalence6 <- function() {
  
  bednetstimesteps <- c(1, 1489)
  
  slice_start <- 1096
  
  prev_mean <- output_bednet6$n_detect_lm_182.5_5109 / output_bednet6$n_age_182.5_5109
  prev_lower <- pmax(0, prev_mean * 0.90)
  prev_upper <- pmin(1, prev_mean * 1.05)
  
  # 🔧 FIX: limit to 2022
  keep <- output_bednet6$timestep >= slice_start & 
    output_bednet6$timestep <= 2922
  
  timesteps <- output_bednet6$timestep[keep]
  prev_mean <- prev_mean[keep]
  prev_lower <- prev_lower[keep]
  prev_upper <- prev_upper[keep]
  
  plot(timesteps, prev_mean,
       type = "n",
       xlab = "Time (Years)", ylab = "Prevalence",
       xaxt = "n",
       ylim = c(0, 0.7),
       xaxs = "i", yaxs = "i",
       # 🔧 FIX: end at 2022
       xlim = c(min(timesteps), 2922),
       main = "Estimated Prevalence (June 2019 onwards)",
       cex.main = 1.4, cex.lab = 1.2)
  
  polygon(
    c(timesteps, rev(timesteps)),
    c(prev_lower, rev(prev_upper)),
    col = adjustcolor(cols[1], alpha.f = 0.15),
    border = NA
  )
  
  lines(timesteps, prev_mean, col = cols[1], lwd = 2)
  
  years <- 2015:2022
  year_days <- seq(1, sim_length, by = 365)
  year_days_visible <- year_days[year_days >= slice_start & year_days <= 2922]
  years_visible <- years[year_days >= slice_start & year_days <= 2922]
  
  axis(1, at = year_days_visible, labels = years_visible)
  
  bednet_visible <- bednetstimesteps[bednetstimesteps >= slice_start]
  if(length(bednet_visible) > 0){
    abline(v = bednet_visible, col = "darkgray", lty = "dashed", lwd = 1.5)
    text(bednet_visible + 30, 0.65, "Bed Net Dist.", pos = 4, cex = 0.8, srt = 90)
  }
  
  obs_x <- c(1400, 1854, 2034, 2220, 2398)
  obs_y <- c(0.423, 0.19, 0.433, 0.407, 0.492)
  
  obs_lower <- c(0.391, 0.170, 0.405, 0.380, 0.485)
  obs_upper <- c(0.450, 0.214, 0.461, 0.434, 0.518)
  
  # 🔧 FIX: keep points within 2022
  obs_keep <- obs_x >= slice_start & obs_x <= 2922
  
  if(any(obs_keep)){
    points(obs_x[obs_keep], obs_y[obs_keep], col = "aquamarine3", pch = 15, cex = 1)
    
    arrows(
      x0 = obs_x[obs_keep],
      y0 = obs_lower[obs_keep],
      x1 = obs_x[obs_keep],
      y1 = obs_upper[obs_keep],
      code = 3, angle = 90, length = 0.05, col = "aquamarine3", lwd = 1.5
    )
  }
  
  #legend(x = min(timesteps) + 200, y = 0.65,
  #    legend = c("Model prevalence ±30% uncertainty", "Observed Chlorfenapyr prevalence ±95% CI"),
  #    col = c(cols[1], "darkblue"),
  #    lty = c(1, NA),
  #   lwd = c(2, 1.5),
  #   pch = c(NA, 15),
  #   pt.cex = 1,
  #   box.lty = 0)
}

plot_combined_prevalence6()





##########great now lets do RG for Tanzania 
#################################################################################################################################

############extract prevalence

# -----------------------------
# Observed days
# -----------------------------
obs_x <- c( 1854, 2034, 2220, 2398)  # the days you want to extract

# -----------------------------
# Function to extract prevalence
# -----------------------------
extract_pfpr <- function(sim_output, days) {
  pfpr <- sim_output$n_detect_lm_182.5_5109 / sim_output$n_age_182.5_5109
  
  data.frame(
    day = days,
    prevalence = sapply(days, function(d) {
      pfpr[which.min(abs(sim_output$timestep - d))]
    })
  )
}

# -----------------------------
# Extract for each model
# -----------------------------
df4 <- extract_pfpr(output_bednet4, obs_x)
df4$model <- "Interceptor G2"

df5 <- extract_pfpr(output_bednet5, obs_x)
df5$model <- "Interceptor "

df6 <- extract_pfpr(output_bednet6, obs_x)
df6$model <- "Olyset plus"

#df7 <- extract_pfpr(output_bednet7, obs_x)
#df7$model <- "Royal guard"

# -----------------------------
# Combine model outputs
# -----------------------------
df_models <- dplyr::bind_rows(df4, df5, df6)

# -----------------------------
# View results
# -----------------------------
print(df_models)

# -----------------------------
# Save to Excel
# -----------------------------
library(openxlsx)
write.xlsx(df_models, "4_model.xlsx", rowNames = FALSE)






################33now lets do RMSE for Tanzania




# -----------------------------
# 1. Load libraries
# -----------------------------
library(readxl)
library(dplyr)
library(openxlsx)




df <- read_excel("Wezzie-Tom/EHT DATA/Benin Trial/benin trial data results/New folder/feeding outcomes/model 9/redo models1-9.xlsx", 
                             sheet = "model4_real")

# 3. Remove Royal Guard
# -----------------------------
df_filtered <- df %>%
  filter(Treatment != "Royal guard", Treatment != "Royal Guard", Site != "Benin") #, Treatment != "Olyset plus")

# Check
unique(df_filtered$Treatment)

# -----------------------------
# 4. Compute RMSE
# -----------------------------
rmse_results <- df_filtered %>%
  group_by(Treatment) %>%
  summarise(
    RMSE = sqrt(mean((predicted - observed)^2, na.rm = TRUE)),
    n_points = n()
  )

# View results
print(rmse_results)


# 4. Compute OVERALL RMSE
# -----------------------------
rmse_overall <- sqrt(mean((df_filtered$predicted - df_filtered$observed)^2, 
                          na.rm = TRUE))

# Print result
rmse_overall




########real RMSE


# Load packages
library(dplyr)
library(openxlsx)

# -----------------------------

df <- read_excel("Wezzie-Tom/EHT DATA/Benin Trial/benin trial data results/New folder/feeding outcomes/model 9/redo models1-9.xlsx", 
                 sheet = "Model1_real_again")


# -----------------------------
# 2. Remove "Royal guard" treatment
# -----------------------------
df_filtered <- df %>%
  filter(Treatment != "Royal guard", Treatment != "Royal Guard")

# -----------------------------
# 3. Compute squared error per observation
# -----------------------------
df_filtered <- df_filtered %>%
  mutate(squared_error = (predicted - observed)^2)

# -----------------------------
# 4. Compute RMSE per treatment
# -----------------------------
rmse_per_treatment <- df_filtered %>%
  group_by(Treatment) %>%
  summarise(
    RMSE = sqrt(mean(squared_error, na.rm = TRUE)),
    n_points = n()
  )

print("RMSE per treatment:")
print(rmse_per_treatment)

# -----------------------------
# 5. Compute overall RMSE (all rows combined)
# -----------------------------
rmse_overall <- sqrt(mean(df_filtered$squared_error, na.rm = TRUE))
print(paste("Overall RMSE:", rmse_overall))




##########################3LOSS


library(dplyr)
library(readxl)





df <- read_excel("Wezzie-Tom/EHT DATA/Benin Trial/benin trial data results/New folder/feeding outcomes/model 9/redo models1-9.xlsx", 
                 sheet = "model3_real")


# -----------------------------
# 2. Remove Royal Guard
# -----------------------------
df_filtered <- df %>%
  filter(Treatment != "Royal guard", Treatment != "Royal Guard")

# -----------------------------
# 3. Compute LOSS (squared error)
# -----------------------------
df_filtered <- df_filtered %>%
  mutate(loss = (predicted - observed)^2)

# -----------------------------
# 4. Loss per treatment
# -----------------------------
loss_per_treatment <- df_filtered %>%
  group_by(Treatment) %>%
  summarise(
    mean_loss = mean(loss, na.rm = TRUE),
    RMSE = sqrt(mean_loss),
    n_points = n()
  )

print("Loss (and RMSE) per treatment:")
print(loss_per_treatment)

# -----------------------------
# 5. Overall loss for this model
# -----------------------------
overall_loss <- mean(df_filtered$loss, na.rm = TRUE)
overall_rmse <- sqrt(overall_loss)

print(paste("Overall Loss (MSE):", overall_loss))
print(paste("Overall RMSE:", overall_rmse))






#######################NOW ESTIMATING DIFFERENCE INMODELS


library(dplyr)
library(tidyr)
library(MCS)


df <- read_excel("Wezzie-Tom/EHT DATA/Benin Trial/benin trial data results/New folder/feeding outcomes/model 9/redo models1-9.xlsx", 
                 sheet = "LOSS_combined")


df <- df %>%
  mutate(loss = (predicted - observed)^2)

df_loss <- df %>%
  select(day, Model, loss)


loss_matrix <- df_loss %>%
  pivot_wider(names_from = Model, values_from = loss)


loss_matrix <- as.matrix(loss_matrix[,-1])  # remove 'day' column


loss_matrix <- na.omit(loss_matrix)

# =========================
# 8. Run MCS procedure
# =========================
mcs_result <- MCSprocedure(
  Loss = loss_matrix,
  alpha = 0.05
)

# =========================
# 9. View results
# =========================
print(mcs_result)

# Optional: see selected models
print(mcs_result$Included)



# Aggregate first
df_loss <- df %>%
  mutate(loss = (predicted - observed)^2) %>%
  group_by(day, Model) %>%
  summarise(loss = mean(loss, na.rm = TRUE), .groups = "drop")

# Convert to wide format (clean)
loss_matrix <- df_loss %>%
  pivot_wider(names_from = Model, values_from = loss)

# Convert to numeric matrix
loss_matrix <- as.data.frame(loss_matrix)

# Remove day column
loss_matrix_only <- loss_matrix[,-1]

# Ensure numeric
loss_matrix_only <- as.matrix(sapply(loss_matrix_only, as.numeric))

# Remove NA rows
loss_matrix_only <- na.omit(loss_matrix_only)

library(MCS)

mcs_result <- MCSprocedure(
  Loss = loss_matrix_only,
  alpha = 0.05
)

print(mcs_result)


##########################extract prevalence on specific days


# -----------------------------
# Load libraries
# -----------------------------
library(dplyr)
library(openxlsx)

# -----------------------------
# 1. Define days
# -----------------------------
net_day <- 1489

# Your follow-up days (already correct)
obs_days <- c(1854, 2034, 2220, 2398)

# -----------------------------
# 2. Function to extract prevalence
# -----------------------------
extract_pfpr_custom <- function(sim_output, obs_days, net_day) {
  
  pfpr <- sim_output$n_detect_lm_182.5_5109 / sim_output$n_age_182.5_5109
  time <- sim_output$timestep
  
  # -----------------------------
  # Point estimates (12,18,24,30 months)
  # -----------------------------
  point_estimates <- sapply(obs_days, function(d) {
    pfpr[which.min(abs(time - d))]
  })
  
  # -----------------------------
  # Interval averages
  # -----------------------------
  month_days <- 365 / 12
  
  get_avg <- function(end_m) {
    start_day <- net_day
    end_day   <- net_day + end_m * month_days
    
    idx <- time >= start_day & time <= end_day
    mean(pfpr[idx], na.rm = TRUE)
  }
  
  interval_estimates <- c(
    get_avg(24),
    get_avg(36)
  )
  
  # -----------------------------
  # Combine results
  # -----------------------------
  data.frame(
    time_period = c("12 months", "18 months", "24 months", "30 months",
                    "0-24 months", "0-36 months"),
    prevalence = c(point_estimates, interval_estimates)
  )
}

# -----------------------------
# 3. Apply to models
# -----------------------------
df1 <- extract_pfpr_custom(output_bednet4, obs_days, net_day)
df1$model <- "Interceptor G2"

df2 <- extract_pfpr_custom(output_bednet5, obs_days, net_day)
df2$model <- "Interceptor"

df3 <- extract_pfpr_custom(output_bednet6, obs_days, net_day)
df3$model <- "Olyset Plus"

# -----------------------------
# 4. Combine
# -----------------------------
df_all <- bind_rows(df1, df2, df3)

# -----------------------------
# 5. View
# -----------------------------
print(df_all)

# -----------------------------
# 6. Save
# -----------------------------
write.xlsx(df_all, "Model8_prevalence_Tanzania_combined.xlsx", rowNames = FALSE)












########let metry model 9 and see if it is better 


