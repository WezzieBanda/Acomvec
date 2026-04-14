

###############best model 8 for Benin ##############################################################################
###############best model 8 for Benin ##############################################################################
###############best model 8 for Benin ##############################################################################
###############best model 8 for Benin ##############################################################################
###############best model 8 for Benin ##############################################################################

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



#target_days <- c(1 , 2190)  # mid-year of years 1, 2, 3

target_days <- c(1 , 2190)  # mid-year of years 1, 2, 3

#target_pfpr <- c(0.42, 0.407)  MODEL 1
#target_pfpr <- c(0.42, 0.407)  # desired PfPR at these days #model 2


#target_pfpr <- c(0.42, 0.407)  # desired PfPR at these days , model 4

#target_pfpr <- c(0.42, 0.407)  # desired PfPR at these days , model 4
#model 4(0.410, 0.407)
#target_pfpr <- c(0.425, 0.407)  # desired PfPR at these days model 5
#target_pfpr <- c(0.42, 0.407)  # desired PfPR at these days model 6

#target_pfpr <- c(0.430, 0.407)  # desired PfPR at these days model 7

#target_pfpr <- c(0.430, 0.407)  # desired PfPR at these days model 8

#target_pfpr <- c(0.430, 0.407)  # desired PfPR at these days model 9

target_pfpr <- c(0.425, 0.407)  # desired PfPR at these days model 8 real
#target_pfpr <- c(0.435, 0.407)  # desired PfPR at these days  real real model 9


point_pfpr_summary <- function(x, days = target_days) {
  pfpr <- x$n_detect_lm_0_36499 / x$n_age_0_36499
  sapply(days, function(d) {
    pfpr[which.min(abs(x$timestep - d))]
  })
}



#available_sites()
BEN <- fetch_site(iso3c = "BEN", admin_level = 1, urban_rural = TRUE)
BEN_Zou <- BEN |> subset_site(site_filter = BEN$sites |> filter(name_1 == "Zou")|> filter(urban_rural =="rural"))




# Convert site information to malariasimulation parameters
site_par <- site_parameters(
  interventions = BEN_Zou$interventions,
  demography = BEN_Zou$demography,
  vectors = BEN_Zou$vectors$vector_species,
  seasonality = BEN_Zou$seasonality$seasonality_parameters,
  #eir = BEN_Zou$eir$eir,
  overrides = list(
    human_population = human_population
  )
) |> set_epi_outputs(
  clinical_incidence = c(0.5, 10) * 365,
  prevalence = c(0, 100)*365 # add prevalence here
)




bednetstimesteps <- c(1, 2312)

parameters <- get_parameters()

parameters <- set_bednets(
  site_par,
  timesteps = bednetstimesteps,
  coverages = c(0.40, 0.58),  # Each round is distributed to 50% of the population.
  retention = 10^10, # Nets are kept on average 5 years
  dn0 = matrix(c(0.34, 0.34, 0.34,
                 0.54, 0.54, 0.54),nrow = 2,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.63, 0.63,0.63,
                0.45, 0.45, 0.45),nrow = 2,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24,
                 0.24, 0.24, 0.24),nrow = 2,ncol = 3,byrow = TRUE),
  gamman = rep(c(1.90 ,2.69)*365) # Vector of bed net half-lives for each distribution timestep
)





parameters$timesteps <- 365 * 10  # total simulation length

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

raw$pfpr <- raw$n_detect_lm_0_36499 / raw$n_age_0_36499

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

output_control1 <- run_simulation(timesteps = parameters$timesteps, parameters = simparams)


output_bednet1 <- run_simulation(
  timesteps = parameters$timesteps,
  parameters = parameters
)



sim_length<- 10*365



sim_length <- 10 * 365  # 10 years
cols <- c("darkgreen", "aquamarine3", "#E69F00")  # colors for line


plot_combined_prevalence1 <- function() {
  # -----------------------------
  # Bed net distribution parameters
  # -----------------------------
  bednetstimesteps <- c(1, 2310)
  
  # -----------------------------
  # Slice start day: June 1, 2019
  # -----------------------------
  slice_start <- 1977
  
  # -----------------------------
  # Model prevalence ±30% uncertainty
  # -----------------------------
  prev_mean <- output_bednet1$n_detect_lm_0_36499 / output_bednet1$n_age_0_36499
  prev_lower <- pmax(0, prev_mean * 0.80)
  prev_upper <- pmin(1, prev_mean * 1.2)
  
  # -----------------------------
  # Slice to show only from slice_start
  # -----------------------------
  keep <- output_bednet1$timestep >= slice_start
  timesteps <- output_bednet1$timestep[keep]
  prev_mean <- prev_mean[keep]
  prev_lower <- prev_lower[keep]
  prev_upper <- prev_upper[keep]
  
  # -----------------------------
  # Plot setup
  # -----------------------------
  plot(timesteps, prev_mean,
       type = "n",
       xlab = "Time (Years)", ylab = "Prevalence",
       xaxt = "n",
       ylim = c(0, 0.7),
       xaxs = "i", yaxs = "i",
       xlim = c(min(timesteps), max(timesteps)),
       main = "Estimated Prevalence (June 2019 onwards)",
       cex.main = 1.4, cex.lab = 1.2)
  
  # -----------------------------
  # Shaded model uncertainty band
  # -----------------------------
  polygon(
    c(timesteps, rev(timesteps)),
    c(prev_lower, rev(prev_upper)),
    col = adjustcolor(cols[1], alpha.f = 0.15),
    border = NA
  )
  
  # -----------------------------
  # Mean prevalence line
  # -----------------------------
  lines(timesteps, prev_mean, col = cols[1], lwd = 2)
  
  # -----------------------------
  # Custom x-axis: yearly ticks
  # -----------------------------
  years <- 2014:2023
  year_days <- seq(1, sim_length, by = 365)
  year_days_visible <- year_days[year_days >= slice_start]
  years_visible <- years[year_days >= slice_start]
  axis(1, at = year_days_visible, labels = years_visible)
  
  # -----------------------------
  # Bed net distribution lines
  # -----------------------------
  bednet_visible <- bednetstimesteps[bednetstimesteps >= slice_start]
  if(length(bednet_visible) > 0){
    abline(v = bednet_visible, col = "darkgray", lty = "dashed", lwd = 1.5)
    text(bednet_visible + 30, 0.65, "Bed Net Dist.", pos = 4, cex = 0.8, srt = 90)
  }
  
  # -----------------------------
  # Observed points + 95% CIs
  # -----------------------------
  # Chlorfenapyr–pyrethroid observed prevalence
  obs_x <- c(2190, 2496, 2859, 3195)            # 6 months and 18 months days
  obs_y <- c(0.407, 0.157, 0.279, 0.255)    # prevalence
  # Calculate Wilson 95% CI manually (rounded)
  obs_lower <- c(0.382, 0.139, 0.257, 0.235)      # 95% CI lower
  obs_upper <- c(0.433, 0.177, 0.302, 0.28 )      # 95% CI upper
  
  #c(2190, 2496, 2859, 3195)  
  #interceptor G2 (0.407, 0.157, 0.279, 0.255)
  #obs_x <- c(2190, 2496, 2859, 3195)
  #obs_y <- c(0.407, 0.157, 0.279, 0.225)
  obs_keep <- obs_x >= slice_start
  if(any(obs_keep)){
    # Points
    points(obs_x[obs_keep], obs_y[obs_keep], col = "darkgreen", pch = 15, cex = 1)
    # Error bars
    arrows(
      x0 = obs_x[obs_keep],
      y0 = obs_lower[obs_keep],
      x1 = obs_x[obs_keep],
      y1 = obs_upper[obs_keep],
      code = 3, angle = 90, length = 0.05, col = "darkgreen", lwd = 1.5
    )
  }
  
  # -----------------------------
  # Legend
  # -----------------------------
  legend(x = min(timesteps) + 200, y = 0.65,
         legend = c("Model prevalence ±30% uncertainty", "Observed Chlorfenapyr prevalence ±95% CI"),
         col = c(cols[1], "darkgreen"),
         lty = c(1, NA),
         lwd = c(2, 1.5),
         pch = c(NA, 15),
         pt.cex = 1,
         box.lty = 0)
}

# Run the plot
plot_combined_prevalence1()









######################################Now for IG1 






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



target_days <- c(1 , 2190)  # mid-year of years 1, 2, 3
#target_pfpr <- c(0.48, 0.465)  # desired PfPR at these days

#target_pfpr <- c(0.475, 0.465) model 1
#0.468/
#target_pfpr <- c(0.47, 0.460)  # desired PfPR at these days model 4

#target_pfpr <- c(0.485, 0.465) ####realmodel 4

#target_pfpr <- c(0.485, 0.460) #model 5

#target_pfpr <- c(0.480, 0.465)#model 6

#target_pfpr <- c(0.480, 0.465)#model 6

#target_pfpr <- c(0.480, 0.465) model 7

target_pfpr <- c(0.450, 0.47) #model 8

#target_pfpr <- c(0.450, 0.465) model 9
#target_pfpr <- c(0.475, 0.465) #model 2

#target_pfpr <- c(0.475, 0.455) #model 7 real
#target_pfpr <- c(0.465, 0.46) #model 9 rel real



point_pfpr_summary <- function(x, days = target_days) {
  pfpr <- x$n_detect_lm_0_36499 / x$n_age_0_36499
  sapply(days, function(d) {
    pfpr[which.min(abs(x$timestep - d))]
  })
}



#available_sites()
#BEN <- fetch_site(iso3c = "BEN", admin_level = 1, urban_rural = TRUE)
#BEN_Zou <- BEN |> subset_site(site_filter = BEN$sites |> filter(name_1 == "Zou")|> filter(urban_rural =="rural"))


# Convert site information to malariasimulation parameters
site_par <- site_parameters(
  interventions = BEN_Zou$interventions,
  demography = BEN_Zou$demography,
  vectors = BEN_Zou$vectors$vector_species,
  seasonality = BEN_Zou$seasonality$seasonality_parameters,
  #eir = BEN_Zou$eir$eir,
  overrides = list(
    human_population = human_population
  )
) |> set_epi_outputs(
  clinical_incidence = c(0.5, 10) * 365,
  prevalence = c(0, 100)*365 # add prevalence here
)




bednetstimesteps <- c(1, 2312)

parameters <- get_parameters()

parameters <- set_bednets(
  site_par,
  timesteps = bednetstimesteps,
  coverages = c(0.58, 0.58),  # Each round is distributed to 50% of the population.
  retention = 10^10, # Nets are kept on average 5 years
  dn0 = matrix(c(0.34, 0.34, 0.34,
                 0.34, 0.34, 0.34),nrow = 2,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.63, 0.63,0.63,
                0.63, 0.63, 0.63),nrow = 2,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24,
                 0.24, 0.24, 0.24),nrow = 2,ncol = 3,byrow = TRUE),
  gamman = rep(c(1.90, 1.90) * 365) # Vector of bed net half-lives for each distribution timestep
)




parameters$timesteps <- 365 * 10  # total simulation length

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

raw$pfpr <- raw$n_detect_lm_0_36499 / raw$n_age_0_36499

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

output_control2 <- run_simulation(timesteps = parameters$timesteps, parameters = simparams)


output_bednet2 <- run_simulation(
  timesteps = parameters$timesteps,
  parameters = parameters
)

sim_length<- 10*365




sim_length <- 10 * 365  # 10 years
cols <- c("darkblue", "aquamarine3", "#E69F00")  # colors for line




plot_combined_prevalence2 <- function() {
  # -----------------------------
  # Bed net distribution parameters
  # -----------------------------
  bednetstimesteps <- c(1, 2310)
  
  # -----------------------------
  # Slice start day: June 1, 2019
  # -----------------------------slice_star
  slice_start <- 1977
  
  # -----------------------------
  # Model prevalence ±30% uncertainty
  # -----------------------------
  prev_mean <- output_bednet2$n_detect_lm_0_36499 / output_bednet2$n_age_0_36499
  prev_lower <- pmax(0, prev_mean * 0.90)
  prev_upper <- pmin(1, prev_mean * 1.10)
  
  # -----------------------------
  # Slice to show only from slice_start
  # -----------------------------
  keep <- output_bednet2$timestep >= slice_start
  timesteps <- output_bednet2$timestep[keep]
  prev_mean <- prev_mean[keep]
  prev_lower <- prev_lower[keep]
  prev_upper <- prev_upper[keep]
  
  # -----------------------------
  # Plot setup
  # -----------------------------
  plot(timesteps, prev_mean,
       type = "n",
       xlab = "Time (Years)", ylab = "Prevalence",
       xaxt = "n",
       ylim = c(0, 0.7),
       xaxs = "i", yaxs = "i",
       xlim = c(min(timesteps), max(timesteps)),
       main = "Estimated Prevalence (June 2019 onwards)",
       cex.main = 1.4, cex.lab = 1.2)
  
  # -----------------------------
  # Shaded model uncertainty band
  # -----------------------------
  polygon(
    c(timesteps, rev(timesteps)),
    c(prev_lower, rev(prev_upper)),
    col = adjustcolor(cols[1], alpha.f = 0.15),
    border = NA
  )
  
  # -----------------------------
  # Mean prevalence line
  # -----------------------------
  lines(timesteps, prev_mean, col = cols[1], lwd = 1)
  
  # -----------------------------
  # Custom x-axis: yearly ticks
  # -----------------------------
  years <- 2014:2023
  year_days <- seq(1, sim_length, by = 365)
  year_days_visible <- year_days[year_days >= slice_start]
  years_visible <- years[year_days >= slice_start]
  axis(1, at = year_days_visible, labels = years_visible)
  
  # -----------------------------
  # Bed net distribution lines
  # -----------------------------
  bednet_visible <- bednetstimesteps[bednetstimesteps >= slice_start]
  if(length(bednet_visible) > 0){
    abline(v = bednet_visible, col = "darkgray", lty = "dashed", lwd = 1)
    text(bednet_visible + 30, 0.65, "Bed Net Dist.", pos = 4, cex = 0.8, srt = 90)
  }
  
  # -----------------------------
  # Observed points + 95% CIs
  
  obs_x <- c(2190, 2496, 2859)
  obs_y <- c(0.465, 0.28, 0.387)
  obs_keep <- obs_x >= slice_start
  # -----------------------------
  # Chlorfenapyr–pyrethroid observed prevalence
  obs_x <- c(2190, 2496, 2859, 3195)            # 6 months and 18 months days
  obs_y <- c(0.465, 0.28, 0.387, 0.262)    # prevalence
  # Calculate Wilson 95% CI manually (rounded)
  obs_lower <- c(0.439, 0.257, 0.362, 0.235)      # 95% CI lower
  obs_upper <- c(0.490, 0.303, 0.412, 0.28 )      # 95% CI upper
  
 
  #obs_x <- c(2190, 2496, 2859, 3195)
  #obs_y <- c(0.407, 0.157, 0.279, 0.225)
  obs_keep <- obs_x >= slice_start
  if(any(obs_keep)){
    # Points
    points(obs_x[obs_keep], obs_y[obs_keep], col = "darkblue", pch = 15, cex = 1)
    # Error bars
    arrows(
      x0 = obs_x[obs_keep],
      y0 = obs_lower[obs_keep],
      x1 = obs_x[obs_keep],
      y1 = obs_upper[obs_keep],
      code = 3, angle = 90, length = 0.05, col = "darkblue", lwd = 1
    )
  }
  
  # -----------------------------
  # Legend
  # -----------------------------
  legend(x = min(timesteps) + 200, y = 0.65,
         legend = c("Model prevalence ±30% uncertainty", "Observed Chlorfenapyr prevalence ±95% CI"),
         col = c(cols[1], "darkblue"),
         lty = c(1, NA),
         lwd = c(2, 1.5),
         pch = c(NA, 15),
         pt.cex = 1,
         box.lty = 0)
}

# Run the plot
plot_combined_prevalence2()



###############best model 8 for Benin ##############################################################################
###############best model 8 for Benin ##############################################################################
###############best model 8 for Benin ##############################################################################
###############best model 8 for Benin ##############################################################################
###############best model 8 for Benin ##############################################################################
############################################################# now for RG










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



target_days <- c(1 , 2190)  # mid-year of years 1, 2, 3
#target_pfpr <- c(0.465, 0.431)  # desired PfPR at these days

#0.468

#target_pfpr <- c(0.445, 0.435)  # desired PfPR at these days model 1
#target_pfpr <- c(0.465, 0.42)  # desired PfPR at these days #model 4
#target_pfpr <- c(0.46, 0.422) model 5

#target_pfpr <- c(0.46, 0.43) model 6

target_pfpr <- c(0.45, 0.435) 

point_pfpr_summary <- function(x, days = target_days) {
  pfpr <- x$n_detect_lm_0_36499 / x$n_age_0_36499
  sapply(days, function(d) {
    pfpr[which.min(abs(x$timestep - d))]
  })
}



#available_sites()
#BEN <- fetch_site(iso3c = "BEN", admin_level = 1, urban_rural = TRUE)
#BEN_Zou <- BEN |> subset_site(site_filter = BEN$sites |> filter(name_1 == "Zou")|> filter(urban_rural =="rural"))


# Convert site information to malariasimulation parameters
site_par <- site_parameters(
  interventions = BEN_Zou$interventions,
  demography = BEN_Zou$demography,
  vectors = BEN_Zou$vectors$vector_species,
  seasonality = BEN_Zou$seasonality$seasonality_parameters,
  #eir = BEN_Zou$eir$eir,
  overrides = list(
    human_population = human_population
  )
) |> set_epi_outputs(
  clinical_incidence = c(0.5, 10) * 365,
  prevalence = c(0, 100)*365 # add prevalence here
)

#0.492684101	0.207918568	0.299397331



bednetstimesteps <- c(1, 2312)

parameters <- get_parameters()

parameters <- set_bednets(
  site_par,
  timesteps = bednetstimesteps,
  coverages = c(0.58, 0.58),  # Each round is distributed to 50% of the population.
  retention = 10^10, # Nets are kept on average 5 years
  dn0 = matrix(c(0.11, 0.11, 0.11,
                 0.30, 0.30, 0.30),nrow = 2,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.64, 0.64,0.64,
                0.492684101, 0.492684101, 0.492684101),nrow = 2,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24,
                 0.24, 0.24, 0.24),nrow = 2,ncol = 3,byrow = TRUE),
  gamman = rep(c(2.64, 2.64)*365) # Vector of bed net half-lives for each distribution timestep
)




parameters$timesteps <- 365 * 10  # total simulation length

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

raw$pfpr <- raw$n_detect_lm_0_36499 / raw$n_age_0_36499

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

output_control3 <- run_simulation(timesteps = parameters$timesteps, parameters = simparams)


output_bednet3 <- run_simulation(
  timesteps = parameters$timesteps,
  parameters = parameters
)

sim_length<- 10*365







sim_length <- 10 * 365  # 10 years
cols <- c("darkred", "aquamarine3", "#E69F00")  # colors for line




plot_combined_prevalence3 <- function() {
  # -----------------------------
  # Bed net distribution parameters
  # -----------------------------
  bednetstimesteps <- c(1, 2310)
  
  # -----------------------------
  # Slice start day: June 1, 2019
  # -----------------------------slice_star
  slice_start <- 1977
  
  # -----------------------------
  # Model prevalence ±30% uncertainty
  # -----------------------------
  prev_mean <- output_bednet3$n_detect_lm_0_36499 / output_bednet3$n_age_0_36499
  prev_lower <- pmax(0, prev_mean * 0.90)
  prev_upper <- pmin(1, prev_mean * 1.05)
  
  # -----------------------------
  # Slice to show only from slice_start
  # -----------------------------
  keep <- output_bednet3$timestep >= slice_start
  timesteps <- output_bednet3$timestep[keep]
  prev_mean <- prev_mean[keep]
  prev_lower <- prev_lower[keep]
  prev_upper <- prev_upper[keep]
  
  # -----------------------------
  # Plot setup
  # -----------------------------
  plot(timesteps, prev_mean,
       type = "n",
       xlab = "Time (Years)", ylab = "Prevalence",
       xaxt = "n",
       ylim = c(0, 0.7),
       xaxs = "i", yaxs = "i",
       xlim = c(min(timesteps), max(timesteps)),
       main = "Estimated Prevalence (June 2019 onwards)",
       cex.main = 1.4, cex.lab = 1.2)
  
  # -----------------------------
  # Shaded model uncertainty band
  # -----------------------------
  polygon(
    c(timesteps, rev(timesteps)),
    c(prev_lower, rev(prev_upper)),
    col = adjustcolor(cols[1], alpha.f = 0.15),
    border = NA
  )
  
  # -----------------------------
  # Mean prevalence line
  # -----------------------------
  lines(timesteps, prev_mean, col = cols[1], lwd = 2)
  
  # -----------------------------
  # Custom x-axis: yearly ticks
  # -----------------------------
  years <- 2014:2023
  year_days <- seq(1, sim_length, by = 365)
  year_days_visible <- year_days[year_days >= slice_start]
  years_visible <- years[year_days >= slice_start]
  axis(1, at = year_days_visible, labels = years_visible)
  
  # -----------------------------
  # Bed net distribution lines
  # -----------------------------
  bednet_visible <- bednetstimesteps[bednetstimesteps >= slice_start]
  if(length(bednet_visible) > 0){
    abline(v = bednet_visible, col = "darkgray", lty = "dashed", lwd = 1.5)
    text(bednet_visible + 30, 0.65, "Bed Net Dist.", pos = 4, cex = 0.8, srt = 90)
  }
  
  # -----------------------------
  # Observed points + 95% CIs
  obs_x <- c(2190, 2496, 2859, 3195)            # 6 months and 18 months days
  obs_y <- c(0.431, 0.269, 0.382, 0.284)  
  #obs_x <- c(2190, 2496, 2859)
  #obs_y <- c(0.465, 0.28, 0.387)
  obs_keep <- obs_x >= slice_start
  # -----------------------------
  # Chlorfenapyr–pyrethroid observed prevalence
  # prevalence
  # Calculate Wilson 95% CI manually (rounded)
  obs_lower <- c(0.406, 0.247, 0.357, 0.266)      # 95% CI lower
  obs_upper <- c(0.456, 0.292, 0.407, 0.32 )      # 95% CI upper
  
  
  #obs_x <- c(2190, 2496, 2859, 3195)
  #obs_y <- c(0.407, 0.157, 0.279, 0.225)
  obs_keep <- obs_x >= slice_start
  if(any(obs_keep)){
    # Points
    points(obs_x[obs_keep], obs_y[obs_keep], col = "darkred", pch = 15, cex = 1)
    # Error bars
    arrows(
      x0 = obs_x[obs_keep],
      y0 = obs_lower[obs_keep],
      x1 = obs_x[obs_keep],
      y1 = obs_upper[obs_keep],
      code = 3, angle = 90, length = 0.05, col = "darkred", lwd = 1.5
    )
  }
  
  # -----------------------------
  # Legend
  # -----------------------------
  legend(x = min(timesteps) + 200, y = 0.65,
         legend = c("Model prevalence ±30% uncertainty", "Observed Chlorfenapyr prevalence ±95% CI"),
         col = c(cols[1], "darkred"),
         lty = c(1, NA),
         lwd = c(2, 1.5),
         pch = c(NA, 15),
         pt.cex = 1,
         box.lty = 0)
}

# Run the plot
plot_combined_prevalence3()

























#################################################plot these three plese 
########333anotherone








plot_combined_prevalence_all <- function() {
  # -----------------------------
  slice_start <- 1977
  cols <- c("darkgreen", "darkblue",  "black", "darkred")
  
  # -----------------------------
  # Extract timesteps and prevalence for all three simulations
  timesteps1 <- output_bednet1$timestep[output_bednet1$timestep >= slice_start]
  prev_mean1 <- (output_bednet1$n_detect_lm_0_36499 / output_bednet1$n_age_0_36499)[output_bednet1$timestep >= slice_start]
  prev_lower1 <- pmax(0, prev_mean1 * 0.80)
  prev_upper1 <- pmin(1, prev_mean1 * 1.2)
  
  timesteps2 <- output_bednet2$timestep[output_bednet2$timestep >= slice_start]
  prev_mean2 <- (output_bednet2$n_detect_lm_0_36499 / output_bednet2$n_age_0_36499)[output_bednet2$timestep >= slice_start]
  prev_lower2 <- pmax(0, prev_mean2 * 0.90)
  prev_upper2 <- pmin(1, prev_mean2 * 1.1)
  
  #timesteps3 <- output_bednet3$timestep[output_bednet3$timestep >= slice_start]
  #prev_mean3 <- (output_bednet3$n_detect_lm_0_36499 / output_bednet3$n_age_0_36499)[output_bednet3$timestep >= slice_start]
  #prev_lower3 <- pmax(0, prev_mean3 * 0.90)
  #prev_upper3 <- pmin(1, prev_mean3 * 1.05)
  
  # -----------------------------
  # Plot setup
  plot(timesteps1, prev_mean1, type = "n",
       xlab = "Time (Years)", ylab = " Malaria Prevalence for all ages",
       xaxt = "n", ylim = c(0, 0.7),
       xaxs = "i", yaxs = "i",
       xlim = c(min(timesteps1), max(timesteps1)),
       main = "Benin Trials",
       cex.main = 1, cex.lab = 1)
  
  # -----------------------------
  # Plot uncertainty polygons
  polygon(c(timesteps1, rev(timesteps1)), c(prev_lower1, rev(prev_upper1)),
          col = adjustcolor(cols[1], alpha.f = 0.15), border = NA)
  polygon(c(timesteps2, rev(timesteps2)), c(prev_lower2, rev(prev_upper2)),
          col = adjustcolor(cols[2], alpha.f = 0.15), border = NA)
  #polygon(c(timesteps3, rev(timesteps3)), c(prev_lower3, rev(prev_upper3)),
   #       col = adjustcolor(cols[3], alpha.f = 0.15), border = NA)
  
  # -----------------------------
  # Plot mean prevalence lines
  lines(timesteps1, prev_mean1, col = cols[1], lwd = 1)
  lines(timesteps2, prev_mean2, col = cols[2], lwd = 1)
  #lines(timesteps3, prev_mean3, col = cols[3], lwd = 1)
  
  # -----------------------------
  # Add observed points with 95% CI
  obs_list <- list(
    list(x = c(2190, 2496, 2859, 3195), y = c(0.407, 0.157, 0.279, 0.225), lower = c(0.382, 0.139, 0.257, 0.20), upper = c(0.433, 0.177, 0.302, 0.28), col = cols[1]),
    list(x = c(2190, 2496, 2859, 3195), y = c(0.465, 0.28, 0.387, 0.262), lower = c(0.439, 0.257, 0.362, 0.235), upper = c(0.490, 0.303, 0.412, 0.28), col = cols[2])
    #list(x = c(2190, 2496, 2859, 3195), y = c(0.431, 0.269, 0.382, 0.284), lower = c(0.406, 0.247, 0.357, 0.266), upper = c(0.456, 0.292, 0.407, 0.30), col = cols[3])
  )
  
  for(obs in obs_list){
    keep <- obs$x >= slice_start
    if(any(keep)){
      points(obs$x[keep], obs$y[keep], col = obs$col, pch = 15, cex = 1)
      arrows(x0 = obs$x[keep], y0 = obs$lower[keep],
             x1 = obs$x[keep], y1 = obs$upper[keep],
             code = 3, angle = 90, length = 0.05, col = obs$col, lwd = 1)
    }
  }
  
  # -----------------------------
  # Custom x-axis: yearly ticks
  years <- 2014:2023
  year_days <- seq(1, sim_length, by = 365)
  year_days_visible <- year_days[year_days >= slice_start]
  years_visible <- years[year_days >= slice_start]
  axis(1, at = year_days_visible, labels = years_visible)
  
  # -----------------------------
  # Bed net distribution lines
  bednetstimesteps <- c(1, 2310)
  bednet_visible <- bednetstimesteps[bednetstimesteps >= slice_start]
  if(length(bednet_visible) > 0){
    abline(v = bednet_visible, col = "darkgray", lty = "dashed", lwd = 1)
    text(bednet_visible + 30, 0.65, "Bed Net Distribution", pos = 4, cex = 0.8, srt = 0)
    # Add short arrow pointing to the dashed line
    arrows(
      x0 = bednet_visible + 25,  # slightly left of the text
      y0 = 0.65,                 # same y-level as text
      x1 = bednet_visible,       # point to dashed line
      y1 = 0.65,                 # same y-level
      length = 0.1,              # arrowhead size
      angle = 20,
      col = "black",
      lwd = 1.2
    )
  }
  
  
  
  
  # -----------------------------
  # Legend
  legend(x = 2800, y = 0.69,
         legend = c("pyrethroid-pyrrole", "pyrethroid only", "Observed prevalence"),
         col = c(cols, "black"),   # line colors for models, black for observed points
         lty = c(1, 1, NA),     # lines for models, no line for observed
         lwd = c(3, 3, NA),     # line widths
         pch = c(NA, NA, 15),  # square point only for observed data
         pt.cex = 0.8,
         box.lty = 0)
  
}

# Run the combined plot with observed points
plot_combined_prevalence_all()




###################plot



# -----------------------------
# 1. Observed days
# -----------------------------
obs_x <- c( 2496, 2859, 3195)

# -----------------------------
# 2. Function to extract prevalence
# -----------------------------
extract_pfpr <- function(sim_output, days) {
  pfpr <- sim_output$n_detect_lm_0_36499 / sim_output$n_age_0_36499
  
  data.frame(
    day = days,
    prevalence = sapply(days, function(d) {
      pfpr[which.min(abs(sim_output$timestep - d))]
    })
  )
}

# -----------------------------
# 3. Extract for each model
# -----------------------------
df1 <- extract_pfpr(output_bednet1, obs_x)
df1$model <- "Interceptor G2"

df2 <- extract_pfpr(output_bednet2, obs_x)
df2$model <- "Interceptor"

df3 <- extract_pfpr(output_bednet3, obs_x)
df3$model <- "Royal Guard"

# -----------------------------
# 4. Combine model outputs
# -----------------------------
df_models <- dplyr::bind_rows(df1, df2 ) #, df3)

# -----------------------------
# 5. Add observed data (your values)
# -----------------------------


df_all <- dplyr::bind_rows(df_models) #, df_obs)

# -----------------------------
# 6. View results
# -----------------------------
print(df_all)

# -----------------------------
# 7. Save
# -----------------------------

library(openxlsx)
write.xlsx(df_all, "model8_2_prevalence_TZ.xlsx") #, row.names = FALSE)


















##################now estimating RMSE


# -----------------------------
# 1. Load libraries
# -----------------------------
library(readxl)
library(dplyr)
library(openxlsx)

# -----------------------------
# 2. Import data
# -----------------------------
#df <- read_excel("your_file.xlsx", sheet = "Sheet1")

df <- read_excel("~/model2_prevalence_benin.xlsx", 
                                      sheet = "model9_real")
# -----------------------------
# 3. Remove Royal Guard
# -----------------------------
df_filtered <- df %>%
  filter(Treatment != "Royal Guard", Site != "Tanzania")

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









#realreal

library(scales)

ggplot(combined_df, aes(x = Observed, y = Predicted, color = Net)) +
  
  # 1:1 reference line
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 1, colour = "black") +
  
  # Points: shape by country
  geom_point(aes(shape = Country), size = 2.5) +
  
  # Horizontal CI lines for the actual data (plot only)
  geom_segment(aes(x = ObsLower, xend = ObsUpper, y = Predicted, yend = Predicted),
               linewidth = 1, alpha = 0.5) +
  
  # Colours and shapes with legend override
  scale_color_manual(values = cols,
                     guide = guide_legend(
                       override.aes = list(shape = 16, size = 2.5, linetype = 0)
                     )) +
  scale_shape_manual(values = c("Benin" = 15, "Tanzania" = 16)) +
  
  # 🔑 Convert axes to percentages (INSERTED HERE)
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     breaks = seq(0, 0.7, 0.1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     breaks = seq(0, 0.7, 0.1)) +
  
  # Axes limits
  coord_equal(xlim = c(0, 0.7), ylim = c(0, 0.7)) +
  
  # Labels
  labs(x = "Observed Malaria prevalence (%)",
       y = "Predicted Malaria prevalence (%)",
       color = "Net type",
       shape = "Country") +
  
  # Theme + legend inside plot
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.25, 0.7),
    legend.background = element_rect(fill = "white", color = "black")
  )














ggplot(combined_df, aes(x = Observed * 100, y = Predicted * 100, color = Net)) +
  
  # 1:1 reference line
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 1, colour = "black") +
  
  # Points
  geom_point(aes(shape = Country), size = 2.5) +
  
  # Horizontal CI lines
  geom_segment(aes(x = ObsLower * 100, xend = ObsUpper * 100, 
                   y = Predicted * 100, yend = Predicted * 100),
               linewidth = 1, alpha = 0.5) +
  
  # Colours and shapes
  scale_color_manual(values = cols,
                     guide = guide_legend(
                       override.aes = list(shape = 16, size = 2.5, linetype = 0)
                     )) +
  scale_shape_manual(values = c("Benin" = 15, "Tanzania" = 16)) +
  
  # 🔑 Axis now in 0–70 (not % symbols)
  scale_x_continuous(breaks = seq(0, 70, 10), limits = c(0, 70)) +
  scale_y_continuous(breaks = seq(0, 70, 10), limits = c(0, 70)) +
  
  coord_equal() +
  
  # Labels show %
  labs(x = "Observed Malaria prevalence (%)",
       y = "Predicted Malaria prevalence (%)",
       color = "Net type",
       shape = "Country") +
  
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.25, 0.7),
    legend.background = element_rect(fill = "white", color = "black")
  )




ggplot(combined_df, aes(x = Observed * 100, y = Predicted * 100, color = Net)) +
  
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 1, colour = "black") +
  
  geom_point(aes(shape = Country), size = 2.5) +
  
  geom_segment(aes(x = ObsLower * 100, xend = ObsUpper * 100, 
                   y = Predicted * 100, yend = Predicted * 100),
               linewidth = 1, alpha = 0.5) +
  
  scale_color_manual(values = cols,
                     guide = guide_legend(
                       override.aes = list(shape = 16, size = 2.5, linetype = 0)
                     )) +
  scale_shape_manual(values = c("Benin" = 15, "Tanzania" = 16)) +
  
  scale_x_continuous(breaks = seq(0, 70, 10), limits = c(0, 70)) +
  scale_y_continuous(breaks = seq(0, 70, 10), limits = c(0, 70)) +
  
  coord_equal() +
  
  labs(x = "Observed Malaria prevalence (%)",
       y = "Predicted Malaria prevalence (%)",
       color = "Net type",
       shape = "Country") +
  
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.25, 0.7),
    legend.background = element_rect(fill = "white", color = "black"),
    
    # 🔑 Ensure same font for both axes
    axis.title.x = element_text(family = "sans", face = "plain"),
    axis.title.y = element_text(family = "sans", face = "plain")
  )




########extract at the specifictimepoints 


# -----------------------------
# Load libraries
# -----------------------------
library(dplyr)
library(openxlsx)

# -----------------------------
# 1. Your observed time points
# -----------------------------
obs_days <- c(2496, 2859, 3195)

# Define net distribution day (needed for intervals)
net_day <- 2310

# -----------------------------
# 2. Function
# -----------------------------
extract_pfpr_custom <- function(sim_output, obs_days, net_day) {
  
  pfpr <- sim_output$n_detect_lm_0_36499 / sim_output$n_age_0_36499
  time <- sim_output$timestep
  
  # -----------------------------
  # Point estimates (use your exact days)
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
    get_avg(18),
    get_avg(24),
    get_avg(36)
  )
  
  # -----------------------------
  # Combine
  # -----------------------------
  data.frame(
    time_period = c("6 months", "18 months", "30 months",
                    "0-18 months", "0-24 months", "0-36 months"),
    prevalence = c(point_estimates, interval_estimates)
  )
}

# -----------------------------
# 3. Apply to models
# -----------------------------
df1 <- extract_pfpr_custom(output_bednet1, obs_days, net_day)
df1$model <- "Interceptor G2"

df2 <- extract_pfpr_custom(output_bednet2, obs_days, net_day)
df2$model <- "Interceptor"

df3 <- extract_pfpr_custom(output_bednet3, obs_days, net_day)
df3$model <- "Royal Guard"

# -----------------------------
# 4. Combine
# -----------------------------
df_all <- bind_rows(df1, df2 ) #, df3)

# -----------------------------
# 5. View
# -----------------------------
print(df_all)

# -----------------------------
# 6. Save
# -----------------------------
write.xlsx(df_all, "bestmodelvalues_prevalence_benin_post_net.xlsx", rowNames = FALSE)

