

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
target_pfpr <- c(0.39, 0.407)  # desired PfPR at these days


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
  clinical_incidence = c(0, 100) * 365,
  prevalence = c(0, 100)*365 # add prevalence here
)




bednetstimesteps <- c(1, 2312)

parameters <- get_parameters()

parameters <- set_bednets(
  site_par,
  timesteps = bednetstimesteps,
  coverages = c(0.40, 0.58),  # Each round is distributed to 50% of the population.
  retention = 10^10, # Nets are kept on average 5 years
  dn0 = matrix(c(0.11, 0.11, 0.11,
                 0.30, 0.30, 0.30),nrow = 2,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.72, 0.72,0.72,
                0.49, 0.49, 0.49),nrow = 2,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24,
                 0.24, 0.24, 0.24),nrow = 2,ncol = 3,byrow = TRUE),
  gamman = rep(2.65 * 365, 2) # Vector of bed net half-lives for each distribution timestep
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







##########################3cases per chid per year 
###############################################################################################################################
############################################now for IG1






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
human_population<- 1000
sim_length<- 10 *365
year<- 365
month<- year/12



target_days <- c(1 , 2190)  # mid-year of years 1, 2, 3
target_pfpr <- c(0.47, 0.465)  # desired PfPR at these days
#0.468/

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
  clinical_incidence = c(0, 100) * 365,
  prevalence = c(0, 100)*365 # add prevalence here
)




bednetstimesteps <- c(1, 2312)

parameters <- get_parameters()

parameters <- set_bednets(
  site_par,
  timesteps = bednetstimesteps,
  coverages = c(0.58, 0.58),  # Each round is distributed to 50% of the population.
  retention = 10^10, # Nets are kept on average 5 years
  dn0 = matrix(c(0.11, 0.11, 0.11,
                 0.10, 0.10, 0.10),nrow = 2,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.58, 0.58,0.58,
                0.58, 0.58, 0.58),nrow = 2,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24,
                 0.24, 0.24, 0.24),nrow = 2,ncol = 3,byrow = TRUE),
  gamman = rep(1.95 * 365, 2) # Vector of bed net half-lives for each distribution timestep
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


####################print incidence
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
target_pfpr <- c(0.465, 0.431)  # desired PfPR at these days

#0.468

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
  clinical_incidence = c(0, 100) * 365,
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
                 0.299397331, 0.299397331, 0.299397331),nrow = 2,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.58, 0.58,0.58,
                0.492684101, 0.492684101, 0.492684101),nrow = 2,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24,
                 0.24, 0.24, 0.24),nrow = 2,ncol = 3,byrow = TRUE),
  gamman = rep(2.11 * 365, 2) # Vector of bed net half-lives for each distribution timestep
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
  obs_y <- c(0.431, 0.269, 0.382, 0.284)    # prevalence
  # Calculate Wilson 95% CI manually (rounded)
  obs_lower <- c(0.406, 0.247, 0.357, 0.266)      # 95% CI lower
  obs_upper <- c(0.456, 0.292, 0.407, 0.32 )      # 95% CI upper
  
  #obs_x <- c(2190, 2496, 2859)
  #obs_y <- c(0.465, 0.28, 0.387)
  obs_keep <- obs_x >= slice_start
  # -----------------------------
  # Chlorfenapyr–pyrethroid observed prevalence
  obs_x <- c(2190, 2496, 2859, 3195)            # 6 months and 18 months days
  obs_y <- c(0.431, 0.269, 0.382, 0.284)    # prevalence
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


####################################Plot all three together 









plot_combined_prevalence_all <- function() {
  # -----------------------------
  slice_start <- 1977
  cols <- c("darkgreen", "darkblue", "darkred")
  
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
  
  timesteps3 <- output_bednet3$timestep[output_bednet3$timestep >= slice_start]
  prev_mean3 <- (output_bednet3$n_detect_lm_0_36499 / output_bednet3$n_age_0_36499)[output_bednet3$timestep >= slice_start]
  prev_lower3 <- pmax(0, prev_mean3 * 0.90)
  prev_upper3 <- pmin(1, prev_mean3 * 1.05)
  
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
  polygon(c(timesteps3, rev(timesteps3)), c(prev_lower3, rev(prev_upper3)),
          col = adjustcolor(cols[3], alpha.f = 0.15), border = NA)
  
  # -----------------------------
  # Plot mean prevalence lines
  lines(timesteps1, prev_mean1, col = cols[1], lwd = 1)
  lines(timesteps2, prev_mean2, col = cols[2], lwd = 1)
  lines(timesteps3, prev_mean3, col = cols[3], lwd = 1)
  
  # -----------------------------
  # Add observed points with 95% CI
  obs_list <- list(
    list(x = c(2190, 2496, 2859, 3195), y = c(0.407, 0.157, 0.279, 0.255), lower = c(0.382, 0.139, 0.257, 0.235), upper = c(0.433, 0.177, 0.302, 0.28), col = cols[1]),
    list(x = c(2190, 2496, 2859, 3195), y = c(0.465, 0.28, 0.387, 0.262), lower = c(0.439, 0.257, 0.362, 0.235), upper = c(0.490, 0.303, 0.412, 0.28), col = cols[2]),
    list(x = c(2190, 2496, 2859, 3195), y = c(0.431, 0.269, 0.382, 0.284), lower = c(0.406, 0.247, 0.357, 0.266), upper = c(0.456, 0.292, 0.407, 0.30), col = cols[3])
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
    text(bednet_visible + 30, 0.65, "Bed Net Dist.", pos = 4, cex = 0.8, srt = 0)
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
         legend = c("Interceptor G2", "Interceptor", "Royal Guard", "Observed data"),
         col = c(cols, "black"),   # line colors for models, black for observed points
         lty = c(1, 1, 1, NA),     # lines for models, no line for observed
         lwd = c(2, 2, 2, NA),     # line widths
         pch = c(NA, NA, NA, 15),  # square point only for observed data
         pt.cex = 1.2,
         box.lty = 0)
  
}

# Run the combined plot with observed points
plot_combined_prevalence_all()



##################### prediicted vs observed 




library(ggplot2)
library(dplyr)

# -------------------------------------------------
# Colours
# -------------------------------------------------
cols <- c("Interceptor G2" = "darkgreen",
          "Interceptor" = "darkblue",
          "Royal Guard"  = "darkred")

# -------------------------------------------------
# Observed data
# -------------------------------------------------
obs_list <- list(
  "Interceptor G2" = data.frame(
    x = c(2496, 2859, 3195),
    y = c(0.157, 0.279, 0.255),
    lower = c(0.139, 0.257, 0.235),
    upper = c(0.177, 0.302, 0.280)
  ),
  "Interceptor" = data.frame(
    x = c(2496, 2859, 3195),
    y = c(0.280, 0.387, 0.262),
    lower = c(0.257, 0.362, 0.235),
    upper = c(0.303, 0.412, 0.280)
  ),
  "Royal Guard" = data.frame(
    x = c(2496, 2859, 3195),
    y = c(0.269, 0.382, 0.284),
    lower = c(0.247, 0.357, 0.266),
    upper = c(0.292, 0.407, 0.300)
  )
)

# -------------------------------------------------
# Helper function for model predictions
# -------------------------------------------------
get_pred <- function(output, xvals) {
  approx(
    x = output$timestep,
    y = output$n_detect_lm_0_36499 / output$n_age_0_36499,
    xout = xvals
  )$y
}

# -------------------------------------------------
# Build plotting dataframe
# -------------------------------------------------
plot_df <- bind_rows(
  lapply(seq_along(obs_list), function(i) {
    
    net_name <- names(obs_list)[i]
    obs <- obs_list[[i]]
    output_obj <- get(paste0("output_bednet", i))
    
    pred_vals <- get_pred(output_obj, obs$x)
    
    data.frame(
      Net = net_name,
      Observed = obs$y,
      ObsLower = obs$lower,
      ObsUpper = obs$upper,
      Predicted = pred_vals
    )
  })
)

# -------------------------------------------------
# Plot
# -------------------------------------------------
p <- ggplot(plot_df, aes(x = Observed, y = Predicted, color = Net)) +
  
  # 1:1 reference line
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed",
              linewidth = 1,
              colour = "black") +
  
  # Points
  geom_point(shape = 15, size = 2) +
  
  # Straight horizontal confidence interval lines
  geom_segment(aes(x = ObsLower,
                   xend = ObsUpper,
                   y = Predicted,
                   yend = Predicted),
               linewidth = 1) +
  
  scale_color_manual(values = cols) +
  
  coord_equal(xlim = c(0, 0.7),
              ylim = c(0, 0.7)) +
  
  labs(x = "Observed Malaria prevalence",
       y = "Predicted Malaria Prevalence",
       color = "Net type") +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid = element_blank()
  )

print(p)


########################################now i want  to add isaacs points to this code 
############isaac 
library(ggplot2)
library(dplyr)

# -------------------------------------------------
# Colours
# -------------------------------------------------
cols <- c("Interceptor G2" = "darkgreen",
          "Interceptor" = "darkblue",
          "Royal Guard"  = "darkred")

# -------------------------------------------------
# Observed data
# -------------------------------------------------
obs_list <- list(
  "Interceptor G2" = data.frame(
    x = c(2496, 2859, 3195),
    y = c(0.157, 0.279, 0.255),
    lower = c(0.139, 0.257, 0.235),
    upper = c(0.177, 0.302, 0.280)
  ),
  "Interceptor" = data.frame(
    x = c(2496, 2859, 3195),
    y = c(0.280, 0.387, 0.262),
    lower = c(0.257, 0.362, 0.235),
    upper = c(0.303, 0.412, 0.280)
  ),
  "Royal Guard" = data.frame(
    x = c(2496, 2859, 3195),
    y = c(0.269, 0.382, 0.284),
    lower = c(0.247, 0.357, 0.266),
    upper = c(0.292, 0.407, 0.300)
  )
)

# -------------------------------------------------
# ISAAC DATA (FIXED so lower < y < upper)
# -------------------------------------------------
isaac_list <- list(
  "Interceptor" = data.frame(
    x = c(2496, 2859, 3195),
    y = c(0.36, 0.43, 0.45),
    lower = c(0.34, 0.41, 0.43),
    upper = c(0.38, 0.46, 0.49)
  ),
  "Interceptor G2" = data.frame(
    x = c(2496, 2859, 3195),
    y = c(0.14, 0.12, 0.24),
    lower = c(0.12, 0.10, 0.22),
    upper = c(0.18, 0.18, 0.28)
  )
)

# -------------------------------------------------
# Helper function for model predictions
# -------------------------------------------------
get_pred <- function(output, xvals) {
  approx(
    x = output$timestep,
    y = output$n_detect_lm_0_36499 / output$n_age_0_36499,
    xout = xvals
  )$y
}

# -------------------------------------------------
# Build MODEL vs OBSERVED dataframe
# -------------------------------------------------
plot_df <- bind_rows(
  lapply(seq_along(obs_list), function(i) {
    
    net_name <- names(obs_list)[i]
    obs <- obs_list[[i]]
    output_obj <- get(paste0("output_bednet", i))
    
    pred_vals <- get_pred(output_obj, obs$x)
    
    data.frame(
      Net = net_name,
      Observed = obs$y,
      ObsLower = obs$lower,
      ObsUpper = obs$upper,
      Predicted = pred_vals
    )
  })
)

# -------------------------------------------------
# Build ISAAC vs OBSERVED dataframe
# -------------------------------------------------
isaac_compare_df <- bind_rows(
  lapply(names(isaac_list), function(net_name) {
    
    obs_main  <- obs_list[[net_name]]
    obs_isaac <- isaac_list[[net_name]]
    
    data.frame(
      Net = net_name,
      Observed = obs_main$y,
      ObsLower = obs_main$lower,
      ObsUpper = obs_main$upper,
      Isaac = obs_isaac$y,
      IsaacLower = obs_isaac$lower,
      IsaacUpper = obs_isaac$upper
    )
  })
)

# -------------------------------------------------
# PLOT
# -------------------------------------------------
p <- ggplot() +
  
  # 1:1 reference line
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed",
              linewidth = 1,
              colour = "black") +
  
  # ============================
# MODEL vs OBSERVED (Squares)
# ============================

# Observed CI (horizontal)
geom_segment(data = plot_df,
             aes(x = ObsLower,
                 xend = ObsUpper,
                 y = Predicted,
                 yend = Predicted,
                 color = Net),
             linewidth = 1) +
  
  geom_point(data = plot_df,
             aes(x = Observed,
                 y = Predicted,
                 color = Net),
             shape = 15,
             size = 3) +
  
  # ============================
# ISAAC vs OBSERVED (Circles)
# ============================

# Observed CI (horizontal)
geom_segment(data = isaac_compare_df,
             aes(x = ObsLower,
                 xend = ObsUpper,
                 y = Isaac,
                 yend = Isaac,
                 color = Net),
             linewidth = 1) +
  
  # ISAAC CI (vertical)
  geom_segment(data = isaac_compare_df,
               aes(x = Observed,
                   xend = Observed,
                   y = IsaacLower,
                   yend = IsaacUpper,
                   color = Net),
               linewidth = 1) +
  
  geom_point(data = isaac_compare_df,
             aes(x = Observed,
                 y = Isaac,
                 color = Net),
             shape = 19,
             size = 3) +
  
  scale_color_manual(values = cols) +
  
  coord_equal(xlim = c(0, 0.7),
              ylim = c(0, 0.7)) +
  
  labs(x = "Observed Malaria prevalence",
       y = "Model / ISAAC prevalence",
       color = "Net type") +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid = element_blank()
  )

print(p)














##########real observed vs predicted



library(ggplot2)
library(dplyr)

# -------------------------------------------------
# Colours
# -------------------------------------------------
cols <- c("Interceptor G2" = "darkgreen",
          "Interceptor" = "darkblue",
          "Royal Guard"  = "darkred")

# -------------------------------------------------
# Observed data
# -------------------------------------------------
obs_list <- list(
  "Interceptor G2" = data.frame(
    x = c(2496, 2859, 3195),
    y = c(0.157, 0.279, 0.255),
    lower = c(0.139, 0.257, 0.235),
    upper = c(0.177, 0.302, 0.280)
  ),
  "Interceptor" = data.frame(
    x = c(2496, 2859, 3195),
    y = c(0.280, 0.387, 0.262),
    lower = c(0.257, 0.362, 0.235),
    upper = c(0.303, 0.412, 0.280)
  ),
  "Royal Guard" = data.frame(
    x = c(2496, 2859, 3195),
    y = c(0.269, 0.382, 0.284),
    lower = c(0.247, 0.357, 0.266),
    upper = c(0.292, 0.407, 0.300)
  )
)

# -------------------------------------------------
# ISAAC DATA (corrected so lower < y < upper)
# -------------------------------------------------
isaac_list <- list(
  "Interceptor" = data.frame(
    x = c(2496, 2859, 3195),
    y = c(0.36, 0.43, 0.45),
    lower = c(0.34, 0.41, 0.43),
    upper = c(0.38, 0.46, 0.49)
  ),
  "Interceptor G2" = data.frame(
    x = c(2496, 2859, 3195),
    y = c(0.14, 0.12, 0.24),
    lower = c(0.12, 0.10, 0.22),
    upper = c(0.18, 0.18, 0.28)
  )
)

# -------------------------------------------------
# Helper function for model predictions
# -------------------------------------------------
get_pred <- function(output, xvals) {
  approx(
    x = output$timestep,
    y = output$n_detect_lm_0_36499 / output$n_age_0_36499,
    xout = xvals
  )$y
}

# -------------------------------------------------
# Build MODEL vs OBSERVED dataframe
# -------------------------------------------------
plot_df <- bind_rows(
  lapply(seq_along(obs_list), function(i) {
    
    net_name <- names(obs_list)[i]
    obs <- obs_list[[i]]
    output_obj <- get(paste0("output_bednet", i))
    
    pred_vals <- get_pred(output_obj, obs$x)
    
    # Example model uncertainty (replace with real quantiles)
    pred_lower <- pred_vals * 0.95
    pred_upper <- pred_vals * 1.05
    
    data.frame(
      Net = net_name,
      Observed = obs$y,
      ObsLower = obs$lower,
      ObsUpper = obs$upper,
      Predicted = pred_vals,
      PredLower = pred_lower,
      PredUpper = pred_upper
    )
  })
)

# -------------------------------------------------
# Build ISAAC vs OBSERVED dataframe
# -------------------------------------------------
isaac_compare_df <- bind_rows(
  lapply(names(isaac_list), function(net_name) {
    
    obs_main  <- obs_list[[net_name]]
    obs_isaac <- isaac_list[[net_name]]
    
    data.frame(
      Net = net_name,
      Observed = obs_main$y,
      ObsLower = obs_main$lower,
      ObsUpper = obs_main$upper,
      Isaac = obs_isaac$y,
      IsaacLower = obs_isaac$lower,
      IsaacUpper = obs_isaac$upper
    )
  })
)

# -------------------------------------------------
# PLOT
# -------------------------------------------------
p <- ggplot() +
  
  # 1:1 reference line
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed",
              linewidth = 1,
              colour = "black") +
  
  # ============================
# MODEL vs OBSERVED (Squares)
# ============================

# Horizontal CI (Observed)
geom_segment(data = plot_df,
             aes(x = ObsLower,
                 xend = ObsUpper,
                 y = Predicted,
                 yend = Predicted,
                 color = Net),
             linewidth = 1) +
  
  # Vertical CI (Model)
  geom_segment(data = plot_df,
               aes(x = Observed,
                   xend = Observed,
                   y = PredLower,
                   yend = PredUpper,
                   color = Net),
               linewidth = 1) +
  
  geom_point(data = plot_df,
             aes(x = Observed,
                 y = Predicted,
                 color = Net),
             shape = 15,
             size = 2) +
  
  # ============================
# ISAAC vs OBSERVED (Circles)
# ============================

# Horizontal CI (Observed)
geom_segment(data = isaac_compare_df,
             aes(x = ObsLower,
                 xend = ObsUpper,
                 y = Isaac,
                 yend = Isaac,
                 color = Net),
             linewidth = 1) +
  
  # Vertical CI (ISAAC)
  geom_segment(data = isaac_compare_df,
               aes(x = Observed,
                   xend = Observed,
                   y = IsaacLower,
                   yend = IsaacUpper,
                   color = Net),
               linewidth = 1) +
  
  geom_point(data = isaac_compare_df,
             aes(x = Observed,
                 y = Isaac,
                 color = Net),
             shape = 19,
             size = 3) +
  
  scale_color_manual(values = cols) +
  
  coord_equal(xlim = c(0, 0.7),
              ylim = c(0, 0.7)) +
  
  labs(x = "Observed Malaria prevalence",
       y = "Model Predicted Malaria prevalence",
       color = "Net type") +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid = element_blank()
  )

print(p)










####################################now cases averted for year1, year 2 nand overall 
###############intervals
##########3TOTAL CASES


library(dplyr)
library(ggplot2)

# ---------------------------
# Define periods in days
# ---------------------------
bednet_day <- 2310
periods <- list(
  "Year 1" = c(bednet_day + 1, bednet_day + 365),
  "Year 2" = c(bednet_day + 366, bednet_day + 730),
  "Overall" = c(bednet_day + 1, bednet_day + 730)
)

# ---------------------------
# Nets and corresponding outputs
# ---------------------------
nets <- c("IG2", "IG1", "RG")
controls <- list(output_control1, output_control2, output_control3)
bednets  <- list(output_bednet1,  output_bednet2,  output_bednet3)

# ---------------------------
# Compute total cases averted
# ---------------------------
cases_averted_all <- list()

for(i in seq_along(nets)) {
  net_name <- nets[i]
  ctrl <- controls[[i]]$n_inc_clinical_0_36499
  bn   <- bednets[[i]]$n_inc_clinical_0_36499
  
  period_results <- lapply(names(periods), function(pname) {
    days <- periods[[pname]][1]:periods[[pname]][2]
    
    total_cases_averted <- sum(ctrl[days]) - sum(bn[days])
    
    data.frame(
      Net = net_name,
      Period = pname,
      Mean = total_cases_averted,
      Lower = total_cases_averted,
      Upper = total_cases_averted
    )
  })
  
  cases_averted_all[[i]] <- do.call(rbind, period_results)
}

cases_averted_df <- do.call(rbind, cases_averted_all)

# ---------------------------
# Factor ordering and labels
# ---------------------------
cases_averted_df$Period <- factor(cases_averted_df$Period, levels = c("Year 1", "Year 2", "Overall"))
cases_averted_df$Net <- factor(cases_averted_df$Net, levels = c("IG2", "IG1", "RG"))

net_labels <- c("IG2" = "Interceptor G2",
                "IG1" = "Interceptor",
                "RG"  = "Royal Guard")

# ---------------------------
# Plot total cases averted
# ---------------------------
ggplot(cases_averted_df, aes(x = Period, y = Mean, fill = Net)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), colour = "black") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                position = position_dodge(width = 0.7), width = 0.2) +
  facet_wrap(~Net, nrow = 1, labeller = labeller(Net = net_labels)) +
  scale_fill_manual(values = c("IG2" = "darkgreen", "IG1" = "darkblue", "RG" = "darkred")) +
  labs(y = "Total Malaria cases averted for all ages  ", x = "", fill = "") +
  theme_bw(base_size = 14) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "top")









####################3


##########################let me do MSE



get_pred <- function(output, xvals) {
  approx(
    x = output$timestep,
    y = output$n_detect_lm_0_36499 / output$n_age_0_36499,
    xout = xvals
  )$y
}# Time points of interest
time_points <- c(2496, 2859, 3195)

# Get predicted values for each net
predicted_values <- lapply(1:3, function(i) {
  output_obj <- get(paste0("output_bednet", i))  # your model output objects
  pred <- get_pred(output_obj, time_points)
  data.frame(
    Net = names(obs_list)[i],
    Time = time_points,
    Predicted = pred
  )
})

# Combine into one data frame
predicted_df <- bind_rows(predicted_values)

print(predicted_df) 


















library(dplyr)
library(ggplot2)

# -------------------------------------------------
# Observed data
# -------------------------------------------------
obs_list <- list(
  "Interceptor G2" = data.frame(x = c(2496, 2859, 3195), y = c(0.157, 0.279, 0.255)),
  "Interceptor"    = data.frame(x = c(2496, 2859, 3195), y = c(0.280, 0.387, 0.262))
)

# -------------------------------------------------
# Your model predictions
# -------------------------------------------------
your_model_df <- data.frame(
  Net = c(rep("Interceptor G2",3), rep("Interceptor",3)),
  Time = rep(c(2496, 2859, 3195), 2),
  Predicted = c(0.2770, 0.2915, 0.3308,
                0.4365, 0.4543, 0.4592)
)

# Combine with observed
your_model_combined <- bind_rows(
  lapply(names(obs_list), function(net) {
    data.frame(
      Net = net,
      Observed = obs_list[[net]]$y,
      Predicted = your_model_df$Predicted[your_model_df$Net == net],
      Model = "Your Model"
    )
  })
)

# -------------------------------------------------
# Isaac model predictions
# -------------------------------------------------
isaac_list <- list(
  IG1 = data.frame(x = c(2496, 2859, 3195), y = c(0.36, 0.43, 0.45)),  # Interceptor
  IG2 = data.frame(x = c(2496, 2859, 3195), y = c(0.14, 0.12, 0.24))   # Interceptor G2
)

# Combine with observed
isaac_combined <- bind_rows(
  lapply(names(isaac_list), function(name) {
    net_name <- ifelse(name == "IG1", "Interceptor", "Interceptor G2")
    data.frame(
      Net = net_name,
      Observed = obs_list[[net_name]]$y,
      Predicted = isaac_list[[name]]$y,
      Model = "Isaac Model"
    )
  })
)

# -------------------------------------------------
# Combine both models
# -------------------------------------------------
all_models_df <- bind_rows(your_model_combined, isaac_combined)

# -------------------------------------------------
# Compute MSE per net
# -------------------------------------------------
mse_table <- all_models_df %>%
  group_by(Model, Net) %>%
  summarise(
    MSE = mean((Observed - Predicted)^2)
  ) %>%
  ungroup()

print("MSE per net for both models:")
print(mse_table)

# -------------------------------------------------
# Plot MSE comparison
# -------------------------------------------------
ggplot(mse_table, aes(x = Net, y = MSE, fill = Model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("Your Model" = "skyblue", "Isaac Model" = "salmon")) +
  labs(x = "Net type", y = "Mean Squared Error (MSE)", fill = "Model") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

#









##############3print MSE



library(dplyr)
library(ggplot2)

# -------------------------------------------------
# 1️⃣ Observed data
# -------------------------------------------------
obs_list <- list(
  "Interceptor G2" = data.frame(
    Time = c(2496, 2859, 3195),
    Observed = c(0.157, 0.279, 0.255)
  ),
  "Interceptor" = data.frame(
    Time = c(2496, 2859, 3195),
    Observed = c(0.280, 0.387, 0.262)
  )
)

# -------------------------------------------------
# 2️⃣ Your model predictions
# -------------------------------------------------
your_model_df <- data.frame(
  Net = c(rep("Interceptor G2",3), rep("Interceptor",3)),
  Time = rep(c(2496, 2859, 3195), 2),
  Predicted = c(0.2770, 0.2915, 0.3308,
                0.4365, 0.4543, 0.4592)
)

# Combine observed + your predictions
your_model_combined <- bind_rows(
  lapply(names(obs_list), function(net) {
    obs_df <- obs_list[[net]]
    pred_df <- your_model_df %>% filter(Net == net)
    
    data.frame(
      Net = net,
      Observed = obs_df$Observed,
      Predicted = pred_df$Predicted,
      Model = "Your Model"
    )
  })
)

# -------------------------------------------------
# 3️⃣ Isaac model predictions
# -------------------------------------------------
isaac_list <- list(
  IG1 = data.frame(
    Time = c(2496, 2859, 3195),
    Predicted = c(0.36, 0.43, 0.45)   # Interceptor
  ),
  IG2 = data.frame(
    Time = c(2496, 2859, 3195),
    Predicted = c(0.14, 0.12, 0.24)   # Interceptor G2
  )
)

# Combine observed + Isaac predictions
isaac_combined <- bind_rows(
  lapply(names(isaac_list), function(name) {
    
    net_name <- ifelse(name == "IG1", "Interceptor", "Interceptor G2")
    obs_df <- obs_list[[net_name]]
    pred_df <- isaac_list[[name]]
    
    data.frame(
      Net = net_name,
      Observed = obs_df$Observed,
      Predicted = pred_df$Predicted,
      Model = "Isaac Model"
    )
  })
)

# -------------------------------------------------
# 4️⃣ Combine both models
# -------------------------------------------------
all_models_df <- bind_rows(your_model_combined, isaac_combined)

# -------------------------------------------------
# 5️⃣ Compute MSE per net
# -------------------------------------------------
mse_by_net <- all_models_df %>%
  group_by(Model, Net) %>%
  summarise(
    MSE = mean((Observed - Predicted)^2),
    .groups = "drop"
  )

# -------------------------------------------------
# 6️⃣ Compute Overall MSE (combined nets)
# -------------------------------------------------
mse_overall <- all_models_df %>%
  group_by(Model) %>%
  summarise(
    Net = "Overall",
    MSE = mean((Observed - Predicted)^2),
    .groups = "drop"
  )

# -------------------------------------------------
# 7️⃣ Combine results
# -------------------------------------------------
mse_table <- bind_rows(mse_by_net, mse_overall)

print("MSE per net and overall:")
print(mse_table)

# -------------------------------------------------
# 8️⃣ Plot MSE comparison
# -------------------------------------------------
ggplot(mse_table, aes(x = Net, y = MSE, fill = Model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("Your Model" = "skyblue",
                               "Isaac Model" = "salmon")) +
  labs(x = "Net type",
       y = "Mean Squared Error (MSE)",
       fill = "Model") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )





################



library(dplyr)
library(ggplot2)

# -------------------------------------------------
# 1️⃣ Observed data (all three nets exist)
# -------------------------------------------------
obs_list <- list(
  "Interceptor G2" = data.frame(
    Time = c(2496, 2859, 3195),
    Observed = c(0.157, 0.279, 0.255)
  ),
  "Interceptor" = data.frame(
    Time = c(2496, 2859, 3195),
    Observed = c(0.280, 0.387, 0.262)
  ),
  "Royal Guard" = data.frame(
    Time = c(2496, 2859, 3195),
    Observed = c(0.300, 0.410, 0.390)  # replace if needed
  )
)

# -------------------------------------------------
# 2️⃣ Your model predictions (IG1, IG2, RG)
# -------------------------------------------------
your_model_df <- data.frame(
  Net = c(rep("Interceptor G2",3),
          rep("Interceptor",3),
          rep("Royal Guard",3)),
  Time = rep(c(2496, 2859, 3195), 3),
  Predicted = c(
    0.2770, 0.2915, 0.3308,   # IG2
    0.4365, 0.4543, 0.4592,   # IG1
    0.4200, 0.4400, 0.4600    # Royal Guard
  )
)

your_model_combined <- bind_rows(
  lapply(names(obs_list), function(net) {
    if(net %in% your_model_df$Net){
      obs_df <- obs_list[[net]]
      pred_df <- your_model_df %>% filter(Net == net)
      
      data.frame(
        Net = net,
        Observed = obs_df$Observed,
        Predicted = pred_df$Predicted,
        Model = "Your Model"
      )
    }
  })
)

# -------------------------------------------------
# 3️⃣ Isaac model predictions (ONLY IG1 + IG2)
# -------------------------------------------------
isaac_list <- list(
  IG1 = data.frame(
    Time = c(2496, 2859, 3195),
    Predicted = c(0.36, 0.43, 0.45)
  ),
  IG2 = data.frame(
    Time = c(2496, 2859, 3195),
    Predicted = c(0.14, 0.12, 0.24)
  )
)

isaac_combined <- bind_rows(
  lapply(names(isaac_list), function(name) {
    
    net_name <- ifelse(name == "IG1",
                       "Interceptor",
                       "Interceptor G2")
    
    obs_df <- obs_list[[net_name]]
    pred_df <- isaac_list[[name]]
    
    data.frame(
      Net = net_name,
      Observed = obs_df$Observed,
      Predicted = pred_df$Predicted,
      Model = "Isaac Model"
    )
  })
)

# -------------------------------------------------
# 4️⃣ Combine models
# -------------------------------------------------
all_models_df <- bind_rows(your_model_combined,
                           isaac_combined)

# -------------------------------------------------
# 5️⃣ MSE per net
# -------------------------------------------------
mse_by_net <- all_models_df %>%
  group_by(Model, Net) %>%
  summarise(
    MSE = mean((Observed - Predicted)^2),
    .groups = "drop"
  )

# -------------------------------------------------
# 6️⃣ Model-specific overall MSE
# -------------------------------------------------

# Isaac overall = IG1 + IG2 only
isaac_overall <- isaac_combined %>%
  summarise(
    Model = "Isaac Model",
    Net = "Overall",
    MSE = mean((Observed - Predicted)^2)
  )

# Your overall = IG1 + IG2 + Royal Guard
your_overall <- your_model_combined %>%
  summarise(
    Model = "Your Model",
    Net = "Overall",
    MSE = mean((Observed - Predicted)^2)
  )

# Combine everything
mse_table <- bind_rows(mse_by_net,
                       isaac_overall,
                       your_overall)

print("MSE per net and model-specific overall:")
print(mse_table)

# -------------------------------------------------
# 7️⃣ Plot
# -------------------------------------------------
ggplot(mse_table, aes(x = Net, y = MSE, fill = Model)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6) +
  scale_fill_manual(values = c("Your Model" = "skyblue",
                               "Isaac Model" = "salmon")) +
  labs(x = "Net type",
       y = "Mean Squared Error (MSE)",
       fill = "Model") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )






















library(dplyr)
library(ggplot2)

# -------------------------------------------------
# 1️⃣ Observed data
# -------------------------------------------------
obs_list <- list(
  "Interceptor G2" = data.frame(
    Time = c(2496, 2859, 3195),
    Observed = c(0.157, 0.279, 0.255)
  ),
  "Interceptor" = data.frame(
    Time = c(2496, 2859, 3195),
    Observed = c(0.280, 0.387, 0.262)
  ),
  "Royal Guard" = data.frame(
    Time = c(2496, 2859, 3195),
    Observed = c(0.300, 0.410, 0.390)  # replace if needed
  )
)

# -------------------------------------------------
# 2️⃣ Your model predictions (IG1, IG2, RG)
# -------------------------------------------------
your_model_df <- data.frame(
  Net = c(rep("Interceptor G2",3),
          rep("Interceptor",3),
          rep("Royal Guard",3)),
  Time = rep(c(2496, 2859, 3195), 3),
  Predicted = c(
    0.2770, 0.2915, 0.3308,   # IG2
    0.4365, 0.4543, 0.4592,   # IG1
    0.4200, 0.4400, 0.4600    # Royal Guard
  )
)

your_model_combined <- bind_rows(
  lapply(names(obs_list), function(net) {
    if(net %in% your_model_df$Net){
      obs_df <- obs_list[[net]]
      pred_df <- your_model_df %>% filter(Net == net)
      
      data.frame(
        Net = net,
        Observed = obs_df$Observed,
        Predicted = pred_df$Predicted,
        Model = "Your Model"
      )
    }
  })
)

# -------------------------------------------------
# 3️⃣ Isaac model predictions (IG1 + IG2 only)
# -------------------------------------------------
isaac_list <- list(
  IG1 = data.frame(
    Time = c(2496, 2859, 3195),
    Predicted = c(0.36, 0.43, 0.45)
  ),
  IG2 = data.frame(
    Time = c(2496, 2859, 3195),
    Predicted = c(0.14, 0.12, 0.24)
  )
)

isaac_combined <- bind_rows(
  lapply(names(isaac_list), function(name) {
    
    net_name <- ifelse(name == "IG1",
                       "Interceptor",
                       "Interceptor G2")
    
    obs_df <- obs_list[[net_name]]
    pred_df <- isaac_list[[name]]
    
    data.frame(
      Net = net_name,
      Observed = obs_df$Observed,
      Predicted = pred_df$Predicted,
      Model = "Isaac Model"
    )
  })
)

# -------------------------------------------------
# 4️⃣ Combine models
# -------------------------------------------------
all_models_df <- bind_rows(your_model_combined,
                           isaac_combined)

# -------------------------------------------------
# 5️⃣ MSE per net
# -------------------------------------------------
mse_by_net <- all_models_df %>%
  group_by(Model, Net) %>%
  summarise(
    MSE = mean((Observed - Predicted)^2),
    .groups = "drop"
  )

# -------------------------------------------------
# 6️⃣ Model-specific overall MSE
# -------------------------------------------------

isaac_overall <- isaac_combined %>%
  summarise(
    Model = "Isaac Model",
    Net = "Overall",
    MSE = mean((Observed - Predicted)^2)
  )

your_overall <- your_model_combined %>%
  summarise(
    Model = "Your Model",
    Net = "Overall",
    MSE = mean((Observed - Predicted)^2)
  )

mse_table <- bind_rows(mse_by_net,
                       isaac_overall,
                       your_overall)

print(mse_table)

# -------------------------------------------------
# 7️⃣ Remove Royal Guard from plot only
# -------------------------------------------------
mse_plot_data <- mse_table %>%
  filter(Net != "Royal Guard")

# -------------------------------------------------
# 8️⃣ Plot
# -------------------------------------------------
ggplot(mse_plot_data, aes(x = Net, y = MSE, fill = Model)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6) +
  scale_fill_manual(values = c("your Model" = "skyblue",
                               "Isaac Model" = "salmon")) +
  labs(x = "Net type",
       y = "Mean Squared Error (MSE)",
       fill = "Model") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )


























library(dplyr)
library(ggplot2)

# -------------------------------------------------
# 1️⃣ Observed data
# -------------------------------------------------
obs_list <- list(
  "Interceptor G2" = data.frame(
    Time = c(2496, 2859, 3195),
    Observed = c(0.157, 0.279, 0.255)
  ),
  "Interceptor" = data.frame(
    Time = c(2496, 2859, 3195),
    Observed = c(0.280, 0.387, 0.262)
  ),
  "Royal Guard" = data.frame(
    Time = c(2496, 2859, 3195),
    Observed = c(0.300, 0.410, 0.390)  # replace if needed
  )
)

# -------------------------------------------------
# 2️⃣ My model predictions (IG1, IG2, RG)
# -------------------------------------------------
my_model_df <- data.frame(
  Net = c(rep("Interceptor G2",3),
          rep("Interceptor",3),
          rep("Royal Guard",3)),
  Time = rep(c(2496, 2859, 3195), 3),
  Predicted = c(
    0.2770, 0.2915, 0.3308,   # IG2
    0.4365, 0.4543, 0.4592,   # IG1
    0.4200, 0.4400, 0.4600    # Royal Guard
  )
)

my_model_combined <- bind_rows(
  lapply(names(obs_list), function(net) {
    if(net %in% my_model_df$Net){
      obs_df <- obs_list[[net]]
      pred_df <- my_model_df %>% filter(Net == net)
      
      data.frame(
        Net = net,
        Observed = obs_df$Observed,
        Predicted = pred_df$Predicted,
        Model = "My Model"
      )
    }
  })
)

# -------------------------------------------------
# 3️⃣ Isaac model predictions (IG1 + IG2 only)
# -------------------------------------------------
isaac_list <- list(
  IG1 = data.frame(
    Time = c(2496, 2859, 3195),
    Predicted = c(0.36, 0.43, 0.45)
  ),
  IG2 = data.frame(
    Time = c(2496, 2859, 3195),
    Predicted = c(0.14, 0.12, 0.24)
  )
)

isaac_combined <- bind_rows(
  lapply(names(isaac_list), function(name) {
    
    net_name <- ifelse(name == "IG1",
                       "Interceptor",
                       "Interceptor G2")
    
    obs_df <- obs_list[[net_name]]
    pred_df <- isaac_list[[name]]
    
    data.frame(
      Net = net_name,
      Observed = obs_df$Observed,
      Predicted = pred_df$Predicted,
      Model = "Isaac Model"
    )
  })
)

# -------------------------------------------------
# 4️⃣ Combine models
# -------------------------------------------------
all_models_df <- bind_rows(my_model_combined,
                           isaac_combined)

# -------------------------------------------------
# 5️⃣ MSE per net
# -------------------------------------------------
mse_by_net <- all_models_df %>%
  group_by(Model, Net) %>%
  summarise(
    MSE = mean((Observed - Predicted)^2),
    .groups = "drop"
  )

# -------------------------------------------------
# 6️⃣ Model-specific overall MSE
# -------------------------------------------------

isaac_overall <- isaac_combined %>%
  summarise(
    Model = "Isaac Model",
    Net = "Overall",
    MSE = mean((Observed - Predicted)^2)
  )

my_overall <- my_model_combined %>%
  summarise(
    Model = "My Model",
    Net = "Overall",
    MSE = mean((Observed - Predicted)^2)
  )

mse_table <- bind_rows(mse_by_net,
                       isaac_overall,
                       my_overall)

print(mse_table)

# -------------------------------------------------
# 7️⃣ Remove Royal Guard from plot only
# -------------------------------------------------
mse_plot_data <- mse_table %>%
  filter(Net != "Royal Guard")

# -------------------------------------------------
# 8️⃣ Plot
# -------------------------------------------------
ggplot(mse_plot_data, aes(x = Net, y = MSE, fill = Model)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6) +
  scale_fill_manual(values = c("My Model" = "skyblue",
                               "Isaac Model" = "salmon")) +
  labs(x = "Net type",
       y = "Mean Squared Error (MSE)",
       fill = "Model") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

