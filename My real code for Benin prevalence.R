

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
target_pfpr <- c(0.42, 0.407)  # desired PfPR at these days


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









########incidence

total_cases <- sum(output_bednet1$n_inc_clinical_182.5_3284)
total_child_years <- sum(output_bednet1$n_age_182.5_3284 / 365)

incidence <- total_cases / total_child_years
cat("Incidence per child-year:", round(incidence, 3), "\n")






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
#0.468/
target_pfpr <- c(0.48, 0.465)  # desired PfPR at these days

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
target_pfpr <- c(0.450, 0.435)  # desired PfPR at these days

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
    list(x = c(2190, 2496, 2859, 3195), y = c(0.407, 0.157, 0.279, 0.225), lower = c(0.382, 0.139, 0.257, 0.20), upper = c(0.433, 0.177, 0.302, 0.28), col = cols[1]),
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
         legend = c("Interceptor G2", "Interceptor", "Royal Guard", "Observed prevalence"),
         col = c(cols, "black"),   # line colors for models, black for observed points
         lty = c(1, 1, 1, NA),     # lines for models, no line for observed
         lwd = c(2, 2, 2, NA),     # line widths
         pch = c(NA, NA, NA, 15),  # square point only for observed data
         pt.cex = 1.2,
         box.lty = 0)
  
}

# Run the combined plot with observed points
plot_combined_prevalence_all()




#############lets do estimated vs observed

###################################################################





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
    y = c(0.157, 0.279, 0.225),
    lower = c(0.139, 0.257, 0.20),
    upper = c(0.177, 0.302, 0.250)
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







######################Toms prevalence




library(ggplot2)
library(dplyr)

# -------------------------------------------------
# Colours
# -------------------------------------------------
cols <- c("IG2" = "darkgreen",
          "IG1" = "darkblue",
          "RG"  = "darkred")

# -------------------------------------------------
# Observed data (your main dataset)
# -------------------------------------------------
obs_list <- list(
  IG2 = data.frame(
    x = c(2496, 2859, 3195),
    y = c(0.157, 0.279, 0.255),
    lower = c(0.139, 0.257, 0.235),
    upper = c(0.177, 0.302, 0.280)
  ),
  IG1 = data.frame(
    x = c(2496, 2859, 3195),
    y = c(0.280, 0.387, 0.262),
    lower = c(0.257, 0.362, 0.235),
    upper = c(0.303, 0.412, 0.280)
  ),
  RG = data.frame(
    x = c(2496, 2859, 3195),
    y = c(0.269, 0.382, 0.284),
    lower = c(0.247, 0.357, 0.266),
    upper = c(0.292, 0.407, 0.300)
  )
)

# -------------------------------------------------
# Helper function for predictions
# -------------------------------------------------
get_pred <- function(output, xvals) {
  approx(
    x = output$timestep,
    y = output$n_detect_lm_0_36499 / output$n_age_0_36499,
    xout = xvals
  )$y
}

# -------------------------------------------------
# Build main plotting dataframe
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
# Isaac study data (circles)
# -------------------------------------------------
isaac_list <- list(
  IG1 = data.frame(
    x = c(2496, 2859, 3195),
    y = c(0.36, 0.43, 0.45),
    lower = c(0.40, 0.433, 0.46),
    upper = c(0.42, 0.460, 0.489)
    
  ),
  IG2 = data.frame(
    x = c(2496, 2859, 3195),
    y = c(0.14, 0.12, 0.24),
    lower = c(0.15, 0.15, 0.235),
    upper = c(0.18, 0.177, 0.280)
  )
)




isaac_df <- bind_rows(
  lapply(names(isaac_list), function(net_name) {
    
    obs <- isaac_list[[net_name]]
    
    output_obj <- get(paste0("output_bednet",
                             which(names(obs_list) == net_name)))
    
    pred_vals <- get_pred(output_obj, obs$x)
    
    data.frame(
      Net = net_name,
      Observed = obs$y,
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
  
  # ■ Your main dataset (squares)
  geom_point(shape = 15, size = 2.5) +
  
  # Straight horizontal CI lines
  geom_segment(aes(x = ObsLower,
                   xend = ObsUpper,
                   y = Predicted,
                   yend = Predicted),
               linewidth = 1) +
  
  # ● Isaac study (circles)
  geom_point(data = isaac_df,
             aes(x = Observed, y = Predicted, color = Net),
             shape = 16,
             size = 2.5) +
  
  scale_color_manual(values = cols) +
  
  coord_equal(xlim = c(0, 0.7),
              ylim = c(0, 0.7)) +
  
  labs(x = "Observed prevalence",
       y = "Predicted prevalence",
       color = "Net type") +
  
  theme_bw(base_size = 15) +
  theme(
    legend.position = "top",
    panel.grid = element_blank()
  )

print(p)












#######################3



library(ggplot2)
library(dplyr)

# -------------------------------------------------
# Colours
# -------------------------------------------------
cols <- c("Interceptor G2" = "darkgreen",
          "Interceptor" = "darkblue",
          "Royal Guard"  = "darkred")

# -------------------------------------------------
# Original observed data (your data)
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
# Isaacs observed data (circle points)
# -------------------------------------------------
isaacs_list <- list(
  "Interceptor G2" = data.frame(
    x = c(2496, 2859, 3195),
    y = c(0.14, 0.12, 0.24)   # replace with Isaacs' actual values
  ),
  "Interceptor" = data.frame(
    x = c(2496, 2859, 3195),
    y = c(0.36, 0.43, 0.45)   # replace with Isaacs' actual values
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
# Build plotting dataframe for predictions
# -------------------------------------------------
plot_df <- bind_rows(
  lapply(seq_along(obs_list), function(i) {
    
    net_name <- names(obs_list)[i]
    obs <- obs_list[[i]]
    output_obj <- get(paste0("output_bednet", i))
    
    pred_vals <- get_pred(output_obj, obs$x)
    
    data.frame(
      Net = net_name,
      x = obs$x,
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
  
  # Your original observed points (squares)
  geom_point(aes(x = Observed, y = Predicted),
             shape = 15, size = 3) +
  
  # Horizontal confidence intervals
  geom_segment(aes(x = ObsLower,
                   xend = ObsUpper,
                   y = Predicted,
                   yend = Predicted),
               linewidth = 1) +
  
  # Overlay Isaacs points (circles) at same time points
  lapply(names(isaacs_list), function(net) {
    geom_point(
      data = isaacs_list[[net]] %>% mutate(Net = net),
      aes(x = y, y = get_pred(get(paste0("output_bednet", which(names(obs_list)==net))), x), color = Net),
      shape = 16,   # circle
      size = 3
    )
  }) +
  
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
################mean square error 






#########extract values at the time poits 


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





#######MSE


library(ggplot2)
library(dplyr)

# -------------------------------------------------
# Colours for nets
# -------------------------------------------------
cols <- c("Interceptor G2" = "darkgreen",
          "Interceptor" = "darkblue",
          "Royal Guard"  = "darkred")

# -------------------------------------------------
# Observed data
# -------------------------------------------------
obs_list <- list(
  "Interceptor G2" = data.frame(x = c(2496, 2859, 3195), y = c(0.157, 0.279, 0.255)),
  "Interceptor"    = data.frame(x = c(2496, 2859, 3195), y = c(0.280, 0.387, 0.262)),
  "Royal Guard"    = data.frame(x = c(2496, 2859, 3195), y = c(0.269, 0.382, 0.284))
)

# -------------------------------------------------
# Predicted data (your model)
# -------------------------------------------------
predicted_df <- data.frame(
  Net = c(rep("Interceptor G2",3), rep("Interceptor",3), rep("Royal Guard",3)),
  Time = rep(c(2496, 2859, 3195), 3),
  Predicted = c(0.2770, 0.2915, 0.3308,
                0.4365, 0.4543, 0.4592,
                0.3303, 0.3701, 0.4025)
)

# -------------------------------------------------
# Combine observed and predicted
# -------------------------------------------------
combined_df <- bind_rows(
  lapply(names(obs_list), function(net) {
    data.frame(
      Net = net,
      Time = obs_list[[net]]$x,
      Observed = obs_list[[net]]$y,
      Predicted = predicted_df$Predicted[predicted_df$Net == net]
    )
  })
)

# -------------------------------------------------
# Compute ±5% interval around predicted values
# -------------------------------------------------
combined_df <- combined_df %>%
  mutate(
    PredLower = Predicted - 0.05,
    PredUpper = Predicted + 0.05
  )

# -------------------------------------------------
# Calculate MSE and RMSE per net
# -------------------------------------------------
mse_rmse_by_net <- combined_df %>%
  group_by(Net) %>%
  summarise(
    n = n(),
    MSE = mean((Observed - Predicted)^2),
    RMSE = sqrt(MSE)
  )

print("MSE/RMSE per net:")
print(mse_rmse_by_net)

# Overall MSE/RMSE
mse_overall <- mean((combined_df$Observed - combined_df$Predicted)^2)
rmse_overall <- sqrt(mse_overall)
cat("Overall MSE:", round(mse_overall, 4), "RMSE:", round(rmse_overall, 4), "\n")

# -------------------------------------------------
# Plot Observed vs Predicted
# -------------------------------------------------
ggplot(combined_df, aes(x = Observed, y = Predicted, color = Net)) +
  
  # 1:1 reference line
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  
  # Observed points (squares)
  geom_point(shape = 15, size = 3) +
  
  # Horizontal ±5% interval for predicted values
  geom_segment(aes(x = PredLower, xend = PredUpper, y = Predicted, yend = Predicted),
               linetype = "dotted", linewidth = 1) +
  
  scale_color_manual(values = cols) +
  
  coord_equal(xlim = c(0, 0.7), ylim = c(0, 0.7)) +
  
  labs(x = "Observed Prevalence",
       y = "Predicted Prevalence",
       color = "Net type") +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid = element_blank()
  )





################share code isaacs mse


library(dplyr)

# -------------------------------------------------
# Observed data for Interceptor G2 and Interceptor
# -------------------------------------------------
obs_list <- list(
  "Interceptor G2" = data.frame(x = c(2496, 2859, 3195), y = c(0.157, 0.279, 0.255)),
  "Interceptor"    = data.frame(x = c(2496, 2859, 3195), y = c(0.280, 0.387, 0.262))
)

# -------------------------------------------------
# Isaac model predictions
# -------------------------------------------------
isaac_list <- list(
  IG1 = data.frame(x = c(2496, 2859, 3195), y = c(0.36, 0.43, 0.45)),  # Interceptor
  IG2 = data.frame(x = c(2496, 2859, 3195), y = c(0.14, 0.12, 0.24))   # Interceptor G2
)

# -------------------------------------------------
# Combine observed and predicted for Isaac model
# -------------------------------------------------
combined_isaac <- bind_rows(
  lapply(names(isaac_list), function(name) {
    # Map IG1 → Interceptor, IG2 → Interceptor G2
    net_name <- ifelse(name == "IG1", "Interceptor", "Interceptor G2")
    data.frame(
      Net = net_name,
      Time = isaac_list[[name]]$x,
      Observed = obs_list[[net_name]]$y,
      Predicted = isaac_list[[name]]$y
    )
  })
)

# -------------------------------------------------
# Compute ±5% interval for plotting
# -------------------------------------------------
combined_isaac <- combined_isaac %>%
  mutate(
    PredLower = Predicted - 0.05,
    PredUpper = Predicted + 0.05
  )

# -------------------------------------------------
# Compute MSE and RMSE per net
# -------------------------------------------------
mse_rmse_by_net <- combined_isaac %>%
  group_by(Net) %>%
  summarise(
    n = n(),
    MSE = mean((Observed - Predicted)^2),
    RMSE = sqrt(MSE)
  )

print("Isaac model MSE/RMSE per net:")
print(mse_rmse_by_net)

# Overall MSE/RMSE
mse_overall <- mean((combined_isaac$Observed - combined_isaac$Predicted)^2)
rmse_overall <- sqrt(mse_overall)
cat("Isaac model overall MSE:", round(mse_overall, 4), "RMSE:", round(rmse_overall, 4), "\n")

























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







####################################Tanzania Trials 
###########################################################################################################################################################################
############3let me try tz agains
#################################now let me do for children only 



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




#target_days <- c( 1, 1369)  # mid-year of years 1, 2, 3
#target_pfpr <- c( 0.50, 0.427)  # desired

target_days <- c(  1369)  # mid-year of years 1, 2, 3

#target_pfpr <- c(  0.50)  # desired
target_pfpr <- c(  0.42)  # desired



point_pfpr_summary <- function(x, days = target_days) {
  pfpr <- x$n_detect_lm_182.5_5109  / x$n_age_182.5_5109
  sapply(days, function(d) {
    pfpr[which.min(abs(x$timestep - d))]
  })
}


TZA <- fetch_site(iso3c = "TZA", admin_level = 1, urban_rural = TRUE)
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
  gamman = rep(c(2.25, 2.25, 2.65) *365)  # Vector of bed net half-lives for each distribution timestep
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
  
  obs_x <- c(1369, 1854, 2034, 2220, 2398)
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




#target_days <- c( 1, 1369)  # mid-year of years 1, 2, 3
#target_pfpr <- c( 0.50, 0.427)  # desired

target_days <- c(  1369)  # mid-year of years 1, 2, 3

#target_pfpr <- c(  0.48)  # desired
target_pfpr <- c(  0.475)  # desired



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
  gamman = rep(c(2.25, 2.25, 2.25) *365)  # Vector of bed net half-lives for each distribution timestep
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
  
  obs_x <- c(1369, 1854, 2034, 2220, 2398)
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




#target_days <- c( 1, 1369)  # mid-year of years 1, 2, 3
#target_pfpr <- c( 0.50, 0.427)  # desired

target_days <- c(  1369)  # mid-year of years 1, 2, 3

#target_pfpr <- c(  0.48)  # desired
target_pfpr <- c(  0.425)  # desired



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
  gamman = rep(c(2.25, 2.25, 2.50) *365)  # Vector of bed net half-lives for each distribution timestep
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
  
  obs_x <- c(1369, 1854, 2034, 2220, 2398)
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




#target_days <- c( 1, 1369)  # mid-year of years 1, 2, 3
#target_pfpr <- c( 0.50, 0.427)  # desired

target_days <- c(   1369)  # mid-year of years 1, 2, 3

#target_pfpr <- c(  0.48)  # desired
target_pfpr <- c(  0.464)  # desired



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
                 0.25, 0.25, 0.25),nrow = 3,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.72, 0.72,0.72,
                0.72, 0.72,0.72,
                0.58, 0.58,0.58 ),nrow = 3,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24,
                 0.24, 0.24, 0.24,
                 0.24, 0.24, 0.24),nrow = 3,ncol = 3,byrow = TRUE),
  gamman = rep(c(2.25, 2.25, 2.44) *365)  # Vector of bed net half-lives for each distribution timestep
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

output_control7 <- run_simulation(timesteps = parameters$timesteps, parameters = simparams)


output_bednet7 <- run_simulation(
  timesteps = parameters$timesteps,
  parameters = parameters
)



sim_length <- 8 * 365  # 10 years
cols <- c( "darkred", "#E69F00")  # colors for line

plot_combined_prevalence7 <- function() {
  
  bednetstimesteps <- c(1, 1489)
  
  slice_start <- 1096
  
  prev_mean <- output_bednet7$n_detect_lm_182.5_5109 / output_bednet7$n_age_182.5_5109
  prev_lower <- pmax(0, prev_mean * 0.90)
  prev_upper <- pmin(1, prev_mean * 1.05)
  
  # 🔧 FIX: limit to 2022
  keep <- output_bednet7$timestep >= slice_start & 
    output_bednet7$timestep <= 2922
  
  timesteps <- output_bednet7$timestep[keep]
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
  
  obs_x <- c(1369, 1854, 2034, 2220, 2398)
  obs_y <- c(0.462, 0.217, 0.506, 0.375, 0.424)
  
  obs_lower <- c(0.433, 0.194, 0.478, 0.349, 0.411)
  obs_upper <- c(0.49, 0.243, 0.535, 0.402, 0.435)
  
  # 🔧 FIX: keep points within 2022
  obs_keep <- obs_x >= slice_start & obs_x <= 2922
  
  if(any(obs_keep)){
    points(obs_x[obs_keep], obs_y[obs_keep], col = "darkred", pch = 15, cex = 1)
    
    arrows(
      x0 = obs_x[obs_keep],
      y0 = obs_lower[obs_keep],
      x1 = obs_x[obs_keep],
      y1 = obs_upper[obs_keep],
      code = 3, angle = 90, length = 0.05, col = "darkred", lwd = 1.5
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

plot_combined_prevalence7()












#################################################plot these three plese 
########333anotherone




plot_combined_prevalence_all2 <- function() {
  # -----------------------------
  slice_start <- 1096
  slice_end <- 2920  # end of 2022
  cols <- c("darkgreen", "darkblue","aquamarine3", "black", "darkred")
  
  # -----------------------------
  # Extract timesteps and prevalence for all four simulations
  timesteps1 <- output_bednet4$timestep[output_bednet4$timestep >= slice_start & output_bednet4$timestep <= slice_end]
  prev_mean1 <- (output_bednet4$n_detect_lm_182.5_5109 / output_bednet4$n_age_182.5_5109)[output_bednet4$timestep >= slice_start & output_bednet4$timestep <= slice_end]
  prev_lower1 <- pmax(0, prev_mean1 * 0.90)
  prev_upper1 <- pmin(1, prev_mean1 * 1.15)
  
  timesteps2 <- output_bednet5$timestep[output_bednet5$timestep >= slice_start & output_bednet5$timestep <= slice_end]
  prev_mean2 <- (output_bednet5$n_detect_lm_182.5_5109 / output_bednet5$n_age_182.5_5109)[output_bednet5$timestep >= slice_start & output_bednet5$timestep <= slice_end]
  prev_lower2 <- pmax(0, prev_mean2 * 0.89)
  prev_upper2 <- pmin(1, prev_mean2 * 1.03)
  
  timesteps3 <- output_bednet6$timestep[output_bednet6$timestep >= slice_start & output_bednet6$timestep <= slice_end]
  prev_mean3 <- (output_bednet6$n_detect_lm_182.5_5109 / output_bednet6$n_age_182.5_5109)[output_bednet6$timestep >= slice_start & output_bednet6$timestep <= slice_end]
  prev_lower3 <- pmax(0, prev_mean3 * 0.90)
  prev_upper3 <- pmin(1, prev_mean3 * 1.05)
  
  #timesteps4 <- output_bednet7$timestep[output_bednet7$timestep >= slice_start & output_bednet7$timestep <= slice_end]
  #prev_mean4 <- (output_bednet7$n_detect_lm_182.5_5109 / output_bednet7$n_age_182.5_5109)[output_bednet7$timestep >= slice_start & output_bednet7$timestep <= slice_end]
  #prev_lower4 <- pmax(0, prev_mean4 * 0.98)
  #prev_upper4 <- pmin(1, prev_mean4 * 1.12)
  
  # -----------------------------
  # Plot setup
  plot(timesteps1, prev_mean1, type = "n",
       xlab = "Time (Years)", ylab = " Malaria Prevalence for children aged 0.5-14 years",
       xaxt = "n", ylim = c(0, 0.7),
       xaxs = "i", yaxs = "i",
       xlim = c(min(timesteps1), max(timesteps1)),
       main = "Tanzania Trials",
       cex.main = 1, cex.lab = 1)
  
  # -----------------------------
  # Plot uncertainty polygons
  polygon(c(timesteps1, rev(timesteps1)), c(prev_lower1, rev(prev_upper1)),
          col = adjustcolor(cols[1], alpha.f = 0.15), border = NA)
  polygon(c(timesteps2, rev(timesteps2)), c(prev_lower2, rev(prev_upper2)),
          col = adjustcolor(cols[2], alpha.f = 0.15), border = NA)
  polygon(c(timesteps3, rev(timesteps3)), c(prev_lower3, rev(prev_upper3)),
          col = adjustcolor(cols[3], alpha.f = 0.15), border = NA)
  #polygon(c(timesteps4, rev(timesteps4)), c(prev_lower4, rev(prev_upper4)),
       #   col = adjustcolor(cols[4], alpha.f = 0.15), border = NA)
  
  # -----------------------------
  # Plot mean prevalence lines
  lines(timesteps1, prev_mean1, col = cols[1], lwd = 1)
  lines(timesteps2, prev_mean2, col = cols[2], lwd = 1)
  lines(timesteps3, prev_mean3, col = cols[3], lwd = 1)
  #lines(timesteps4, prev_mean4, col = cols[4], lwd = 1)
  
  # -----------------------------
  # Add observed points with 95% CI
  obs_list <- list(
    list(x = c(1400, 1854, 2034, 2220, 2398), 
         y =c(0.420, 0.156, 0.4092, 0.256, 0.492), 
         lower = c(0.398, 0.136, 0.382, 0.223, 0.475), 
         upper = c(0.456, 0.178, 0.437, 0.281, 0.51), col = cols[1]),
    list(x = c(1400, 1854, 2034, 2220, 2398), 
         y = c(0.459, 0.312, 0.523, 0.458, 0.530), 
         lower = c(0.450, 0.285, 0.495, 0.430, 0.520), 
         upper = c(0.488, 0.340, 0.551, 0.486, 0.547), col = cols[2]),
    list(x = c(1400, 1854, 2034, 2220, 2398), 
         y = c(0.423, 0.19, 0.433, 0.407, 0.492), 
         lower = c(0.391, 0.170, 0.405, 0.380, 0.485), 
         upper = c(0.450, 0.214, 0.461, 0.434, 0.518), col = cols[3])
    #list(x = c(1400, 1854, 2034, 2220, 2398), 
     #    y = c(0.462, 0.217, 0.506, 0.375, 0.424), 
      #   lower = c(0.433, 0.194, 0.478, 0.349, 0.411), 
       #  upper = c(0.49, 0.243, 0.535, 0.402, 0.435), col = cols[4])
  )
  
  for(obs in obs_list){
    keep <- obs$x >= slice_start & obs$x <= slice_end
    if(any(keep)){
      points(obs$x[keep], obs$y[keep], col = obs$col, pch = 19, cex = 1)
      arrows(x0 = obs$x[keep], y0 = obs$lower[keep],
             x1 = obs$x[keep], y1 = obs$upper[keep],
             code = 3, angle = 90, length = 0.05, col = obs$col, lwd = 1)
    }
  }
  
  # -----------------------------
  # Custom x-axis: yearly ticks
  years <- 2015:2022
  year_days <- seq(1, sim_length, by = 365)
  year_days_visible <- year_days[year_days >= slice_start & year_days <= slice_end]
  years_visible <- years[year_days >= slice_start & year_days <= slice_end]
  axis(1, at = year_days_visible, labels = years_visible)
  
  # -----------------------------
  # Bed net distribution lines
  bednetstimesteps <- c(1489)
  bednet_visible <- bednetstimesteps[bednetstimesteps >= slice_start & bednetstimesteps <= slice_end]
  if(length(bednet_visible) > 0){
    abline(v = bednet_visible, col = "darkgray", lty = "dashed", lwd = 1)
    text(bednet_visible + 30, 0.65, "Bed Net Distribution.", pos = 4, cex = 0.8, srt = 0)
    arrows(
      x0 = bednet_visible + 25,
      y0 = 0.65,
      x1 = bednet_visible,
      y1 = 0.65,
      length = 0.1,
      angle = 20,
      col = "black",
      lwd = 1.2
    )
  }
  
  # -----------------------------
  # Legend
  legend(x = 2200, y = 0.20,
         legend = c("Pyrethroid-pyrrole", "pyrethroid only",  "pyrethroid-PBO",
                     "Observed prevalence"),
         col = c(cols, "black"),
         lty = c(1, 1, 1, NA),
         lwd = c(2, 2, 2, NA),
         pch = c(NA, NA, NA, 19),
         pt.cex = 1.2,
         cex = 0.8,
         box.lty = 0)
}

# Run the combined plot
plot_combined_prevalence_all2()



#########################predicted vs estimated 


library(ggplot2)
library(dplyr)

# -------------------------------------------------
# Net names + colours
# -------------------------------------------------
net_names <- c("Pyrethroid-chlorfenapyr", "pyrethroid only",  "pyrethroid-PBO")

cols <- c("darkgreen", "darkblue", "aquamarine3")
names(cols) <- net_names

# -------------------------------------------------
# Observed data
# -------------------------------------------------
obs_list <- list(
  list(x = c( 1854, 2034, 2220, 2398), 
       y = c( 0.156, 0.4092, 0.256, 0.492), 
       lower = c( 0.136, 0.382, 0.223, 0.470), 
       upper = c( 0.178, 0.437, 0.281, 0.53)),
  
  list(x = c( 1854, 2034, 2220, 2398), 
       y = c(0.312, 0.523, 0.458, 0.530), 
       lower = c( 0.285, 0.495, 0.430, 0.520), 
       upper = c( 0.340, 0.551, 0.486, 0.547)),
  
  list(x = c( 1854, 2034, 2220, 2398), 
       y = c( 0.19, 0.433, 0.407, 0.492), 
       lower = c(0.170, 0.405, 0.380, 0.470), 
       upper = c( 0.214, 0.461, 0.434, 0.525))
)

# -------------------------------------------------
# Helper function (UPDATED prevalence)
# -------------------------------------------------
get_pred <- function(output, xvals) {
  approx(
    x = output$timestep,
    y = output$n_detect_lm_182.5_5109 / output$n_age_182.5_5109,
    xout = xvals
  )$y
}

# -------------------------------------------------
# Build plotting dataframe
# -------------------------------------------------
plot_df <- bind_rows(
  lapply(seq_along(obs_list), function(i) {
    
    obs <- obs_list[[i]]
    
    # map to output_bednet4, 5, 6
    output_obj <- get(paste0("output_bednet", i + 3))
    
    pred_vals <- get_pred(output_obj, obs$x)
    
    data.frame(
      Net = net_names[i],
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
p1 <- ggplot(plot_df, aes(x = Observed, y = Predicted, color = Net)) +
  
  # 1:1 reference line
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed",
              linewidth = 1,
              colour = "black") +
  
  # Points
  geom_point(shape = 19, size = 2) +
  
  # Horizontal confidence intervals
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

print(p1)



















#benin + Z



library(ggplot2)
library(dplyr)

# -----------------------------
# Colour mapping (consistent across countries)
# -----------------------------
cols <- c(
  "Pyrethroid-chlorfenapyr" = "darkgreen",
  "pyrethroid only" = "darkblue",
  "pyrethroid-PBO" = "aquamarine3",
  "pyrethroid-Pyriproxyfen" = "darkred"
)

# -----------------------------
# Tanzania data
# -----------------------------
tanzania_names <- c("Pyrethroid-chlorfenapyr", "pyrethroid only", "pyrethroid-PBO")
tanzania_obs <- list(
  list(x = c(1854, 2034, 2220, 2398),
       y = c(0.156, 0.4092, 0.256, 0.492),
       lower = c(0.136, 0.382, 0.223, 0.470),
       upper = c(0.178, 0.437, 0.281, 0.53)),
  list(x = c(1854, 2034, 2220, 2398),
       y = c(0.312, 0.523, 0.458, 0.530),
       lower = c(0.285, 0.495, 0.430, 0.520),
       upper = c(0.340, 0.551, 0.486, 0.547)),
  list(x = c(1854, 2034, 2220, 2398),
       y = c(0.19, 0.433, 0.407, 0.492),
       lower = c(0.170, 0.405, 0.380, 0.470),
       upper = c(0.214, 0.461, 0.434, 0.525))
)

tanzania_df <- bind_rows(
  lapply(seq_along(tanzania_obs), function(i) {
    obs <- tanzania_obs[[i]]
    output_obj <- get(paste0("output_bednet", i + 3)) # output_bednet4,5,6
    pred_vals <- approx(
      x = output_obj$timestep,
      y = output_obj$n_detect_lm_182.5_5109 / output_obj$n_age_182.5_5109,
      xout = obs$x
    )$y
    data.frame(
      Country = "Tanzania",
      Net = tanzania_names[i],
      Observed = obs$y,
      ObsLower = obs$lower,
      ObsUpper = obs$upper,
      Predicted = pred_vals
    )
  })
)

# -----------------------------
# Benin data
# -----------------------------
benin_names <- c("Pyrethroid-chlorfenapyr", "pyrethroid only", "pyrethroid-Pyriproxyfen")
benin_obs <- list(
  list(x = c(2496, 2859, 3195),
       y = c(0.157, 0.279, 0.225),
       lower = c(0.139, 0.257, 0.20),
       upper = c(0.177, 0.302, 0.250)),
  list(x = c(2496, 2859, 3195),
       y = c(0.280, 0.387, 0.262),
       lower = c(0.257, 0.362, 0.235),
       upper = c(0.303, 0.412, 0.280))
  #list(x = c(2496, 2859, 3195),
   #    y = c(0.269, 0.382, 0.284),
    #   lower = c(0.247, 0.357, 0.266),
     #  upper = c(0.292, 0.407, 0.300))
)

benin_df <- bind_rows(
  lapply(seq_along(benin_obs), function(i) {
    obs <- benin_obs[[i]]
    output_obj <- get(paste0("output_bednet", i)) # output_bednet1,2,3
    pred_vals <- approx(
      x = output_obj$timestep,
      y = output_obj$n_detect_lm_0_36499 / output_obj$n_age_0_36499,
      xout = obs$x
    )$y
    data.frame(
      Country = "Benin",
      Net = benin_names[i],
      Observed = obs$y,
      ObsLower = obs$lower,
      ObsUpper = obs$upper,
      Predicted = pred_vals
    )
  })
)

# -----------------------------
# Combine both datasets
# -----------------------------
combined_df <- bind_rows(tanzania_df, benin_df)

# -----------------------------
# Plot
# -----------------------------
ggplot(combined_df, aes(x = Observed, y = Predicted, color = Net)) +
  
  # 1:1 reference line
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 1, colour = "black") +
  
  # Points: shape by country
  geom_point(aes(shape = Country), size = 3) +
  
  # Horizontal confidence intervals
  geom_segment(aes(x = ObsLower, xend = ObsUpper, y = Predicted, yend = Predicted), linewidth = 1) +
  
  # Colours and shapes
  scale_color_manual(values = cols) +
  scale_shape_manual(values = c("Benin" = 15, "Tanzania" = 16)) +
  
  # Axes
  coord_equal(xlim = c(0, 0.7), ylim = c(0, 0.7)) +
  labs(x = "Observed Malaria prevalence",
       y = "Predicted Malaria prevalence",
       color = "Net type",
       shape = "Country") +
  
  theme_bw(base_size = 14) +
  theme(legend.position = "top", panel.grid = element_blank())














ggplot(combined_df, aes(x = Observed, y = Predicted, color = Net)) +
  
  # 1:1 reference line
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 1, colour = "black") +
  
  # Points: shape by country
  geom_point(aes(shape = Country), size = 3) +
  
  # Optional: horizontal CI lines for actual points
  geom_segment(aes(x = ObsLower, xend = ObsUpper, y = Predicted, yend = Predicted),
               linewidth = 1, alpha = 0.5) +
  
  # Colours and shapes
  scale_color_manual(values = cols) +
  scale_shape_manual(values = c("Benin" = 15, "Tanzania" = 16)) +
  
  # Axes limits
  coord_equal(xlim = c(0, 0.7), ylim = c(0, 0.7)) +
  
  # Labels
  labs(x = "Observed Malaria prevalence",
       y = "Predicted Malaria prevalence",
       color = "Net type",
       shape = "Country") +
  
  # Theme + legend inside plot
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.75, 0.25), # inside plot: x,y as fraction
    legend.background = element_rect(fill = "white", color = "black")
  )







#####finalcombined plot 


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
  
  # Axes limits
  coord_equal(xlim = c(0, 0.7), ylim = c(0, 0.7)) +
  
  # Labels
  labs(x = "Observed Malaria prevalence",
       y = "Predicted Malaria prevalence",
       color = "Net type",
       shape = "Country") +
  
  # Theme + legend inside plot
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.25, 0.7),
    legend.background = element_rect(fill = "white", color = "black")
  )















################################MSE##############################################################################################
###############################################################################################################################################3
################mean square error 
#########extract values at the time poits 




#########REAL MSE 




library(dplyr)
library(ggplot2)

# -----------------------------
# 1. Net names + colors
# -----------------------------
tanzania_names <- c("Pyrethroid-chlorfenapyr", "pyrethroid only", "pyrethroid-PBO", "pyrethroid-Pyriproxyfen")
benin_names <- c("Pyrethroid-chlorfenapyr", "pyrethroid only", "pyrethroid-Pyriproxyfen")

cols <- c(
  "Pyrethroid-chlorfenapyr" = "darkgreen",
  "pyrethroid only" = "darkblue",
  "pyrethroid-PBO" = "aquamarine3",
  "bednet7" = "darkred",
  "pyrethroid-Pyriproxyfen" = "orange"
)

# -----------------------------
# 2. Observed data
# -----------------------------
tanzania_obs <- list(
  list(x = c( 1854, 2034, 2220, 2398), 
       #y = c( 0.217, 0.506, 0.375, 0.424)
       y = c( 0.312, 0.523, 0.458, 0.53),
       lower = c( 0.194, 0.478, 0.349, 0.411),
       upper = c( 0.243, 0.535, 0.402, 0.435)),
  list(x = c( 1854, 2034, 2220, 2398), 
       y = c( 0.156, 0.4092, 0.256, 0.492),
       lower = c( 0.136, 0.382, 0.223, 0.475),
       upper = c( 0.178, 0.437, 0.281, 0.51)),
  list(x = c( 1854, 2034, 2220, 2398),
       y = c( 0.19, 0.433, 0.407, 0.492),
       lower = c( 0.170, 0.405, 0.380, 0.485),
       upper = c( 0.214, 0.461, 0.434, 0.518)),
  list(x = c( 1854, 2034, 2220, 2398),
       y = c( 0.217, 0.506, 0.375, 0.424),
       lower = c( 0.194, 0.478, 0.349, 0.411),
       upper = c( 0.243, 0.535, 0.402, 0.435))
)

benin_obs <- list(
  list(x = c(2496, 2859, 3195),
       y = c(0.157, 0.279, 0.225),
       lower = c(0.139, 0.257, 0.20),
       upper = c(0.177, 0.302, 0.250)),
  list(x = c(2496, 2859, 3195),
       y = c(0.280, 0.387, 0.262),
       lower = c(0.257, 0.362, 0.235),
       upper = c(0.303, 0.412, 0.280)),
  list(x = c(2496, 2859, 3195),
       y = c(0.269, 0.382, 0.284),
       lower = c(0.247, 0.357, 0.266),
       upper = c(0.292, 0.407, 0.300))
)

# -----------------------------
# 3. Helper: get predicted values
# -----------------------------
get_pred <- function(output, xvals){
  approx(
    x = output$timestep,
    y = output$n_detect_lm_182.5_5109 / output$n_age_182.5_5109,
    xout = xvals
  )$y
}

# -----------------------------
# 4. Build Tanzania dataframe
# -----------------------------
tanzania_df <- bind_rows(
  lapply(seq_along(tanzania_obs), function(i){
    obs <- tanzania_obs[[i]]
    output_obj <- get(paste0("output_bednet", i+3))  # output_bednet4-7
    pred <- get_pred(output_obj, obs$x)
    data.frame(
      Country = "Tanzania",
      Net = tanzania_names[i],
      Observed = obs$y,
      ObsLower = obs$lower,
      ObsUpper = obs$upper,
      Predicted = pred
    )
  })
)

# -----------------------------
# 5. Build Benin dataframe
# -----------------------------
benin_df <- bind_rows(
  lapply(seq_along(benin_obs), function(i){
    obs <- benin_obs[[i]]
    output_obj <- get(paste0("output_bednet", i))  # output_bednet1-3
    pred <- approx(
      x = output_obj$timestep,
      y = output_obj$n_detect_lm_0_36499 / output_obj$n_age_0_36499,
      xout = obs$x
    )$y
    data.frame(
      Country = "Benin",
      Net = benin_names[i],
      Observed = obs$y,
      ObsLower = obs$lower,
      ObsUpper = obs$upper,
      Predicted = pred
    )
  })
)

# -----------------------------
# 6. Combine both
# -----------------------------
combined_df <- bind_rows(tanzania_df, benin_df)

# -----------------------------
# 7. Compute residuals, MSE, RMSE
# -----------------------------
combined_df <- combined_df %>%
  mutate(Residual = Observed - Predicted)

mse_country <- combined_df %>%
  group_by(Country) %>%
  summarise(
    MSE = mean(Residual^2),
    RMSE = sqrt(mean(Residual^2))
  )

print(mse_country)

# -----------------------------
# 8. Plot residuals per net
# -----------------------------
ggplot(combined_df, aes(x = Net, y = Residual, fill = Net)) +
  geom_boxplot(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~Country) +
  scale_fill_manual(values = cols) +
  labs(y = "Residual (Observed - Predicted)",
       title = "Residuals per Net by Country") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

# -----------------------------
# 9. Optional: Scatter residuals vs Observed
# -----------------------------
ggplot(combined_df, aes(x = Observed, y = Residual, color = Net)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~Country) +
  scale_color_manual(values = cols) +
  labs(y = "Residual", x = "Observed Prevalence",
       title = "Residuals vs Observed Prevalence") +
  theme_bw(base_size = 12)





library(ggplot2)
library(dplyr)
library(tidyr)

# -----------------------------
# Compute residuals
# -----------------------------
residuals_df <- combined_df %>%
  mutate(Residual = Predicted - Observed)

# RMSE per net
rmse_net <- residuals_df %>%
  group_by(Country, Net) %>%
  summarise(RMSE = sqrt(mean(Residual^2, na.rm = TRUE)), .groups = "drop") %>%
  mutate(Type = "Per Net")

# Overall RMSE per



library(ggplot2)
library(dplyr)

# -----------------------------
# Compute residuals
# -----------------------------
residuals_df <- combined_df %>%
  mutate(Residual = Predicted - Observed)

# RMSE per net
rmse_net <- residuals_df %>%
  group_by(Country, Net) %>%
  summarise(RMSE = sqrt(mean(Residual^2, na.rm = TRUE)), .groups = "drop") %>%
  mutate(Type = "Per Net")

# Overall RMSE per country
rmse_country <- residuals_df %>%
  group_by(Country) %>%
  summarise(RMSE = sqrt(mean(Residual^2, na.rm = TRUE)), .groups = "drop") %>%
  mutate(Net = "Overall", Type = "Overall")

# Combine for plotting
rmse_plot_df <- bind_rows(rmse_net, rmse_country)

# -----------------------------
# Publication-ready RMSE plot with full border
# -----------------------------
ggplot(rmse_plot_df, aes(x = Net, y = RMSE, fill = Country)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "black") + # bars with borders
  geom_text(aes(label = round(RMSE, 3)), 
            position = position_dodge(width = 0.7), 
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Tanzania" = "darkgreen", "Benin" = "darkblue")) +
  labs(
    title = "RMSE",
    x = "Net Type",
    y = "RMSE",
    fill = "Country"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, size = 1.2) # full border
  ) +
  ylim(0, 0.25) # set y-axis limit

##################################################END OF MSE


##################COMPARE MY OWN PREDICTIONS WITH ISAACS PREDICTED TO FIND WHICH MODEL HAS BETTER MSE
library(ggplot2)
library(dplyr)

# -------------------------------------------------
# Colours for nets
# -------------------------------------------------
cols <- c("Interceptor G2" = "darkgreen",
          "Interceptor" = "darkblue",
          "Royal Guard"  = "darkred")

# -------------------------------------------------
# Observed data
# -------------------------------------------------
obs_list <- list(
  "Interceptor G2" = data.frame(x = c(2496, 2859, 3195), y = c(0.157, 0.279, 0.255)),
  "Interceptor"    = data.frame(x = c(2496, 2859, 3195), y = c(0.280, 0.387, 0.262)),
  "Royal Guard"    = data.frame(x = c(2496, 2859, 3195), y = c(0.269, 0.382, 0.284))
)

# -------------------------------------------------
# Predicted data (your model)
# -------------------------------------------------
predicted_df <- data.frame(
  Net = c(rep("Interceptor G2",3), rep("Interceptor",3), rep("Royal Guard",3)),
  Time = rep(c(2496, 2859, 3195), 3),
  Predicted = c(0.2770, 0.2915, 0.3308,
                0.4365, 0.4543, 0.4592,
                0.3303, 0.3701, 0.4025)
)

# -------------------------------------------------
# Combine observed and predicted
# -------------------------------------------------
combined_df <- bind_rows(
  lapply(names(obs_list), function(net) {
    data.frame(
      Net = net,
      Time = obs_list[[net]]$x,
      Observed = obs_list[[net]]$y,
      Predicted = predicted_df$Predicted[predicted_df$Net == net]
    )
  })
)

# -------------------------------------------------
# Compute ±5% interval around predicted values
# -------------------------------------------------
combined_df <- combined_df %>%
  mutate(
    PredLower = Predicted - 0.05,
    PredUpper = Predicted + 0.05
  )

# -------------------------------------------------
# Calculate MSE and RMSE per net
# -------------------------------------------------
mse_rmse_by_net <- combined_df %>%
  group_by(Net) %>%
  summarise(
    n = n(),
    MSE = mean((Observed - Predicted)^2),
    RMSE = sqrt(MSE)
  )

print("MSE/RMSE per net:")
print(mse_rmse_by_net)

# Overall MSE/RMSE
mse_overall <- mean((combined_df$Observed - combined_df$Predicted)^2)
rmse_overall <- sqrt(mse_overall)
cat("Overall MSE:", round(mse_overall, 4), "RMSE:", round(rmse_overall, 4), "\n")

# -------------------------------------------------
# Plot Observed vs Predicted
# -------------------------------------------------
ggplot(combined_df, aes(x = Observed, y = Predicted, color = Net)) +
  
  # 1:1 reference line
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  
  # Observed points (squares)
  geom_point(shape = 15, size = 3) +
  
  # Horizontal ±5% interval for predicted values
  geom_segment(aes(x = PredLower, xend = PredUpper, y = Predicted, yend = Predicted),
               linetype = "dotted", linewidth = 1) +
  
  scale_color_manual(values = cols) +
  
  coord_equal(xlim = c(0, 0.7), ylim = c(0, 0.7)) +
  
  labs(x = "Observed Prevalence",
       y = "Predicted Prevalence",
       color = "Net type") +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid = element_blank()
  )





################share code isaacs mse


library(dplyr)

# -------------------------------------------------
# Observed data for Interceptor G2 and Interceptor
# -------------------------------------------------
obs_list <- list(
  "Interceptor G2" = data.frame(x = c(2496, 2859, 3195), y = c(0.157, 0.279, 0.255)),
  "Interceptor"    = data.frame(x = c(2496, 2859, 3195), y = c(0.280, 0.387, 0.262))
)

# -------------------------------------------------
# Isaac model predictions
# -------------------------------------------------
isaac_list <- list(
  IG1 = data.frame(x = c(2496, 2859, 3195), y = c(0.36, 0.43, 0.45)),  # Interceptor
  IG2 = data.frame(x = c(2496, 2859, 3195), y = c(0.14, 0.12, 0.24))   # Interceptor G2
)

# -------------------------------------------------
# Combine observed and predicted for Isaac model
# -------------------------------------------------
combined_isaac <- bind_rows(
  lapply(names(isaac_list), function(name) {
    # Map IG1 → Interceptor, IG2 → Interceptor G2
    net_name <- ifelse(name == "IG1", "Interceptor", "Interceptor G2")
    data.frame(
      Net = net_name,
      Time = isaac_list[[name]]$x,
      Observed = obs_list[[net_name]]$y,
      Predicted = isaac_list[[name]]$y
    )
  })
)

# -------------------------------------------------
# Compute ±5% interval for plotting
# -------------------------------------------------
combined_isaac <- combined_isaac %>%
  mutate(
    PredLower = Predicted - 0.05,
    PredUpper = Predicted + 0.05
  )

# -------------------------------------------------
# Compute MSE and RMSE per net
# -------------------------------------------------
mse_rmse_by_net <- combined_isaac %>%
  group_by(Net) %>%
  summarise(
    n = n(),
    MSE = mean((Observed - Predicted)^2),
    RMSE = sqrt(MSE)
  )

print("Isaac model MSE/RMSE per net:")
print(mse_rmse_by_net)

# Overall MSE/RMSE
mse_overall <- mean((combined_isaac$Observed - combined_isaac$Predicted)^2)
rmse_overall <- sqrt(mse_overall)
cat("Isaac model overall MSE:", round(mse_overall, 4), "RMSE:", round(rmse_overall, 4), "\n")

























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









#################to be continued 



library(dplyr)

# -----------------------------
# 1. Observed data
# -----------------------------
tanzania_obs <- list(
  list(x = c(1854, 2034, 2220, 2398), y = c(0.40, 0.43, 0.50, 0.47), lower = c(0.35, 0.38, 0.48, 0.45), upper = c(0.45, 0.48, 0.52, 0.50)),
  list(x = c(1854, 2034, 2220, 2398), y = c(0.26, 0.28, 0.34, 0.36), lower = c(0.22, 0.24, 0.31, 0.33), upper = c(0.30, 0.32, 0.37, 0.39)),
  list(x = c(1854, 2034, 2220, 2398), y = c(0.29, 0.32, 0.40, 0.42), lower = c(0.25, 0.28, 0.36, 0.39), upper = c(0.33, 0.36, 0.44, 0.45))
)

benin_obs <- list(
  list(x = c(2496, 2859, 3195), y = c(0.15, 0.13, 0.25), lower = c(0.12, 0.10, 0.20), upper = c(0.18, 0.16, 0.30)),
  list(x = c(2496, 2859, 3195), y = c(0.37, 0.42, 0.46), lower = c(0.33, 0.38, 0.42), upper = c(0.41, 0.46, 0.50)),
  list(x = c(2496, 2859, 3195), y = c(0.31, 0.34, 0.42), lower = c(0.28, 0.30, 0.38), upper = c(0.34, 0.38, 0.46))
)

# -----------------------------
# 2. Isaac predictions
# -----------------------------
isaac_tanzania <- list(
  list(x = c(1854, 2034, 2220, 2398), y = c(0.42, 0.433, 0.51, 0.48)),
  list(x = c(1854, 2034, 2220, 2398), y = c(0.25, 0.27, 0.35, 0.37)),
  list(x = c(1854, 2034, 2220, 2398), y = c(0.30, 0.33, 0.41, 0.41))
)

isaac_benin <- list(
  list(x = c(2496, 2859, 3195), y = c(0.14, 0.12, 0.24)),
  list(x = c(2496, 2859, 3195), y = c(0.36, 0.43, 0.45)),
  list(x = c(2496, 2859, 3195), y = c(0.30, 0.33, 0.41))
)

# -----------------------------
# 3. Helper: get predicted values from your model output
# -----------------------------
get_pred <- function(output, xvals){
  approx(
    x = output$timestep,
    y = output$n_detect_lm_182.5_5109 / output$n_age_182.5_5109,
    xout = xvals
  )$y
}

# -----------------------------
# 4. Build Tanzania dataframe
# -----------------------------
tanzania_df <- bind_rows(
  lapply(seq_along(tanzania_obs), function(i){
    obs <- tanzania_obs[[i]]
    output_obj <- get(paste0("output_bednet", i+3))  # output_bednet4-7
    pred <- get_pred(output_obj, obs$x)
    data.frame(
      Country = "Tanzania",
      Net = paste0("Net", i),
      Observed = obs$y,
      Predicted = pred,
      Isaac = isaac_tanzania[[i]]$y
    )
  })
)

# -----------------------------
# 5. Build Benin dataframe
# -----------------------------
benin_df <- bind_rows(
  lapply(seq_along(benin_obs), function(i){
    obs <- benin_obs[[i]]
    output_obj <- get(paste0("output_bednet", i))  # output_bednet1-3
    pred <- approx(
      x = output_obj$timestep,
      y = output_obj$n_detect_lm_0_36499 / output_obj$n_age_0_36499,
      xout = obs$x
    )$y
    data.frame(
      Country = "Benin",
      Net = paste0("Net", i),
      Observed = obs$y,
      Predicted = pred,
      Isaac = isaac_benin[[i]]$y
    )
  })
)

# -----------------------------
# 6. Combine and compute RMSE
# -----------------------------
all_df <- bind_rows(tanzania_df, benin_df)

rmse <- function(obs, pred) sqrt(mean((obs - pred)^2))

rmse_summary <- all_df %>%
  summarise(
    RMSE_YourModel = rmse(Observed, Predicted),
    RMSE_Isaac = rmse(Observed, Isaac)
  )

print(all_df)
print(rmse_summary)










###########Tanzania












#############cases averted 

library(dplyr)

# -----------------------------
# Define day ranges for Year 1 and Year 2
# -----------------------------
year1_start <- 1489
year1_end   <- 1489 + 365 - 1   # Year 1: day 1489 to 1853
year2_start <- 1489 + 365       # Year 2: day 1854 to 2218
year2_end   <- year2_start + 365 - 1

# Column with incidence
colname <- "n_inc_clinical_182.5_5109"

# Bednets and names
bednets <- 4:7
net_names <- c("Pyrethroid-Chlorfenapyr", "Pyrethroid-Only", "Pyrethroid-PBO", "pyrethroid-pyriproxyfen")

# -----------------------------
# Function to calculate cases averted over a day range
# -----------------------------
cases_averted_days <- function(control_df, bednet_df, start_day, end_day, colname) {
  control_sum <- control_df %>% filter(timestep >= start_day & timestep <= end_day) %>% summarise(total = sum(.data[[colname]], na.rm = TRUE)) %>% pull(total)
  bednet_sum  <- bednet_df  %>% filter(timestep >= start_day & timestep <= end_day) %>% summarise(total = sum(.data[[colname]], na.rm = TRUE)) %>% pull(total)
  
  cases_averted <- control_sum - bednet_sum
  pct_averted <- 100 * cases_averted / control_sum
  
  return(c(cases_averted = cases_averted, pct_averted = pct_averted))
}

# -----------------------------
# Loop over bednets
# -----------------------------
cases_summary <- lapply(seq_along(bednets), function(i) {
  bed <- bednets[i]
  control_df <- get(paste0("output_control", bed))
  bednet_df  <- get(paste0("output_bednet", bed))
  
  y1 <- cases_averted_days(control_df, bednet_df, year1_start, year1_end, colname)
  y2 <- cases_averted_days(control_df, bednet_df, year2_start, year2_end, colname)
  overall <- cases_averted_days(control_df, bednet_df, year1_start, year2_end, colname)
  
  data.frame(
    Net = net_names[i],
    Year1_Cases = y1["cases_averted"],
    Year1_Pct   = y1["pct_averted"],
    Year2_Cases = y2["cases_averted"],
    Year2_Pct   = y2["pct_averted"],
    Overall_Cases = overall["cases_averted"],
    Overall_Pct   = overall["pct_averted"]
  )
}) %>% bind_rows()

print(cases_summary)



library(ggplot2)
library(tidyr)
library(dplyr)

# -----------------------------
# Pivot data for plotting
# -----------------------------
plot_df <- cases_summary %>%
  select(Net, Year1_Pct, Year2_Pct, Overall_Pct) %>%
  pivot_longer(cols = -Net, names_to = "Period", values_to = "Pct_Averted") %>%
  mutate(Period = recode(Period, 
                         "Year1_Pct" = "Year 1",
                         "Year2_Pct" = "Year 2",
                         "Overall_Pct" = "Overall"))

# -----------------------------
# Plot % Cases Averted
# -----------------------------
ggplot(plot_df, aes(x = Net, y = Pct_Averted, fill = Period)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("Year 1" = "skyblue", "Year 2" = "salmon", "Overall" = "seagreen")) +
  labs(
    x = "Bednet Type",
    y = "Cases Averted (%)",
    fill = "Period",
    title = "Tanzania"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )




#################Benin


library(dplyr)
library(tidyr)
library(ggplot2)

# -----------------------------
# Define day ranges for Year 1 and Year 2 (Benin)
# -----------------------------
year1_start <- 2310
year1_end   <- 2310 + 365 - 1     # Year 1: day 1977 to 2341
year2_start <- 2130 + 365         # Year 2: day 2342 to 2706
year2_end   <- year2_start + 365 - 1

# Column with incidence
colname <- "n_inc_clinical_182.5_3649"

# Bednets and names
bednets <- 1:3
net_names <- c("Pyrethroid-Chlorfenapyr", "Pyrethroid-Only", "Pyrethroid-Pyriproxyfen")

# -----------------------------
# Function to calculate cases averted over a day range
# -----------------------------
cases_averted_days <- function(control_df, bednet_df, start_day, end_day, colname) {
  control_sum <- control_df %>% 
    filter(timestep >= start_day & timestep <= end_day) %>% 
    summarise(total = sum(.data[[colname]], na.rm = TRUE)) %>% pull(total)
  
  bednet_sum  <- bednet_df %>% 
    filter(timestep >= start_day & timestep <= end_day) %>% 
    summarise(total = sum(.data[[colname]], na.rm = TRUE)) %>% pull(total)
  
  cases_averted <- control_sum - bednet_sum
  pct_averted <- ifelse(control_sum == 0, NA, 100 * cases_averted / control_sum)
  
  return(c(cases_averted = cases_averted, pct_averted = pct_averted))
}

# -----------------------------
# Loop over bednets
# -----------------------------
cases_summary_benin <- lapply(seq_along(bednets), function(i) {
  bed <- bednets[i]
  control_df <- get(paste0("output_control", bed))
  bednet_df  <- get(paste0("output_bednet", bed))
  
  y1 <- cases_averted_days(control_df, bednet_df, year1_start, year1_end, colname)
  y2 <- cases_averted_days(control_df, bednet_df, year2_start, year2_end, colname)
  overall <- cases_averted_days(control_df, bednet_df, year1_start, year2_end, colname)
  
  data.frame(
    Net = net_names[i],
    Year1_Cases = y1["cases_averted"],
    Year1_Pct   = y1["pct_averted"],
    Year2_Cases = y2["cases_averted"],
    Year2_Pct   = y2["pct_averted"],
    Overall_Cases = overall["cases_averted"],
    Overall_Pct   = overall["pct_averted"]
  )
}) %>% bind_rows()

# -----------------------------
# Convert all columns to numeric
# -----------------------------
cases_summary_benin_clean <- cases_summary_benin %>%
  mutate(
    Year1_Cases   = as.numeric(Year1_Cases),
    Year1_Pct     = as.numeric(Year1_Pct),
    Year2_Cases   = as.numeric(Year2_Cases),
    Year2_Pct     = as.numeric(Year2_Pct),
    Overall_Cases = as.numeric(Overall_Cases),
    Overall_Pct   = as.numeric(Overall_Pct)
  )

# -----------------------------
# Replace Inf/-Inf with NA
# -----------------------------
cases_summary_benin_clean[cases_summary_benin_clean == Inf | 
                            cases_summary_benin_clean == -Inf] <- NA

# -----------------------------
# Pivot data for plotting
# -----------------------------
plot_df_benin <- cases_summary_benin_clean %>%
  select(Net, Year1_Pct, Year2_Pct, Overall_Pct) %>%
  pivot_longer(cols = -Net, names_to = "Period", values_to = "Pct_Averted") %>%
  mutate(Period = recode(Period, 
                         "Year1_Pct" = "Year 1",
                         "Year2_Pct" = "Year 2",
                         "Overall_Pct" = "Overall")) %>%
  filter(!is.na(Pct_Averted))  # Skip invalid values

# -----------------------------
# Plot % Cases Averted
# -----------------------------
ggplot(plot_df_benin, aes(x = Net, y = Pct_Averted, fill = Period)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("Year 1" = "skyblue", "Year 2" = "salmon", "Overall" = "seagreen")) +
  labs(
    x = "Bednet Type",
    y = "Cases Averted (%)",
    fill = "Period",
    title = "Benin"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )








#######################################################################################3

########combined in one plot

