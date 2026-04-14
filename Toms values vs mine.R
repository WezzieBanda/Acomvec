

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
human_population<- 10000
sim_length<- 10 *365
year<- 365
month<- year/12



target_days <- c(1 , 2190)  # mid-year of years 1, 2, 3
#target_pfpr <- c(0.47, 0.465)  # desired PfPR at these days
target_pfpr <- c(0.39, 0.407)  # desired PfPR at these days

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
  coverages = c(0.40, 0.90),  # Each round is distributed to 50% of the population.
  retention = 10^10, # Nets are kept on average 5 years
  dn0 = matrix(c(0.392, 0.392 , 0.392,
                 0.57, 0.57, 0.57),nrow = 2,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.59, 0.59,0.59,
                0.42, 0.42, 0.42),nrow = 2,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24,
                 0.24, 0.24, 0.24),nrow = 2,ncol = 3,byrow = TRUE),
  gamman = rep(3.19 * 365, 2) # Vector of bed net half-lives for each distribution timestep
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
      # 95% CI upper
  
  obs_x <- c(2190, 2496, 2859, 3195)            # 6 months and 18 months days
  obs_y <- c(0.407, 0.157, 0.279, 0.225)    # prevalence
  # Calculate Wilson 95% CI manually (rounded)
  obs_lower <- c(0.382, 0.139, 0.257, 0.235)      # 95% CI lower
  obs_upper <- c(0.433, 0.177, 0.302, 0.28 ) 
  
  
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
####################################Plot all three together 









plot_combined_prevalence_all <- function() {
  # -----------------------------
  slice_start <- 1977
  cols <- c("darkgreen", "lightgreen", "darkred")
  
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
 # lines(timesteps3, prev_mean3, col = cols[3], lwd = 1)
  
  # -----------------------------
  # Add observed points with 95% CI
  obs_list <- list(
    list(x = c(2190, 2496, 2859, 3195), y = c(0.407, 0.157, 0.279, 0.255), lower = c(0.382, 0.139, 0.257, 0.235), upper = c(0.433, 0.177, 0.302, 0.28), col = cols[1])
    #list(x = c(2190, 2496, 2859, 3195), y = c(0.465, 0.28, 0.387, 0.262), lower = c(0.439, 0.257, 0.362, 0.235), upper = c(0.490, 0.303, 0.412, 0.28), col = cols[2])
   # list(x = c(2190, 2496, 2859, 3195), y = c(0.431, 0.269, 0.382, 0.284), lower = c(0.406, 0.247, 0.357, 0.266), upper = c(0.456, 0.292, 0.407, 0.30), col = cols[3])
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


#(2190, 2496, 2859, 3195)   
target_days <- c(  2190, 3195) # 3195 )  # mid-year of years 1, 2, 3
#target_pfpr <- c(0.47, 0.465)  # desired PfPR at these days
#target_pfpr <- c(0.28, 0.407, 0.14, 0.12, 0.24)  # desired PfPR at these days
#target_pfpr <- c( 0.407, 0.14, 0.18, 0.24)  # desired PfPR at these days


#target_pfpr <- c( 0.402, 0.14 ) #  0.24)  # desired PfPR at these days

target_pfpr <- c( 0.4,0.24 ) #  0.24)  # desired PfPR at these days


#target_pfpr <- c( 0.30, 0.14 )
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
  coverages = c(0.40, 0.90),  # Each round is distributed to 50% of the population.
  retention = 10^10, # Nets are kept on average 5 years
  dn0 = matrix(c(0.392, 0.392 , 0.392,
                 0.57, 0.57, 0.57),nrow = 2,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.59, 0.59,0.59,
                0.42, 0.42, 0.42),nrow = 2,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24,
                 0.24, 0.24, 0.24),nrow = 2,ncol = 3,byrow = TRUE),
  gamman = rep(3.19 * 365, 2) # Vector of bed net half-lives for each distribution timestep
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
  # 95% CI upper
  
  obs_x <- c(2190, 2496, 2859, 3195)            # 6 months and 18 months days
  obs_y <- c(0.407, 0.157, 0.279, 0.225)    # prevalence
  # Calculate Wilson 95% CI manually (rounded)
  obs_lower <- c(0.382, 0.139, 0.257, 0.235)      # 95% CI lower
  obs_upper <- c(0.433, 0.177, 0.302, 0.28 ) 
  
  
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

length(target_days)
target_days

