


#################################################################################################################################
##################################now i want to add previous distribution in 2014
####IG2



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


output_bednet1 <- run_simulation(
  timesteps = parameters$timesteps,
  parameters = parameters
)

sim_length<- 10*365






sim_length <- 10 * 365  # 10 years
cols <- c("darkgreen", "aquamarine3", "#E69F00")  # colors for lines

plot_combined_prevalence1 <- function() {
  # -----------------------------
  # Bed net distribution parameters (original days)
  # -----------------------------
  bednetstimesteps <- c(1, 2310)  # actual distribution days
  
  # -----------------------------
  # Slice start day: June 1, 2019
  # -----------------------------
  slice_start <- 1977  # 1 June 2019
  
  # -----------------------------
  # Prevalence and ±10% uncertainty
  # -----------------------------
  prev_mean <- output_bednet1$n_detect_lm_0_36499 / output_bednet1$n_age_0_36499
  prev_lower <- pmax(0, prev_mean * 0.9)
  prev_upper <- pmin(1, prev_mean * 1.1)
  
  # -----------------------------
  # Slice the vectors
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
  # Shaded uncertainty band
  # -----------------------------
  polygon(
    c(timesteps, rev(timesteps)),
    c(prev_lower, rev(prev_upper)),
    col = adjustcolor(cols[1], alpha.f = 0.25),
    border = NA
  )
  
  # -----------------------------
  # Mean prevalence line
  # -----------------------------
  lines(timesteps, prev_mean, col = cols[1], lwd = 2)
  
  # -----------------------------
  # Custom x-axis: yearly ticks (only visible years)
  # -----------------------------
  years <- 2014:2023
  year_days <- seq(1, sim_length, by = 365)
  year_days_visible <- year_days[year_days >= slice_start]
  years_visible <- years[year_days >= slice_start]
  axis(1, at = year_days_visible, labels = years_visible)
  
  # -----------------------------
  # Bed net distribution lines (only visible in slice)
  # -----------------------------
  bednet_visible <- bednetstimesteps[bednetstimesteps >= slice_start]
  if(length(bednet_visible) > 0){
    abline(v = bednet_visible, col = "darkgray", lty = "dashed", lwd = 1.5)
    text(bednet_visible + 30, 0.65, "Bed Net Dist.", pos = 4, cex = 0.8, srt = 90)
  }
  
  # -----------------------------
  # Legend
  # -----------------------------
  legend(x = min(timesteps) + 200, y = 0.65,
         legend = c("Model prevalence - Interceptor G2", "Observed data prevalence"),
         col = c(cols[1], adjustcolor(cols[1], alpha.f = 0.25)),
         lty = c(1, NA),
         lwd = c(2, NA),
         pch = c(NA, 15),
         pt.cex = 1,
         box.lty = 0)
  
  # -----------------------------
  # Observed points (only if visible)
  # -----------------------------
  obs_x <- c(2190, 2496, 2859, 3195)
  obs_y <- c(0.407, 0.157, 0.279, 0.225)
  obs_keep <- obs_x >= slice_start
  if(any(obs_keep)){
    points(obs_x[obs_keep], obs_y[obs_keep], col = "darkgreen", pch = 15, cex = 1)
  }
}

# Run the sliced plot
plot_combined_prevalence1()










#######################################################################################################################
##############################################################################################################################
####IG1



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
target_pfpr <- c(0.47, 0.465)  # desired PfPR at these days


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
  coverages = c(0.58, 0.58),  # Each round is distributed to 50% of the population.
  retention = 10^10, # Nets are kept on average 5 years
  dn0 = matrix(c(0.11, 0.11, 0.11,
                 0.11, 0.11, 0.11),nrow = 2,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.58, 0.58,0.58,
                0.58, 0.58, 0.58),nrow = 2,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24,
                 0.24, 0.24, 0.24),nrow = 2,ncol = 3,byrow = TRUE),
 # gamman = rep(1.95 * 365, 2) # Vector of bed net half-lives for each distribution timestep
  gamman = c(3.17, 1.95)  # * 365, 2) # Vector of bed net half-lives for each distribution timestep
  
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


output_bednet2 <- run_simulation(
  timesteps = parameters$timesteps,
  parameters = parameters
)

sim_length<- 10*365






sim_length <- 10 * 365  # 10 years
cols <- c("darkblue" ,"darkgreen", "aquamarine3", "#E69F00")  # colors for lines

plot_combined_prevalence2 <- function() {
  # -----------------------------
  # Bed net distribution parameters (original days)
  # -----------------------------
  bednetstimesteps <- c(1, 2310)  # actual distribution days
  
  # -----------------------------
  # Slice start day: June 1, 2019
  # -----------------------------
  slice_start <- 1977  # 1 June 2019
  
  # -----------------------------
  # Prevalence and ±10% uncertainty
  # -----------------------------
  prev_mean <- output_bednet2$n_detect_lm_0_36499 / output_bednet2$n_age_0_36499
  prev_lower <- pmax(0, prev_mean * 0.9)
  prev_upper <- pmin(1, prev_mean * 1.1)
  
  # -----------------------------
  # Slice the vectors
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
  # Shaded uncertainty band
  # -----------------------------
  polygon(
    c(timesteps, rev(timesteps)),
    c(prev_lower, rev(prev_upper)),
    col = adjustcolor(cols[1], alpha.f = 0.25),
    border = NA
  )
  
  # -----------------------------
  # Mean prevalence line
  # -----------------------------
  lines(timesteps, prev_mean, col = cols[1], lwd = 2)
  
  # -----------------------------
  # Custom x-axis: yearly ticks (only visible years)
  # -----------------------------
  years <- 2014:2023
  year_days <- seq(1, sim_length, by = 365)
  year_days_visible <- year_days[year_days >= slice_start]
  years_visible <- years[year_days >= slice_start]
  axis(1, at = year_days_visible, labels = years_visible)
  
  # -----------------------------
  # Bed net distribution lines (only visible in slice)
  # -----------------------------
  bednet_visible <- bednetstimesteps[bednetstimesteps >= slice_start]
  if(length(bednet_visible) > 0){
    abline(v = bednet_visible, col = "darkgray", lty = "dashed", lwd = 1.5)
    text(bednet_visible + 30, 0.65, "Bed Net Dist.", pos = 4, cex = 0.8, srt = 90)
  }
  
  # -----------------------------
  # Legend
  # -----------------------------
  legend(x = min(timesteps) + 200, y = 0.65,
         legend = c("Model prevalence - Interceptor G2", "Observed data prevalence"),
         col = c(cols[1], adjustcolor(cols[1], alpha.f = 0.25)),
         lty = c(1, NA),
         lwd = c(2, NA),
         pch = c(NA, 15),
         pt.cex = 1,
         box.lty = 0)
  
  # -----------------------------
  # Observed points (only if visible)
  # -----------------------------
  obs_x <- c(2190, 2496, 2859, 3196)
  obs_y <- c(0.465, 0.28, 0.387, 0.262)
  obs_keep <- obs_x >= slice_start
  if(any(obs_keep)){
    points(obs_x[obs_keep], obs_y[obs_keep], col = "darkblue", pch = 15, cex = 1)
  }
}

# Run the sliced plot
plot_combined_prevalence2()




#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################

#RG









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
target_pfpr <- c(0.47, 0.431)  # desired PfPR at these days


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


output_bednet3 <- run_simulation(
  timesteps = parameters$timesteps,
  parameters = parameters
)

sim_length<- 10*365






sim_length <- 10 * 365  # 10 years
cols <- c( "darkorange", "darkblue" ,"darkgreen", "aquamarine3", "#E69F00")  # colors for lines



plot_combined_prevalence3 <- function() {
  # -----------------------------
  # Bed net distribution parameters (original days)
  # -----------------------------
  bednetstimesteps <- c(1, 2310)  # actual distribution days
  
  # -----------------------------
  # Slice start day: June 1, 2019
  # -----------------------------
  slice_start <- 1977  # 1 June 2019
  
  # -----------------------------
  # Prevalence and ±10% uncertainty
  # -----------------------------
  prev_mean <- output_bednet3$n_detect_lm_0_36499 / output_bednet3$n_age_0_36499
  prev_lower <- pmax(0, prev_mean * 0.9)
  prev_upper <- pmin(1, prev_mean * 1.1)
  
  # -----------------------------
  # Slice the vectors
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
  # Shaded uncertainty band
  # -----------------------------
  polygon(
    c(timesteps, rev(timesteps)),
    c(prev_lower, rev(prev_upper)),
    col = adjustcolor(cols[1], alpha.f = 0.25),
    border = NA
  )
  
  # -----------------------------
  # Mean prevalence line
  # -----------------------------
  lines(timesteps, prev_mean, col = cols[1], lwd = 2)
  
  # -----------------------------
  # Custom x-axis: yearly ticks (only visible years)
  # -----------------------------
  years <- 2014:2023
  year_days <- seq(1, sim_length, by = 365)
  year_days_visible <- year_days[year_days >= slice_start]
  years_visible <- years[year_days >= slice_start]
  axis(1, at = year_days_visible, labels = years_visible)
  
  # -----------------------------
  # Bed net distribution lines (only visible in slice)
  # -----------------------------
  bednet_visible <- bednetstimesteps[bednetstimesteps >= slice_start]
  if(length(bednet_visible) > 0){
    abline(v = bednet_visible, col = "darkgray", lty = "dashed", lwd = 1.5)
    text(bednet_visible + 30, 0.65, "Bed Net Dist.", pos = 4, cex = 0.8, srt = 90)
  }
  
  # -----------------------------
  # Legend
  # -----------------------------
  legend(x = min(timesteps) + 200, y = 0.65,
         legend = c("Model prevalence - Interceptor G2", "Observed data prevalence"),
         col = c(cols[1], adjustcolor(cols[1], alpha.f = 0.25)),
         lty = c(1, NA),
         lwd = c(2, NA),
         pch = c(NA, 15),
         pt.cex = 1,
         box.lty = 0)
  
  # -----------------------------
  # Observed points (only if visible)
  # -----------------------------
  obs_x <- c(2190, 2496, 2859, 3196)
  obs_y <- c(0.431, 0.269, 0.382, 0.284)
  obs_keep <- obs_x >= slice_start
  if(any(obs_keep)){
    points(obs_x[obs_keep], obs_y[obs_keep], col = "darkorange", pch = 15, cex = 1)
  }
}

# Run the sliced plot
plot_combined_prevalence3()
















#################################################plot these three plese 
########333anotherone





plot_combined_all_original <- function() {
  
  slice_start <- 1977  # 1 June 2019
  sim_length <- 10*365
  
  # -----------------------------
  # Extract timesteps & prevalence
  # -----------------------------
  timesteps <- output_bednet1$timestep[output_bednet1$timestep >= slice_start]
  
  # IG2
  prev_IG2_mean <- output_bednet1$n_detect_lm_0_36499[output_bednet1$timestep >= slice_start] /
    output_bednet1$n_age_0_36499[output_bednet1$timestep >= slice_start]
  prev_IG2_lower <- pmax(0, prev_IG2_mean * 0.9)
  prev_IG2_upper <- pmin(1, prev_IG2_mean * 1.1)
  
  # IG1
  prev_IG1_mean <- output_bednet2$n_detect_lm_0_36499[output_bednet2$timestep >= slice_start] /
    output_bednet2$n_age_0_36499[output_bednet2$timestep >= slice_start]
  prev_IG1_lower <- pmax(0, prev_IG1_mean * 0.9)
  prev_IG1_upper <- pmin(1, prev_IG1_mean * 1.1)
  
  # RG
  prev_RG_mean <- output_bednet3$n_detect_lm_0_36499[output_bednet3$timestep >= slice_start] /
    output_bednet3$n_age_0_36499[output_bednet3$timestep >= slice_start]
  prev_RG_lower <- pmax(0, prev_RG_mean * 0.9)
  prev_RG_upper <- pmin(1, prev_RG_mean * 1.1)
  
  cols <- c("darkgreen", "darkblue", "darkorange")
  
  # -----------------------------
  # Plot setup
  # -----------------------------
  plot(timesteps, prev_IG2_mean, type="n",
       xlab="Years", ylab="Malaria Prevalence for all ages",
       ylim=c(0,0.8), xaxs="i", yaxs="i",
       main="Benin Trials",
       cex.main=1.4, cex.lab=1.2, xaxt="n")
  
  # -----------------------------
  # Shaded uncertainty bands
  # -----------------------------
  polygon(c(timesteps, rev(timesteps)), c(prev_IG2_lower, rev(prev_IG2_upper)),
          col=adjustcolor(cols[1], alpha.f=0.25), border=NA)
  polygon(c(timesteps, rev(timesteps)), c(prev_IG1_lower, rev(prev_IG1_upper)),
          col=adjustcolor(cols[2], alpha.f=0.25), border=NA)
  polygon(c(timesteps, rev(timesteps)), c(prev_RG_lower, rev(prev_RG_upper)),
          col=adjustcolor(cols[3], alpha.f=0.25), border=NA)
  
  # -----------------------------
  # Mean prevalence lines
  # -----------------------------
  lines(timesteps, prev_IG2_mean, col=cols[1], lwd=1)
  lines(timesteps, prev_IG1_mean, col=cols[2], lwd=1)
  lines(timesteps, prev_RG_mean,  col=cols[3], lwd=1)
  
  # -----------------------------
  # Custom x-axis: yearly ticks
  # -----------------------------
  years <- 2014:2023
  year_days <- seq(1, sim_length, by=365)
  year_days_visible <- year_days[year_days >= slice_start]
  years_visible <- years[year_days >= slice_start]
  axis(1, at=year_days_visible, labels=years_visible)
  
  # -----------------------------
  # Bed net distribution lines
  # -----------------------------
  bednet_days <- c(1, 2310)
  bednet_visible <- bednet_days[bednet_days >= slice_start]
  if(length(bednet_visible) > 0){
    abline(v = bednet_visible, col="darkgray", lty="dashed", lwd=1.5)
    text(bednet_visible + 30, 0.65, "Bed Net Dist.", srt=0, pos=4, cex=0.8)
  }
  
  # -----------------------------
  # Observed points
  # -----------------------------
  
  
  points(c(2190, 2496, 2859), c(0.407, 0.157, 0.279), col=cols[1], pch=15)
  points(c(2190, 2496, 2859), c(0.465, 0.28, 0.387), col=cols[2], pch=15)
  points(c(2190, 2496, 2859), c(0.431, 0.269, 0.382), col=cols[3], pch=15)
  

  
  #legend(x = 2800, y = 0.74,  # adjust these to shift horizontally and vertically
  #      legend = c("Interceptor G2", "Interceptor", "Royal Guard"),
  #     col = cols,            # only the colors of your lines
  #    lty = 1,               # solid lines
  #   lwd = 1.5,             # line width
  #  pch =c(NA,NA,NA,15),              # no points
  # pt.cex = 1,
  #box.lty = 0) 
  
  
  
  legend(x = 2800, y = 0.74,
         legend = c("Interceptor G2","Interceptor","Royal Guard","Observed prevalence"),
         col = c(cols, "black"),        # colors for lines and points
         lty = c(1,1,1,NA),             # NA for points
         lwd = c(1,1,1,NA),             # NA for points
         pch = c(NA,NA,NA,15),          # points for observed only
         pt.cex = 1,
         box.lty = 0)
}

# Run the combined plot
plot_combined_all_original()











################################plot predicted vs actual



# Your observed timesteps
obs_times <- c(2190, 2496, 2859)

# Extract predicted prevalence at these exact timesteps
pred_IG2 <- approx(
  x = output_bednet1$timestep,
  y = output_bednet1$n_detect_lm_0_36499 / output_bednet1$n_age_0_36499,
  xout = obs_times
)$y

pred_IG1 <- approx(
  x = output_bednet2$timestep,
  y = output_bednet2$n_detect_lm_0_36499 / output_bednet2$n_age_0_36499,
  xout = obs_times
)$y

pred_RG <- approx(
  x = output_bednet3$timestep,
  y = output_bednet3$n_detect_lm_0_36499 / output_bednet3$n_age_0_36499,
  xout = obs_times
)$y

# Combine into a clean table
predicted_values <- data.frame(
  time = obs_times,
  IG2_pred = pred_IG2,
  IG1_pred = pred_IG1,
  RG_pred  = pred_RG
)

print(predicted_values)













##########################################################################################################################################################################################
#################################################################################TANZANIA TRIAL LETS SEE #########################################################################################################
##########################################################################################################################################################################################




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
sim_length<- 7*365
year<- 365
month<- year/12



target_days <- c( 1)  # mid-year of years 1, 2, 3
target_pfpr <- c( 0.427)  # desired PfPR at these days


point_pfpr_summary <- function(x, days = target_days) {
  pfpr <- x$n_detect_lm_0_36499 / x$n_age_0_36499
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
site_par <- get_parameters(
  list(
    human_population = human_population,
    clinical_incidence_rendering_min_ages = 182.5,
    clinical_incidence_rendering_max_ages = 14 * year,
    prevalence_rendering_min_ages = 182.5,
    prevalence_rendering_max_ages = 14 * year,
    age_group_rendering_min_ages = 182.5,
    age_group_rendering_max_ages = 14* year,
    model_seasonality = TRUE,
    g0 = 89.02751,
    g  = c(60.86729, -27.62949, 5.56094),
    h  = c(30.55385, -27.37595 , -28.53617),
    Q0 = c(0.92, 0.71, 0.94)
    
  )
)



## -------------------------------
## Mosquito species (ORDER MATTERS)
## -------------------------------
site_par <- set_species(
  parameters = site_par,
  species = list(fun_params, gamb_params, arab_params),
  proportions = c(0.1055746 , 0.5411252, 0.3533002 )
)





bednetstimesteps <- c(1,1096, 1369)

#parameters <- get_parameters()

parameters <- set_bednets(
  site_par,
  timesteps = bednetstimesteps,
  coverages = c(0.40, 0.58, 0.58),  # Each round is distributed to 50% of the population.
  retention = 10^10, # Nets are kept on average 5 years
  dn0 = matrix(c(0.11, 0.11, 0.11,
                 0.30, 0.30, 0.30,
                 0.30, 0.30, 0.30),nrow = 3,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.72, 0.72,0.72,
                0.49, 0.49, 0.49,
                0.30, 0.30, 0.30 ),nrow = 3,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24,
                 0.24, 0.24, 0.24,
                 0.24, 0.24, 0.24),nrow = 3,ncol = 3,byrow = TRUE),
  gamman = rep(2.65 * 365, 3) # Vector of bed net half-lives for each distribution timestep
)




parameters$timesteps <- 365 * 7  # total simulation length

#  Calibration 

out <- calibrate(
  parameters = parameters,
  target = target_pfpr,
  summary_function = point_pfpr_summary,
  eq_prevalence = target_pfpr[1]   # baseline starting PfPR
)


# Run simulation using calibrated EIR


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


output_bednet4 <- run_simulation(
  timesteps = parameters$timesteps,
  parameters = parameters
)

sim_length<- 7*365






sim_length <- 7 * 365  # 10 years
cols <- c("darkgreen", "aquamarine3", "#E69F00")  # colors for lines

plot_combined_prevalence4 <- function() {
  # -----------------------------
  # Bed net distribution parameters (original days)
  # -----------------------------
  bednetstimesteps <- c(1, 1096, 1396)  # actual distribution days
  
  # -----------------------------
  # Slice start day: June 1, 2019
  # -----------------------------
  slice_start <- 1096  # 1 June 2019
  
  # -----------------------------
  # Prevalence and ±10% uncertainty
  # -----------------------------
  prev_mean <- output_bednet4$n_detect_lm_0_36499 / output_bednet4$n_age_0_36499
  prev_lower <- pmax(0, prev_mean * 0.9)
  prev_upper <- pmin(1, prev_mean * 1.1)
  
  # -----------------------------
  # Slice the vectors
  # -----------------------------
  keep <- output_bednet4$timestep >= slice_start
  timesteps <- output_bednet4$timestep[keep]
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
  # Shaded uncertainty band
  # -----------------------------
  polygon(
    c(timesteps, rev(timesteps)),
    c(prev_lower, rev(prev_upper)),
    col = adjustcolor(cols[1], alpha.f = 0.25),
    border = NA
  )
  
  # -----------------------------
  # Mean prevalence line
  # -----------------------------
  lines(timesteps, prev_mean, col = cols[1], lwd = 2)
  
  # -----------------------------
  # Custom x-axis: yearly ticks (only visible years)
  # -----------------------------
  years <- 2015:2023
  year_days <- seq(1, sim_length, by = 365)
  year_days_visible <- year_days[year_days >= slice_start]
  years_visible <- years[year_days >= slice_start]
  axis(1, at = year_days_visible, labels = years_visible)
  
  # -----------------------------
  # Bed net distribution lines (only visible in slice)
  # -----------------------------
  bednet_visible <- bednetstimesteps[bednetstimesteps >= slice_start]
  if(length(bednet_visible) > 0){
    abline(v = bednet_visible, col = "darkgray", lty = "dashed", lwd = 1.5)
    text(bednet_visible + 30, 0.65, "Bed Net Dist.", pos = 4, cex = 0.8, srt = 90)
  }
  
  # -----------------------------
  # Legend
  # -----------------------------
  legend(x = min(timesteps) + 200, y = 0.65,
         legend = c("Model prevalence - Interceptor G2", "Observed data prevalence"),
         col = c(cols[1], adjustcolor(cols[1], alpha.f = 0.25)),
         lty = c(1, NA),
         lwd = c(2, NA),
         pch = c(NA, 15),
         pt.cex = 1,
         box.lty = 0)
  

  # -----------------------------
  # Observed points (only if visible)
  # -----------------------------
  obs_x <- c(1369, 1826, 2009, 2196)
  obs_y <- c(0.427, 0.156, 0.409, 0.256)
  obs_keep <- obs_x >= slice_start
  if(any(obs_keep)){
    points(obs_x[obs_keep], obs_y[obs_keep], col = "darkgreen", pch = 15, cex = 1)
  }
}

# Run the sliced plot
plot_combined_prevalence4()






###########let me try with benin data



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
  prev_lower <- pmax(0, prev_mean * 0.7)
  prev_upper <- pmin(1, prev_mean * 1.3)
  
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
  obs_y <- c(0.40, 0.157, 0.279, 0.225)    # prevalence
  # Calculate Wilson 95% CI manually (rounded)
  obs_lower <- c(0.38, 0.139, 0.257, 0.235)      # 95% CI lower
  obs_upper <- c(0.429, 0.177, 0.302, 0.28 )      # 95% CI upper
  
  
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































