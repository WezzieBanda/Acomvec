
library(malariasimulation)
library(curl)
library(site)
library(dplyr)
library(ggplot2)

# Simulation parameters



year <- 365
month<-year/12
sim_length <- 4 * year
human_population <- 10000
starting_EIR <- 65

#available_sites()
BEN <- fetch_site(iso3c = "BEN", admin_level = 1, urban_rural = TRUE)
BEN_Zou <- BEN |> subset_site(site_filter = BEN$sites |> filter(name_1 == "Zou")|> filter(urban_rural =="rural"))


# Convert site information to malariasimulation parameters
site_par <- site_parameters(
  interventions = BEN_Zou$interventions,
  demography = BEN_Zou$demography,
  vectors = BEN_Zou$vectors$vector_species,
  seasonality = BEN_Zou$seasonality$seasonality_parameters,
  eir = BEN_Zou$eir$eir,
  overrides = list(
    human_population =human_population
  )
) |> set_epi_outputs(clinical_incidence = c(2,10)*365)




#site_par <- site_parameters(
 # interventions = BEN_Zou$interventions,
  #demography = BEN_Zou$demography,
  #vectors = BEN_Zou$vectors$vector_species,
  #seasonality = BEN_Zou$seasonality$seasonality_parameters,
  #eir = BEN_Zou$eir$eir,
  #overrides = list(
  #  human_population = human_population
  #)
#)


# Set drug parameters with seasonality
drug_params <- set_drugs(site_par,  list(AL_params, SP_AQ_params, DHA_PQP_params))



# Set treatment parameters with seasonality
treatment_params1 <- set_clinical_treatment(
  parameters = drug_params,
  drug = 1,
  timesteps =  c(100),
  coverages =  c(0.4)
)

treatment_params2 <- set_clinical_treatment(
  parameters = drug_params,
  drug = 2,
  timesteps =  c(100),
  coverages =  c(0.4)
)
treatment_params<- set_clinical_treatment(
  parameters = treatment_params1,
  drug = 3,
  timesteps =  c(100),
  coverages =  c(0.2)
)
# Use set_equilibrium to update the parameter set for a given initial EIR
treatment_params <- set_equilibrium(parameters=treatment_params, BEN_Zou$eir$eir)

###set bednets

# Set equilibrium
#species_params <- set_equilibrium(parameters = species_params, init_EIR = starting_EIR)

# Run control simulation (without bed nets)
output_control <- run_simulation(timesteps = sim_length, parameters = treatment_params)



bednetstimesteps <- c(1*365)

IG1 <- set_bednets(
  treatment_params,
  timesteps = bednetstimesteps,
  coverages = c(0.362),  # Each round is distributed to 50% of the population.
  retention = 10^10, # Nets are kept on average 5 years
  dn0 = matrix(c(0.062242872, 0.062242872, 0.062242872),nrow = 1,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.667323298, 0.667323298, 0.667323298),nrow = 1,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24),nrow = 1,ncol = 3,byrow = TRUE),
  gamman = rep(1.90* 365) # Vector of bed net half-lives for each distribution timestep
)



IG2 <- set_bednets(
  treatment_params,
  timesteps = bednetstimesteps,
  coverages = c(.9),  # Each round is distributed to 50% of the population.
  retention = 10^10, # Nets are kept on average 5 years
  dn0 = matrix(c(0.220157756, 0.220157756, 0.220157756),nrow = 1,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.584550395, 0.584550395, 0.584550395),nrow = 1,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24),nrow = 1,ncol = 3,byrow = TRUE),
  gamman = rep(2.65 * 365) # Vector of bed net half-lives for each distribution timestep
)



P3 <- set_bednets(
  treatment_params,
  timesteps = bednetstimesteps,
  coverages = c(.7),  # Each round is distributed to 50% of the population.
  retention = 10^10, # Nets are kept on average 5 years
  dn0 = matrix(c(0.143356157, 0.143356157, 0.143356157),nrow = 1,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.679186299, 0.679186299, 0.679186299),nrow = 1,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24),nrow = 1,ncol = 3,byrow = TRUE),
  gamman = rep(2.16* 365) # Vector of bed net half-lives for each distribution timestep
)




RG <- set_bednets(
  treatment_params,
  timesteps = bednetstimesteps,
  coverages = c(.7),  # Each round is distributed to 50% of the population.
  retention = 10^10, # Nets are kept on average 5 years
  dn0 = matrix(c(0.162232052, 0.162232052, 0.162232052),nrow = 1,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.680353824, 0.680353824, 0.680353824),nrow = 1,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24),nrow = 1,ncol = 3,byrow = TRUE),
  gamman = rep(2.11* 365) # Vector of bed net half-lives for each distribution timestep
)




output_bednet1 <- run_simulation( timesteps = sim_length, parameters = IG1)
output_bednet2 <- run_simulation( timesteps = sim_length, parameters = IG2)
output_bednet3 <- run_simulation( timesteps = sim_length, parameters = P3)
output_bednet4 <- run_simulation( timesteps = sim_length, parameters = RG)




cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")




plot_combined_prevalence <- function() {
  
  bednetstimesteps <- c(1 * 365)
  bednetstimesteps_shifted <- bednetstimesteps + 25
  
  x_positions <- seq(from = bednetstimesteps, by = 365, length.out = 4)
  x_label_text <- 1:4
  
  
  par(mar = c(5, 6, 4, 2))  # bottom, left, top, right
  
  plot(
    output_control$timestep,
    output_control$n_detect_lm_730_3649 / output_control$n_age_730_3649,
    type = "l", col = cols[1], lwd = 2,
    xlab = "Time in Years after Bednet Distribution",
    ylab = "% of Malaria Prevalence in children aged \n 2-10 years (per 10,000 population)",
    xaxt = "n",
    ylim = c(0, 1),
    xaxs = "i", yaxs = "i",
    xlim = c(bednetstimesteps, bednetstimesteps + 3 * 365)
  )
  
  axis(1, at = x_positions, labels = x_label_text)
  
  lines(output_bednet1$timestep,
        output_bednet1$n_detect_lm_730_3649 / output_bednet1$n_age_730_3649,
        col = cols[2], lwd = 2)
  
  lines(output_bednet2$timestep,
        output_bednet2$n_detect_lm_730_3649 / output_bednet2$n_age_730_3649,
        col = cols[3], lwd = 2)
  
  lines(output_bednet3$timestep,
        output_bednet3$n_detect_lm_730_3649 / output_bednet3$n_age_730_3649,
        col = cols[4], lwd = 2)
  
  lines(output_bednet4$timestep,
        output_bednet4$n_detect_lm_730_3649 / output_bednet4$n_age_730_3649,
        col = cols[5], lwd = 2)
  
  abline(v = bednetstimesteps_shifted, col = "darkgray", lty = "dashed")
  text(bednetstimesteps_shifted + 5, 0.9, "36M Aged Nets", pos = 4, cex = 0.9)
  
  legend(
    "topright",
    legend = c("No bednets", "IG1", "IG2", "P3", "RG"),
    col = cols[1:5],
    lty = 1, lwd = 2, box.lty = 0
  )
}


plot_combined_prevalence()









##################################################################cases avrted Benin################################################################
#################################################################################################################################################################
#################################################################################################################################################################


plot_combined_prevalence <- function() {
  # Time of actual bed net distribution
  bednetstimesteps <- 365
  bednetstimesteps_shifted <- bednetstimesteps + 25
  
  # Create numeric x-axis labels for years since distribution: 1 to 4
  x_positions <- seq(from = bednetstimesteps, by = 365, length.out = 4)
  x_label_text <- 1:4  # Label as "1", "2", "3", "4" years after intervention
  
  # Set up the plot
  plot(output_control$timestep, output_control$n_detect_lm_730_3649/ output_control$n_age_730_3649,
       type = "l", col = cols[1], lwd = 2,
       xlab = "Time in Years", ylab = "Prevalence",
       xaxt = "n",  # Remove default x-axis
       ylim = c(0, 1), xaxs = "i", yaxs = "i",
       xlim = c(bednetstimesteps, bednetstimesteps + 3 * 365),  # Shift x-axis to end at x = 4
       main = "Estimated Prevalence")
  
  # Add custom x-axis with year labels
  axis(1, at = x_positions, labels = x_label_text)
  
  # Add lines for each intervention scenario
  lines(output_bednet1$timestep, output_bednet1$n_detect_lm_730_3649 / output_bednet1$n_age_730_3649, col = cols[2], lwd = 2)
  lines(output_bednet2$timestep, output_bednet2$n_detect_lm_730_3649 / output_bednet2$n_age_730_3649, col = cols[3], lwd = 2)
  lines(output_bednet3$timestep, output_bednet3$n_detect_lm_730_3649 / output_bednet3$n_age_730_3649, col = cols[4], lwd = 2)
  lines(output_bednet4$timestep, output_bednet4$n_detect_lm_730_3649 / output_bednet4$n_age_730_3649, col = cols[5], lwd = 2)
  # lines(output_bednet5$timestep, output_bednet5$n_detect_91_25550 / output_bednet5$n_91_25550, col = cols[6], lwd = 2)
  
  # Add vertical line for bed net distribution (shifted)
  abline(v = bednetstimesteps_shifted, col = "darkgray", lty = "dashed")
  
  # Add label near the vertical line
  text(x = bednetstimesteps_shifted + 5, y = 0.9, "Bed Net Distribution", pos = 4, cex = 0.9)
  
  # Add legend
  legend("topright", legend = c("Control", "IG1", "IG2", "P3", "RG"),
         col = cols[1:6], lty = 1, lwd = 2, box.lty = 0)
}

plot_combined_prevalence()

########################################################################################
################realistic cases averted in %


# -----------------------------
# Colours and net labels
# -----------------------------
cols <- c("#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
net_names <- c("IG1", "IG2", "OP", "RG")

# -----------------------------
# Function: cumulative cases averted
# -----------------------------
cases_averted_ts <- function(out_int, out_ctrl, start_day) {
  idx <- out_ctrl$timestep >= start_day
  cumsum(out_ctrl$n_inc_clinical_730_3649[idx]) -
    cumsum(out_int$n_inc_clinical_730_3649[idx])
}

# -----------------------------
# Time axis (years after distribution)
# -----------------------------
t_days <- output_control$timestep[output_control$timestep >= bednetstimesteps]
t_years <- (t_days - bednetstimesteps) / 365
t_years_plot <- t_years[t_years <= 3]
t_years_plot <- c(0, t_years_plot)  # include 0 at start

# -----------------------------
# Calculate cumulative cases averted
# -----------------------------
cases_list <- lapply(
  list(output_bednet1, output_bednet2, output_bednet3, output_bednet4),
  function(net) {
    tmp <- cases_averted_ts(net, output_control, bednetstimesteps)
    tmp <- tmp[t_years <= 3]
    c(0, tmp)  # include 0 at start
  }
)

# -----------------------------
# Convert to % of total control cumulative cases
# -----------------------------
total_control <- sum(
  output_control$n_inc_clinical_730_3649[
    output_control$timestep >= bednetstimesteps
  ][t_years <= 3]
)

cases_list_percent <- lapply(cases_list, function(x) {
  (x / total_control) * 100
})

# -----------------------------
# Plot limits
# -----------------------------
ymin <- 0
ymax <- 100  # percentage

# -----------------------------
# Plot setup
# -----------------------------
par(mar = c(4, 6, 2, 6) + 0.1)

plot(
  NA, NA,
  xlim = c(0, 3),
  ylim = c(ymin, ymax),
  xlab = "Time in Years after bednet distribution",
  ylab = "Cumulative Malaria Cases Averted in
          \nchildren aged 2–10 years (%)",
  xaxs = "i", yaxs = "i",
  xaxt = "n", yaxt = "n"
)

title("24M Aged Nets")

# -----------------------------
# Axes
# -----------------------------
axis(side = 1, at = seq(0, 3, by = 0.5), tcl = -0.3)
axis(side = 2, at = seq(0, 100, by = 20), las = 1, tcl = -0.3)

# -----------------------------
# Add slim lines and points
# -----------------------------
for (i in 1:4) {
  lines(
    t_years_plot,
    cases_list_percent[[i]],
    col = cols[i],
    lwd = 0.6     # slimmer, smooth lines
  )
  points(
    t_years_plot,
    cases_list_percent[[i]],
    col = cols[i],
    pch = 16,
    cex = 0.4     # smaller points
  )
}

# -----------------------------
# Legend outside plot (matched line width)
# -----------------------------
par(xpd = TRUE)
legend(
  x = 3.05, y = ymax,
  legend = net_names,
  col = cols[1:4],
  lwd = 0.6,
  pch = 16,
  bty = "n",
  xjust = 0,
  yjust = 1
)




#############################TANZANIA TRIAL########################################################################################################################
###########################################################################################################################################################################
#############################TANZANIA TRIAL########################################################################################################################
###########################################################################################################################################################################
#############################TANZANIA TRIAL########################################################################################################################
###########################################################################################################################################################################
#############################TANZANIA TRIAL########################################################################################################################
###########################################################################################################################################################################
################################################################################################
library(malariasimulation)
library(curl)
library(site)
library(dplyr)
library(ggplot2)

# Simulation parameters
year <- 365
month<-year/12
sim_length <- 4 * year
human_population <- 10000
starting_EIR <- 65

available_sites()
TZA <- fetch_site(iso3c = "TZA", admin_level = 1, urban_rural = TRUE)
TZA_Mwanza <- TZA|> subset_site(site_filter = TZA$sites |> filter(name_1 == "Mwanza")|> filter(urban_rural =="rural"))


# Convert site information to malariasimulation parameters
site_par <- site_parameters(
  interventions = TZA_Mwanza$interventions,
  demography = TZA_Mwanza$demography,
  vectors = TZA_Mwanza$vectors$vector_species,
  seasonality = TZA_Mwanza$seasonality$seasonality_parameters,
  eir = TZA_Mwanza$eir$eir,
  overrides = list(
    human_population =human_population
  )
) |> set_epi_outputs(clinical_incidence = c(2,10)*365)



site_par <- set_equilibrium(parameters=site_par, TZA_Mwanza$eir$eir)
###set bednets

# Set equilibrium
#species_params <- set_equilibrium(parameters = species_params, init_EIR = starting_EIR)

# Run control simulation (without bed nets)
output_control <- run_simulation(timesteps = sim_length, parameters = site_par)



bednetstimesteps <- c(1*365)

IG1 <- set_bednets(
  site_par,
  timesteps = bednetstimesteps,
  coverages = c(.66),  # Each round is distributed to 50% of the population.
  retention = 10^10, # Nets are kept on average 5 years
  dn0 = matrix(c(0.081154132, 0.081154132, 0.081154132),nrow = 1,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.650762568, 0.650762568, 0.650762568),nrow = 1,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24),nrow = 1,ncol = 3,byrow = TRUE),
  gamman = rep(2.25 * 365) # Vector of bed net half-lives for each distribution timestep
)



IG2 <- set_bednets(
  site_par,
  timesteps = bednetstimesteps,
  coverages = c(.9),  # Each round is distributed to 50% of the population.
  retention = 10^10, # Nets are kept on average 5 years
  dn0 = matrix(c(0.164971208, 0.164971208, 0.164971208),nrow = 1,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.564740171, 0.564740171, 0.564740171),nrow = 1,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24),nrow = 1,ncol = 3,byrow = TRUE),
  gamman = rep(2.65 * 365) # Vector of bed net half-lives for each distribution timestep
)





OP <- set_bednets(
  site_par,
  timesteps = bednetstimesteps,
  coverages = c(.66),  # Each round is distributed to 50% of the population.
  retention = 10^10, # Nets are kept on average 5 years
  dn0 = matrix(c(0.170526109, 0.170526109, 0.170526109),nrow = 1,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.583537302, 0.583537302, 0.583537302),nrow = 1,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24),nrow = 1,ncol = 3,byrow = TRUE),
  gamman = rep(2.50* 365) # Vector of bed net half-lives for each distribution timestep
)



RG <- set_bednets(
  site_par,
  timesteps = bednetstimesteps,
  coverages = c(.66),  # Each round is distributed to 50% of the population.
  retention = 10^10, # Nets are kept on average 5 years
  dn0 = matrix(c(0.152538966, 0.152538966, 0.152538966),nrow = 1,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.566177051,0.566177051, 0.566177051),nrow = 1,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24),nrow = 1,ncol = 3,byrow = TRUE),
  gamman = rep(2.44* 365) # Vector of bed net half-lives for each distribution timestep
)




output_bednet1 <- run_simulation( timesteps = sim_length, parameters = IG1)
output_bednet2 <- run_simulation( timesteps = sim_length, parameters = IG2)
output_bednet3 <- run_simulation( timesteps = sim_length, parameters = OP)
output_bednet4 <- run_simulation( timesteps = sim_length, parameters = RG)




cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")




plot_combined_prevalence <- function() {
  
  bednetstimesteps <- c(1 * 365)
  bednetstimesteps_shifted <- bednetstimesteps + 25
  
  x_positions <- seq(from = bednetstimesteps, by = 365, length.out = 4)
  x_label_text <- 1:4
  
  
  par(mar = c(5, 6, 4, 2))  # bottom, left, top, right
  
  plot(
    output_control$timestep,
    output_control$n_detect_lm_730_3649 / output_control$n_age_730_3649,
    type = "l", col = cols[1], lwd = 2,
    xlab = "Time in Years after Bednet Distribution",
    ylab = "% of Malaria Prevalence in children aged \n 2-10 years (per 10,000 population)",
    xaxt = "n",
    ylim = c(0, 1),
    xaxs = "i", yaxs = "i",
    xlim = c(bednetstimesteps, bednetstimesteps + 3 * 365)
  )
  
  axis(1, at = x_positions, labels = x_label_text)
  
  lines(output_bednet1$timestep,
        output_bednet1$n_detect_lm_730_3649 / output_bednet1$n_age_730_3649,
        col = cols[2], lwd = 2)
  
  lines(output_bednet2$timestep,
        output_bednet2$n_detect_lm_730_3649 / output_bednet2$n_age_730_3649,
        col = cols[3], lwd = 2)
  
  lines(output_bednet3$timestep,
        output_bednet3$n_detect_lm_730_3649 / output_bednet3$n_age_730_3649,
        col = cols[4], lwd = 2)
  
  lines(output_bednet4$timestep,
        output_bednet4$n_detect_lm_730_3649 / output_bednet4$n_age_730_3649,
        col = cols[5], lwd = 2)
  
  abline(v = bednetstimesteps_shifted, col = "darkgray", lty = "dashed")
  text(bednetstimesteps_shifted + 5, 0.9, "12M Aged Nets", pos = 4, cex = 0.9)
  
  legend(
    "topright",
    legend = c("No bednets", "IG1", "IG2", "OP", "RG"),
    col = cols[1:5],
    lty = 1, lwd = 2, box.lty = 0
  )
}







plot_combined_prevalence()




##################################################################cases avrted tanzania################################################################
#################################################################################################################################################################
#################################################################################################################################################################


plot_combined_prevalence <- function() {
  # Time of actual bed net distribution
  bednetstimesteps <- 365
  bednetstimesteps_shifted <- bednetstimesteps + 25
  
  # Create numeric x-axis labels for years since distribution: 1 to 4
  x_positions <- seq(from = bednetstimesteps, by = 365, length.out = 4)
  x_label_text <- 1:4  # Label as "1", "2", "3", "4" years after intervention
  
  # Set up the plot
  plot(output_control$timestep, output_control$n_detect_lm_730_3649/ output_control$n_age_730_3649,
       type = "l", col = cols[1], lwd = 2,
       xlab = "Time in Years", ylab = "Prevalence",
       xaxt = "n",  # Remove default x-axis
       ylim = c(0, 1), xaxs = "i", yaxs = "i",
       xlim = c(bednetstimesteps, bednetstimesteps + 3 * 365),  # Shift x-axis to end at x = 4
       main = "Estimated Prevalence")
  
  # Add custom x-axis with year labels
  axis(1, at = x_positions, labels = x_label_text)
  
  # Add lines for each intervention scenario
  lines(output_bednet1$timestep, output_bednet1$n_detect_lm_730_3649 / output_bednet1$n_age_730_3649, col = cols[2], lwd = 2)
  lines(output_bednet2$timestep, output_bednet2$n_detect_lm_730_3649 / output_bednet2$n_age_730_3649, col = cols[3], lwd = 2)
  lines(output_bednet3$timestep, output_bednet3$n_detect_lm_730_3649 / output_bednet3$n_age_730_3649, col = cols[4], lwd = 2)
  lines(output_bednet4$timestep, output_bednet4$n_detect_lm_730_3649 / output_bednet4$n_age_730_3649, col = cols[5], lwd = 2)
  # lines(output_bednet5$timestep, output_bednet5$n_detect_91_25550 / output_bednet5$n_91_25550, col = cols[6], lwd = 2)
  
  # Add vertical line for bed net distribution (shifted)
  abline(v = bednetstimesteps_shifted, col = "darkgray", lty = "dashed")
  
  # Add label near the vertical line
  text(x = bednetstimesteps_shifted + 5, y = 0.9, "Bed Net Distribution", pos = 4, cex = 0.9)
  
  # Add legend
  legend("topright", legend = c("Control", "IG1", "IG2", "P3", "RG"),
         col = cols[1:6], lty = 1, lwd = 2, box.lty = 0)
}

plot_combined_prevalence()

########################################################################################
################realistic cases averted in %


# -----------------------------
# Colours and net labels
# -----------------------------
cols <- c("#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
net_names <- c("IG1", "IG2", "OP", "RG")

# -----------------------------
# Function: cumulative cases averted
# -----------------------------
cases_averted_ts <- function(out_int, out_ctrl, start_day) {
  idx <- out_ctrl$timestep >= start_day
  cumsum(out_ctrl$n_inc_clinical_730_3649[idx]) -
    cumsum(out_int$n_inc_clinical_730_3649[idx])
}

# -----------------------------
# Time axis (years after distribution)
# -----------------------------
t_days <- output_control$timestep[output_control$timestep >= bednetstimesteps]
t_years <- (t_days - bednetstimesteps) / 365
t_years_plot <- t_years[t_years <= 3]
t_years_plot <- c(0, t_years_plot)  # include 0 at start

# -----------------------------
# Calculate cumulative cases averted
# -----------------------------
cases_list <- lapply(
  list(output_bednet1, output_bednet2, output_bednet3, output_bednet4),
  function(net) {
    tmp <- cases_averted_ts(net, output_control, bednetstimesteps)
    tmp <- tmp[t_years <= 3]
    c(0, tmp)  # include 0 at start
  }
)

# -----------------------------
# Convert to % of total control cumulative cases
# -----------------------------
total_control <- sum(
  output_control$n_inc_clinical_730_3649[
    output_control$timestep >= bednetstimesteps
  ][t_years <= 3]
)

cases_list_percent <- lapply(cases_list, function(x) {
  (x / total_control) * 100
})

# -----------------------------
# Plot limits
# -----------------------------
ymin <- 0
ymax <- 100  # percentage

# -----------------------------
# Plot setup
# -----------------------------
par(mar = c(4, 6, 2, 6) + 0.1)

plot(
  NA, NA,
  xlim = c(0, 3),
  ylim = c(ymin, ymax),
  xlab = "Time in Years after bednet distribution",
  ylab = "Cumulative Malaria Cases Averted in
          \nchildren aged 2–10 years (%)",
  xaxs = "i", yaxs = "i",
  xaxt = "n", yaxt = "n"
)

title("24M Aged Nets")

# -----------------------------
# Axes
# -----------------------------
axis(side = 1, at = seq(0, 3, by = 0.5), tcl = -0.3)
axis(side = 2, at = seq(0, 100, by = 20), las = 1, tcl = -0.3)

# -----------------------------
# Add slim lines and points
# -----------------------------
for (i in 1:4) {
  lines(
    t_years_plot,
    cases_list_percent[[i]],
    col = cols[i],
    lwd = 0.6     # slimmer, smooth lines
  )
  points(
    t_years_plot,
    cases_list_percent[[i]],
    col = cols[i],
    pch = 16,
    cex = 0.4     # smaller points
  )
}

# -----------------------------
# Legend outside plot (matched line width)
# -----------------------------
par(xpd = TRUE)
legend(
  x = 3.05, y = ymax,
  legend = net_names,
  col = cols[1:4],
  lwd = 0.6,
  pch = 16,
  bty = "n",
  xjust = 0,
  yjust = 1
)



########################################################################################
################realistic cases averted in %




# Define colors for each scenario
plot_colors <- c(
  "Control" = "#56B4E9",        # blue
  "Ellie_Pyrrole" = "#009E73",  # green
  "IG2" = "#F0E442",            # yellow
  "Callum_RG" = "#0072B2",      # darker blue
  "RG" = "#D55E00"              # red
)

plot_combined_prevalence <- function() {
  bednetstimesteps <- 365
  bednetstimesteps_shifted <- bednetstimesteps + 25
  
  x_positions <- seq(from = bednetstimesteps, by = 365, length.out = 4)
  x_label_text <- 1:4
  
  # Start the plot with the control line
  plot(
    output_control$timestep,
    output_control$n_detect_lm_730_3649 / output_control$n_age_730_3649,
    type = "l", col = plot_colors["Control"], lwd = 2,
    xlab = "Time in Years",
    ylab = "% of Malaria Prevalence in children aged \n 2-10 years (per 10,000 population)",
    xaxt = "n",
    ylim = c(0, 1),
    xaxs = "i", yaxs = "i",
    xlim = c(bednetstimesteps, bednetstimesteps + 3 * 365),
    main = "Estimated Prevalence"
  )
  
  axis(1, at = x_positions, labels = x_label_text)
  
  # Plot each bednet scenario
  lines(output_bednet1$timestep, output_bednet1$n_detect_lm_730_3649 / output_bednet1$n_age_730_3649,
        col = plot_colors["Ellie_Pyrrole"], lwd = 2)
  lines(output_bednet2$timestep, output_bednet2$n_detect_lm_730_3649 / output_bednet2$n_age_730_3649,
        col = plot_colors["IG2"], lwd = 2)
  lines(output_bednet3$timestep, output_bednet3$n_detect_lm_730_3649 / output_bednet3$n_age_730_3649,
        col = plot_colors["Callum_RG"], lwd = 2)
  lines(output_bednet4$timestep, output_bednet4$n_detect_lm_730_3649 / output_bednet4$n_age_730_3649,
        col = plot_colors["RG"], lwd = 2)
  
  # Bednet distribution line
  abline(v = bednetstimesteps_shifted, col = "darkgray", lty = "dashed")
  text(x = bednetstimesteps_shifted + 5, y = 0.9, "Bed Net Distribution", pos = 4, cex = 0.9)
  
  # Legend
  legend(
    "topright",
    legend = names(plot_colors),
    col = plot_colors,
    lty = 1,
    lwd = 2,
    box.lty = 0
  )
}

plot_combined_prevalence()


######################################################################################################################################################################################################
######################################################################################################################################################################################################
######################################################################################################################################################################################################
#########################prevalence real interceptor G2
######################################################################################################################################################################################################





library(malariasimulation)
library(curl)
library(site)
library(dplyr)
library(ggplot2)

# Simulation parameters
year <- 365
month<-year/12
sim_length <- 4 * year
human_population <- 10000
starting_EIR <- 34

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





# Set drug parameters with seasonality
drug_params <- set_drugs(site_par,  list(AL_params, SP_AQ_params, DHA_PQP_params))



# Set treatment parameters with seasonality
treatment_params1 <- set_clinical_treatment(
  parameters = drug_params,
  drug = 1,
  timesteps =  c(1),
  coverages =  c(0.3)
)

treatment_params2 <- set_clinical_treatment(
  parameters = drug_params,
  drug = 2,
  timesteps =  c(1),
  coverages =  c(0.1)
)
treatment_params<- set_clinical_treatment(
  parameters = treatment_params1,
  drug = 3,
  timesteps =  c(1),
  coverages =  c(0)
)
# Use set_equilibrium to update the parameter set for a given initial EIR
treatment_params <- set_equilibrium(parameters=treatment_params, starting_EIR)

###set bednets

# Set equilibrium
#species_params <- set_equilibrium(parameters = species_params, init_EIR = starting_EIR)

# Run control simulation (without bed nets)
output_control <- run_simulation(timesteps = sim_length, parameters = treatment_params)


bednetstimesteps <- 120

IG2 <- set_bednets(
  treatment_params,
  timesteps = bednetstimesteps,
  coverages = c(.70),  # Each round is distributed to 50% of the population.
  retention = 10^10, # Nets are kept on average 5 years
  dn0 = matrix(c(0.30, 0.30, 0.30),nrow = 1,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.49, 0.49, 0.49),nrow = 1,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24),nrow = 1,ncol = 3,byrow = TRUE),
  gamman = rep(2.65 * 365) # Vector of bed net half-lives for each distribution timestep
)




output_bednet2 <- run_simulation(
  timesteps = sim_length,
  parameters = IG2
)




cols <- c("darkgreen", "aquamarine3", "#E69F00")  # colors for lines

plot_combined_prevalence <- function() {
  # -----------------------------
  # Bed net distribution parameters
  # -----------------------------
  bednetstimesteps <- 1
  bednetstimesteps_shifted <- bednetstimesteps + 120
  
  # x-axis positions and labels
  x_positions <- seq(from = bednetstimesteps, by = 365, length.out = 4)
  x_label_text <- c("2020", "2021", "2022", "2023")
  
  # -----------------------------
  # Prevalence and ±10% uncertainty
  # -----------------------------
  prev_mean <- output_bednet2$n_detect_lm_0_36499 / output_bednet2$n_age_0_36499
  prev_lower <- pmax(0, prev_mean * 0.9)
  prev_upper <- pmin(1, prev_mean * 1.1)
  
  # -----------------------------
  # Plot setup
  # -----------------------------
  plot(output_bednet2$timestep, prev_mean,
       type = "n",  # draw axes first
       xlab = "Time in Years", ylab = "Prevalence",
       xaxt = "n",
       ylim = c(0, 0.7),
       xaxs = "i", yaxs = "i",
       xlim = c(bednetstimesteps, bednetstimesteps + 3 * 365),
       main = "Estimated Prevalence")
  
  # -----------------------------
  # Uncertainty band (shaded)
  # -----------------------------
  polygon(
    c(output_bednet2$timestep, rev(output_bednet2$timestep)),
    c(prev_lower, rev(prev_upper)),
    col = adjustcolor(cols[1], alpha.f = 0.25),
    border = NA
  )
  
  # -----------------------------
  # Mean prevalence line
  # -----------------------------
  lines(output_bednet2$timestep, prev_mean, col = cols[1], lwd = 2)
  
  # -----------------------------
  # Custom x-axis
  # -----------------------------
  axis(1, at = x_positions, labels = x_label_text)
  
  # -----------------------------
  # Bed net distribution line
  # -----------------------------
  abline(v = bednetstimesteps_shifted, col = "darkgray", lty = "dashed")
  text(bednetstimesteps_shifted + 5, 0.65,
       "Bed Net Distribution", pos = 4, cex = 0.9)
  
  # -----------------------------
  # Legend inside plot
  # -----------------------------
  legend(x = 500, y = 0.65,
         legend = c("model prevalance-Interceptor G2", "observed data prevalence"),
         col = c(cols[1], adjustcolor(cols[1], alpha.f = 0.25)),
         lty = c(1, NA),
         lwd = c(2, NA),
         pch = c(NA, 15),
         pt.cex = 2,
         box.lty = 0)
  
  # -----------------------------
  # Observed points
  # -----------------------------
  points(x = c(5, 300, 668),
         y = c(0.40, 0.157, 0.279),
         col = "darkgreen", pch = 15, cex = 1.6)
}

plot_combined_prevalence()






######################################################################################################################################################################################################
#######################let me nowdointerceptor 
#####################################################################################################################################################################################################



library(malariasimulation)
library(curl)
library(site)
library(dplyr)
library(ggplot2)

# Simulation parameters
year <- 365
month<-year/12
sim_length <- 4 * year
human_population <- 10000
starting_EIR <- 66

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





# Set drug parameters with seasonality
drug_params <- set_drugs(site_par,  list(AL_params, SP_AQ_params, DHA_PQP_params))



# Set treatment parameters with seasonality
treatment_params1 <- set_clinical_treatment(
  parameters = drug_params,
  drug = 1,
  timesteps =  c(1),
  coverages =  c(0.3)
)

treatment_params2 <- set_clinical_treatment(
  parameters = drug_params,
  drug = 2,
  timesteps =  c(1),
  coverages =  c(0.1)
)
treatment_params<- set_clinical_treatment(
  parameters = treatment_params1,
  drug = 3,
  timesteps =  c(1),
  coverages =  c(0)
)

treatment_params <- set_equilibrium(parameters=treatment_params, starting_EIR)


output_control <- run_simulation(timesteps = sim_length, parameters = treatment_params)



bednetstimesteps <- 120

IG1U <- set_bednets(
  treatment_params,
  timesteps = bednetstimesteps,
  coverages = c(.70),  # Each round is distributed to 50% of the population.
  retention = 10^10, # Nets are kept on average 5 years
  dn0 = matrix(c(0.11, 0.11, 0.11),nrow = 1,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.58, 0.58, 0.58),nrow = 1,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24),nrow = 1,ncol = 3,byrow = TRUE),
  gamman = rep(1.90 * 365) # Vector of bed net half-lives for each distribution timestep
)




output_bednet1 <- run_simulation(
  timesteps = sim_length,
  parameters = IG1U
)




cols <- c("darkblue"  , "lightblue","darkgreen", "aquamarine3", "#E69F00")  # colors for lines

plot_combined_prevalence1 <- function() {
  # -----------------------------
  # Bed net distribution parameters
  # -----------------------------
  bednetstimesteps <- 1
  bednetstimesteps_shifted <- bednetstimesteps + 120
  
  # x-axis positions and labels
  x_positions <- seq(from = bednetstimesteps, by = 365, length.out = 4)
  x_label_text <- c("2020", "2021", "2022", "2023")
  
  # -----------------------------
  # Prevalence and ±10% uncertainty
  # -----------------------------
  prev_mean <- output_bednet1$n_detect_lm_0_36499 / output_bednet1$n_age_0_36499
  prev_lower <- pmax(0, prev_mean * 0.9)
  prev_upper <- pmin(1, prev_mean * 1.1)
  
  # -----------------------------
  # Plot setup
  # -----------------------------
  plot(output_bednet1$timestep, prev_mean,
       type = "n",  # draw axes first
       xlab = "Time in Years", ylab = "Prevalence",
       xaxt = "n",
       ylim = c(0, 0.7),
       xaxs = "i", yaxs = "i",
       xlim = c(bednetstimesteps, bednetstimesteps + 3 * 365),
       main = "Estimated Prevalence")
  
  # -----------------------------
  # Uncertainty band (shaded)
  # -----------------------------
  polygon(
    c(output_bednet1$timestep, rev(output_bednet1$timestep)),
    c(prev_lower, rev(prev_upper)),
    col = adjustcolor(cols[1], alpha.f = 0.25),
    border = NA
  )
  
  # -----------------------------
  # Mean prevalence line
  # -----------------------------
  lines(output_bednet1$timestep, prev_mean, col = cols[1], lwd = 2)
  
  # -----------------------------
  # Custom x-axis
  # -----------------------------
  axis(1, at = x_positions, labels = x_label_text)
  
  # -----------------------------
  # Bed net distribution line
  # -----------------------------
  abline(v = bednetstimesteps_shifted, col = "darkgray", lty = "dashed")
  text(bednetstimesteps_shifted + 5, 0.65,
       "Bed Net Distribution", pos = 4, cex = 0.9)
  
  # -----------------------------
  # Legend inside plot
  # -----------------------------
  legend(x = 500, y = 0.65,
         legend = c("model prevalance-Interceptor", "observed data prevalence"),
         col = c(cols[1], adjustcolor(cols[1], alpha.f = 0.25)),
         lty = c(1, NA),
         lwd = c(2, NA),
         pch = c(NA, 15),
         pt.cex = 2,
         box.lty = 0)
  
  # -----------------------------
  # Observed points
  # -----------------------------
  points(x = c(4, 300, 668),
         y = c(0.45, 0.28, 0.387),
         col = "darkblue", pch = 19, cex = 1.6)
}

plot_combined_prevalence1()









######################################################################################################################################################################################################
###########################royoal guard
######################################################################################################################################################################################################




library(malariasimulation)
library(curl)
library(site)
library(dplyr)
library(ggplot2)

# Simulation parameters
year <- 365
month<-year/12
sim_length <- 4 * year
human_population <- 10000
starting_EIR <- 54

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





# Set drug parameters with seasonality
drug_params <- set_drugs(site_par,  list(AL_params, SP_AQ_params, DHA_PQP_params))



# Set treatment parameters with seasonality
treatment_params1 <- set_clinical_treatment(
  parameters = drug_params,
  drug = 1,
  timesteps =  c(1),
  coverages =  c(0.3)
)

treatment_params2 <- set_clinical_treatment(
  parameters = drug_params,
  drug = 2,
  timesteps =  c(1),
  coverages =  c(0.1)
)
treatment_params<- set_clinical_treatment(
  parameters = treatment_params1,
  drug = 3,
  timesteps =  c(1),
  coverages =  c(0)
)
# Use set_equilibrium to update the parameter set for a given initial EIR
treatment_params <- set_equilibrium(parameters=treatment_params, starting_EIR)


# Run control simulation (without bed nets)
output_control <- run_simulation(timesteps = sim_length, parameters = treatment_params)


bednetstimesteps <- 120

RG <- set_bednets(
  treatment_params,
  timesteps = bednetstimesteps,
  coverages = c(.70),  # Each round is distributed to 50% of the population.
  retention = 10^10, # Nets are kept on average 5 years
  dn0 = matrix(c(0.29, 0.29, 0.29),nrow = 1,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.49, 0.49, 0.49),nrow = 1,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24),nrow = 1,ncol = 3,byrow = TRUE),
  gamman = rep(2.11 * 365) # Vector of bed net half-lives for each distribution timestep
)



output_bednet3 <- run_simulation(
  timesteps = sim_length,
  parameters = RG
)



cols <- c("#E69F00", "orange", "darkblue"  , "lightblue","darkgreen", "aquamarine3", "#E69F00")  # colors for lines

plot_combined_prevalence2 <- function() {
  # -----------------------------
  # Bed net distribution parameters
  # -----------------------------
  bednetstimesteps <- 1
  bednetstimesteps_shifted <- bednetstimesteps + 120
  
  # x-axis positions and labels
  x_positions <- seq(from = bednetstimesteps, by = 365, length.out = 4)
  x_label_text <- c("2020", "2021", "2022", "2023")
  
  # -----------------------------
  # Prevalence and ±10% uncertainty
  # -----------------------------
  prev_mean <- output_bednet3$n_detect_lm_0_36499 / output_bednet3$n_age_0_36499
  prev_lower <- pmax(0, prev_mean * 0.9)
  prev_upper <- pmin(1, prev_mean * 1.1)
  
  # -----------------------------
  # Plot setup
  # -----------------------------
  plot(output_bednet1$timestep, prev_mean,
       type = "n",  # draw axes first
       xlab = "Time in Years", ylab = "Prevalence",
       xaxt = "n",
       ylim = c(0, 0.7),
       xaxs = "i", yaxs = "i",
       xlim = c(bednetstimesteps, bednetstimesteps + 3 * 365),
       main = "Estimated Prevalence")
  
  # -----------------------------
  # Uncertainty band (shaded)
  # -----------------------------
  polygon(
    c(output_bednet3$timestep, rev(output_bednet1$timestep)),
    c(prev_lower, rev(prev_upper)),
    col = adjustcolor(cols[1], alpha.f = 0.25),
    border = NA
  )
  
  # -----------------------------
  # Mean prevalence line
  # -----------------------------
  lines(output_bednet3$timestep, prev_mean, col = cols[1], lwd = 2)
  
  # -----------------------------
  # Custom x-axis
  # -----------------------------
  axis(1, at = x_positions, labels = x_label_text)
  
  # -----------------------------
  # Bed net distribution line
  # -----------------------------
  abline(v = bednetstimesteps_shifted, col = "darkgray", lty = "dashed")
  text(bednetstimesteps_shifted + 5, 0.65,
       "Bed Net Distribution", pos = 4, cex = 0.9)
  
  # -----------------------------
  # Legend inside plot
  # -----------------------------
  legend(x = 500, y = 0.65,
         legend = c("model prevalance-Interceptor", "observed data prevalence"),
         col = c(cols[1], adjustcolor(cols[1], alpha.f = 0.25)),
         lty = c(1, NA),
         lwd = c(2, NA),
         pch = c(NA, 15),
         pt.cex = 2,
         box.lty = 0)
  
  # -----------------------------
  # Observed points
  # -----------------------------
  points(x = c(4, 300, 668),
         y = c(0.42, 0.269, 0.382),
         col = "#E69F00", pch = 17, cex = 1.6)
}

plot_combined_prevalence2()







######################################################################################################################################################################################################
#############3nw P3
######################################################################################################################################################################################################





library(malariasimulation)
library(curl)
library(site)
library(dplyr)
library(ggplot2)

# Simulation parameters
year <- 365
month<-year/12
sim_length <- 4 * year
human_population <- 10000
starting_EIR <- 54

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



# Set drug parameters with seasonality
drug_params <- set_drugs(site_par,  list(AL_params, SP_AQ_params, DHA_PQP_params))



# Set treatment parameters with seasonality
treatment_params1 <- set_clinical_treatment(
  parameters = drug_params,
  drug = 1,
  timesteps =  c(1),
  coverages =  c(0.3)
)

treatment_params2 <- set_clinical_treatment(
  parameters = drug_params,
  drug = 2,
  timesteps =  c(1),
  coverages =  c(0.1)
)
treatment_params<- set_clinical_treatment(
  parameters = treatment_params1,
  drug = 3,
  timesteps =  c(1),
  coverages =  c(0)
)
# Use set_equilibrium to update the parameter set for a given initial EIR
treatment_params <- set_equilibrium(parameters=treatment_params, starting_EIR)


output_control <- run_simulation(timesteps = sim_length, parameters = treatment_params)



bednetstimesteps <- 120

P3 <- set_bednets(
  treatment_params,
  timesteps = bednetstimesteps,
  coverages = c(.70),  # Each round is distributed to 50% of the population.
  retention = 10^10, # Nets are kept on average 5 years
  dn0 = matrix(c(0.189, 0.189, 0.189),nrow = 1,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.64, 0.64, 0.64),nrow = 1,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24),nrow = 1,ncol = 3,byrow = TRUE),
  gamman = rep(2.16 * 365) # Vector of bed net half-lives for each distribution timestep
)



output_bednet4 <- run_simulation(
  timesteps = sim_length,
  parameters = P3
)




cols <- c("aquamarine3", "darkgreen", "aquamarine3", "#E69F00")  # colors for lines

plot_combined_prevalence4 <- function() {
  # -----------------------------
  # Bed net distribution parameters
  # -----------------------------
  bednetstimesteps <- 1
  bednetstimesteps_shifted <- bednetstimesteps + 120
  
  # x-axis positions and labels
  x_positions <- seq(from = bednetstimesteps, by = 365, length.out = 4)
  x_label_text <- c("2020", "2021", "2022", "2023")
  
  # -----------------------------
  # Prevalence and ±10% uncertainty
  # -----------------------------
  prev_mean <- output_bednet4$n_detect_lm_0_36499 / output_bednet4$n_age_0_36499
  prev_lower <- pmax(0, prev_mean * 0.9)
  prev_upper <- pmin(1, prev_mean * 1.1)
  
  # -----------------------------
  # Plot setup
  # -----------------------------
  plot(output_bednet4$timestep, prev_mean,
       type = "n",  # draw axes first
       xlab = "Time in Years", ylab = "Prevalence",
       xaxt = "n",
       ylim = c(0, 0.7),
       xaxs = "i", yaxs = "i",
       xlim = c(bednetstimesteps, bednetstimesteps + 3 * 365),
       main = "Estimated Prevalence")
  
  # -----------------------------
  # Uncertainty band (shaded)
  # -----------------------------
  polygon(
    c(output_bednet4$timestep, rev(output_bednet4$timestep)),
    c(prev_lower, rev(prev_upper)),
    col = adjustcolor(cols[1], alpha.f = 0.25),
    border = NA
  )
  
  # -----------------------------
  # Mean prevalence line
  # -----------------------------
  lines(output_bednet4$timestep, prev_mean, col = cols[1], lwd = 1)
  
  # -----------------------------
  # Custom x-axis
  # -----------------------------
  axis(1, at = x_positions, labels = x_label_text)
  
  # -----------------------------
  # Bed net distribution line
  # -----------------------------
  abline(v = bednetstimesteps_shifted, col = "darkgray", lty = "dashed")
  text(bednetstimesteps_shifted + 5, 0.65,
       "Bed Net Distribution", pos = 4, cex = 0.9)
  
  # -----------------------------
  # Legend inside plot
  # -----------------------------
  legend(x = 500, y = 0.65,
         legend = c("model prevalance-Interceptor Permanet 3.0" ), #, "observed data prevalence"),
         col = c(cols[1], adjustcolor(cols[1], alpha.f = 0.25)),
         lty = c(1, NA),
         lwd = c(2, NA),
         pch = c(NA, 15),
         pt.cex = 2,
         box.lty = 0)
  
  
}

plot_combined_prevalence4()











######################################################################################################################################################################################################
######################################################################################################################################################################################################
#########################plot together 
######################################################################################################################################################################################################

plot_all_prevalence <- function() {
  # -----------------------------
  # Bed net distribution parameters
  # -----------------------------
  bednetstimesteps <- 1
  bednetstimesteps_shifted <- bednetstimesteps + 120
  x_positions <- seq(from = bednetstimesteps, by = 365, length.out = 4)
  x_label_text <- c("2020", "2021", "2022", "2023")
  
  # -----------------------------
  # Prevalence for each intervention
  # -----------------------------
  prev_IG1U <- output_bednet1$n_detect_lm_0_36499 / output_bednet1$n_age_0_36499
  prev_IG2  <- output_bednet2$n_detect_lm_0_36499 / output_bednet2$n_age_0_36499
  prev_RG   <- output_bednet3$n_detect_lm_0_36499 / output_bednet3$n_age_0_36499
  prev_P3   <- output_bednet4$n_detect_lm_0_36499 / output_bednet4$n_age_0_36499
  
  # Uncertainty bands ±10%
  lower_IG1U <- pmax(0, prev_IG1U * 0.9)
  upper_IG1U <- pmin(1, prev_IG1U * 1.1)
  lower_IG2 <- pmax(0, prev_IG2 * 0.9)
  upper_IG2 <- pmin(1, prev_IG2 * 1.1)
  lower_RG <- pmax(0, prev_RG * 0.9)
  upper_RG <- pmin(1, prev_RG * 1.1)
  lower_P3 <- pmax(0, prev_P3 * 0.9)
  upper_P3 <- pmin(1, prev_P3 * 1.1)
  
  # -----------------------------
  # Colors and symbols
  # -----------------------------
  cols <- c("darkblue", "darkgreen", "#E69F00", "aquamarine3")  # IG1U, IG2, RG, P3
  pchs <- c(15, 15, 15, 15)                                     # symbols for observed points
  
  obs_points <- list(
    IG1U = list(x = c(4, 300, 668), y = c(0.45, 0.28, 0.387)),
    IG2  = list(x = c(5, 300, 668), y = c(0.40, 0.157, 0.279)),
    RG   = list(x = c(4, 300, 668), y = c(0.42, 0.269, 0.382))
    #P3   = list(x = c(4, 300, 668), y = c(0.38, 0.21, 0.32))
  )
  
  # -----------------------------
  # Base plot (no borders)
  # -----------------------------
  plot(output_bednet1$timestep, prev_IG1U,
       type = "n",
       xlab = "Time in Years", ylab = "Malaria Prevalence for all ages (%)",
       xaxt = "n",
       yaxt = "n",
       ylim = c(0, 0.8),
       xaxs = "i", yaxs = "i",
       xlim = c(bednetstimesteps, bednetstimesteps + 3*365),
       main = "Benin Trial Arms",
       bty = "n")  # removes all borders
  
  # Draw left and bottom borders manually
  segments(x0 = bednetstimesteps, y0 = 0, 
           x1 = bednetstimesteps, y1 = 0.7, lwd = 1)
  segments(x0 = bednetstimesteps, y0 = 0, 
           x1 = bednetstimesteps + 3*365, y1 = 0, lwd = 1)
  
  # Custom axes
  axis(1, at = x_positions, labels = x_label_text)
  axis(2, las = 1)
  
  # Bed net distribution line
  abline(v = bednetstimesteps_shifted, col = "darkgray", lty = "dashed")
  text(bednetstimesteps_shifted + 5, 0.65,
       "Bed Net Distribution", pos = 4, cex = 0.9)
  
  # -----------------------------
  # Draw uncertainty bands
  # -----------------------------
  polygon(c(output_bednet1$timestep, rev(output_bednet1$timestep)),
          c(lower_IG1U, rev(upper_IG1U)),
          col = adjustcolor(cols[1], alpha.f = 0.10), border = NA)
  polygon(c(output_bednet2$timestep, rev(output_bednet2$timestep)),
          c(lower_IG2, rev(upper_IG2)),
          col = adjustcolor(cols[2], alpha.f = 0.25), border = NA)
  polygon(c(output_bednet3$timestep, rev(output_bednet3$timestep)),
          c(lower_RG, rev(upper_RG)),
          col = adjustcolor(cols[3], alpha.f = 0.20), border = NA)
  polygon(c(output_bednet4$timestep, rev(output_bednet4$timestep)),
          c(lower_P3, rev(upper_P3)),
          col = adjustcolor(cols[4], alpha.f = 0.25), border = NA)
  
  # -----------------------------
  # Draw mean prevalence lines
  # -----------------------------
  lines(output_bednet1$timestep, prev_IG1U, col = cols[1], lwd = 0.5)
  lines(output_bednet2$timestep, prev_IG2, col = cols[2], lwd = 0.5)
  lines(output_bednet3$timestep, prev_RG, col = cols[3], lwd = 0.5)
  lines(output_bednet4$timestep, prev_P3, col = cols[4], lwd = 0.5)
  
  # -----------------------------
  # Draw observed points
  # -----------------------------
  points(obs_points$IG1U$x, obs_points$IG1U$y, col = cols[1], pch = pchs[1], cex = 1)
  points(obs_points$IG2$x, obs_points$IG2$y, col = cols[2], pch = pchs[2], cex = 1)
  points(obs_points$RG$x, obs_points$RG$y, col = cols[3], pch = pchs[3], cex = 1)
  #points(obs_points$P3$x, obs_points$P3$y, col = cols[4], pch = pchs[4], cex = 1.6)
  
  # -----------------------------
  # Legend inside the plot
  # -----------------------------
  #legend(x = 500, y = 0.65,
  #    legend = c("Interceptor",
  #       "Interceptor G2",
  #     "Royal Guard",
  #      "Permanet 3.0"),
  #col = cols,
  # lty = c(1,1,1,1),
  # lwd = c(1,1,1,1),
  #   pch = pchs,
  #   pt.cex = 1,
  #  box.lty = 0)
  
  legend(x = 500, y = 0.8,
         legend = c("Interceptor","Interceptor G2","Royal Guard","Olyset Plus","Observed prevalence"),
         col = c(cols,"black"),
         lty = c(1,1,1,1,NA),
         lwd = c(1, 1, 1, 1, NA), #,2,2,2,NA),
         pch = c(NA,NA,NA,NA,15),
         pt.cex = 1,
         box.lty = 0)
}

# Plot all four together
plot_all_prevalence()




######################################################################################################################################################################################################
################Only three plots now IG2,IG1 nad RG
######################################################################################################################################################################################################




plot_all_prevalence <- function() {
  # -----------------------------
  # Bed net distribution parameters
  # -----------------------------
  bednetstimesteps <- 1
  bednetstimesteps_shifted <- bednetstimesteps + 120
  x_positions <- seq(from = bednetstimesteps, by = 365, length.out = 4)
  x_label_text <- c("2020", "2021", "2022", "2023")
  
  # -----------------------------
  # Prevalence for each intervention
  # -----------------------------
  prev_IG1U <- output_bednet1$n_detect_lm_0_36499 / output_bednet1$n_age_0_36499
  prev_IG2  <- output_bednet2$n_detect_lm_0_36499 / output_bednet2$n_age_0_36499
  prev_RG   <- output_bednet3$n_detect_lm_0_36499 / output_bednet3$n_age_0_36499
  
  # Uncertainty bands ±10%
  lower_IG1U <- pmax(0, prev_IG1U * 0.9)
  upper_IG1U <- pmin(1, prev_IG1U * 1.1)
  lower_IG2 <- pmax(0, prev_IG2 * 0.9)
  upper_IG2 <- pmin(1, prev_IG2 * 1.1)
  lower_RG <- pmax(0, prev_RG * 0.9)
  upper_RG <- pmin(1, prev_RG * 1.1)
  
  # -----------------------------
  # Colors and symbols
  # -----------------------------
  cols <- c("darkblue", "darkgreen", "#E69F00")  # IG1U, IG2, RG
  pchs <- c(19, 15, 17)                            # IG1U=circle, IG2=box, RG=triangle
  
  obs_points <- list(
    IG1U = list(x = c(4, 300, 668), y = c(0.45, 0.28, 0.387)),
    IG2  = list(x = c(5, 300, 668), y = c(0.40, 0.157, 0.279)),
    RG   = list(x = c(4, 300, 668), y = c(0.42, 0.269, 0.382))
  )
  
  # -----------------------------
  # Base plot
  # -----------------------------
  plot(output_bednet1$timestep, prev_IG1U,
       type = "n",
       xlab = "Time in Years", ylab = "Prevalence",
       xaxt = "n",
       ylim = c(0, 0.7),
       xaxs = "i", yaxs = "i",
       xlim = c(bednetstimesteps, bednetstimesteps + 3*365),
       main = "Estimated Prevalence")
  
  # Custom x-axis
  axis(1, at = x_positions, labels = x_label_text)
  
  # Bed net distribution line
  abline(v = bednetstimesteps_shifted, col = "darkgray", lty = "dashed")
  text(bednetstimesteps_shifted + 5, 0.65,
       "Bed Net Distribution", pos = 4, cex = 0.9)
  
  # -----------------------------
  # Draw uncertainty bands
  # -----------------------------
  polygon(c(output_bednet1$timestep, rev(output_bednet1$timestep)),
          c(lower_IG1U, rev(upper_IG1U)),
          col = adjustcolor(cols[1], alpha.f = 0.25), border = NA)
  polygon(c(output_bednet2$timestep, rev(output_bednet2$timestep)),
          c(lower_IG2, rev(upper_IG2)),
          col = adjustcolor(cols[2], alpha.f = 0.25), border = NA)
  polygon(c(output_bednet3$timestep, rev(output_bednet3$timestep)),
          c(lower_RG, rev(upper_RG)),
          col = adjustcolor(cols[3], alpha.f = 0.25), border = NA)
  
  # -----------------------------
  # Draw mean prevalence lines
  # -----------------------------
  lines(output_bednet1$timestep, prev_IG1U, col = cols[1], lwd = 1)
  lines(output_bednet2$timestep, prev_IG2, col = cols[2], lwd = 1)
  lines(output_bednet3$timestep, prev_RG, col = cols[3], lwd = 1)
  
  # -----------------------------
  # Draw observed points
  # -----------------------------
  points(obs_points$IG1U$x, obs_points$IG1U$y, col = cols[1], pch = pchs[1], cex = 1.6)
  points(obs_points$IG2$x, obs_points$IG2$y, col = cols[2], pch = pchs[2], cex = 1.6)  # now solid boxes
  points(obs_points$RG$x, obs_points$RG$y, col = cols[3], pch = pchs[3], cex = 1.6)
  
  # -----------------------------
  # Legend
  # -----------------------------
  legend("topright",
         legend = c("IG1U - Interceptor", "IG2 - Interceptor G2", "RG - Royal Guard"),
         col = cols,
         lty = 1, lwd = 1,
         pch = pchs,
         pt.cex = 1.5,
         box.lty = 0)
}

plot_all_prevalence()





######################################################################################################################################################################################################
###################model predicted efficacy and real eficacy 
######################################################################################################################################################################################################



# -----------------------------
# Observed data
# -----------------------------
obs_times <- c(4, 300, 668)  # timesteps of observed data

obs_data <- data.frame(
  time = obs_times,
  IG1U_obs = c(0.45, 0.28, 0.387),
  IG2_obs  = c(0.40, 0.157, 0.279),
  RG_obs   = c(0.42, 0.269, 0.382)
)

# -----------------------------
# Extract model-predicted prevalence at observed times
# -----------------------------
pred_IG1U <- approx(
  x = output_bednet1$timestep,
  y = output_bednet1$n_detect_lm_0_36499 / output_bednet1$n_age_0_36499,
  xout = obs_times
)$y

pred_IG2 <- approx(
  x = output_bednet2$timestep,
  y = output_bednet2$n_detect_lm_0_36499 / output_bednet2$n_age_0_36499,
  xout = obs_times
)$y

pred_RG <- approx(
  x = output_bednet3$timestep,
  y = output_bednet3$n_detect_lm_0_36499 / output_bednet3$n_age_0_36499,
  xout = obs_times
)$y

# Combine observed and predicted prevalence
obs_data$IG1U_pred <- pred_IG1U
obs_data$IG2_pred  <- pred_IG2
obs_data$RG_pred   <- pred_RG

print(obs_data)

# -----------------------------
# Plot Observed vs Predicted Prevalence
# -----------------------------
plot(obs_data$IG1U_obs, obs_data$IG1U_pred,
     xlim = c(0, 0.5), ylim = c(0, 0.5),
     xlab = "Observed Prevalence",
     ylab = "Model-predicted Prevalence",
     main = "Observed vs Predicted Prevalence",
     pch = 15, col = "darkblue", cex = 1)

points(obs_data$IG2_obs, obs_data$IG2_pred, pch = 15, col = "darkgreen", cex = 1)
points(obs_data$RG_obs, obs_data$RG_pred, pch = 15, col = "#E69F00", cex = 1)

abline(0,1, lty = 2, col = "gray")  # 1:1 reference line

legend("topleft",
       legend = c("Interceptor","Interceptor G2","Royal Guard"),
       pch = c(15,15,15),
       col = c("darkblue","darkgreen","#E69F00"),
       pt.cex = 1)

# -----------------------------
# Compute correlations and RMSE
# -----------------------------
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2))

cor_IG1U <- cor(obs_data$IG1U_obs, obs_data$IG1U_pred)
cor_IG2  <- cor(obs_data$IG2_obs, obs_data$IG2_pred)
cor_RG   <- cor(obs_data$RG_obs, obs_data$RG_pred)

rmse_IG1U <- rmse(obs_data$IG1U_obs, obs_data$IG1U_pred)
rmse_IG2  <- rmse(obs_data$IG2_obs, obs_data$IG2_pred)
rmse_RG   <- rmse(obs_data$RG_obs, obs_data$RG_pred)

cat("\nCorrelations:\n")
cat("IG1U:", cor_IG1U, " IG2:", cor_IG2, " RG:", cor_RG, "\n")

cat("\nRMSE:\n")
cat("IG1U:", rmse_IG1U, " IG2:", rmse_IG2, " RG:", rmse_RG, "\n")











########################Bothbenin and tanzania
########################################################NOW FOR TANZANIA#################################################################################################################################
###############################################################################################################################################################################################
########################################################################################################################################################################################
###############################################################################################################################################################################################
########################################################################################################################################################################################
###############################################################################################################################################################################################





##################################################################################################################################################################33
#IG2
##################################################################################################################################################################




library(malariasimulation)
library(curl)
library(site)
library(dplyr)
library(ggplot2)

# -----------------------------
# Simulation parameters
# -----------------------------
year <- 365
month <- year/12
sim_length <- 4 * year
human_population <- 10000
starting_EIR <- 14.5




## -------------------------------
simparams <- get_parameters(
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
specie_params <- set_species(
  parameters = simparams,
  species = list(fun_params, gamb_params, arab_params),
  proportions = c(0.1055746 , 0.5411252, 0.3533002 )
)




# Set drug parameters with seasonality
drug_params <- set_drugs(specie_params,  list(AL_params, SP_AQ_params, DHA_PQP_params))



# Set treatment parameters with seasonality
treatment_params1 <- set_clinical_treatment(
  parameters = drug_params,
  drug = 1,
  timesteps =  c(1),
  coverages =  c(0.3)
)

treatment_params2 <- set_clinical_treatment(
  parameters = drug_params,
  drug = 2,
  timesteps =  c(1),
  coverages =  c(0.1)
)
treatment_params<- set_clinical_treatment(
  parameters = treatment_params1,
  drug = 3,
  timesteps =  c(1),
  coverages =  c(0)
)
# Use set_equilibrium to update the parameter set for a given initial EIR
treatment_params <- set_equilibrium(parameters=treatment_params, starting_EIR)


output_control <- run_simulation(timesteps = sim_length, parameters = treatment_params)




bednet_day <- 123  # Jan 1, 2019 relative to Oct 1, 2018 = day 93

IG2 <- set_bednets(
  treatment_params,
  timesteps = bednet_day,
  coverages = c(0.80),
  retention = 10^10,
  dn0 = matrix(c(0.28, 0.28, 0.28), nrow = 1, ncol = 3, byrow = TRUE),
  rn = matrix(c(0.59, 0.59, 0.59), nrow = 1, ncol = 3, byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24), nrow = 1, ncol = 3, byrow = TRUE),
  gamman = rep(2.65 * 365)
)

output_bednet2TZ <- run_simulation(
  timesteps = sim_length,
  parameters = IG2
)

cols <- c("darkgreen", "aquamarine3")  # line colors

# -----------------------------
# Plot function
# -----------------------------
plot_combined_prevalence <- function() {
  
  start_day <- 1  # October 1, 2018
  bednet_distribution <- bednet_day
  
  # x-axis labels: show only 2019-2022
  x_positions <- c(93 + 365*(0:3))  # Jan 2019, 2020, 2021, 2022
  x_label_text <- c("2019", "2020", "2021", "2022")
  
  # Prevalence ±10% uncertainty
  prev_mean <- output_bednet2TZ$n_detect_lm_182.5_5110/ output_bednet2TZ$n_age_182.5_5110
  prev_lower <- pmax(0, prev_mean * 0.9)
  prev_upper <- pmin(1, prev_mean * 1.1)
  
  # -----------------------------
  # Base plot
  # -----------------------------
  plot(output_bednet2TZ$timestep, prev_mean,
       type = "n",
       xlab = "Time in Years", ylab = "Malaria Prevalence for all ages",
       xaxt = "n",
       ylim = c(0, 1),
       xaxs = "i", yaxs = "i",
       xlim = c(start_day, start_day + 4*365),
       bty = "l")  # only bottom and left borders
  
  # Custom x-axis
  axis(1, at = x_positions, labels = x_label_text)
  
  # Bed net distribution line
  abline(v = bednet_distribution, col = "darkgray", lty = "dashed")
  text(bednet_distribution + 10, 0.9, "Bed Net Distribution", pos = 4, cex = 0.9)
  
  # Uncertainty band
  polygon(
    c(output_bednet2TZ$timestep, rev(output_bednet2TZ$timestep)),
    c(prev_lower, rev(prev_upper)),
    col = adjustcolor(cols[1], alpha.f = 0.25),
    border = NA
  )
  
  # Mean prevalence line
  lines(output_bednet2TZ$timestep, prev_mean, col = cols[1], lwd = 2)
  
  # Observed points
  points(x = c(1, 495, 712, 867), 
         y = c(0.462, 0.156, 0.409, 0.256), 
         col = cols[1], pch = 19, cex = 1)
  
  # Legend inside plot
  legend(x = 500, y = 0.9,
         legend = c("Model prevalence - Interceptor G2", "Observed prevalence"),
         col = c(cols[1], cols[1]),
         lty = c(1, NA),
         lwd = c(2, NA),
         pch = c(NA, 19),
         pt.cex = 1,
         box.lty = 0)
}

# -----------------------------
# Plot
# -----------------------------
plot_combined_prevalence()








##########################################################################################################################################################################
#OP
##############################################################################################################################################################################




library(malariasimulation)
library(curl)
library(site)
library(dplyr)
library(ggplot2)

# -----------------------------
# Simulation parameters
# -----------------------------
year <- 365
month <- year/12
sim_length <- 4 * year
human_population <- 10000
starting_EIR <- 11




## -------------------------------
simparams <- get_parameters(
  list(
    human_population = human_population,
    clinical_incidence_rendering_min_ages = 182.5,
    clinical_incidence_rendering_max_ages = 10 * year,
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
specie_params <- set_species(
  parameters = simparams,
  species = list(fun_params, gamb_params, arab_params),
  proportions = c(0.1055746 , 0.5411252, 0.3533002 )
)


# Set drug parameters with seasonality
drug_params <- set_drugs(specie_params,  list(AL_params, SP_AQ_params, DHA_PQP_params))



# Set treatment parameters with seasonality
treatment_params1 <- set_clinical_treatment(
  parameters = drug_params,
  drug = 1,
  timesteps =  c(1),
  coverages =  c(0.3)
)

treatment_params2 <- set_clinical_treatment(
  parameters = drug_params,
  drug = 2,
  timesteps =  c(1),
  coverages =  c(0.1)
)
treatment_params<- set_clinical_treatment(
  parameters = treatment_params1,
  drug = 3,
  timesteps =  c(1),
  coverages =  c(0)
)
# Use set_equilibrium to update the parameter set for a given initial EIR
treatment_params <- set_equilibrium(parameters=treatment_params, starting_EIR)


output_control <- run_simulation(timesteps = sim_length, parameters = treatment_params)



# Use set_equilibrium to update the parameter set for a given initial EIR
bednet_day <- 123  # Jan 1, 2019 relative to Oct 1, 2018 = day 93

OP <- set_bednets(
  treatment_params,
  timesteps = bednet_day,
  coverages = c(0.48),
  retention = 10^10,
  dn0 = matrix(c(0.27, 0.27, 0.27), nrow = 1, ncol = 3, byrow = TRUE),
  rn = matrix(c(0.58, 0.58, 0.58), nrow = 1, ncol = 3, byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24), nrow = 1, ncol = 3, byrow = TRUE),
  gamman = rep(2.11* 365)
)

output_bednet4TZ <- run_simulation(
  timesteps = sim_length,
  parameters = OP
)

cols <- c("darkred", "aquamarine3","darkblue", "darkgreen", "aquamarine3")  # line colors

# -----------------------------
# Plot function
# -----------------------------
plot_combined_prevalence2 <- function() {
  
  start_day <- 1  # October 1, 2018
  bednet_distribution <- bednet_day
  
  # x-axis labels: show only 2019-2022
  x_positions <- c(93 + 365*(0:3))  # Jan 2019, 2020, 2021, 2022
  x_label_text <- c("2019", "2020", "2021", "2022")
  
  # Prevalence ±10% uncertainty
  prev_mean <- output_bednet4TZ$n_detect_lm_182.5_5110 / output_bednet4TZ$n_age_182.5_5110
  prev_lower <- pmax(0, prev_mean * 0.9)
  prev_upper <- pmin(1, prev_mean * 1.1)
  
  # -----------------------------
  # Base plot
  # -----------------------------
  plot(output_bednet4TZ$timestep, prev_mean,
       type = "n",
       xlab = "Time in Years", ylab = "Malaria Prevalence for all ages",
       xaxt = "n",
       ylim = c(0, 1),
       xaxs = "i", yaxs = "i",
       xlim = c(start_day, start_day + 4*365),
       bty = "l")  # only bottom and left borders
  
  # Custom x-axis
  axis(1, at = x_positions, labels = x_label_text)
  
  # Bed net distribution line
  abline(v = bednet_distribution, col = "darkgray", lty = "dashed")
  text(bednet_distribution + 10, 0.9, "Bed Net Distribution", pos = 4, cex = 0.9)
  
  # Uncertainty band
  polygon(
    c(output_bednet4TZ$timestep, rev(output_bednet4TZ$timestep)),
    c(prev_lower, rev(prev_upper)),
    col = adjustcolor(cols[1], alpha.f = 0.25),
    border = NA
  )
  
  # Mean prevalence line
  lines(output_bednet4TZ$timestep, prev_mean, col = cols[1], lwd = 1)
  
  # Observed points
  points(x = c(1, 495, 712, 867), 
         y = c(0.42, 0.192, 0.433, 0.407), 
         col = cols[1], pch = 19, cex = 1)
  
  # Legend inside plot
  legend(x = 500, y = 0.9,
         legend = c(" Olyeset Plus", "Observed prevalence"),
         col = c(cols[1], cols[1]),
         lty = c(1, NA),
         lwd = c(2, NA),
         pch = c(NA, 19),
         pt.cex = 1,
         box.lty = 0)
}

# -----------------------------
# Plot
# -----------------------------
plot_combined_prevalence2()


#######################################################################################################################royal guard#################################################################
#RG
##############################################################################################################################################################################################


# -----------------------------
year <- 365
month <- year/12
sim_length <- 4 * year
human_population <- 10000
starting_EIR <- 11.2




## -------------------------------
simparams <- get_parameters(
  list(
    human_population = human_population,
    clinical_incidence_rendering_min_ages = 182.5,
    clinical_incidence_rendering_max_ages = 10 * year,
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
specie_params <- set_species(
  parameters = simparams,
  species = list(fun_params, gamb_params, arab_params),
  proportions = c(0.1055746 , 0.5411252, 0.3533002 )
)

# Set drug parameters with seasonality
drug_params <- set_drugs(specie_params,  list(AL_params, SP_AQ_params, DHA_PQP_params))



# Set treatment parameters with seasonality
treatment_params1 <- set_clinical_treatment(
  parameters = drug_params,
  drug = 1,
  timesteps =  c(1),
  coverages =  c(0.3)
)

treatment_params2 <- set_clinical_treatment(
  parameters = drug_params,
  drug = 2,
  timesteps =  c(1),
  coverages =  c(0.1)
)
treatment_params<- set_clinical_treatment(
  parameters = treatment_params1,
  drug = 3,
  timesteps =  c(1),
  coverages =  c(0)
)
# Use set_equilibrium to update the parameter set for a given initial EIR
treatment_params <- set_equilibrium(parameters=treatment_params, starting_EIR)
# Use set_equilibrium to update the parameter set for a given initial EIR
treatment_params <- set_equilibrium(parameters=specie_params, starting_EIR)


# -----------------------------
# Simulation parameters

# Bed net distribution
# -----------------------------
bednet_day <- 123  # Jan 1, 2019 relative to Oct 1, 2018 = day 93

RG <- set_bednets(
  treatment_params,
  timesteps = bednet_day,
  coverages = c(0.58),
  retention = 10^10,
  dn0 = matrix(c(0.25, 0.25, 0.25), nrow = 1, ncol = 3, byrow = TRUE),
  rn = matrix(c(0.58, 0.58, 0.58), nrow = 1, ncol = 3, byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24), nrow = 1, ncol = 3, byrow = TRUE),
  gamman = rep(2.11 * 365)
)

output_bednet5TZ <- run_simulation(
  timesteps = sim_length,
  parameters = RG
)

cols <- c( "darkorange","darkblue", "darkgreen", "aquamarine3")  # line colors

# -----------------------------
# Plot function
# -----------------------------
plot_combined_prevalence3 <- function() {
  
  start_day <- 1  # October 1, 2018
  bednet_distribution <- bednet_day
  
  # x-axis labels: show only 2019-2022
  x_positions <- c(93 + 365*(0:3))  # Jan 2019, 2020, 2021, 2022
  x_label_text <- c("2019", "2020", "2021", "2022")
  
  # Prevalence ±10% uncertainty
  prev_mean <- output_bednet5TZ$n_detect_lm_182.5_5110 / output_bednet5TZ$n_age_182.5_5110
  prev_lower <- pmax(0, prev_mean * 0.9)
  prev_upper <- pmin(1, prev_mean * 1.1)
  
  # -----------------------------
  # Base plot
  # -----------------------------
  plot(output_bednet5TZ$timestep, prev_mean,
       type = "n",
       xlab = "Time in Years", ylab = "Malaria Prevalence for all ages",
       xaxt = "n",
       ylim = c(0, 1),
       xaxs = "i", yaxs = "i",
       xlim = c(start_day, start_day + 4*365),
       bty = "l")  # only bottom and left borders
  
  # Custom x-axis
  axis(1, at = x_positions, labels = x_label_text)
  
  # Bed net distribution line
  abline(v = bednet_distribution, col = "darkgray", lty = "dashed")
  text(bednet_distribution + 10, 0.9, "Bed Net Distribution", pos = 4, cex = 0.9)
  
  # Uncertainty band
  polygon(
    c(output_bednet5TZ$timestep, rev(output_bednet5TZ$timestep)),
    c(prev_lower, rev(prev_upper)),
    col = adjustcolor(cols[1], alpha.f = 0.25),
    border = NA
  )
  
  # Mean prevalence line
  lines(output_bednet5TZ$timestep, prev_mean, col = cols[1], lwd = 2)
  
  # Observed points
  points(x = c(1, 495, 712, 867), 
         y = c(0.427, 0.156, 0.409, 0.256), 
         col = cols[1], pch = 19, cex = 1)
  
  # Legend inside plot
  legend(x = 500, y = 0.9,
         legend = c(" Royal guard", "Observed prevalence"),
         col = c(cols[1], cols[1]),
         lty = c(1, NA),
         lwd = c(2, NA),
         pch = c(NA, 19),
         pt.cex = 1,
         box.lty = 0)
}

# -----------------------------
# Plot
# -----------------------------
plot_combined_prevalence3()




##################################################################################################################################################################33
#IG1
##################################################################################################################################################################




library(malariasimulation)
library(curl)
library(site)
library(dplyr)
library(ggplot2)


year <- 365
month <- year/12
sim_length <- 4 * year
human_population <- 10000
starting_EIR <- 13.5



## -------------------------------
simparams <- get_parameters(
  list(
    human_population = human_population,
    clinical_incidence_rendering_min_ages = 182.5,
    clinical_incidence_rendering_max_ages = 10 * year,
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
specie_params <- set_species(
  parameters = simparams,
  species = list(fun_params, gamb_params, arab_params),
  proportions = c(0.1055746 , 0.5411252, 0.3533002 )
)



# Set drug parameters with seasonality
drug_params <- set_drugs(specie_params,  list(AL_params, SP_AQ_params, DHA_PQP_params))



# Set treatment parameters with seasonality
treatment_params1 <- set_clinical_treatment(
  parameters = drug_params,
  drug = 1,
  timesteps =  c(1),
  coverages =  c(0.3)
)

treatment_params2 <- set_clinical_treatment(
  parameters = drug_params,
  drug = 2,
  timesteps =  c(1),
  coverages =  c(0.1)
)
treatment_params<- set_clinical_treatment(
  parameters = treatment_params1,
  drug = 3,
  timesteps =  c(1),
  coverages =  c(0)
)
# Use set_equilibrium to update the parameter set for a given initial EIR
treatment_params <- set_equilibrium(parameters=treatment_params, starting_EIR)
# Use set_equilibrium to update the parameter set for a given initial EIR
treatment_params <- set_equilibrium(parameters=specie_params, starting_EIR)


# -----------------------------
# Bed net distribution
# -----------------------------
bednet_day <- 123  # Jan 1, 2019 relative to Oct 1, 2018 = day 93

IG1 <- set_bednets(
  treatment_params,
  timesteps = bednet_day,
  coverages = c(0.68),
  retention = 10^10,
  dn0 = matrix(c(0.13, 0.13, 0.13), nrow = 1, ncol = 3, byrow = TRUE),
  rn = matrix(c(0.72, 0.72, 0.72), nrow = 1, ncol = 3, byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24), nrow = 1, ncol = 3, byrow = TRUE),
  gamman = rep(1.90 * 365)
)

output_bednet3TZ <- run_simulation(
  timesteps = sim_length,
  parameters = IG1
)

cols <- c("darkblue", "darkgreen", "aquamarine3")  # line colors

# -----------------------------
# Plot function
# -----------------------------
plot_combined_prevalence1 <- function() {
  
  start_day <- 1  # October 1, 2018
  bednet_distribution <- bednet_day
  
  # x-axis labels: show only 2019-2022
  x_positions <- c(93 + 365*(0:3))  # Jan 2019, 2020, 2021, 2022
  x_label_text <- c("2019", "2020", "2021", "2022")
  
  # Prevalence ±10% uncertainty
  prev_mean <- output_bednet3TZ$n_detect_lm_182.5_5110 / output_bednet3TZ$n_age_182.5_5110
  prev_lower <- pmax(0, prev_mean * 0.9)
  prev_upper <- pmin(1, prev_mean * 1.1)
  
  # -----------------------------
  # Base plot
  # -----------------------------
  plot(output_bednet3TZ$timestep, prev_mean,
       type = "n",
       xlab = "Time in Years", ylab = "Malaria Prevalence for children aged 6 months -14 years (%)",
       xaxt = "n",
       ylim = c(0, 0.7),
       xaxs = "i", yaxs = "i",
       xlim = c(start_day, start_day + 4*365),
       main = "Tanzania Trials",
       bty = "l")  # only bottom and left borders
  
  # Custom x-axis
  axis(1, at = x_positions, labels = x_label_text)
  
  # Bed net distribution line
  abline(v = bednet_distribution, col = "darkgray", lty = "dashed")
  text(bednet_distribution + 10, 0.9, "Bed Net Distribution", pos = 4, cex = 0.9)
  
  # Uncertainty band
  polygon(
    c(output_bednet3TZ$timestep, rev(output_bednet3TZ$timestep)),
    c(prev_lower, rev(prev_upper)),
    col = adjustcolor(cols[1], alpha.f = 0.25),
    border = NA
  )
  
  # Mean prevalence line
  lines(output_bednet3TZ$timestep, prev_mean, col = cols[1], lwd = 2)
  
  # Observed points
  points(x = c(1, 495, 712, 867), 
         y = c(0.459, 0.312, 0.523, 0.458), 
         col = cols[1], pch = 19, cex = 1.6)
  
  # Legend inside plot
  legend(x = 500, y = 0.9,
         legend = c(" Interceptor", "Observed prevalence"),
         col = c(cols[1], cols[1]),
         lty = c(1, NA),
         lwd = c(2, NA),
         pch = c(NA, 19),
         pt.cex = 2,
         box.lty = 0)
}

# -----------------------------
# Plot
# -----------------------------
plot_combined_prevalence1()




















###########################################combined plot#######################################################################################################################
#############################################################################################################################################################################################



# -----------------------------
# Combined plot for 4 simulations
# -----------------------------
plot_combined_prevalence_all <- function() {
  
  start_day <- 1            # Simulation start: Oct 1, 2018
  bednet_distribution <- 123  # 31 Jan 2019 relative to start
  
  # X-axis labels for years (approx Jan each year)
  x_positions <- 93 + 365*(0:3)  # Jan 2019, 2020, 2021, 2022
  x_label_text <- c("2019", "2020", "2021", "2022")
  
  # Colors for each simulation
  sim_cols <- c("darkblue", "darkgreen", "darkorange", "darkred")
  
  # List of simulation outputs
  sims <- list(output_bednet3TZ, output_bednet2TZ, output_bednet5TZ, output_bednet4TZ)
  
  # Observed points for each simulation
  obs_points <- list(
    c(0.459, 0.312, 0.523, 0.458),  # IG1
    c(0.462, 0.156, 0.409, 0.256),  # IG2
    c(0.427, 0.217, 0.506, 0.375),  # RG
    c(0.42, 0.192, 0.433, 0.407)    # OP
  )
  
  par(mar = c(5, 6, 4, 2) + 0.1)  # bottom, left, top, right
  
  # -----------------------------
  # Base plot (empty)
  # -----------------------------
  plot(sims[[1]]$timestep, sims[[1]]$n_detect_lm_182.5_5110 / sims[[1]]$n_age_182.5_5110,
       type = "n",
       xlab = "Time in Years", ylab = "Malaria Prevalence for  children \n aged 6 months- 14 years (%)",
       xaxt = "n",
       ylim = c(0, 0.8),
       xaxs = "i", yaxs = "i",
       xlim = c(start_day, start_day + 4*365),
       main="Tanzania Trial",
       bty = "l")
  
  # Custom x-axis
  axis(1, at = x_positions, labels = x_label_text)
  
  # Bed net distribution line
  abline(v = bednet_distribution, col = "darkgray", lty = "dashed")
  text(bednet_distribution + 10, 0.75, "Bed Net Distribution", pos = 4, cex = 0.75)
  
  # -----------------------------
  # Overlay all simulations
  # -----------------------------
  for(i in 1:4) {
    sim <- sims[[i]]
    
    # Prevalence ±10% uncertainty
    prev_mean <- sim$n_detect_lm_182.5_5110 / sim$n_age_182.5_5110
    prev_lower <- pmax(0, prev_mean * 0.9)
    prev_upper <- pmin(1, prev_mean * 1.1)
    
    # Uncertainty polygon
    polygon(
      c(sim$timestep, rev(sim$timestep)),
      c(prev_lower, rev(prev_upper)),
      col = adjustcolor(sim_cols[i], alpha.f = 0.25),
      border = NA
    )
    
    # Mean prevalence line
    lines(sim$timestep, prev_mean, col = sim_cols[i], lwd = 1)
    
    # Observed points
    points(x = c(1, 495, 712, 867),
           y = obs_points[[i]],
           col = sim_cols[i], pch = 19, cex = 1)
  }
  
  # -----------------------------
  # Legend
  # -----------------------------
  legend(x = 700, y = 0.8,
         legend = c("Interceptor","Interceptor G2","Royal Guard","Olyset Plus","Observed prevalence"),
         col = c(sim_cols,"black"),
         lty = c(1,1,1,1,NA),
         lwd = c(1, 1, 1, 1, NA), #,2,2,2,NA),
         pch = c(NA,NA,NA,NA,19),
         pt.cex = 1,
         box.lty = 0)
}

# -----------------------------
# Plot
# -----------------------------
plot_combined_prevalence_all()





#########################now lets plot the predicted vs model







# -----------------------------
# Observed time points
# -----------------------------
obs_days <- c(1, 495, 712, 867)

# -----------------------------
# Function to extract predicted prevalence
# -----------------------------
get_predicted_prev <- function(output, obs_days) {
  prev <- output$n_detect_lm_182.5_5110 / output$n_age_182.5_5110
  prev[match(obs_days, output$timestep)]
}

# -----------------------------
# Build combined dataset
# -----------------------------
pred_obs <- data.frame(
  Treatment = rep(c("IG1", "IG2", "RG", "OP"), each = length(obs_days)),
  
  Observed = c(
    0.459, 0.312, 0.523, 0.458,   # IG1
    0.462, 0.156, 0.409, 0.256,   # IG2
    0.427, 0.156, 0.409, 0.256,   # RG
    0.420, 0.192, 0.433, 0.407    # OP
  ),
  
  Predicted = c(
    get_predicted_prev(output_bednet3TZ, obs_days),
    get_predicted_prev(output_bednet2TZ, obs_days),
    get_predicted_prev(output_bednet5TZ, obs_days),
    get_predicted_prev(output_bednet4TZ, obs_days)
  )
)

# Colors
cols <- c("IG1" = "darkblue",
          "IG2" = "darkgreen",
          "RG"  = "darkorange",
          "OP"  = "darkred")

# -----------------------------
# Combined Predicted vs Observed plot
# -----------------------------
plot(pred_obs$Observed, pred_obs$Predicted,
     col = cols[pred_obs$Treatment],
     pch = 16,
     xlab = "Observed prevalence",
     ylab = "Predicted prevalence",
     xlim = c(0, 1),
     ylim = c(0, 1),
     bty = "l")

# 1:1 reference line
abline(0, 1, lty = 2, col = "gray40")

# Legend
legend("topleft",
       legend = c("IG1 Interceptor", "IG2 Interceptor G2",
                  "Royal Guard", "Olyset Plus"),
       col = cols,
       pch = 16,
       bty = "n")






##################LETS PLOT POINTS FR IG2 and the other one 








# -----------------------------
# Combine Benin and Tanzania data
# -----------------------------

# Benin
benin_obs <- c(obs_data$IG1U_obs, obs_data$IG2_obs, obs_data$RG_obs)
benin_pred <- c(obs_data$IG1U_pred, obs_data$IG2_pred, obs_data$RG_pred)
benin_treat <- rep(c("IG1","IG2","RG"), each = length(obs_times))
benin_country <- rep("Benin", length(benin_obs))

# Tanzania
tanz_obs <- pred_obs$Observed
tanz_pred <- pred_obs$Predicted
tanz_treat <- pred_obs$Treatment
tanz_country <- rep("Tanzania", length(tanz_obs))

# Combine into one data frame
combined <- data.frame(
  Observed = c(benin_obs, tanz_obs),
  Predicted = c(benin_pred, tanz_pred),
  Treatment = c(benin_treat, tanz_treat),
  Country = c(benin_country, tanz_country)
)

# -----------------------------
# Colors and shapes
# -----------------------------
cols <- c("IG1U" = "darkblue", "IG2" = "darkgreen", "RG" = "darkorange",
          "IG1" = "darkblue", "IG2" = "darkgreen", "RG" = "darkorange", "OP" = "darkred")

shapes <- c("Benin" = 15, "Tanzania" = 16)  # squares vs circles

# -----------------------------
# Plot
# -----------------------------
plot(combined$Observed, combined$Predicted,
     col = cols[combined$Treatment],
     pch = shapes[combined$Country],
     xlab = "Observed prevalence",
     ylab = "Predicted prevalence",
     xlim = c(0,1), ylim = c(0,0.7),
     main = "Observed vs Predicted Prevalence (Benin & Tanzania)")

# 1:1 reference line
abline(0,1,lty=2,col="gray40")

# -----------------------------
# Legend
# -----------------------------
legend("topleft",
       legend = c("Benin","Tanzania",
                  "Interceptor","Interceptor G2","Royal Guard","Olyset Plus"),
       col = c("black","black",cols[c("IG1U","IG2","RG","OP")]),
       pch = c(15,16,NA,NA,NA,NA),
       pt.cex = c(1,1,1,1,1,1),
       lty = c(NA,NA,1,1,1,1),
       bty = "n")


legend("topleft",
       legend = c("Benin", "Tanzania",
                  "Interceptor", "Interceptor G2",
                  "Royal Guard", "Olyset Plus"),
       col = c("black", "black",
               cols[c("IG1U","IG2","RG","OP")]),
       pch = c(15, 16, 16, 16, 16, 16),  # squares & circles only
       pt.cex = 1,
       bty = "n")






















################################################My EIR calibrartion models ##################################################################################################################################
#########################################################################################################################################################################################################################







library(malariasimulation)
library(curl)
library(site)
rm(list = ls())
set.seed(123)



library(cali)
library(malariasimulation)
library(ggplot2)
library(dplyr)
library(ggplot2)
parameters$human_population <- 5000
parameters$timesteps <- 365 * 4  # total simulation length
human_population<- 5000
sim_length<- 4 *365

# Simulation parameters
#year <- 365
#month<-year/12
#sim_length <- 4 * year
#starting_EIR <- 34


target_days <- c( 1)  # mid-year of years 1, 2, 3
target_pfpr <- c( 0.4)  # desired PfPR at these days


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





# Set drug parameters with seasonality
drug_params <- set_drugs(site_par,  list(AL_params, SP_AQ_params, DHA_PQP_params))



# Set treatment parameters with seasonality
treatment_params1 <- set_clinical_treatment(
  parameters = drug_params,
  drug = 1,
  timesteps =  c(1),
  coverages =  c(0.3)
)

treatment_params2 <- set_clinical_treatment(
  parameters = drug_params,
  drug = 2,
  timesteps =  c(1),
  coverages =  c(0.1)
)
treatment_params<- set_clinical_treatment(
  parameters = treatment_params1,
  drug = 3,
  timesteps =  c(1),
  coverages =  c(0)
)
# Use set_equilibrium to update the parameter set for a given initial EIR
#treatment_params <- set_equilibrium(parameters=treatment_params, starting_EIR)

###set bednets

# Set equilibrium
#species_params <- set_equilibrium(parameters = species_params, init_EIR = starting_EIR)

# Run control simulation (without bed nets)
#output_control <- run_simulation(timesteps = sim_length, parameters = treatment_params)


bednetstimesteps <- 120

parameters <- get_parameters()

parameters <- set_bednets(
  treatment_params,
  timesteps = bednetstimesteps,
  coverages = c(.70),  # Each round is distributed to 50% of the population.
  retention = 10^10, # Nets are kept on average 5 years
  dn0 = matrix(c(0.30, 0.30, 0.30),nrow = 1,ncol = 3,byrow = TRUE),
  rn = matrix(c(0.49, 0.49, 0.49),nrow = 1,ncol = 3,byrow = TRUE),
  rnm = matrix(c(0.24, 0.24, 0.24),nrow = 1,ncol = 3,byrow = TRUE),
  gamman = rep(2.65 * 365) # Vector of bed net half-lives for each distribution timestep
)




parameters$timesteps <- 365 * 4  # total simulation length

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

parameters <- set_equilibrium(
  parameters,
  init_EIR = out
)

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







parameters$timesteps <- 365 * 4  # total simulation length


output_bednet2 <- run_simulation(
  timesteps = parameters$timesteps,
  parameters = parameters
)

sim_length<- 4*365


cols <- c("darkgreen", "aquamarine3", "#E69F00")  # colors for lines

plot_combined_prevalence <- function() {
  # -----------------------------
  # Bed net distribution parameters
  # -----------------------------
  bednetstimesteps <- 1
  bednetstimesteps_shifted <- bednetstimesteps + 120
  
  # x-axis positions and labels
  x_positions <- seq(from = bednetstimesteps, by = 365, length.out = 4)
  x_label_text <- c("2020", "2021", "2022", "2023")
  
  # -----------------------------
  # Prevalence and ±10% uncertainty
  # -----------------------------
  prev_mean <- output_bednet2$n_detect_lm_0_36499 / output_bednet2$n_age_0_36499
  prev_lower <- pmax(0, prev_mean * 0.9)
  prev_upper <- pmin(1, prev_mean * 1.1)
  
  # -----------------------------
  # Plot setup
  # -----------------------------
  plot(output_bednet2$timestep, prev_mean,
       type = "n",  # draw axes first
       xlab = "Time in Years", ylab = "Prevalence",
       xaxt = "n",
       ylim = c(0, 0.7),
       xaxs = "i", yaxs = "i",
       xlim = c(bednetstimesteps, bednetstimesteps + 3 * 365),
       main = "Estimated Prevalence")
  
  # -----------------------------
  # Uncertainty band (shaded)
  # -----------------------------
  polygon(
    c(output_bednet2$timestep, rev(output_bednet2$timestep)),
    c(prev_lower, rev(prev_upper)),
    col = adjustcolor(cols[1], alpha.f = 0.25),
    border = NA
  )
  
  # -----------------------------
  # Mean prevalence line
  # -----------------------------
  lines(output_bednet2$timestep, prev_mean, col = cols[1], lwd = 2)
  
  # -----------------------------
  # Custom x-axis
  # -----------------------------
  axis(1, at = x_positions, labels = x_label_text)
  
  # -----------------------------
  # Bed net distribution line
  # -----------------------------
  abline(v = bednetstimesteps_shifted, col = "darkgray", lty = "dashed")
  text(bednetstimesteps_shifted + 5, 0.65,
       "Bed Net Distribution", pos = 4, cex = 0.9)
  
  # -----------------------------
  # Legend inside plot
  # -----------------------------
  legend(x = 500, y = 0.65,
         legend = c("model prevalance-Interceptor G2", "observed data prevalence"),
         col = c(cols[1], adjustcolor(cols[1], alpha.f = 0.25)),
         lty = c(1, NA),
         lwd = c(2, NA),
         pch = c(NA, 15),
         pt.cex = 2,
         box.lty = 0)
  
  # -----------------------------
  # Observed points
  # -----------------------------
  points(x = c(5, 300, 668),
         y = c(0.40, 0.157, 0.279),
         col = "darkgreen", pch = 15, cex = 1.6)
}

plot_combined_prevalence()









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
  obs_x <- c(2190, 2496, 2859)
  obs_y <- c(0.407, 0.157, 0.279)
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
  obs_x <- c(2190, 2496, 2859)
  obs_y <- c(0.465, 0.28, 0.387)
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
  obs_x <- c(2190, 2496, 2859)
  obs_y <- c(0.431, 0.269, 0.382)
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




































##########forkeeps 





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
  obs_y <- c(0.40, 0.157, 0.279, 0.255)    # prevalence
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




#######################################################################################################################################################################
##############################################################################################################################################################################################################################################################################################################################################
#######################################################################################################################################################################
#################now for IG1





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
  obs_x <- c(2190, 2496, 2859)
  obs_y <- c(0.465, 0.28, 0.387)
  obs_keep <- obs_x >= slice_start
  if(any(obs_keep)){
    points(obs_x[obs_keep], obs_y[obs_keep], col = "darkblue", pch = 15, cex = 1)
  }
}

# Run the sliced plot
plot_combined_prevalence2()







