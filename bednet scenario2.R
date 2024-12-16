library(malariasimulation)

# Simulation parameters
year <- 365
sim_length <- 5 * year
human_population <- 8000
starting_EIR <- 15


# Initialize simulation parameters
simparams <- get_parameters(
  list(
    human_population = human_population, 
    clinical_incidence_rendering_min_ages = 91,
    clinical_incidence_rendering_max_ages = 70 * 365, 
    prevalence_rendering_min_ages = 91, 
    prevalence_rendering_max_ages = 70 * 365
  )
)

# Set species parameters
species_params <- set_species(
  parameters = simparams,
  species = list(fun_params, gamb_params),
  proportions = c(0.408, 0.592)
)

# Set equilibrium
species_params <- set_equilibrium(parameters = species_params, init_EIR = starting_EIR)

# Run control simulation (without bed nets)
output_control <- run_simulation(timesteps = sim_length, parameters = species_params)

# Define bednet parameters for five different scenarios
IG1U <- set_bednets(
  parameters = species_params,
  timesteps = c(1, 3) * year,
  coverages = c(0.9, 0.9),
  retention = 2 * 365,
  dn0 = matrix(cbind(c(0.354621849, 0.354621849), c(0.155488951, 0.155488951)), nrow = 2, ncol = 2),
  rn = matrix(cbind(c(0.565286611, 0.565286611), c(0.526669837, 0.526669837)), nrow = 2, ncol = 2),
  rnm = matrix(cbind(c(0.22, 0.22), c(0.22, 0.22)), nrow = 2, ncol = 2),
  gamman = rep(2.895257* 365, 2)
)

IG2U <- set_bednets(
  parameters = species_params,
  timesteps = c(1, 3) * year,
  coverages = c(0.9, 0.9),
  retention = 2 * 365,
  dn0 = matrix(cbind(c(0.37838552, 0.37838552), c(0.165655393, 0.165655393)), nrow = 2, ncol = 2),
  rn = matrix(cbind(c(0.541522693	, 0.541522693	), c(0.516502589, 0.516502589)), nrow = 2, ncol = 2),
  rnm = matrix(cbind(c(0.22, 0.22), c(0.22, 0.22)), nrow = 2, ncol = 2),
  gamman = rep(2.920989* 365, 2)
)


OPU <- set_bednets(
  parameters = species_params,
  timesteps = c(1, 3) * year,
  coverages = c(0.9, 0.9),
  retention = 2 * 365,
  dn0 = matrix(cbind(c(0.465751648, 0.465751648), c(0.210416336, 0.210416336)), nrow = 2, ncol = 2),
  rn = matrix(cbind(c(0.450628788, 0.450628788), c(0.4602951, 0.4602951)), nrow = 2, ncol = 2),
  rnm = matrix(cbind(c(0.22, 0.22), c(0.22, 0.22)), nrow = 2, ncol = 2),
  gamman = rep(3.194394 * 365, 2)
)



P2U <- set_bednets(
  parameters = species_params,
  timesteps = c(1, 3) * year,
  coverages = c(0.9, 0.9),
  retention = 2 * 365,
  dn0 = matrix(cbind(c(0.203309951, 0.203309951), c(	0.164190551, 	0.164190551)), nrow = 2, ncol = 2),
  rn = matrix(cbind(c(.735943401, .735943401), c(.587085102, .587085102)), nrow = 2, ncol = 2),
  rnm = matrix(cbind(c(0.22, 0.22), c(0.22, 0.22)), nrow = 2, ncol = 2),
  #gamman = matrix(rep(2.65196 * 365, 2), nrow = 2, ncol = 1) # 2x1 matrix for gamman
  gamman = rep(2.938527 * 365, 2)
)


RGU <- set_bednets(
  parameters = species_params,
  timesteps = c(1, 3) * year,
  coverages = c(0.9, 0.9),
  retention = 2 * 365,
  dn0 = matrix(cbind(c(0.169164837, 0.169164837), c(0.132083126, 0.132083126)), nrow = 2, ncol = 2),
  rn = matrix(cbind(c(0.754104889, 0.754104889), c(0.565991711, 0.565991711)), nrow = 2, ncol = 2),
  rnm = matrix(cbind(c(0.22, 0.22), c(0.22, 0.22)), nrow = 2, ncol = 2),
  gamman = rep(2.874849  * 365, 2)
)

# Run simulations for each bednet parameter
output_bednet1 <- run_simulation(timesteps = sim_length, parameters = IG1U)
output_bednet2 <- run_simulation(timesteps = sim_length, parameters = IG2U)
output_bednet3 <- run_simulation(timesteps = sim_length, parameters = OPU)
output_bednet4 <- run_simulation(timesteps = sim_length, parameters = P2U)
output_bednet5 <- run_simulation(timesteps = sim_length, parameters = RGU)

# Combine prevalence plots for all scenarios
cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Function to plot combined prevalence for all scenarios
plot_combined_prevalence <- function() {
  par(mar = c(5, 5, 2, 2)) # Set margins
  
  # Plot control scenario
  plot(output_control$timestep, output_control$n_detect_91_25550 / output_control$n_91_25550,
       type = "l", col = cols[1], lwd = 2,
       xlab = "Time (days)", ylab = "Prevalence",
       ylim = c(0, 1), xaxs = "i", yaxs = "i",
       main = "Prevalence Under Different Bed Net Scenarios")
  
  # Add bednet scenarios
  lines(output_bednet1$timestep, output_bednet1$n_detect_91_25550 / output_bednet1$n_91_25550, col = cols[2], lwd = 2)
  lines(output_bednet2$timestep, output_bednet2$n_detect_91_25550 / output_bednet2$n_91_25550, col = cols[3], lwd = 2)
  lines(output_bednet3$timestep, output_bednet3$n_detect_91_25550 / output_bednet3$n_91_25550, col = cols[4], lwd = 2)
  lines(output_bednet4$timestep, output_bednet4$n_detect_91_25550 / output_bednet4$n_91_25550, col = cols[5], lwd = 2)
  lines(output_bednet5$timestep, output_bednet5$n_detect_91_25550 / output_bednet5$n_91_25550, col = cols[6], lwd = 2)
  
  # Add bednet distribution vertical lines
  abline(v = c(1, 3) * year, col = "gray", lty = 1)
  
  # Add legend
  legend("topright", legend = c("Control", "IG1U", "IG2U", "OPU", "P2U", "RGW"),
         col = cols[1:6], lty = 1, lwd = 2, box.lty = 0)
}

# Plot combined prevalence
plot_combined_prevalence()

# Function to plot combined incidence for all scenarios
plot_combined_incidence <- function() {
  par(mar = c(5, 5, 2, 2)) # Set margins
  
  # Plot control scenario
  plot(output_control$timestep, output_control$n_inc_91_25550 / output_control$n_91_25550,
       type = "l", col = cols[1], lwd = 2,
       xlab = "Time (days)", ylab = "Incidence",
       ylim = c(0, 1), xaxs = "i", yaxs = "i",
       main = "Incidence Under Different Bed Net Scenarios")
  
  # Add bednet scenarios
  lines(output_bednet1$timestep, output_bednet1$n_inc_91_25550 / output_bednet1$n_91_25550, col = cols[2], lwd = 2)
  lines(output_bednet2$timestep, output_bednet2$n_inc_91_25550 / output_bednet2$n_91_25550, col = cols[3], lwd = 2)
  lines(output_bednet3$timestep, output_bednet3$n_inc_91_25550 / output_bednet3$n_91_25550, col = cols[4], lwd = 2)
  lines(output_bednet4$timestep, output_bednet4$n_inc_91_25550 / output_bednet4$n_91_25550, col = cols[5], lwd = 2)
  lines(output_bednet5$timestep, output_bednet5$n_inc_91_25550 / output_bednet5$n_91_25550, col = cols[6], lwd = 2)
  
  # Add bednet distribution vertical lines
  abline(v = c(1,3) * year, col = "gray", lty = 0)
  
  # Add legend
  legend("topright", legend = c("Control", "IG1U", "IG2U", "OPU", "P2U", "RGW"),
         col = cols[1:6], lty = 1, lwd = 2, box.lty = 0)
}

# Plot combined incidence
plot_combined_incidence()






#########################################################################################let us see for washed
#########################################################################################let us see for washed
#########################################################################################let us see for washed



# Extract the output for year 1 (from day 0 to day 365)
year_1_timestep <- output_control$timestep[output_control$timestep <= 365]
year_1_prevalence_control <- output_control$n_detect_91_25550[output_control$timestep <= 365] / output_control$n_91_25550[output_control$timestep <= 365]

year_1_prevalence_bednet1 <- output_bednet1$n_detect_91_25550[output_bednet1$timestep <= 365] / output_bednet1$n_91_25550[output_bednet1$timestep <= 365]
year_1_prevalence_bednet2 <- output_bednet2$n_detect_91_25550[output_bednet2$timestep <= 365] / output_bednet2$n_91_25550[output_bednet2$timestep <= 365]
year_1_prevalence_bednet3 <- output_bednet3$n_detect_91_25550[output_bednet3$timestep <= 365] / output_bednet3$n_91_25550[output_bednet3$timestep <= 365]
year_1_prevalence_bednet4 <- output_bednet4$n_detect_91_25550[output_bednet4$timestep <= 365] / output_bednet4$n_91_25550[output_bednet4$timestep <= 365]
year_1_prevalence_bednet5 <- output_bednet5$n_detect_91_25550[output_bednet5$timestep <= 365] / output_bednet5$n_91_25550[output_bednet5$timestep <= 365]

# Plot the prevalence for year 1
plot(year_1_timestep, year_1_prevalence_control,
     type = "l", col = "#E69F00", lwd = 2,
     xlab = "Time (days)", ylab = "Prevalence",
     ylim = c(0, 1), xaxs = "i", yaxs = "i",
     main = "Prevalence for Year 1 Under Different Bed Net Scenarios")

lines(year_1_timestep, year_1_prevalence_bednet1, col = "#56B4E9", lwd = 2)
lines(year_1_timestep, year_1_prevalence_bednet2, col = "#009E73", lwd = 2)
lines(year_1_timestep, year_1_prevalence_bednet3, col = "#F0E442", lwd = 2)
lines(year_1_timestep, year_1_prevalence_bednet4, col = "#0072B2", lwd = 2)
lines(year_1_timestep, year_1_prevalence_bednet5, col = "#D55E00", lwd = 2)

# Add legend
legend("topright", legend = c("Control", "IG1U", "IG2U", "OPU", "P2U", "RGW"),
       col = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00"),
       lty = 1, lwd = 2, box.lty = 0)

# Ensure all prevalence vectors have the same length, trimming if necessary
n_steps <- length(year_1_timestep)  # Length of time steps
year_1_prevalence_control <- year_1_prevalence_control[1:n_steps]
year_1_prevalence_bednet1 <- year_1_prevalence_bednet1[1:n_steps]
year_1_prevalence_bednet2 <- year_1_prevalence_bednet2[1:n_steps]
year_1_prevalence_bednet3 <- year_1_prevalence_bednet3[1:n_steps]
year_1_prevalence_bednet4 <- year_1_prevalence_bednet4[1:n_steps]
year_1_prevalence_bednet5 <- year_1_prevalence_bednet5[1:n_steps]

# Create a data frame for year 1 prevalence data
year_1_prevalence <- data.frame(
  time = rep(year_1_timestep, 6),  # Adjusted to include 6 scenarios
  prevalence = c(year_1_prevalence_control, 
                 year_1_prevalence_bednet1, 
                 year_1_prevalence_bednet2, 
                 year_1_prevalence_bednet3, 
                 year_1_prevalence_bednet4, 
                 year_1_prevalence_bednet5),
  scenario = rep(c("Control", "IG1U", "IG2U", "OPU", "P2U", "RGW"), each = n_steps)
)

# Box plot for prevalence across different scenarios
library(ggplot2)

ggplot(year_1_prevalence, aes(x = scenario, y = prevalence, fill = scenario)) +
  geom_boxplot() +
  labs(title = "Prevalence for Year 1 Across Bed Net Scenarios", 
       x = "Scenario", 
       y = "Prevalence") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +  # Add a nice color palette
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels for better readability








#ear 2

# Extract the output for year 2 (from day 366 to day 730)
year_2_timestep <- output_control$timestep[output_control$timestep > 365 & output_control$timestep <= 730]
year_2_prevalence_control <- output_control$n_detect_91_25550[output_control$timestep > 365 & output_control$timestep <= 730] / output_control$n_91_25550[output_control$timestep > 365 & output_control$timestep <= 730]

year_2_prevalence_bednet1 <- output_bednet1$n_detect_91_25550[output_bednet1$timestep > 365 & output_bednet1$timestep <= 730] / output_bednet1$n_91_25550[output_bednet1$timestep > 365 & output_bednet1$timestep <= 730]
year_2_prevalence_bednet2 <- output_bednet2$n_detect_91_25550[output_bednet2$timestep > 365 & output_bednet2$timestep <= 730] / output_bednet2$n_91_25550[output_bednet2$timestep > 365 & output_bednet2$timestep <= 730]
year_2_prevalence_bednet3 <- output_bednet3$n_detect_91_25550[output_bednet3$timestep > 365 & output_bednet3$timestep <= 730] / output_bednet3$n_91_25550[output_bednet3$timestep > 365 & output_bednet3$timestep <= 730]
year_2_prevalence_bednet4 <- output_bednet4$n_detect_91_25550[output_bednet4$timestep > 365 & output_bednet4$timestep <= 730] / output_bednet4$n_91_25550[output_bednet4$timestep > 365 & output_bednet4$timestep <= 730]
year_2_prevalence_bednet5 <- output_bednet5$n_detect_91_25550[output_bednet5$timestep > 365 & output_bednet5$timestep <= 730] / output_bednet5$n_91_25550[output_bednet5$timestep > 365 & output_bednet5$timestep <= 730]

# Plot the prevalence for year 2
plot(year_2_timestep, year_2_prevalence_control,
     type = "l", col = "#E69F00", lwd = 2,
     xlab = "Time (days)", ylab = "Prevalence",
     ylim = c(0, 1), xaxs = "i", yaxs = "i",
     main = "Prevalence for Year 2 Under Different Bed Net Scenarios")

lines(year_2_timestep, year_2_prevalence_bednet1, col = "#56B4E9", lwd = 2)
lines(year_2_timestep, year_2_prevalence_bednet2, col = "#009E73", lwd = 2)
lines(year_2_timestep, year_2_prevalence_bednet3, col = "#F0E442", lwd = 2)
lines(year_2_timestep, year_2_prevalence_bednet4, col = "#0072B2", lwd = 2)
lines(year_2_timestep, year_2_prevalence_bednet5, col = "#D55E00", lwd = 2)

# Add legend
legend("topright", legend = c("Control", "IG1U", "IG2U", "OPU", "P2U", "RGW"),
       col = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00"),
       lty = 1, lwd = 2, box.lty = 0)

# Ensure all prevalence vectors have the same length, trimming if necessary
n_steps <- length(year_2_timestep)  # Length of time steps
year_2_prevalence_control <- year_2_prevalence_control[1:n_steps]
year_2_prevalence_bednet1 <- year_2_prevalence_bednet1[1:n_steps]
year_2_prevalence_bednet2 <- year_2_prevalence_bednet2[1:n_steps]
year_2_prevalence_bednet3 <- year_2_prevalence_bednet3[1:n_steps]
year_2_prevalence_bednet4 <- year_2_prevalence_bednet4[1:n_steps]
year_2_prevalence_bednet5 <- year_2_prevalence_bednet5[1:n_steps]

# Create a data frame for year 2 prevalence data
year_2_prevalence <- data.frame(
  time = rep(year_2_timestep, 6),  # Adjusted to include 6 scenarios
  prevalence = c(year_2_prevalence_control, 
                 year_2_prevalence_bednet1, 
                 year_2_prevalence_bednet2, 
                 year_2_prevalence_bednet3, 
                 year_2_prevalence_bednet4, 
                 year_2_prevalence_bednet5),
  scenario = rep(c("Control", "IG1U", "IG2U", "OPU", "P2U", "RGW"), each = n_steps)
)

# Box plot for prevalence across different scenarios
ggplot(year_2_prevalence, aes(x = scenario, y = prevalence, fill = scenario)) +
  geom_boxplot() +
  labs(title = "Prevalence for Year 2 Across Bed Net Scenarios", 
       x = "Scenario", 
       y = "Prevalence") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +  # Add a nice color palette
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels for better readability







###################year 3


# Extract the output for year 3 (from day 731 to day 1095)
year_3_timestep <- output_control$timestep[output_control$timestep > 730 & output_control$timestep <= 1095]
year_3_prevalence_control <- output_control$n_detect_91_25550[output_control$timestep > 730 & output_control$timestep <= 1095] / output_control$n_91_25550[output_control$timestep > 730 & output_control$timestep <= 1095]

year_3_prevalence_bednet1 <- output_bednet1$n_detect_91_25550[output_bednet1$timestep > 730 & output_bednet1$timestep <= 1095] / output_bednet1$n_91_25550[output_bednet1$timestep > 730 & output_bednet1$timestep <= 1095]
year_3_prevalence_bednet2 <- output_bednet2$n_detect_91_25550[output_bednet2$timestep > 730 & output_bednet2$timestep <= 1095] / output_bednet2$n_91_25550[output_bednet2$timestep > 730 & output_bednet2$timestep <= 1095]
year_3_prevalence_bednet3 <- output_bednet3$n_detect_91_25550[output_bednet3$timestep > 730 & output_bednet3$timestep <= 1095] / output_bednet3$n_91_25550[output_bednet3$timestep > 730 & output_bednet3$timestep <= 1095]
year_3_prevalence_bednet4 <- output_bednet4$n_detect_91_25550[output_bednet4$timestep > 730 & output_bednet4$timestep <= 1095] / output_bednet4$n_91_25550[output_bednet4$timestep > 730 & output_bednet4$timestep <= 1095]
year_3_prevalence_bednet5 <- output_bednet5$n_detect_91_25550[output_bednet5$timestep > 730 & output_bednet5$timestep <= 1095] / output_bednet5$n_91_25550[output_bednet5$timestep > 730 & output_bednet5$timestep <= 1095]

# Plot the prevalence for year 3
plot(year_3_timestep, year_3_prevalence_control,
     type = "l", col = "#E69F00", lwd = 2,
     xlab = "Time (days)", ylab = "Prevalence",
     ylim = c(0, 1), xaxs = "i", yaxs = "i",
     main = "Prevalence for Year 3 Under Different Bed Net Scenarios")

lines(year_3_timestep, year_3_prevalence_bednet1, col = "#56B4E9", lwd = 2)
lines(year_3_timestep, year_3_prevalence_bednet2, col = "#009E73", lwd = 2)
lines(year_3_timestep, year_3_prevalence_bednet3, col = "#F0E442", lwd = 2)
lines(year_3_timestep, year_3_prevalence_bednet4, col = "#0072B2", lwd = 2)
lines(year_3_timestep, year_3_prevalence_bednet5, col = "#D55E00", lwd = 2)

# Add legend
legend("topright", legend = c("Control", "IG1U", "IG2U", "OPU", "P2U", "RGW"),
       col = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00"),
       lty = 1, lwd = 2, box.lty = 0)

# Ensure all prevalence vectors have the same length, trimming if necessary
n_steps <- length(year_3_timestep)  # Length of time steps
year_3_prevalence_control <- year_3_prevalence_control[1:n_steps]
year_3_prevalence_bednet1 <- year_3_prevalence_bednet1[1:n_steps]
year_3_prevalence_bednet2 <- year_3_prevalence_bednet2[1:n_steps]
year_3_prevalence_bednet3 <- year_3_prevalence_bednet3[1:n_steps]
year_3_prevalence_bednet4 <- year_3_prevalence_bednet4[1:n_steps]
year_3_prevalence_bednet5 <- year_3_prevalence_bednet5[1:n_steps]

# Create a data frame for year 3 prevalence data
year_3_prevalence <- data.frame(
  time = rep(year_3_timestep, 6),  # Adjusted to include 6 scenarios
  prevalence = c(year_3_prevalence_control, 
                 year_3_prevalence_bednet1, 
                 year_3_prevalence_bednet2, 
                 year_3_prevalence_bednet3, 
                 year_3_prevalence_bednet4, 
                 year_3_prevalence_bednet5),
  scenario = rep(c("Control", "IG1U", "IG2U", "OPU", "P2U", "RGW"), each = n_steps)
)

# Box plot for prevalence across different scenarios
ggplot(year_3_prevalence, aes(x = scenario, y = prevalence, fill = scenario)) +
  geom_boxplot() +
  labs(title = "Prevalence for Year 3 Across Bed Net Scenarios", 
       x = "Scenario", 
       y = "Prevalence") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +  # Add a nice color palette
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels for better readability







######year 4
# Extract the output for year 4 (from day 1096 to day 1460)
year_4_timestep <- output_control$timestep[output_control$timestep > 1095 & output_control$timestep <= 1460]
year_4_prevalence_control <- output_control$n_detect_91_25550[output_control$timestep > 1095 & output_control$timestep <= 1460] / output_control$n_91_25550[output_control$timestep > 1095 & output_control$timestep <= 1460]

year_4_prevalence_bednet1 <- output_bednet1$n_detect_91_25550[output_bednet1$timestep > 1095 & output_bednet1$timestep <= 1460] / output_bednet1$n_91_25550[output_bednet1$timestep > 1095 & output_bednet1$timestep <= 1460]
year_4_prevalence_bednet2 <- output_bednet2$n_detect_91_25550[output_bednet2$timestep > 1095 & output_bednet2$timestep <= 1460] / output_bednet2$n_91_25550[output_bednet2$timestep > 1095 & output_bednet2$timestep <= 1460]
year_4_prevalence_bednet3 <- output_bednet3$n_detect_91_25550[output_bednet3$timestep > 1095 & output_bednet3$timestep <= 1460] / output_bednet3$n_91_25550[output_bednet3$timestep > 1095 & output_bednet3$timestep <= 1460]
year_4_prevalence_bednet4 <- output_bednet4$n_detect_91_25550[output_bednet4$timestep > 1095 & output_bednet4$timestep <= 1460] / output_bednet4$n_91_25550[output_bednet4$timestep > 1095 & output_bednet4$timestep <= 1460]
year_4_prevalence_bednet5 <- output_bednet5$n_detect_91_25550[output_bednet5$timestep > 1095 & output_bednet5$timestep <= 1460] / output_bednet5$n_91_25550[output_bednet5$timestep > 1095 & output_bednet5$timestep <= 1460]

# Plot the prevalence for year 4
plot(year_4_timestep, year_4_prevalence_control,
     type = "l", col = "#E69F00", lwd = 2,
     xlab = "Time (days)", ylab = "Prevalence",
     ylim = c(0, 1), xaxs = "i", yaxs = "i",
     main = "Prevalence for Year 4 Under Different Bed Net Scenarios")

lines(year_4_timestep, year_4_prevalence_bednet1, col = "#56B4E9", lwd = 2)
lines(year_4_timestep, year_4_prevalence_bednet2, col = "#009E73", lwd = 2)
lines(year_4_timestep, year_4_prevalence_bednet3, col = "#F0E442", lwd = 2)
lines(year_4_timestep, year_4_prevalence_bednet4, col = "#0072B2", lwd = 2)
lines(year_4_timestep, year_4_prevalence_bednet5, col = "#D55E00", lwd = 2)

# Add legend
legend("topright", legend = c("Control", "IG1U", "IG2U", "OPU", "P2U", "RGW"),
       col = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00"),
       lty = 1, lwd = 2, box.lty = 0)

# Ensure all prevalence vectors have the same length, trimming if necessary
n_steps <- length(year_4_timestep)  # Length of time steps
year_4_prevalence_control <- year_4_prevalence_control[1:n_steps]
year_4_prevalence_bednet1 <- year_4_prevalence_bednet1[1:n_steps]
year_4_prevalence_bednet2 <- year_4_prevalence_bednet2[1:n_steps]
year_4_prevalence_bednet3 <- year_4_prevalence_bednet3[1:n_steps]
year_4_prevalence_bednet4 <- year_4_prevalence_bednet4[1:n_steps]
year_4_prevalence_bednet5 <- year_4_prevalence_bednet5[1:n_steps]

# Create a data frame for year 4 prevalence data
year_4_prevalence <- data.frame(
  time = rep(year_4_timestep, 6),  # Adjusted to include 6 scenarios
  prevalence = c(year_4_prevalence_control, 
                 year_4_prevalence_bednet1, 
                 year_4_prevalence_bednet2, 
                 year_4_prevalence_bednet3, 
                 year_4_prevalence_bednet4, 
                 year_4_prevalence_bednet5),
  scenario = rep(c("Control", "IG1U", "IG2U", "OPU", "P2U", "RGW"), each = n_steps)
)

# Box plot for prevalence across different scenarios
ggplot(year_4_prevalence, aes(x = scenario, y = prevalence, fill = scenario)) +
  geom_boxplot() +
  labs(title = "Prevalence for Year 4 Across Bed Net Scenarios", 
       x = "Scenario", 
       y = "Prevalence") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +  # Add a nice color palette
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate


###year 5

# Extract the output for year 5 (from day 1461 to day 1825)
year_5_timestep <- output_control$timestep[output_control$timestep > 1460 & output_control$timestep <= 1825]
year_5_prevalence_control <- output_control$n_detect_91_25550[output_control$timestep > 1460 & output_control$timestep <= 1825] / output_control$n_91_25550[output_control$timestep > 1460 & output_control$timestep <= 1825]

year_5_prevalence_bednet1 <- output_bednet1$n_detect_91_25550[output_bednet1$timestep > 1460 & output_bednet1$timestep <= 1825] / output_bednet1$n_91_25550[output_bednet1$timestep > 1460 & output_bednet1$timestep <= 1825]
year_5_prevalence_bednet2 <- output_bednet2$n_detect_91_25550[output_bednet2$timestep > 1460 & output_bednet2$timestep <= 1825] / output_bednet2$n_91_25550[output_bednet2$timestep > 1460 & output_bednet2$timestep <= 1825]
year_5_prevalence_bednet3 <- output_bednet3$n_detect_91_25550[output_bednet3$timestep > 1460 & output_bednet3$timestep <= 1825] / output_bednet3$n_91_25550[output_bednet3$timestep > 1460 & output_bednet3$timestep <= 1825]
year_5_prevalence_bednet4 <- output_bednet4$n_detect_91_25550[output_bednet4$timestep > 1460 & output_bednet4$timestep <= 1825] / output_bednet4$n_91_25550[output_bednet4$timestep > 1460 & output_bednet4$timestep <= 1825]
year_5_prevalence_bednet5 <- output_bednet5$n_detect_91_25550[output_bednet5$timestep > 1460 & output_bednet5$timestep <= 1825] / output_bednet5$n_91_25550[output_bednet5$timestep > 1460 & output_bednet5$timestep <= 1825]

# Plot the prevalence for year 5
plot(year_5_timestep, year_5_prevalence_control,
     type = "l", col = "#E69F00", lwd = 2,
     xlab = "Time (days)", ylab = "Prevalence",
     ylim = c(0, 1), xaxs = "i", yaxs = "i",
     main = "Prevalence for Year 5 Under Different Bed Net Scenarios")

lines(year_5_timestep, year_5_prevalence_bednet1, col = "#56B4E9", lwd = 2)
lines(year_5_timestep, year_5_prevalence_bednet2, col = "#009E73", lwd = 2)
lines(year_5_timestep, year_5_prevalence_bednet3, col = "#F0E442", lwd = 2)
lines(year_5_timestep, year_5_prevalence_bednet4, col = "#0072B2", lwd = 2)
lines(year_5_timestep, year_5_prevalence_bednet5, col = "#D55E00", lwd = 2)

# Add legend
legend("topright", legend = c("Control", "IG1U", "IG2U", "OPU", "P2U", "RGW"),
       col = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00"),
       lty = 1, lwd = 2, box.lty = 0)

# Ensure all prevalence vectors have the same length, trimming if necessary
n_steps <- length(year_5_timestep)  # Length of time steps
year_5_prevalence_control <- year_5_prevalence_control[1:n_steps]
year_5_prevalence_bednet1 <- year_5_prevalence_bednet1[1:n_steps]
year_5_prevalence_bednet2 <- year_5_prevalence_bednet2[1:n_steps]
year_5_prevalence_bednet3 <- year_5_prevalence_bednet3[1:n_steps]
year_5_prevalence_bednet4 <- year_5_prevalence_bednet4[1:n_steps]
year_5_prevalence_bednet5 <- year_5_prevalence_bednet5[1:n_steps]

# Create a data frame for year 5 prevalence data
year_5_prevalence <- data.frame(
  time = rep(year_5_timestep, 6),  # Adjusted to include 6 scenarios
  prevalence = c(year_5_prevalence_control, 
                 year_5_prevalence_bednet1, 
                 year_5_prevalence_bednet2, 
                 year_5_prevalence_bednet3, 
                 year_5_prevalence_bednet4, 
                 year_5_prevalence_bednet5),
  scenario = rep(c("Control", "IG1U", "IG2U", "OPU", "P2U", "RGW"), each = n_steps)
)

# Box plot for prevalence across different scenarios
ggplot(year_5_prevalence, aes(x = scenario, y = prevalence, fill = scenario)) +
  geom_boxplot() +
  labs(title = "Prevalence for Year 5 Across Bed Net Scenarios", 
       x = "Scenario", 
       y = "Prevalence") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +  # Add a nice color palette
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels for better readability












#Add a 'year' column to each data frame
year_1_prevalence$year <- "Year 1"
year_2_prevalence$year <- "Year 2"
year_3_prevalence$year <- "Year 3"
year_4_prevalence$year <- "Year 4"
year_5_prevalence$year <- "Year 5"

# Combine all years into one data frame
all_years_prevalence <- rbind(year_1_prevalence, year_2_prevalence, 
                              year_3_prevalence, year_4_prevalence, 
                              year_5_prevalence)




ggplot(all_years_prevalence, aes(x = year, y = prevalence, fill = scenario)) +
  geom_boxplot(width = 1) +  # Increase box plot width
  geom_vline(xintercept = c(1, 3), linetype = "dotted", color = "black") +  # Add dotted vertical lines
  labs(title = "Prevalence Across Years and Bed Net Scenarios", 
       x = "Year", 
       y = "Prevalence") +
  ylim(0, 0.45) +  # Set y-axis limits
  scale_fill_brewer(palette = "Set3") +  # Add a nice color palette
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    panel.grid = element_blank(),  # Remove grid lines
    #panel.border = element_blank(),  # Remove full border
    axis.line = element_line(color = "black", size = 0.5),  # Add x and y axis lines
    plot.margin = unit(c(1, 1, 0.7, 1), "cm")  # Add margin around the graph
  )




ggplot(all_years_prevalence, aes(x = year, y = prevalence, fill = scenario)) +
  geom_boxplot(width = 1) +  # Box plot width
  geom_vline(xintercept = c(1, 3), linetype = "dotted", color = "black") +  # Dotted vertical lines at year 1 and year 3
  labs(title = "Prevalence Across Years and Bed Net Scenarios", 
       x = "Year", 
       y = "Prevalence") +
  ylim(0, 0.45) +  # Set y-axis limits
  scale_fill_brewer(palette = "Set3") +  # Color palette
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black", size = 0.5),  # Add x and y axis lines
    plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "cm"),  # Add margins around the plot
    #panel.border = element_rect(color = "black", size = 1),  # Border around the plot panel
    plot.border = element_rect(color = "black", size = 2)  # Border around the entire plot
  )




