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

# Define bednet parameters for three different scenarios
IG1U <- set_bednets(
  parameters = species_params,
  timesteps = c(1, 3) * year,
  coverages = c(0.9, 0.9),
  retention = 2 * 365,
  dn0 = matrix(cbind(c(0.13, 0.13), c(0.11, 0.11)), nrow = 2, ncol = 2),
  rn = matrix(cbind(c(0.203, 0.203), c(0.54, 0.54)), nrow = 2, ncol = 2),
  rnm = matrix(cbind(c(0.22, 0.22), c(0.22, 0.22)), nrow = 2, ncol = 2),
  gamman = rep(1.97 * 365, 2)
)

IG2U <- set_bednets(
  parameters = species_params,
  timesteps = c(1, 3) * year,
  coverages = c(0.9, 0.9),
  retention = 2 * 365,
  dn0 = matrix(cbind(c(0.32, 0.32), c(0.09, 0.09)), nrow = 2, ncol = 2),
  rn = matrix(cbind(c(0.31, 0.31), c(0.49, 0.49)), nrow = 2, ncol = 2),
  rnm = matrix(cbind(c(0.22, 0.22), c(0.22, 0.22)), nrow = 2, ncol = 2),
  gamman = rep(1.5 * 365, 2)
)

OPU <- set_bednets(
  parameters = species_params,
  timesteps = c(1, 3) * year,
  coverages = c(0.7, 0.7),
  retention = 2 * 365,
  dn0 = matrix(cbind(c(0.19, 0.19), c(0.19, 0.19)), nrow = 2, ncol = 2),
  rn = matrix(cbind(c(0.12, 0.12), c(0.58, 0.58)), nrow = 2, ncol = 2),
  rnm = matrix(cbind(c(0.22, 0.22), c(0.22, 0.22)), nrow = 2, ncol = 2),
  gamman = rep(2 * 365, 2)
)

# Run simulations for each bednet parameter
output_bednet1 <- run_simulation(timesteps = sim_length, parameters = IG1U)
output_bednet2 <- run_simulation(timesteps = sim_length, parameters = IG2U)
output_bednet3 <- run_simulation(timesteps = sim_length, parameters = OPU)

# Combine prevalence plots for all scenarios
cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

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
  
  # Add bednet distribution vertical lines
  abline(v = c(1, 2, 3, 4, 5) * year, col = "gray", lty = 1)
  
  # Add legend
  legend("topright", legend = c("Control", "IGIU", "IG2U", "OPU"),
         col = cols[1:4], lty = 1, lwd = 2, box.lty = 0)
}

# Plot combined prevalence
plot_combined_prevalence()


############now incidence#################################

