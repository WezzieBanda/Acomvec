
library(malariasimulation)

# Simulation parameters
year <- 365
sim_length <- 5 * year
human_population <- 8000
starting_EIR <- 20

#?set_drugs()
#?set_bednets()





# Update simulation parameters with seasonality
simparams <- get_parameters(
  list(
    human_population = human_population, 
    clinical_incidence_rendering_min_ages = 91,
    clinical_incidence_rendering_max_ages = 70 * 365, 
    prevalence_rendering_min_ages = 91, 
    prevalence_rendering_max_ages = 70 * 365,
    model_seasonality = TRUE,
    g0 = 3.19575186,
    g = c(4.685451743, 1.541884506, 0.076556693),
    h = c(2.867258988, 2.483677998, 0.915187495)
  )
)

#3	Malawi	MWI	Chikwawa	2.298702841	3.094974077	1.926434341	0.672370947	1.342555063	1.020904811	0.309695984
#26	Malawi	MWI	Salima	3.19575186	4.685451743	1.541884506	0.076556693	2.867258988	2.483677998	0.915187495


# Set species parameters with seasonality
species_params <- set_species(
  parameters = simparams,
  species = list(fun_params, gamb_params),
  proportions = c(0.408, 0.592)
)

# Set drug parameters with seasonality
drug_params <- set_drugs(species_params, list(AL_params, DHA_PQP_params))




# Set treatment parameters with seasonality
treatment_params <- set_clinical_treatment(
  parameters = drug_params,
  drug = 1,
  timesteps =  c(300,600),
  coverages =  c(0.4, 0)
)


# Use set_equilibrium to update the parameter set for a given initial EIR
treatment_params <- set_equilibrium(parameters=treatment_params, init_EIR=starting_EIR)

###set bednets

# Set equilibrium
#species_params <- set_equilibrium(parameters = species_params, init_EIR = starting_EIR)

# Run control simulation (without bed nets)
output_control <- run_simulation(timesteps = sim_length, parameters = treatment_params)




# Define bednet parameters for five different scenarios
IG1U <- set_bednets(
  parameters = treatment_params,
  timesteps = c(1, 3) * year,
  coverages = c(0.9, 0.9),
  retention = 2 * 365,
  dn0 = matrix(cbind(c(0.354621849, 0.354621849), c(0.155488951, 0.155488951)), nrow = 2, ncol = 2),
  rn = matrix(cbind(c(0.565286611, 0.565286611), c(0.526669837, 0.526669837)), nrow = 2, ncol = 2),
  rnm = matrix(cbind(c(0.22, 0.22), c(0.22, 0.22)), nrow = 2, ncol = 2),
  gamman = rep(2.895257* 365, 2)
)

IG2U <- set_bednets(
  parameters = treatment_params,
  timesteps = c(1, 3) * year,
  coverages = c(0.9, 0.9),
  retention = 2 * 365,
  dn0 = matrix(cbind(c(0.37838552, 0.37838552), c(0.165655393, 0.165655393)), nrow = 2, ncol = 2),
  rn = matrix(cbind(c(0.541522693	, 0.541522693	), c(0.516502589, 0.516502589)), nrow = 2, ncol = 2),
  rnm = matrix(cbind(c(0.22, 0.22), c(0.22, 0.22)), nrow = 2, ncol = 2),
  gamman = rep(2.920989* 365, 2)
)


OPU <- set_bednets(
  parameters = treatment_params,
  timesteps = c(1, 3) * year,
  coverages = c(0.9, 0.9),
  retention = 2 * 365,
  dn0 = matrix(cbind(c(0.465751648, 0.465751648), c(0.210416336, 0.210416336)), nrow = 2, ncol = 2),
  rn = matrix(cbind(c(0.450628788, 0.450628788), c(0.4602951, 0.4602951)), nrow = 2, ncol = 2),
  rnm = matrix(cbind(c(0.22, 0.22), c(0.22, 0.22)), nrow = 2, ncol = 2),
  gamman = rep(3.194394 * 365, 2)
)



P2U <- set_bednets(
  parameters = treatment_params,
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
  parameters = treatment_params,
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
  abline(v = c(1) * year, col = "darkgray", lty = "dashed")
  
  # Add text label near the dashed line
  text(365, 0.9, "Bed Net Distribution", srt = 0, pos = 4, cex = 0.9)
  
  # Add legend
  legend("topright", legend = c("Control", "IG1U", "IG2U", "OPU", "P2U", "RGU"),
         col = cols[1:6], lty = 1, lwd = 2, box.lty = 0)
}

# Plot combined prevalence
plot_combined_prevalence()

