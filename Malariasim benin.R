############################################
## FULL MALARIA SIMULATION – CORRECT DIMS
## 3 species × 4 bednet time steps
############################################

library(malariasimulation)
library(dplyr)
library(ggplot2)

## -------------------------------
## Time definitions
## -------------------------------
year  <- 365
month <- year / 12
sim_length <- 4 * year

## -------------------------------
## Population & transmission
## -------------------------------
human_population <- 10000
starting_EIR <- 180

## -------------------------------
## Core simulation parameters
## -------------------------------
simparams <- get_parameters(
  list(
    human_population = human_population,
    clinical_incidence_rendering_min_ages = 91,
    clinical_incidence_rendering_max_ages = 70 * year,
    prevalence_rendering_min_ages = 91,
    prevalence_rendering_max_ages = 70 * year,
    model_seasonality = TRUE,
    g0 = 2.965573076,
    g  = c(4.403474653, 1.237759796, -0.219862959),
    h  = c(2.567688899, 2.400592056, 0.256744044)
  )
)

## -------------------------------
## Mosquito species (ORDER MATTERS)
## -------------------------------
simparams <- set_species(
  parameters = simparams,
  species = list(fun_params, gamb_params, arab_params),
  proportions = c(0.40, 0.40, 0.20)
)

## -------------------------------
## Drugs and treatment
## -------------------------------
simparams <- set_drugs(
  parameters = simparams,
  drugs = list(AL_params, DHA_PQP_params)
)

simparams <- set_clinical_treatment(
  parameters = simparams,
  drug = 1,
  timesteps = c(300, 600),
  coverages = c(0.4, 0)
)

## -------------------------------
## Set equilibrium
## -------------------------------
simparams <- set_equilibrium(
  parameters = simparams,
  init_EIR = starting_EIR
)

## -------------------------------
## Bednet time steps (4 events)
## -------------------------------
bednet_timesteps <- c(1,4) * year

## -------------------------------
## Bednet parameters (3 × 4)
## rows = species, cols = time
## -------------------------------

## Killing effect
dn0 <- matrix(
  c(
    ## funestus
    0.30, 0.20, 0.10, 0.05,
    ## gambiae
    0.35, 0.25, 0.15, 0.05,
    ## arabiensis
    0.25, 0.18, 0.10, 0.04
  ),
  nrow = 3,
  ncol = 4,
  byrow = TRUE
)

## Repellency
rn <- matrix(
  c(
    0.70, 0.55, 0.35, 0.20,
    0.75, 0.60, 0.40, 0.25,
    0.60, 0.45, 0.30, 0.15
  ),
  nrow = 3,
  ncol = 4,
  byrow = TRUE
)

## Mortality modifier
rnm <- matrix(
  rep(0.22, 3 * 4),
  nrow = 3,
  ncol = 4,
  byrow = TRUE
)

## Net decay rate (species-specific)
gamman <- rep(2.67 * year, 3)

## -------------------------------
## Sanity checks (DO NOT SKIP)
## -------------------------------
stopifnot(
  length(bednet_timesteps) == 4,
  ncol(dn0) == 4,
  ncol(rn) == 4,
  ncol(rnm) == 4,
  nrow(dn0) == 3,
  length(gamman) == 3
)

## -------------------------------
## Set bednets
## -------------------------------
simparams <- set_bednets(
  parameters = simparams,
  timesteps  = bednet_timesteps,
  coverages  = c(0.8, 0.8, 0.8, 0.8),
  retention  = 1e10,
  dn0        = dn0,
  rn         = rn,
  rnm        = rnm,
  gamman     = gamman
)

## -------------------------------
## Run simulation
## -------------------------------
output <- run_simulation(
  simparams,
  timesteps = sim_length
)

## -------------------------------
## Quick diagnostic plot
## -------------------------------
plot(
  output$timestep / year,
  output$PfPR,
  type = "l",
  lwd = 2,
  xlab = "Time (years)",
  ylab = "PfPR"
)






















































library(malariasimulation)

library(dplyr)
library(ggplot2)

# Simulation parameters
year <- 365
sim_length <- 6 * year
human_population <- 10000
starting_EIR <- 20


# Update simulation parameters with seasonality
simparams <- get_parameters(
  list(
    human_population = human_population, 
    clinical_incidence_rendering_min_ages = 91,
    clinical_incidence_rendering_max_ages = 70 * 365, 
    prevalence_rendering_min_ages = 91, 
    prevalence_rendering_max_ages = 70 * 365,
    model_seasonality = TRUE,
    g0 = 2.965573076,
    g = c(4.403474653,	1.237759796,	-0.219862959),
    h = c(2.567688899,	2.400592056,	0.256744044)
  )
)

#3	Malawi	MWI	Chikwawa	2.298702841	3.094974077	1.926434341	0.672370947	1.342555063	1.020904811	0.309695984
#26	Malawi	MWI	Salima	3.19575186	4.685451743	1.541884506	0.076556693	2.867258988	2.483677998	0.915187495
#Kasungu	2.965573076	4.403474653	1.237759796	-0.219862959	2.567688899	2.400592056	0.256744044


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
  timesteps = c(1, 4) * year,
  coverages = c(0.9, 0.9),
  retention = 2 * 365,
  dn0 = matrix(cbind(c(0.0114 , 0.0114 ), c(0.0114 , 0.0114 )), nrow = 2, ncol = 2),
  rn = matrix(cbind(c(0.981, 0.981), c(0.981, 0.981)), nrow = 2, ncol = 2),
  rnm = matrix(cbind(c(0.22, 0.22), c(0.22, 0.22)), nrow = 2, ncol = 2),
  gamman = rep( 3.37* 365, 2)
)

IG2U <- set_bednets(
  parameters = treatment_params,
  timesteps = c(1, 4) * year,
  coverages = c(0.9, 0.9),
  retention = 2 * 365,
  dn0 = matrix(cbind(c(0.0780 , 0.0780 ), c(0.0780 , 0.0780 )), nrow = 2, ncol = 2),
  rn = matrix(cbind(c(0.917	, 0.917	), c(0.917, 0.917)), nrow = 2, ncol = 2),
  rnm = matrix(cbind(c(0.22, 0.22), c(0.22, 0.22)), nrow = 2, ncol = 2),
  gamman = rep(6.93* 365, 2)
)


P3 <- set_bednets(
  parameters = treatment_params,
  timesteps = c(1, 4) * year,
  coverages = c(0.9, 0.9),
  retention = 2 * 365,
  dn0 = matrix(cbind(c(0.0313 , 0.0313 ), c(0.0313 , 0.0313 )), nrow = 2, ncol = 2),
  rn = matrix(cbind(c(0.962, 0.962), c(0.962, 0.962)), nrow = 2, ncol = 2),
  rnm = matrix(cbind(c(0.22, 0.22), c(0.22, 0.22)), nrow = 2, ncol = 2),
  gamman = rep(4.48 * 365, 2)
)



RG <- set_bednets(
  parameters = treatment_params,
  timesteps = c(1, 4) * year,
  coverages = c(0.9, 0.9),
  retention = 2 * 365,
  dn0 = matrix(cbind(c(0.0273, 0.0273), c(	0.0273, 	0.0273)), nrow = 2, ncol = 2),
  rn = matrix(cbind(c(0.966, 0.966), c(0.966, 0.966)), nrow = 2, ncol = 2),
  rnm = matrix(cbind(c(0.22, 0.22), c(0.22, 0.22)), nrow = 2, ncol = 2),
  #gamman = matrix(rep(2.65196 * 365, 2), nrow = 2, ncol = 1) # 2x1 matrix for gamman
  gamman = rep(4.27 * 365, 2)
)


#RGU <- set_bednets(
 # parameters = treatment_params,
  #timesteps = c(1, 4) * year,
  #coverages = c(0.9, 0.9),
  #retention = 2 * 365,
  #dn0 = matrix(cbind(c(0.169164837, 0.169164837), c(0.132083126, 0.132083126)), nrow = 2, ncol = 2),
 # rn = matrix(cbind(c(0.754104889, 0.754104889), c(0.565991711, 0.565991711)), nrow = 2, ncol = 2),
  #rnm = matrix(cbind(c(0.22, 0.22), c(0.22, 0.22)), nrow = 2, ncol = 2),
  #gamman = rep(2.874849  * 365, 2)
#)

# Run simulations for each bednet parameter
output_bednet1 <- run_simulation(timesteps = sim_length, parameters = IG1U)
output_bednet2 <- run_simulation(timesteps = sim_length, parameters = IG2U)
output_bednet3 <- run_simulation(timesteps = sim_length, parameters = P3)
output_bednet4 <- run_simulation(timesteps = sim_length, parameters = RG)
#output_bednet5 <- run_simulation(timesteps = sim_length, parameters = RGU)

# Combine prevalence plots for all scenarios
cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


#############adding november 2024 

# Function to plot combined prevalence with specified x-axis labels
plot_combined_prevalence <- function() {
  # Define the specific x-axis labels
  start_date <- as.Date("2024-11-01")  # Start date: November 2024
  x_labels <- seq(start_date, by = "1 year", length.out = 6)  # Add 1 year intervals for 5 years
  x_label_text <- format(x_labels, "%b %Y")  # Format as "Nov YYYY"
  
  # Calculate positions for these labels (365 days apart from each other)
  x_positions <- seq(365, by = 365, length.out = length(x_labels))
  
  # Set up the plot with customized x-axis
  plot(output_control$timestep, output_control$n_detect_91_25550 / output_control$n_91_25550,
       type = "l", col = cols[1], lwd = 2,
       xlab = "Time", ylab = "Prevalence",
       xaxt = "n", # Disable default x-axis
       ylim = c(0, 1), xaxs = "i", yaxs = "i",
       main = "Prevalence Under unwashed Bednets")
  
  # Add the custom x-axis labels
  axis(1, at = x_positions, labels = x_label_text)
  
  # Add bednet scenarios
  lines(output_bednet1$timestep, output_bednet1$n_detect_91_25550 / output_bednet1$n_91_25550, col = cols[2], lwd = 2)
  lines(output_bednet2$timestep, output_bednet2$n_detect_91_25550 / output_bednet2$n_91_25550, col = cols[3], lwd = 2)
  lines(output_bednet3$timestep, output_bednet3$n_detect_91_25550 / output_bednet3$n_91_25550, col = cols[4], lwd = 2)
  lines(output_bednet4$timestep, output_bednet4$n_detect_91_25550 / output_bednet4$n_91_25550, col = cols[5], lwd = 2)
  #lines(output_bednet5$timestep, output_bednet5$n_detect_91_25550 / output_bednet5$n_91_25550, col = cols[6], lwd = 2)
  
  # Add bednet distribution vertical lines
  abline(v = c(1, 4) * year, col = "darkgray", lty = "dashed")
  
  # Add text label near the dashed line
  text(365, 0.9, "Bed Net Distribution", srt = 0, pos = 4, cex = 0.9)
  text(1095, 0.9, "Bed Net Distribution", srt = 0, pos = 4, cex = 0.9)
  
  
  # Add legend
  legend("topright", legend = c("Control", "IG1U", "IG2U", "P3", "RG"),
         col = cols[1:6], lty = 1, lwd = 2, box.lty = 0)
}

# Plot combined prevalence
plot_combined_prevalence()

