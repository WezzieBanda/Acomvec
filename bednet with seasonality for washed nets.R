
############################################washed bednets####################################
##########################################################################################33

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
    g0 = 2.298702841,
    g = c(3.094974077,	1.926434341,	0.672370947),
    h = c(1.342555063, 1.020904811, 	0.309695984)
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
  timesteps = c(1, 4) * year,
  coverages = c(0.9, 0.9),
  retention = 2 * 365,
  dn0 = matrix(cbind(c(0.216110006	, 0.216110006	), c(0.208768885, 0.208768885)), nrow = 2, ncol = 2),
  rn = matrix(cbind(c(0.413668575, 0.413668575), c(0.497126749, 0.497126749)), nrow = 2, ncol = 2),
  rnm = matrix(cbind(c(0.22, 0.22), c(0.22, 0.22)), nrow = 2, ncol = 2),
  gamman = rep(2.895257* 365, 2)
)

IG2U <- set_bednets(
  parameters = treatment_params,
  timesteps = c(1, 4) * year,
  coverages = c(0.9, 0.9),
  retention = 2 * 365,
  dn0 = matrix(cbind(c(0.226342172, 0.226342172), c(0.218671298, 0.218671298)), nrow = 2, ncol = 2),
  rn = matrix(cbind(c(0.403435514	, 0.403435514	), c(0.487223575, 0.487223575)), nrow = 2, ncol = 2),
  rnm = matrix(cbind(c(0.22, 0.22), c(0.22, 0.22)), nrow = 2, ncol = 2),
  gamman = rep(2.920989* 365, 2)
)

#0.403435514	0.370222314	0.226342172		IG2W	0.487223575	0.294105128	0.218671298

OPU <- set_bednets(
  parameters = treatment_params,
  timesteps = c(1, 4) * year,
  coverages = c(0.9, 0.9),
  retention = 2 * 365,
  dn0 = matrix(cbind(c(0.26432241, 0.26432241), c(0.261759846, 0.261759846)), nrow = 2, ncol = 2),
  rn = matrix(cbind(c(0.35277145, 0.35277145), c(0.433309073, 0.433309073)), nrow = 2, ncol = 2),
  rnm = matrix(cbind(c(0.22, 0.22), c(0.22, 0.22)), nrow = 2, ncol = 2),
  gamman = rep(3.194394 * 365, 2)
)

#0.35277145	0.38290614	0.26432241		OPW	0.433309073	0.304931081	0.261759846


P2U <- set_bednets(
  parameters = treatment_params,
  timesteps = c(1, 4) * year,
  coverages = c(0.9, 0.9),
  retention = 2 * 365,
  dn0 = matrix(cbind(c(0.236059753, 0.236059753), c(0.223662126, 	0.223662126)), nrow = 2, ncol = 2),
  rn = matrix(cbind(c(0.472185623, 0.472185623), c(0.54691576, 0.54691576)), nrow = 2, ncol = 2),
  rnm = matrix(cbind(c(0.22, 0.22), c(0.22, 0.22)), nrow = 2, ncol = 2),
  #gamman = matrix(rep(2.65196 * 365, 2), nrow = 2, ncol = 1) # 2x1 matrix for gamman
  gamman = rep(2.938527 * 365, 2)
)

#0.472185623	0.291754624	0.236059753		P2W	0.54691576	0.229422114	0.223662126
#0.454675457	0.351027914	0.194296629		RGW	0.532837814	0.279613204	0.187548982

RGU <- set_bednets(
  parameters = treatment_params,
  timesteps = c(1, 4) * year,
  coverages = c(0.9, 0.9),
  retention = 2 * 365,
  dn0 = matrix(cbind(c(0.194296629, 0.194296629), c(0.187548982, 0.187548982)), nrow = 2, ncol = 2),
  rn = matrix(cbind(c(0.454675457, 0.454675457), c(0.532837814, 0.532837814)), nrow = 2, ncol = 2),
  rnm = matrix(cbind(c(0.22, 0.22), c(0.22, 0.22)), nrow = 2, ncol = 2),
  gamman = rep(2.874849  * 365, 2)
)

#0.454675457	0.351027914	0.194296629		RGW	0.532837814	0.279613204	0.187548982

# Run simulations for each bednet parameter
output_bednet1 <- run_simulation(timesteps = sim_length, parameters = IG1U)
output_bednet2 <- run_simulation(timesteps = sim_length, parameters = IG2U)
output_bednet3 <- run_simulation(timesteps = sim_length, parameters = OPU)
output_bednet4 <- run_simulation(timesteps = sim_length, parameters = P2U)
output_bednet5 <- run_simulation(timesteps = sim_length, parameters = RGU)

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
       main = "Prevalence Under Different Bed Net Scenarios")
  
  # Add the custom x-axis labels
  axis(1, at = x_positions, labels = x_label_text)
  
  # Add bednet scenarios
  lines(output_bednet1$timestep, output_bednet1$n_detect_91_25550 / output_bednet1$n_91_25550, col = cols[2], lwd = 2)
  lines(output_bednet2$timestep, output_bednet2$n_detect_91_25550 / output_bednet2$n_91_25550, col = cols[3], lwd = 2)
  lines(output_bednet3$timestep, output_bednet3$n_detect_91_25550 / output_bednet3$n_91_25550, col = cols[4], lwd = 2)
  lines(output_bednet4$timestep, output_bednet4$n_detect_91_25550 / output_bednet4$n_91_25550, col = cols[5], lwd = 2)
  lines(output_bednet5$timestep, output_bednet5$n_detect_91_25550 / output_bednet5$n_91_25550, col = cols[6], lwd = 2)
  
  # Add bednet distribution vertical lines
  abline(v = c(1, 4) * year, col = "darkgray", lty = "dashed")
  
  # Add text label near the dashed line
  text(365, 0.9, "Bed Net Distribution", srt = 0, pos = 4, cex = 0.9)
  text(1095, 0.9, "Bed Net Distribution", srt = 0, pos = 4, cex = 0.9)
  
  
  # Add legend
  legend("topright", legend = c("Control", "IG1U", "IG2U", "OPU", "P2U", "RGU"),
         col = cols[1:6], lty = 1, lwd = 2, box.lty = 0)
}

# Plot combined prevalence
plot_combined_prevalence()


############################################################################################################################################
######################################################################################################################################################
#########################total number of cases
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################


# Add a "Day" column to each dataset
output_bednet1 <- output_bednet1 %>% mutate(Day = 1:n())
output_bednet2 <- output_bednet2 %>% mutate(Day = 1:n())
output_bednet3 <- output_bednet3 %>% mutate(Day = 1:n())
output_bednet4 <- output_bednet4 %>% mutate(Day = 1:n())
output_bednet5 <- output_bednet5 %>% mutate(Day = 1:n())
output_control <- output_control %>% mutate(Day = 1:n(), Bednet = "Control")

# Combine all datasets into one, adding a "Bednet" column
all_data <- bind_rows(
  output_bednet1 %>% mutate(Bednet = "IG1"),
  output_bednet2 %>% mutate(Bednet = "IG2"),
  output_bednet3 %>% mutate(Bednet = "OP"),
  output_bednet4 %>% mutate(Bednet = "P2"),
  output_bednet5 %>% mutate(Bednet = "RG"),
  output_control
)

# Define custom year ranges
year_ranges <- list(
  "Nov 2024 - Nov 2025" = 365:730,
  "Nov 2025 - Nov 2026" = 731:1096,
  "Nov 2026 - Nov 2027" = 1097:1462,
  "Nov 2027 - Nov 2028" = 1463:1828,
  "Nov 2028 - Nov 2029" = 1829:2190
)

# Summarize cases by year range and bednet type
cases_per_range <- all_data %>%
  rowwise() %>%
  mutate(Year_Range = case_when(
    Day %in% year_ranges[["Nov 2024 - Nov 2025"]] ~ "Nov 2024 - Nov 2025",
    Day %in% year_ranges[["Nov 2025 - Nov 2026"]] ~ "Nov 2025 - Nov 2026",
    Day %in% year_ranges[["Nov 2026 - Nov 2027"]] ~ "Nov 2026 - Nov 2027",
    Day %in% year_ranges[["Nov 2027 - Nov 2028"]] ~ "Nov 2027 - Nov 2028",
    Day %in% year_ranges[["Nov 2028 - Nov 2029"]] ~ "Nov 2028 - Nov 2029",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Year_Range)) %>%
  group_by(Year_Range, Bednet) %>%
  summarise(Total_Cases = sum(n_inc_clinical_91_25550, na.rm = TRUE), .groups = "drop")

# Custom colors for the plot
custom_colors <- c(
  "IG1" = "#56B4E9",   # Light Blue
  "IG2" = "#009E73",   # Green
  "OP" = "#F0E442",    # Yellow
  "P2" = "#0072B2",    # Blue
  "RG" = "#D55E00",    # Orange
  "Control" = "#E69F00" # Gold
)

# Create the bar plot
ggplot(cases_per_range, aes(x = Year_Range, y = Total_Cases, fill = Bednet)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Total Clinical Cases by Custom Year Range and Bednet",
    x = "Year",
    y = "Total Number of Cases (Per 10,000 Population)",
    fill = "Bednet Type"
  ) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", size = 1.5, fill = NA),  # Add border around the plot
    # plot.title = element_text(hjust = 0.5, size = 14, face = "bold")# Center and style title
  )

############################################################################################################################################
######################################################################################################################################################
#########line graph version of the graph
################################################################################################################
###############################################################################################################

# Summarize the cases by year range and bednet type
cases_per_range <- all_data %>%
  rowwise() %>%
  mutate(Year_Range = case_when(
    Day %in% year_ranges[["Nov 2024 - Nov 2025"]] ~ "Nov 2024 - Nov 2025",
    Day %in% year_ranges[["Nov 2025 - Nov 2026"]] ~ "Nov 2025 - Nov 2026",
    Day %in% year_ranges[["Nov 2026 - Nov 2027"]] ~ "Nov 2026 - Nov 2027",
    Day %in% year_ranges[["Nov 2027 - Nov 2028"]] ~ "Nov 2027 - Nov 2028",
    Day %in% year_ranges[["Nov 2028 - Nov 2029"]] ~ "Nov 2028 - Nov 2029",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Year_Range)) %>%
  group_by(Year_Range, Bednet) %>%
  summarise(Total_Cases = sum(n_inc_clinical_91_25550, na.rm = TRUE), .groups = "drop")

# Convert Year_Range into a factor with the correct order
cases_per_range <- cases_per_range %>%
  mutate(Year_Range = factor(Year_Range, 
                             levels = c("Nov 2024 - Nov 2025", 
                                        "Nov 2025 - Nov 2026", 
                                        "Nov 2026 - Nov 2027", 
                                        "Nov 2027 - Nov 2028", 
                                        "Nov 2028 - Nov 2029")))

# Create the line plot with correct year range labels on the x-axis
ggplot(cases_per_range, aes(x = Year_Range, y = Total_Cases, color = Bednet, group = Bednet)) +
  geom_line(size = 1.2) +  # Line for each bednet
  geom_point(size = 3) +   # Points for each bednet on the line
  labs(
    title = "Total Clinical Cases by Custom Year Range and Bednet", 
    x = "Year Range",
    y = "Total Number of Cases (Per 10,000 Population)",
    color = "Bednet Type"
  ) +
  scale_color_manual(values = custom_colors) +  # Custom colors for each bednet type
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", size = 1.5, fill = NA),  # Add border around the plot
    #plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),      # Center and style title
    legend.position = "bottom"  # Optional: position the legend at the bottom
  )



############################################################################################################################################
######################################################################################################################################################
#########now can we do cases averted
############################################################################################################################################
######################################################################################################################################################


# Summarize total clinical cases for all bednet types and year ranges
cases_per_range <- all_data %>%
  rowwise() %>%
  mutate(Year_Range = case_when(
    Day %in% year_ranges[["Nov 2024 - Nov 2025"]] ~ "Nov 2024 - Nov 2025",
    Day %in% year_ranges[["Nov 2025 - Nov 2026"]] ~ "Nov 2025 - Nov 2026",
    Day %in% year_ranges[["Nov 2026 - Nov 2027"]] ~ "Nov 2026 - Nov 2027",
    Day %in% year_ranges[["Nov 2027 - Nov 2028"]] ~ "Nov 2027 - Nov 2028",
    Day %in% year_ranges[["Nov 2028 - Nov 2029"]] ~ "Nov 2028 - Nov 2029",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Year_Range)) %>%
  group_by(Year_Range, Bednet) %>%
  summarise(Total_Cases = sum(n_inc_clinical_91_25550, na.rm = TRUE), .groups = "drop")

# Extract control cases for each year range
control_cases <- cases_per_range %>%
  filter(Bednet == "Control") %>%
  select(Year_Range, Total_Cases) %>%
  rename(Control_Cases = Total_Cases)

# Merge the control cases back into the treatment data
cases_with_control <- cases_per_range %>%
  left_join(control_cases, by = "Year_Range") %>%
  mutate(Cases_Averted = Control_Cases - Total_Cases)

# Custom colors for the plot
custom_colors <- c(
  "IG1" = "#56B4E9",   # Light Blue
  "IG2" = "#009E73",   # Green
  "OP" = "#F0E442",    # Yellow
  "P2" = "#0072B2",    # Blue
  "RG" = "#D55E00",    # Orange
  "Control" = "#E69F00" # Gold
)

# Create the bar plot for cases averted
ggplot(cases_with_control, aes(x = Year_Range, y = Cases_Averted, fill = Bednet)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    # title = "Cases Averted by Bednet Type and Year Range",
    x = "Year Range",
    y = "Cases Averted per 10,000 population per year",
    fill = "Bednet Type"
  ) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", size = 1.5, fill = NA),  # Add border around the plot
    #plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),      # Center and style title
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),         # Customize x-axis text
    axis.text.y = element_text(size = 10)                                 # Customize y-axis text
  )
