

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
human_population<- 1500
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
  clinical_incidence = c(0.5, 9) * 365,
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




# Output from your simulation
cases <- output_bednet1$n_inc_clinical_182.5_3284    # clinical episodes
population <- output_bednet1$n_age_182.5_3284        # children at risk
timesteps <- output_bednet1$timestep                 # days

# Define years
year_starts <- seq(1, max(timesteps), by = 365)
year_ends <- year_starts + 364
n_years <- length(year_starts)

# Calculate incidence per year
incidence_per_year <- numeric(n_years)

for(i in 1:n_years){
  keep <- timesteps >= year_starts[i] & timesteps <= year_ends[i]
  
  total_cases <- sum(cases[keep])
  total_child_years <- sum(population[keep]) / 365  # convert to child-years
  
  incidence_per_year[i] <- total_cases / total_child_years
}

# Show results
data.frame(
  Year = 1:n_years,
  Incidence_per_child_per_year = incidence_per_year
)




# Simulation output
cases <- output_bednet1$n_inc_clinical_182.5_3284
population <- output_bednet1$n_age_182.5_3284
timesteps <- output_bednet1$timestep

# Bed net distribution day
bednet_day <- 2310

# Define yearly windows from the bednet distribution
years <- 1:3  # number of years you want to calculate
year_length <- 365

incidence_per_year <- numeric(length(years))

for(i in seq_along(years)){
  start_day <- bednet_day + (i - 1) * year_length
  end_day   <- start_day + year_length - 1
  
  keep <- timesteps >= start_day & timesteps <= end_day
  
  total_cases <- sum(cases[keep])
  total_child_years <- sum(population[keep]) / year_length  # child-years
  
  incidence_per_year[i] <- total_cases / total_child_years
}

data.frame(
  Year_since_bednet = years,
  Incidence_per_child_per_year = incidence_per_year
)






#########################overall  incidence


# ======================================
# Extract simulation outputs
# ======================================

cases <- output_bednet1$n_inc_clinical_182.5_3284
population <- output_bednet1$n_age_182.5_3284
timesteps <- output_bednet1$timestep

year_length <- 365
bednet_day <- 2310

# ======================================
# Year 1 (days 2310–2674)
# ======================================

start_y1 <- bednet_day
end_y1   <- bednet_day + year_length - 1

keep_y1 <- timesteps >= start_y1 & timesteps <= end_y1

cases_y1 <- sum(cases[keep_y1])
child_years_y1 <- sum(population[keep_y1]) / year_length

inc_y1 <- cases_y1 / child_years_y1


# ======================================
# Year 2 (days 2675–3039)
# ======================================

start_y2 <- bednet_day + year_length
end_y2   <- start_y2 + year_length - 1

keep_y2 <- timesteps >= start_y2 & timesteps <= end_y2

cases_y2 <- sum(cases[keep_y2])
child_years_y2 <- sum(population[keep_y2]) / year_length

inc_y2 <- cases_y2 / child_years_y2


# ======================================
# Overall (Year 1 + Year 2 combined)
# ======================================

keep_y1_y2 <- timesteps >= start_y1 & timesteps <= end_y2

cases_total <- sum(cases[keep_y1_y2])
child_years_total <- sum(population[keep_y1_y2]) / year_length

inc_overall_2yrs <- cases_total / child_years_total


# ======================================
# Results table
# ======================================

results <- data.frame(
  Measure = c("Year 1 after bednet",
              "Year 2 after bednet",
              "Overall (Year 1 + Year 2)"),
  Incidence_per_child_per_year = c(inc_y1,
                                   inc_y2,
                                   inc_overall_2yrs),
  Incidence_per_1000_children_per_year = c(inc_y1,
                                           inc_y2,
                                           inc_overall_2yrs) * 1000
)

print(results)






###############################################################################################################################





library(ggplot2)

# ======================================
# 1️⃣ Observed trial data
# ======================================

observed <- data.frame(
  Period = c("Year 1", "Year 2", "Overall"),
  Incidence = c(0.36, 0.69, 0.56),
  Lower = c(0.30, 0.62, 0.51),
  Upper = c(0.42, 0.76, 0.61),
  Source = "Observed"
)

# ======================================
# 2️⃣ Model 95% CI (Poisson approximation)
# ======================================

# --- Year 1 ---
se_y1 <- sqrt(cases_y1) / child_years_y1
lower_y1 <- inc_y1 - 1.96 * se_y1
upper_y1 <- inc_y1 + 1.96 * se_y1

# --- Year 2 ---
se_y2 <- sqrt(cases_y2) / child_years_y2
lower_y2 <- inc_y2 - 1.96 * se_y2
upper_y2 <- inc_y2 + 1.96 * se_y2

# --- Overall ---
se_total <- sqrt(cases_total) / child_years_total
lower_total <- inc_overall_2yrs - 1.96 * se_total
upper_total <- inc_overall_2yrs + 1.96 * se_total

# Prevent negative lower bounds
lower_y1 <- max(0, lower_y1)
lower_y2 <- max(0, lower_y2)
lower_total <- max(0, lower_total)

# ======================================
# 3️⃣ Model dataframe
# ======================================

model <- data.frame(
  Period = c("Year 1", "Year 2", "Overall"),
  Incidence = c(inc_y1, inc_y2, inc_overall_2yrs),
  Lower = c(lower_y1, lower_y2, lower_total),
  Upper = c(upper_y1, upper_y2, upper_total),
  Source = "Model"
)

# ======================================
# 4️⃣ Combine observed + model
# ======================================

plot_data <- rbind(observed, model)

plot_data$Period <- factor(plot_data$Period,
                           levels = c("Year 1", "Year 2", "Overall"))

# ======================================
# 5️⃣ Bar Plot with CI
# ======================================

p1<- ggplot(plot_data, aes(x = Period, y = Incidence, fill = Source)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.7),
           width = 0.6) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                position = position_dodge(width = 0.7),
                width = 0.2,
                linewidth = 1) +
  labs(
    y = "Incidence per child per year",
    x = "",
    title = "Model vs Observed Malaria Incidence",
    fill = ""
  ) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c("Observed" = "darkgreen",
                               "Model" = "deeppink")) +
  ylim(0, max(plot_data$Upper) * 1.1)



p1
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
human_population<- 1500
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
  clinical_incidence = c(0.5, 9) * 365,
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



##########cases averted for IG1
#########################overall  incidence


# ======================================
# Extract simulation outputs
# ======================================

cases <- output_bednet2$n_inc_clinical_182.5_3284
population <- output_bednet2$n_age_182.5_3284
timesteps <- output_bednet2$timestep

year_length <- 365
bednet_day <- 2310

# ======================================
# Year 1 (days 2310–2674)
# ======================================

start_y1 <- bednet_day
end_y1   <- bednet_day + year_length - 1

keep_y1 <- timesteps >= start_y1 & timesteps <= end_y1

cases_y1 <- sum(cases[keep_y1])
child_years_y1 <- sum(population[keep_y1]) / year_length

inc_y1 <- cases_y1 / child_years_y1


# ======================================
# Year 2 (days 2675–3039)
# ======================================

start_y2 <- bednet_day + year_length
end_y2   <- start_y2 + year_length - 1

keep_y2 <- timesteps >= start_y2 & timesteps <= end_y2

cases_y2 <- sum(cases[keep_y2])
child_years_y2 <- sum(population[keep_y2]) / year_length

inc_y2 <- cases_y2 / child_years_y2


# ======================================
# Overall (Year 1 + Year 2 combined)
# ======================================

keep_y1_y2 <- timesteps >= start_y1 & timesteps <= end_y2

cases_total <- sum(cases[keep_y1_y2])
child_years_total <- sum(population[keep_y1_y2]) / year_length

inc_overall_2yrs <- cases_total / child_years_total


# ======================================
# Results table
# ======================================

results <- data.frame(
  Measure = c("Year 1 after bednet",
              "Year 2 after bednet",
              "Overall (Year 1 + Year 2)"),
  Incidence_per_child_per_year = c(inc_y1,
                                   inc_y2,
                                   inc_overall_2yrs),
  Incidence_per_1000_children_per_year = c(inc_y1,
                                           inc_y2,
                                           inc_overall_2yrs) * 1000
)

print(results)



##################################################################################plotting 




library(ggplot2)

# ======================================
# 1️⃣ Observed trial data
# ======================================

observed <- data.frame(
  Period = c("Year 1", "Year 2", "Overall"),
  Incidence = c(0.77, 1.19, 1.03),
  Lower = c(0.69, 1.10, 0.96),
  Upper = c(0.87, 1.29, 1.09),
  Source = "Observed"
)

# --- Year 1 ---
se_y1 <- sqrt(cases_y1) / child_years_y1
lower_y1 <- inc_y1 - 1.96 * se_y1
upper_y1 <- inc_y1 + 1.96 * se_y1

# --- Year 2 ---
se_y2 <- sqrt(cases_y2) / child_years_y2
lower_y2 <- inc_y2 - 1.96 * se_y2
upper_y2 <- inc_y2 + 1.96 * se_y2

# --- Overall ---
se_total <- sqrt(cases_total) / child_years_total
lower_total <- inc_overall_2yrs - 1.96 * se_total
upper_total <- inc_overall_2yrs + 1.96 * se_total

# Prevent negative lower bounds
lower_y1 <- max(0, lower_y1)
lower_y2 <- max(0, lower_y2)
lower_total <- max(0, lower_total)

# ======================================
# 3️⃣ Model dataframe
# ======================================

model <- data.frame(
  Period = c("Year 1", "Year 2", "Overall"),
  Incidence = c(inc_y1, inc_y2, inc_overall_2yrs),
  Lower = c(lower_y1, lower_y2, lower_total),
  Upper = c(upper_y1, upper_y2, upper_total),
  Source = "Model"
)

# ======================================
# 4️⃣ Combine observed + model
# ======================================

plot_data <- rbind(observed, model)

plot_data$Period <- factor(plot_data$Period,
                           levels = c("Year 1", "Year 2", "Overall"))

# ======================================
# 5️⃣ Bar Plot with CI
# ======================================

p2 <- ggplot(plot_data, aes(x = Period, y = Incidence, fill = Source)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.7),
           width = 0.6) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                position = position_dodge(width = 0.7),
                width = 0.2,
                linewidth = 1) +
  labs(
    y = "Incidence per child per year",
    x = "",
    title = "Model vs Observed Malaria Incidence",
    fill = ""
  ) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c("Observed" = "darkgreen",
                               "Model" = "deeppink")) +
  ylim(0, max(plot_data$Upper) * 1.1)



p2










#################now for RG



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
target_pfpr <- c(0.47, 0.431)  # desired PfPR at these days

#0.468

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
  clinical_incidence = c(0.5, 9) * 365,
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


#################now incidence per child per something 












# Output from your simulation
cases <- output_bednet3$n_inc_clinical_182.5_3284    # clinical episodes
population <- output_bednet3$n_age_182.5_3284        # children at risk
timesteps <- output_bednet3$timestep                 # days

# Define years
year_starts <- seq(1, max(timesteps), by = 365)
year_ends <- year_starts + 364
n_years <- length(year_starts)

# Calculate incidence per year
incidence_per_year <- numeric(n_years)

for(i in 1:n_years){
  keep <- timesteps >= year_starts[i] & timesteps <= year_ends[i]
  
  total_cases <- sum(cases[keep])
  total_child_years <- sum(population[keep]) / 365  # convert to child-years
  
  incidence_per_year[i] <- total_cases / total_child_years
}

# Show results
data.frame(
  Year = 1:n_years,
  Incidence_per_child_per_year = incidence_per_year
)




# Simulation output
cases <- output_bednet3$n_inc_clinical_182.5_3284
population <- output_bednet3$n_age_182.5_3284
timesteps <- output_bednet3$timestep

# Bed net distribution day
bednet_day <- 2310

# Define yearly windows from the bednet distribution
years <- 1:3  # number of years you want to calculate
year_length <- 365

incidence_per_year <- numeric(length(years))

for(i in seq_along(years)){
  start_day <- bednet_day + (i - 1) * year_length
  end_day   <- start_day + year_length - 1
  
  keep <- timesteps >= start_day & timesteps <= end_day
  
  total_cases <- sum(cases[keep])
  total_child_years <- sum(population[keep]) / year_length  # child-years
  
  incidence_per_year[i] <- total_cases / total_child_years
}

data.frame(
  Year_since_bednet = years,
  Incidence_per_child_per_year = incidence_per_year
)






#########################overall  incidence


# ======================================
# Extract simulation outputs
# ======================================

cases <- output_bednet3$n_inc_clinical_182.5_3284
population <- output_bednet3$n_age_182.5_3284
timesteps <- output_bednet3$timestep

year_length <- 365
bednet_day <- 2310

# ======================================
# Year 1 (days 2310–2674)
# ======================================

start_y1 <- bednet_day
end_y1   <- bednet_day + year_length - 1

keep_y1 <- timesteps >= start_y1 & timesteps <= end_y1

cases_y1 <- sum(cases[keep_y1])
child_years_y1 <- sum(population[keep_y1]) / year_length

inc_y1 <- cases_y1 / child_years_y1


# ======================================
# Year 2 (days 2675–3039)
# ======================================

start_y2 <- bednet_day + year_length
end_y2   <- start_y2 + year_length - 1

keep_y2 <- timesteps >= start_y2 & timesteps <= end_y2

cases_y2 <- sum(cases[keep_y2])
child_years_y2 <- sum(population[keep_y2]) / year_length

inc_y2 <- cases_y2 / child_years_y2


# ======================================
# Overall (Year 1 + Year 2 combined)
# ======================================

keep_y1_y2 <- timesteps >= start_y1 & timesteps <= end_y2

cases_total <- sum(cases[keep_y1_y2])
child_years_total <- sum(population[keep_y1_y2]) / year_length

inc_overall_2yrs <- cases_total / child_years_total


# ======================================
# Results table
# ======================================

results <- data.frame(
  Measure = c("Year 1 after bednet",
              "Year 2 after bednet",
              "Overall (Year 1 + Year 2)"),
  Incidence_per_child_per_year = c(inc_y1,
                                   inc_y2,
                                   inc_overall_2yrs),
  Incidence_per_1000_children_per_year = c(inc_y1,
                                           inc_y2,
                                           inc_overall_2yrs) * 1000
)

print(results)



##################################################################################plotting 




library(ggplot2)

# ======================================
# 1️⃣ Observed trial data
# ======================================

observed <- data.frame(
  Period = c("Year 1", "Year 2", "Overall"),
  Incidence = c(0.62, 0.98, 0.84),
  Lower = c(0.54, 0.90, 0.78),
  Upper = c(0.70, 1.07, 0.90),
  Source = "Observed"
)



# --- Year 1 ---
se_y1 <- sqrt(cases_y1) / child_years_y1
lower_y1 <- inc_y1 - 1.96 * se_y1
upper_y1 <- inc_y1 + 1.96 * se_y1

# --- Year 2 ---
se_y2 <- sqrt(cases_y2) / child_years_y2
lower_y2 <- inc_y2 - 1.96 * se_y2
upper_y2 <- inc_y2 + 1.96 * se_y2

# --- Overall ---
se_total <- sqrt(cases_total) / child_years_total
lower_total <- inc_overall_2yrs - 1.96 * se_total
upper_total <- inc_overall_2yrs + 1.96 * se_total

# Prevent negative lower bounds
lower_y1 <- max(0, lower_y1)
lower_y2 <- max(0, lower_y2)
lower_total <- max(0, lower_total)

# ======================================
# 3️⃣ Model dataframe
# ======================================

model <- data.frame(
  Period = c("Year 1", "Year 2", "Overall"),
  Incidence = c(inc_y1, inc_y2, inc_overall_2yrs),
  Lower = c(lower_y1, lower_y2, lower_total),
  Upper = c(upper_y1, upper_y2, upper_total),
  Source = "Model"
)

# ======================================
# 4️⃣ Combine observed + model
# ======================================

plot_data <- rbind(observed, model)

plot_data$Period <- factor(plot_data$Period,
                           levels = c("Year 1", "Year 2", "Overall"))

# ======================================
# 5️⃣ Bar Plot with CI
# ======================================

p3 <- ggplot(plot_data, aes(x = Period, y = Incidence, fill = Source)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.7),
           width = 0.6) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                position = position_dodge(width = 0.7),
                width = 0.2,
                linewidth = 1) +
  labs(
    y = "Incidence per child per year",
    x = "",
    title = "Model vs Observed Malaria Incidence",
    fill = ""
  ) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c("Observed" = "darkgreen",
                               "Model" = "deeppink")) +
  ylim(0, max(plot_data$Upper) * 1.1)



p3


################################combine all three together 



library(ggplot2)
library(dplyr)

############################################################
# 1️⃣ FUNCTION TO EXTRACT MODEL INCIDENCE
############################################################

extract_incidence <- function(output, net_name){
  
  cases <- output$n_inc_clinical_182.5_3284
  population <- output$n_age_182.5_3284
  timesteps <- output$timestep
  
  bednet_day <- 2310
  year_length <- 365
  
  # Year 1
  keep_y1 <- timesteps >= bednet_day &
    timesteps < bednet_day + year_length
  
  inc_y1 <- sum(cases[keep_y1]) /
    (sum(population[keep_y1]) / year_length)
  
  # Year 2
  keep_y2 <- timesteps >= bednet_day + year_length &
    timesteps < bednet_day + 2*year_length
  
  inc_y2 <- sum(cases[keep_y2]) /
    (sum(population[keep_y2]) / year_length)
  
  # Overall
  keep_total <- timesteps >= bednet_day &
    timesteps < bednet_day + 2*year_length
  
  inc_total <- sum(cases[keep_total]) /
    (sum(population[keep_total]) / year_length)
  
  data.frame(
    Net = net_name,
    Period = c("Year 1","Year 2","Overall"),
    Incidence = c(inc_y1, inc_y2, inc_total),
    Lower = c(inc_y1, inc_y2, inc_total),   # change if you have model CI
    Upper = c(inc_y1, inc_y2, inc_total),
    Source = "Model"
  )
}

############################################################
# 2️⃣ EXTRACT MODEL RESULTS
############################################################

model_all <- rbind(
  extract_incidence(output_bednet1, "IG2"),
  extract_incidence(output_bednet2, "IG1"),
  extract_incidence(output_bednet3, "RG")
)

############################################################
# 3️⃣ ADD OBSERVED TRIAL DATA (EDIT CI IF NEEDED)
############################################################

observed_all <- rbind(
  
  data.frame(
    Net = "IG2",
    Period = c("Year 1","Year 2","Overall"),
    Incidence = c(0.36,0.69,0.56),
    Lower = c(0.30,0.60,0.50),
    Upper = c(0.42,0.78,0.62),
    Source = "Observed"
  ),
  
  data.frame(
    Net = "IG1",
    Period = c("Year 1","Year 2","Overall"),
    Incidence = c(0.77,1.19,1.03),
    Lower = c(0.70,1.10,0.95),
    Upper = c(0.84,1.28,1.11),
    Source = "Observed"
  ),
  
  data.frame(
    Net = "RG",
    Period = c("Year 1","Year 2","Overall"),
    Incidence = c(0.62,0.98,0.84),
    Lower = c(0.55,0.90,0.77),
    Upper = c(0.69,1.06,0.91),
    Source = "Observed"
  )
)

############################################################
# 4️⃣ COMBINE EVERYTHING
############################################################

plot_data <- rbind(model_all, observed_all)

plot_data$Period <- factor(plot_data$Period,
                           levels = c("Year 1","Year 2","Overall"))

plot_data$Net <- factor(plot_data$Net,
                        levels = c("IG2","IG1","RG"))

############################################################
# 5️⃣ FINAL COMBINED PLOT
############################################################

ggplot(plot_data,
       aes(x = Period,
           y = Incidence,
           fill = Net,
           alpha = Source)) +
  
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.9),
           width = 0.7,
           colour = "black") +
  
  geom_errorbar(aes(ymin = Lower,
                    ymax = Upper),
                position = position_dodge(width = 0.9),
                width = 0.25,
                linewidth = 0.8) +
  
  scale_alpha_manual(values = c("Observed" = 1,
                                "Model" = 0.5)) +
  
  scale_fill_manual(values = c(
    "IG2" = "darkgreen",
    "IG1" = "darkblue",
    "RG"  = "darkred"
  )) +
  
  labs(
    y = "Incidence per child per year",
    x = "",
    fill = "Net type",
    alpha = "Source",
    title = "Observed vs Model Malaria Incidence"
  ) +
  
  theme_bw(base_size = 14)








library(ggplot2)
library(dplyr)

############################################################
# Combine model + observed (same as before)
############################################################

plot_data <- rbind(model_all, observed_all)

plot_data$Period <- factor(plot_data$Period,
                           levels = c("Year 1","Year 2","Overall"))

plot_data$Net <- factor(plot_data$Net,
                        levels = c("IG2","IG1","RG"))

############################################################
# Create combined x-axis variable
############################################################

plot_data$Net_Period <- interaction(plot_data$Net,
                                    plot_data$Period,
                                    sep = " - ")

############################################################
# ORDER AXIS EXACTLY HOW YOU WANT
############################################################

plot_data$Net_Period <- factor(
  plot_data$Net_Period,
  levels = c(
    "IG2 - Year 1","IG2 - Year 2","IG2 - Overall",
    "IG1 - Year 1","IG1 - Year 2","IG1 - Overall",
    "RG - Year 1","RG - Year 2","RG - Overall"
  )
)

############################################################
# FINAL PLOT
############################################################

ggplot(plot_data,
       aes(x = Net_Period,
           y = Incidence,
           fill = Source)) +
  
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.7,
           colour = "black") +
  
  geom_errorbar(aes(ymin = Lower,
                    ymax = Upper),
                position = position_dodge(width = 0.8),
                width = 0.25,
                linewidth = 0.8) +
  
  scale_fill_manual(values = c(
    "Observed" = "grey70",
    "Model" = "darkgreen"
  )) +
  
  labs(
    y = "Incidence per child per year",
    x = "",
    fill = "Sample",
    title = "Observed vs Predicted Incidence by Net Type"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )










library(ggplot2)
library(dplyr)

############################################################
# Combine data
############################################################

plot_data <- rbind(model_all, observed_all)

plot_data$Period <- factor(plot_data$Period,
                           levels = c("Year 1","Year 2","Overall"))

plot_data$Net <- factor(plot_data$Net,
                        levels = c("IG2","IG1","RG"))

############################################################
# FINAL PLOT (Grouped by Net using facets)
############################################################

ggplot(plot_data,
       aes(x = Period,
           y = Incidence,
           fill = Source)) +
  
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.7,
           colour = "black") +
  
  geom_errorbar(aes(ymin = Lower,
                    ymax = Upper),
                position = position_dodge(width = 0.8),
                width = 0.25,
                linewidth = 0.8) +
  
  scale_fill_manual(values = c(
    "Observed" = "grey70",
    "Model" = "darkgreen"
  )) +
  
  facet_wrap(~ Net, nrow = 1) +
  
  labs(
    y = "Incidence per child per year",
    x = "",
    fill = "Sample",
    title = "Observed vs Predicted Incidence by Net Type"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    panel.grid.major.x = element_blank()
  )




















library(ggplot2)
library(dplyr)

############################################################
# Combine data
############################################################

plot_data <- rbind(model_all, observed_all)

plot_data$Period <- factor(plot_data$Period,
                           levels = c("Year 1","Year 2","Overall"))

plot_data$Net <- factor(plot_data$Net,
                        levels = c("IG2","IG1","RG"))

############################################################
# Create ordered Net-Period axis
############################################################

plot_data$Net_Period <- interaction(plot_data$Net,
                                    plot_data$Period,
                                    sep = " - ")

plot_data$Net_Period <- factor(
  plot_data$Net_Period,
  levels = c(
    "Year 1","Year 2","Overall",
    "Year 1","Year 2","Overall",
    "Year 1","Year 2","Overall"
  )
)

############################################################
# FINAL SINGLE-PANEL PLOT
############################################################

ggplot(plot_data,
       aes(x = Net_Period,
           y = Incidence,
           fill = Source)) +
  
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.7,
           colour = "black") +
  
  geom_errorbar(aes(ymin = Lower,
                    ymax = Upper),
                position = position_dodge(width = 0.8),
                width = 0.25,
                linewidth = 0.8) +
  
  scale_fill_manual(values = c(
    "Observed" = "grey70",
    "Model" = "darkgreen"
  )) +
  
  labs(
    y = "Incidence per child per year",
    x = "",
    fill = "Sample",
    title = "Observed vs Predicted Incidence by Net Type"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )















library(ggplot2)
library(dplyr)

############################################################
# Combine data
############################################################

plot_data <- rbind(model_all, observed_all)

plot_data$Period <- factor(plot_data$Period,
                           levels = c("Year 1","Year 2","Overall"))

plot_data$Net <- factor(plot_data$Net,
                        levels = c("IG2","IG1","RG"))

############################################################
# Create numeric x positions for grouped layout
############################################################

plot_data <- plot_data %>%
  mutate(
    Period_num = as.numeric(Period),
    Net_num = as.numeric(Net),
    x_position = Period_num + (Net_num - 2) * 0.25
  )

############################################################
# Single panel grouped plot
############################################################

ggplot(plot_data,
       aes(x = x_position,
           y = Incidence,
           fill = Source)) +
  
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.2),
           width = 0.18,
           colour = "black") +
  
  geom_errorbar(aes(ymin = Lower,
                    ymax = Upper),
                position = position_dodge(width = 0.2),
                width = 0.05,
                linewidth = 0.8) +
  
  scale_x_continuous(
    breaks = 1:3,
    labels = c("Year 1","Year 2","Overall")
  ) +
  
  scale_fill_manual(values = c(
    "Observed" = "grey70",
    "Model" = "darkgreen"
  )) +
  
  labs(
    y = "Incidence per child per year",
    x = "",
    fill = "Sample",
    title = "Observed vs Predicted Incidence by Net Type"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank()
  ) +
  
  # Vertical separators between nets
  geom_vline(xintercept = c(1.5, 2.5),
             linetype = "dashed",
             colour = "grey40") +
  
  # Net labels on top
  annotate("text", x = 1, y = max(plot_data$Upper) * 1.08, label = "IG2", size = 5, fontface = "bold") +
  annotate("text", x = 2, y = max(plot_data$Upper) * 1.08, label = "IG1", size = 5, fontface = "bold") +
  annotate("text", x = 3, y = max(plot_data$Upper) * 1.08, label = "RG",  size = 5, fontface = "bold")





































library(ggplot2)
library(dplyr)

plot_data$Period <- factor(plot_data$Period,
                           levels = c("Year 1","Year 2","Overall"))

plot_data$Net <- factor(plot_data$Net,
                        levels = c("IG2","IG1","RG"))

plot_data$Source <- factor(plot_data$Source,
                           levels = c("Observed","Model"))

# Create combined grouping for Year within Net
plot_data <- plot_data %>%
  mutate(Net_Period = interaction(Net, Period, sep = "\n"))

# Plot
ggplot(plot_data,
       aes(x = Net_Period,
           y = Incidence,
           fill = Source)) +
  
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.7),
           width = 0.6,
           colour = "black") +
  
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                position = position_dodge(width = 0.7),
                width = 0.15,
                linewidth = 0.8) +
  
  scale_fill_manual(values = c("Observed" = "grey80",
                               "Model" = "darkgreen")) +
  
  labs(
    y = "Incidence per child per year",
    x = "",
    fill = "Sample",
    title = "Observed vs Predicted Incidence by Net Type"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 12),
    legend.position = "top"
  )







library(ggplot2)
library(dplyr)

plot_data$Period <- factor(plot_data$Period,
                           levels = c("Year 1","Year 2","Overall"))

plot_data$Net <- factor(plot_data$Net,
                        levels = c("IG2","IG1","RG"))

plot_data$Source <- factor(plot_data$Source,
                           levels = c("Observed","Model"))

# Create interaction for Year within Net
plot_data <- plot_data %>%
  mutate(Net_Period = interaction(Net, Period, sep = " - "))

ggplot(plot_data,
       aes(x = Net_Period,
           y = Incidence,
           fill = Source)) +
  
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.7),
           width = 0.6,
           colour = "black") +
  
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                position = position_dodge(width = 0.7),
                width = 0.15,
                linewidth = 0.8) +
  
  scale_fill_manual(values = c("Observed" = "grey80",
                               "Model" = "darkgreen")) +
  
  labs(
    y = "Incidence per child per year",
    x = "",
    fill = "Sample",
    title = "Observed vs Predicted Incidence by Net Type"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )



















library(ggplot2)
library(dplyr)

# -----------------------------
# Prepare factors
# -----------------------------
plot_data$Period <- factor(plot_data$Period, levels = c("Year 1","Year 2","Overall"))
plot_data$Net <- factor(plot_data$Net, levels = c("IG2","IG1","RG"))
plot_data$Source <- factor(plot_data$Source, levels = c("Observed","Model"))

# Keep interaction but reorder by Net first
plot_data <- plot_data %>%
  arrange(Net, Period, Source) %>%
  mutate(Net_Period = interaction(Period, Source, sep = " - "))

# -----------------------------
# Plot
# -----------------------------
ggplot(plot_data,
       aes(x = Period,
           y = Incidence,
           fill = Source)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.7),
           width = 0.6,
           colour = "black") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                position = position_dodge(width = 0.7),
                width = 0.15,
                linewidth = 0.8) +
  scale_fill_manual(values = c("Observed" = "grey80",
                               "Model" = "darkgreen")) +
  labs(
    y = "Incidence cases per child per year\n (children aged 6 months-9 years)",
    x = "",
    fill = "Sample",
    title = "Observed vs Predicted Incidence by Net Type"
  ) +
  facet_wrap(~Net, scales = "free_x", nrow = 1) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "top"
  )
