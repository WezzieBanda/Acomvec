# Load necessary package
library(dplyr)
library(readr)
library(ggplot2)

# Assuming your data is imported as `data`
# Replace "your_file.csv" with your actual CSV file name
data <- read.csv("Wezzie-Tom/EHT DATA/Malawi_EHT_data_all.csv")

# Summarize the data
species_counts <- data %>%
  group_by(compartment, species) %>%
  summarise(
    total_bf_dead = sum(bf_dead, na.rm = TRUE),
    total_bf_live = sum(bf_live, na.rm = TRUE),
    total_gravid_dead = sum(gravid_dead, na.rm = TRUE),
    total_gravid_live = sum(gravid_live, na.rm = TRUE),
    total_unf_dead = sum(unf_dead, na.rm = TRUE),
    total_unf_live = sum(unf_live, na.rm = TRUE),
    total_caught = sum(bf_dead, na.rm = TRUE) + 
      sum(bf_live, na.rm = TRUE) + 
      sum(gravid_dead, na.rm = TRUE) + 
      sum(gravid_live, na.rm = TRUE) + 
      sum(unf_dead, na.rm = TRUE) + 
      sum(unf_live, na.rm = TRUE)
  ) %>%
  ungroup() # Remove grouping for further operations

# View the results
print(species_counts)



library(MASS)
library(tidyr)

# Poisson regression model to predict total_caught based on compartment
glm_poisson <- glm(total_caught ~ compartment + species, family = poisson(), data = species_counts)

# Check the model summary
summary(glm_poisson)





# Extract coefficients from the regression output
intercept <- 6.17531
beta_room <- 0.17778
beta_veranda <- 0.50178
beta_gambiae <- 0.37043

# Calculate predicted log counts for each combination
log_counts <- data.frame(
  compartment = c("Net", "Room", "Veranda", "Net", "Room", "Veranda"),
  species = c("funestus", "funestus", "funestus", "gambiae", "gambiae", "gambiae"),
  log_count = c(
    intercept, 
    intercept + beta_room, 
    intercept + beta_veranda, 
    intercept + beta_gambiae, 
    intercept + beta_room + beta_gambiae, 
    intercept + beta_veranda + beta_gambiae
  )
)

# Exponentiate to get total counts
log_counts$total_count <- exp(log_counts$log_count)

# View the results
print(log_counts)


#######################now plot######################################





ggplot(log_counts, aes(x = interaction(compartment, species), y = total_count, fill = species)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_text(aes(label = round(total_count, 1)), 
            position = position_dodge(width = 0.95), vjust = -0.2, size = 3) +
  labs(
    title = ".",
    x = "Compartments",
    y = "Total Mosquitoes",
    fill = "Species"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_line(color = "gray", linetype = "dotted")
  )

library(readxl)

library(readr)
Book1 <- read_csv("C:/Users/user/Desktop/Book1.csv")

library(ggplot2)

# Assuming your dataset is named `data2` and contains the columns `day_continuous`, `total_caught`, and `compartment`
ggplot(Book1, aes(x = day_continuous, y = total_caught, color = compartment)) +
  geom_point(size = 2, alpha = 0.6) +  # Dot plot with transparency
  labs(
    title = "Dot Plot of Total Caught by Day (Continuous) and Compartment",
    x = "Day (Continuous)",
    y = "Total Caught",
    color = "Compartment"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )


