# Load necessary libraries
library(lme4)
library(MASS)
library(VGAM)
library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)
library(patchwork)

#############################IG1


funestus_unwashed <- read_excel("Malawi_EHT_data_gambiae.xlsx", 
                                sheet = "washed")

# Convert columns to appropriate types
funestus_unwashed$compartment <- as.factor(funestus_unwashed$compartment)
funestus_unwashed$Total <- as.numeric(funestus_unwashed$Total)

# Use dplyr to filter for 'UT' treatment and then fit the model
funestus_unwashed_UT <- funestus_unwashed %>%
  filter(treatment == "IG1")  # Filter only 'UT' treatment

# Fit the multinomial model for UT treatment
fit1.ms <- vglm(compartment ~ Total, multinomial(refLevel = 1), data = funestus_unwashed_UT)

# Check the model summary
summary(fit1.ms)

# Create prediction data based on the 'Total' variable
pred_data <- data.frame(Total = seq(min(funestus_unwashed_UT$Total), max(funestus_unwashed_UT$Total), length.out = 100))

# Generate predicted probabilities for the UT treatment
pred_probs <- predict(fit1.ms, newdata = pred_data, type = "response")

# Convert predictions to long format
pred_data_long <- as.data.frame(pred_probs)
long_pred_data <- pivot_longer(pred_data_long, 
                               cols = everything(), 
                               names_to = "Compartment", 
                               values_to = "Probability")

# Add the 'Total' variable back to the long format data
long_pred_data$Total <- rep(pred_data$Total, each = ncol(pred_data_long))

# Transform the 'Total' values into percentages
max_total <- max(long_pred_data$Total, na.rm = TRUE)  # Get the maximum Total value for scaling
long_pred_data$Total_percentage <- (long_pred_data$Total / max_total)  # Convert to percentage

# Plot the predicted probabilities
plot1 <- ggplot(long_pred_data, aes(x = Total_percentage, y = Probability, fill = Compartment, group = Compartment)) +
  geom_area(alpha = 0.4) +
  scale_fill_manual(values = c("red", "blue", "green")) +
  labs(
    x = "Total Mosquitoes Caught (%)", 
    y = "Predicted Probability",
    fill = "Compartment",
    title = "IG1"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add border around plot area
    axis.text = element_text(size = 12)  # Optional: Customize axis text size
  ) +
  scale_x_continuous(
    expand = c(0, 0), 
    breaks = seq(0, 100, by = 0.5)  # Set x-axis breaks in percentage
  ) +
  scale_y_continuous(expand = c(0, 0))    # Remove gap on y-axis





#############################IG2


funestus_unwashed <- read_excel("Malawi_EHT_data_gambiae.xlsx", 
                                sheet = "washed")

# Convert columns to appropriate types
funestus_unwashed$compartment <- as.factor(funestus_unwashed$compartment)
funestus_unwashed$Total <- as.numeric(funestus_unwashed$Total)

# Use dplyr to filter for 'UT' treatment and then fit the model
funestus_unwashed_UT <- funestus_unwashed %>%
  filter(treatment == "IG2")  # Filter only 'UT' treatment

# Fit the multinomial model for UT treatment
fit2.ms <- vglm(compartment ~ Total, multinomial(refLevel = 1), data = funestus_unwashed_UT)

# Check the model summary
summary(fit2.ms)

# Create prediction data based on the 'Total' variable
pred_data <- data.frame(Total = seq(min(funestus_unwashed_UT$Total), max(funestus_unwashed_UT$Total), length.out = 100))

# Generate predicted probabilities for the UT treatment
pred_probs <- predict(fit2.ms, newdata = pred_data, type = "response")

# Convert predictions to long format
pred_data_long <- as.data.frame(pred_probs)
long_pred_data <- pivot_longer(pred_data_long, 
                               cols = everything(), 
                               names_to = "Compartment", 
                               values_to = "Probability")

# Add the 'Total' variable back to the long format data
long_pred_data$Total <- rep(pred_data$Total, each = ncol(pred_data_long))

# Transform the 'Total' values into percentages
max_total <- max(long_pred_data$Total, na.rm = TRUE)  # Get the maximum Total value for scaling
long_pred_data$Total_percentage <- (long_pred_data$Total / max_total)  # Convert to percentage

# Plot the predicted probabilities
plot2 <- ggplot(long_pred_data, aes(x = Total_percentage, y = Probability, fill = Compartment, group = Compartment)) +
  geom_area(alpha = 0.4) +
  scale_fill_manual(values = c("red", "blue", "green")) +
  labs(
    x = "Total Mosquitoes Caught (%)", 
    y = "Predicted Probability",
    fill = "Compartment",
    title = "IG2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add border around plot area
    axis.text = element_text(size = 12)  # Optional: Customize axis text size
  ) +
  scale_x_continuous(
    expand = c(0, 0), 
    breaks = seq(0, 100, by = 0.5)  # Set x-axis breaks in percentage
  ) +
  scale_y_continuous(expand = c(0, 0))    # Remove gap on y-axis





#############################OP


funestus_unwashed <- read_excel("Malawi_EHT_data_gambiae.xlsx", 
                                sheet = "washed")

# Convert columns to appropriate types
funestus_unwashed$compartment <- as.factor(funestus_unwashed$compartment)
funestus_unwashed$Total <- as.numeric(funestus_unwashed$Total)

# Use dplyr to filter for 'UT' treatment and then fit the model
funestus_unwashed_UT <- funestus_unwashed %>%
  filter(treatment == "OP")  # Filter only 'UT' treatment

# Fit the multinomial model for UT treatment
fit3.ms <- vglm(compartment ~ Total, multinomial(refLevel = 1), data = funestus_unwashed_UT)

# Check the model summary
summary(fit3.ms)

# Create prediction data based on the 'Total' variable
pred_data <- data.frame(Total = seq(min(funestus_unwashed_UT$Total), max(funestus_unwashed_UT$Total), length.out = 100))

# Generate predicted probabilities for the UT treatment
pred_probs <- predict(fit3.ms, newdata = pred_data, type = "response")

# Convert predictions to long format
pred_data_long <- as.data.frame(pred_probs)
long_pred_data <- pivot_longer(pred_data_long, 
                               cols = everything(), 
                               names_to = "Compartment", 
                               values_to = "Probability")

# Add the 'Total' variable back to the long format data
long_pred_data$Total <- rep(pred_data$Total, each = ncol(pred_data_long))

# Transform the 'Total' values into percentages
max_total <- max(long_pred_data$Total, na.rm = TRUE)  # Get the maximum Total value for scaling
long_pred_data$Total_percentage <- (long_pred_data$Total / max_total)  # Convert to percentage

# Plot the predicted probabilities
plot3 <- ggplot(long_pred_data, aes(x = Total_percentage, y = Probability, fill = Compartment, group = Compartment)) +
  geom_area(alpha = 0.4) +
  scale_fill_manual(values = c("red", "blue", "green")) +
  labs(
    x = "Total Mosquitoes Caught (%)", 
    y = "Predicted Probability",
    fill = "Compartment",
    title = "OP"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add border around plot area
    axis.text = element_text(size = 12)  # Optional: Customize axis text size
  ) +
  scale_x_continuous(
    expand = c(0, 0), 
    breaks = seq(0, 100, by = 0.5)  # Set x-axis breaks in percentage
  ) +
  scale_y_continuous(expand = c(0, 0))    # Remove gap on y-axis






#############################P2


funestus_unwashed <- read_excel("Malawi_EHT_data_gambiae.xlsx", 
                                sheet = "washed")

# Convert columns to appropriate types
funestus_unwashed$compartment <- as.factor(funestus_unwashed$compartment)
funestus_unwashed$Total <- as.numeric(funestus_unwashed$Total)

# Use dplyr to filter for 'UT' treatment and then fit the model
funestus_unwashed_UT <- funestus_unwashed %>%
  filter(treatment == "P2")  # Filter only 'UT' treatment

# Fit the multinomial model for UT treatment
fit4.ms <- vglm(compartment ~ Total, multinomial(refLevel = 1), data = funestus_unwashed_UT)

# Check the model summary
summary(fit4.ms)

# Create prediction data based on the 'Total' variable
pred_data <- data.frame(Total = seq(min(funestus_unwashed_UT$Total), max(funestus_unwashed_UT$Total), length.out = 100))

# Generate predicted probabilities for the UT treatment
pred_probs <- predict(fit4.ms, newdata = pred_data, type = "response")

# Convert predictions to long format
pred_data_long <- as.data.frame(pred_probs)
long_pred_data <- pivot_longer(pred_data_long, 
                               cols = everything(), 
                               names_to = "Compartment", 
                               values_to = "Probability")

# Add the 'Total' variable back to the long format data
long_pred_data$Total <- rep(pred_data$Total, each = ncol(pred_data_long))

# Transform the 'Total' values into percentages
max_total <- max(long_pred_data$Total, na.rm = TRUE)  # Get the maximum Total value for scaling
long_pred_data$Total_percentage <- (long_pred_data$Total / max_total)  # Convert to percentage

# Plot the predicted probabilities
plot4 <- ggplot(long_pred_data, aes(x = Total_percentage, y = Probability, fill = Compartment, group = Compartment)) +
  geom_area(alpha = 0.4) +
  scale_fill_manual(values = c("red", "blue", "green")) +
  labs(
    x = "Total Mosquitoes caught (%)", 
    y = "Predicted Probability",
    fill = "Compartment",
    title = "P2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add border around plot area
    axis.text = element_text(size = 12)  # Optional: Customize axis text size
  ) +
  scale_x_continuous(
    expand = c(0, 0), 
    breaks = seq(0, 100, by = 0.5)  # Set x-axis breaks in percentage
  ) +
  scale_y_continuous(expand = c(0, 0))    # Remove gap on y-axis



#############################RG


funestus_unwashed <- read_excel("Malawi_EHT_data_gambiae.xlsx", 
                                sheet = "washed")

# Convert columns to appropriate types
funestus_unwashed$compartment <- as.factor(funestus_unwashed$compartment)
funestus_unwashed$Total <- as.numeric(funestus_unwashed$Total)

# Use dplyr to filter for 'UT' treatment and then fit the model
funestus_unwashed_UT <- funestus_unwashed %>%
  filter(treatment == "RG")  # Filter only 'UT' treatment

# Fit the multinomial model for UT treatment
fit5.ms <- vglm(compartment ~ Total, multinomial(refLevel = 1), data = funestus_unwashed_UT)

# Check the model summary
summary(fit5.ms)

# Create prediction data based on the 'Total' variable
pred_data <- data.frame(Total = seq(min(funestus_unwashed_UT$Total), max(funestus_unwashed_UT$Total), length.out = 100))

# Generate predicted probabilities for the UT treatment
pred_probs <- predict(fit5.ms, newdata = pred_data, type = "response")

# Convert predictions to long format
pred_data_long <- as.data.frame(pred_probs)
long_pred_data <- pivot_longer(pred_data_long, 
                               cols = everything(), 
                               names_to = "Compartment", 
                               values_to = "Probability")

# Add the 'Total' variable back to the long format data
long_pred_data$Total <- rep(pred_data$Total, each = ncol(pred_data_long))

# Transform the 'Total' values into percentages
max_total <- max(long_pred_data$Total, na.rm = TRUE)  # Get the maximum Total value for scaling
long_pred_data$Total_percentage <- (long_pred_data$Total / max_total)  # Convert to percentage

# Plot the predicted probabilities
plot5 <- ggplot(long_pred_data, aes(x = Total_percentage, y = Probability, fill = Compartment, group = Compartment)) +
  geom_area(alpha = 0.4) +
  scale_fill_manual(values = c("red", "blue", "green")) +
  labs(
    x = "Total Mosquitoes Caught (%)", 
    y = "Predicted Probability",
    fill = "Compartment",
    title = "RG"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add border around plot area
    axis.text = element_text(size = 12)  # Optional: Customize axis text size
  ) +
  scale_x_continuous(
    expand = c(0, 0), 
    breaks = seq(0, 100, by = 0.5)  # Set x-axis breaks in percentage
  ) +
  scale_y_continuous(expand = c(0, 0))    # Remove gap on y-axis







####################################UT



funestus_unwashed <- read_excel("Malawi_EHT_data_gambiae.xlsx", 
                                sheet = "washed")

# Convert columns to appropriate types
funestus_unwashed$compartment <- as.factor(funestus_unwashed$compartment)
funestus_unwashed$Total <- as.numeric(funestus_unwashed$Total)

# Use dplyr to filter for 'UT' treatment and then fit the model
funestus_unwashed_UT <- funestus_unwashed %>%
  filter(treatment == "UT")  # Filter only 'UT' treatment

# Fit the multinomial model for UT treatment
fit6.ms <- vglm(compartment ~ Total, multinomial(refLevel = 1), data = funestus_unwashed_UT)

# Check the model summary
summary(fit6.ms)

# Create prediction data based on the 'Total' variable
pred_data <- data.frame(Total = seq(min(funestus_unwashed_UT$Total), max(funestus_unwashed_UT$Total), length.out = 100))

# Generate predicted probabilities for the UT treatment
pred_probs <- predict(fit6.ms, newdata = pred_data, type = "response")

# Convert predictions to long format
pred_data_long <- as.data.frame(pred_probs)
long_pred_data <- pivot_longer(pred_data_long, 
                               cols = everything(), 
                               names_to = "Compartment", 
                               values_to = "Probability")

# Add the 'Total' variable back to the long format data
long_pred_data$Total <- rep(pred_data$Total, each = ncol(pred_data_long))

# Transform the 'Total' values into percentages
max_total <- max(long_pred_data$Total, na.rm = TRUE)  # Get the maximum Total value for scaling
long_pred_data$Total_percentage <- (long_pred_data$Total / max_total)  # Convert to percentage

# Plot the predicted probabilities
plot6 <- ggplot(long_pred_data, aes(x = Total_percentage, y = Probability, fill = Compartment, group = Compartment)) +
  geom_area(alpha = 0.4) +
  scale_fill_manual(values = c("red", "blue", "green")) +
  labs(
    x = "Total Mosquitoes Caught (%)", 
    y = "Predicted Probability",
    fill = "Compartment",
    title = "UT"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add border around plot area
    axis.text = element_text(size = 12)  # Optional: Customize axis text size
  ) +
  scale_x_continuous(
    expand = c(0, 0), 
    breaks = seq(0, 100, by = 0.5)  # Set x-axis breaks in percentage
  ) +
  scale_y_continuous(expand = c(0, 0))    # Remove gap on y-axis

##################combine plots


library(patchwork)# 



combined_plot <- plot1 + plot2 + plot3 + plot4 + plot5 +plot6+
  plot_layout(ncol = 3) + 
  plot_annotation(
    title = "gambiae caught per compartment (washed)",  # Add title for the entire plot
    theme = theme(
      plot.background = element_rect(color = "black", size = 1),  # Add border around the combined plot
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold")  # Customize title appearance
    )
  )

# Display the combined plot
print(combined_plot)


