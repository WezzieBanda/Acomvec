library(lme4)
library(devtools)
library(MASS)
library(RColorBrewer)
library(cowplot)
library(ggplot2)
library(dplyr)
library(reshape2)
library(fitdistrplus)
library(lattice)
library(performance)
library(gtsummary)
library(readxl)
library(patchwork)
library(gridExtra)




########3import data 
#Malawi_EHT_data_funestus <- read_excel("Wezzie-Tom/EHT DATA/Malawi_EHT_data_funestus.xlsx")

##filter data to unwashed  and unwashed

Malawi_EHT_data_funestus <- read_excel("Malawi_EHT_data_funestus.xlsx")



# Filter data for washed and unwashed nets
funestus_washed <- Malawi_EHT_data_funestus %>%
  filter(net_status == "Washed")

funestus_unwashed <- Malawi_EHT_data_funestus %>%
  filter(net_status == "Unwashed")

# Aggregate Total mosquitoes per treatment for washed nets (sum across all rows)
funestus_washed_summary <- funestus_washed %>%
  group_by(treatment) %>%
  summarise(Total_sum = sum(Total, na.rm = TRUE))

# Aggregate Total mosquitoes per treatment for unwashed nets (sum across all rows)
funestus_unwashed_summary <- funestus_unwashed %>%
  group_by(treatment) %>%
  summarise(Total_sum = sum(Total, na.rm = TRUE))

# Plot for Washed Nets (total counts per treatment) with borders around the plot
plot_washed <- ggplot(funestus_washed_summary, aes(x = treatment, y = Total_sum, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8, color = "black") +  # Added border color
  labs(
    title = "Washed Nets",
    x = "Treatment",
    y = "Total Mosquitoes"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, size = 1))  # Add border around the plot

# Plot for Unwashed Nets (total counts per treatment) with borders around the plot
plot_unwashed <- ggplot(funestus_unwashed_summary, aes(x = treatment, y = Total_sum, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8, color = "black") +  # Added border color
  labs(
    title = "Unwashed Nets",
    x = "Treatment",
    y = "Total Mosquitoes"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, size = 1))  # Add border around the plot

# Get the y-axis limits for both plots to use a consistent scale
y_max <- max(c(max(funestus_washed_summary$Total_sum, na.rm = TRUE), max(funestus_unwashed_summary$Total_sum, na.rm = TRUE)))

# Apply consistent y-axis scale to both plots
plot_washed <- plot_washed + ylim(0, y_max)
plot_unwashed <- plot_unwashed + ylim(0, y_max)

# Arrange the plots side by side with one title
grid.arrange(plot_washed, plot_unwashed, ncol = 2, top = "Total Mosquitoes Caught for Washed vs Unwashed Nets")





####################combine to plot one graph of total mosquitoescaught##############3

# Combine washed and unwashed nets into one dataset
funestus_combined <- bind_rows(
  funestus_washed %>%
    mutate(net_status = "Washed"),
  funestus_unwashed %>%
    mutate(net_status = "Unwashed")
)

# Aggregate Total mosquitoes per treatment for each net status (washed/unwashed)
funestus_summary <- funestus_combined %>%
  group_by(treatment, net_status) %>%
  summarise(Total_sum = sum(Total, na.rm = TRUE))

# Plot for both Washed and Unwashed Nets in one graph
plot_combined <- ggplot(funestus_summary, aes(x = treatment, y = Total_sum, fill = net_status)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +  # Dodge for side-by-side bars
  geom_line(aes(group = net_status, color = net_status), 
            size = 1, 
            linetype = "dashed",  # Dotted line (dashed)
            position = position_dodge(width = 0.7)) +  # Line on top of bars
  labs(
    title = "funestus",
    x = "Treatment",
    y = "Total Mosquitoes"
  ) +
  scale_fill_manual(values = c("Washed" = "lightgray", "Unwashed" = "darkgray")) +  # Different shades of gray
  scale_color_manual(values = c("Washed" = "blue", "Unwashed" = "blue")) +  # Line color is blue for both
  theme_minimal() +
  theme(
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Border around the plot
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10),  # Space around the plot (top, right, bottom, left)
    plot.title = element_text(size = 18, face = "bold"),  # Bigger and bold title
    axis.title = element_text(size = 14),  # Bigger axis titles
    axis.text = element_text(size = 12)   # Bigger axis labels
  )

# Display the plot
plot_combined




########unf_live###############
# Combine washed and unwashed nets into one dataset
funestus_combined_unfed <- bind_rows(
  funestus_washed %>%
    mutate(net_status = "Washed"),
  funestus_unwashed %>%
    mutate(net_status = "Unwashed")
)

# Aggregate 'unf_live' (unfed live mosquitoes) for each 'treatment' and 'net_status' group
funestus_summary_unfed <- funestus_combined_unfed %>%
  filter(treatment == "RG") %>%  # Filter to include only 'IG1' treatment
  group_by(treatment, net_status) %>%
  summarise(Unfed_Live_sum = sum(unf_live, na.rm = TRUE))  # Sum 'unf_live' per treatment and net status

# View the summary
print(funestus_summary_unfed)





# Combine washed and unwashed nets into one dataset
funestus_combined_unfed <- bind_rows(
  funestus_washed %>%
    mutate(net_status = "Washed"),
  funestus_unwashed %>%
    mutate(net_status = "Unwashed")
)

# Aggregate `unf_live` per treatment for each net status (washed/unwashed)
funestus_summary_unfed <- funestus_combined_unfed %>%
  group_by(treatment, net_status) %>%  # Group by treatment and net status
  summarise(Unfed_Live_sum = sum(unf_live, na.rm = TRUE))  # Sum of unf_live

# Plot for both Washed and Unwashed Nets in one graph
plot_combined_unfed <- ggplot(funestus_summary_unfed, aes(x = treatment, y = Unfed_Live_sum, fill = net_status)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +  # Dodge for side-by-side bars
  geom_line(aes(group = net_status, color = net_status), 
            size = 1, 
            linetype = "dashed",  # Dotted line (dashed)
            position = position_dodge(width = 0.7)) +  # Line on top of bars
  labs(
    title = "funestus",
    x = "Treatment",
    y = "Unfed Live Mosquitoes"
  ) +
  scale_fill_manual(values = c("Washed" = "lightgray", "Unwashed" = "darkgray")) +  # Different shades of gray
  scale_color_manual(values = c("Washed" = "blue", "Unwashed" = "blue")) +  # Line color is blue for both
  scale_y_continuous(limits = c(0, 200)) +  # Set y-axis limits from 0 to 200
  theme_minimal() +
  theme(
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Border around the plot
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10),  # Space around the plot (top, right, bottom, left)
    plot.title = element_text(size = 18, face = "bold"),  # Bigger and bold title
    axis.title = element_text(size = 14),  # Bigger axis titles
    axis.text = element_text(size = 12)   # Bigger axis labels
  )

# Display the plot
plot_combined_unfed










#################bf_live mosquitoes

# Combine washed and unwashed nets into one dataset
funestus_combined_bf <- bind_rows(
  funestus_washed %>%
    mutate(net_status = "Washed"),
  funestus_unwashed %>%
    mutate(net_status = "Unwashed")
)

# Aggregate `bf_live` per treatment for each net status (washed/unwashed)
funestus_summary_bf <- funestus_combined_bf %>%
  group_by(treatment, net_status) %>%  # Group by treatment and net status
  summarise(Bf_Live_sum = sum(bf_live, na.rm = TRUE))  # Sum of bf_live

# Plot for both Washed and Unwashed Nets in one graph for bf_live
plot_combined_bf <- ggplot(funestus_summary_bf, aes(x = treatment, y = Bf_Live_sum, fill = net_status)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +  # Dodge for side-by-side bars
  geom_line(aes(group = net_status, color = net_status), 
            size = 1, 
            linetype = "dashed",  # Dotted line (dashed)
            position = position_dodge(width = 0.7)) +  # Line on top of bars
  labs(
    title = "funestus",
    x = "Treatment",
    y = "bloof fed live Mosquitoes"
  ) +
  scale_fill_manual(values = c("Washed" = "lightgray", "Unwashed" = "darkgray")) +  # Different shades of gray
  scale_color_manual(values = c("Washed" = "blue", "Unwashed" = "blue")) +  # Line color is blue for both
  scale_y_continuous(limits = c(0, 200)) +  # Set y-axis limits from 0 to 200
  theme_minimal() +
  theme(
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Border around the plot
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10),  # Space around the plot (top, right, bottom, left)
    plot.title = element_text(size = 18, face = "bold"),  # Bigger and bold title
    axis.title = element_text(size = 14),  # Bigger axis titles
    axis.text = element_text(size = 12)   # Bigger axis labels
  )

# Display the plot
plot_combined_bf




##############################mtot_dead#####################################################################3




# Combine washed and unwashed nets into one dataset
funestus_combined_tot_dead <- bind_rows(
  funestus_washed %>%
    mutate(net_status = "Washed"),
  funestus_unwashed %>%
    mutate(net_status = "Unwashed")
)

# Aggregate `tot_dead` per treatment for each net status (washed/unwashed)
funestus_summary_tot_dead <- funestus_combined_tot_dead %>%
  group_by(treatment, net_status) %>%  # Group by treatment and net status
  summarise(Tot_Dead_sum = sum(tot_dead, na.rm = TRUE))  # Sum of tot_dead

# Plot for both Washed and Unwashed Nets in one graph for tot_dead
plot_combined_tot_dead <- ggplot(funestus_summary_tot_dead, aes(x = treatment, y = Tot_Dead_sum, fill = net_status)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +  # Dodge for side-by-side bars
  geom_line(aes(group = net_status, color = net_status), 
            size = 1, 
            linetype = "dashed",  # Dotted line (dashed)
            position = position_dodge(width = 0.7)) +  # Line on top of bars
  labs(
    title = "funestus",
    x = "Treatment",
    y = "Total Dead Mosquitoes"
  ) +
  scale_fill_manual(values = c("Washed" = "lightgray", "Unwashed" = "darkgray")) +  # Different shades of gray
  scale_color_manual(values = c("Washed" = "blue", "Unwashed" = "blue")) +  # Line color is blue for both
  scale_y_continuous(limits = c(0, 200)) +  # Set y-axis limits from 0 to 200
  theme_minimal() +
  theme(
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Border around the plot
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10),  # Space around the plot (top, right, bottom, left)
    plot.title = element_text(size = 18, face = "bold"),  # Bigger and bold title
    axis.title = element_text(size = 14),  # Bigger axis titles
    axis.text = element_text(size = 12)   # Bigger axis labels
  )

# Display the plot
plot_combined_tot_dead




# Aggregate 'unf_live' (unfed live mosquitoes) for each 'treatment' and 'net_status' group
funestus_summary_unfed <- Malawi_EHT_data_funestus %>%
  filter(treatment == "UT") %>%  # Filter to include only 'IG1' treatment
  group_by(treatment, net_status) %>%
  summarise(Unfed_Live_sum = sum(tot_dead, na.rm = TRUE))  # Sum 'unf_live' per treatment and net status

# View the summary
print(funestus_summary_unfed)


