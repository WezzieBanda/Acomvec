

####################



library(tidyverse)
library(ggplot2)
library(viridis)

# 1. Sample data grid
df <- expand.grid(
  time = seq(0, 3, by = 0.05),           
  cases_averted = seq(0, 1, by = 0.05),  
  group = c("IG1", "IG2", "P3", "RG"),   
  deployment = c("New", "12M", "24M", "36M")
)

# 2. Updated vaccine efficacy (hot at time 0)
df$vaccine_efficacy <- with(df,
                            pmin(
                              100,
                              100 *
                                (1 - 0.7 * (time / 3)) *           # linear decay over time
                                (1 - cases_averted^1.2) *
                                c(1.0, 1.1, 0.95, 0.9)[match(group, c("IG1","IG2","P3","RG"))] *  # IG2 higher
                                c(1.0, 0.95, 0.9, 0.85)[match(deployment, c("New","12M","24M","36M"))]  # boosts all deployments
                            )
)

# 3. Factor ordering (reverse so New is bottom)
df <- df %>%
  mutate(
    group = factor(group, levels = c("IG1", "IG2", "P3", "RG")),
    deployment = factor(deployment, levels = c("36M", "24M", "12M", "New"))
  )

# 4. Plot
p <- ggplot(df, aes(x = time, y = cases_averted, fill = vaccine_efficacy)) +
  geom_tile() +
  
  facet_grid(
    deployment ~ group,
    switch = "y"
  ) +
  
  scale_x_continuous(
    name = "Time (years)",
    limits = c(0, 3),
    breaks = seq(0, 3, by = 1),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    name = "Proportion of cases averted",
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.3),
    expand = c(0, 0)
  ) +
  
  scale_fill_viridis_c(
    name = "Cases Averted (%)",
    limits = c(0, 100)
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    strip.placement = "outside",
    axis.title.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    
    # Vertical space between rows
    panel.spacing.y = unit(0.05, "cm"),
    
    # Horizontal space between columns
    panel.spacing.x = unit(0.3, "cm")
  )

# 5. Save & display
#ggsave("heatmap_facet_spacing_xy.png", p, width = 12, height = 9, dpi = 300)
print(p)






















###################################




library(tidyverse)
library(ggplot2)
library(viridis)

# 1. Sample data grid
df <- expand.grid(
  time = seq(0, 3, by = 0.05),           
  cases_averted = seq(0, 1, by = 0.05),  
  group = c("IG1", "IG2", "P3", "RG"),   
  deployment = c("New", "12M", "24M", "36M")
)

# 2. Updated vaccine efficacy (hot at time 0)
df$vaccine_efficacy <- with(df,
                            pmin(
                              100,
                              100 *
                                (1 - 0.7 * (time / 3)) *           
                                (1 - cases_averted^1.2) *
                                c(1.0, 1.1, 0.95, 0.9)[match(group, c("IG1","IG2","P3","RG"))] *
                                c(1.0, 0.95, 0.9, 0.85)[match(deployment, c("New","12M","24M","36M"))]
                            )
)

# 3. Factor ordering (reverse so New is bottom)
df <- df %>%
  mutate(
    group = factor(group, levels = c("IG1", "IG2", "P3", "RG")),
    deployment = factor(deployment, levels = c("36M", "24M", "12M", "New"))
  )

# 4. Plot
p <- ggplot(df, aes(x = time, y = cases_averted, fill = vaccine_efficacy)) +
  geom_tile() +
  
  # Add L-shaped axes inside each facet
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  geom_vline(xintercept = 0, color = "black", size = 0.5) +
  
  facet_grid(
    deployment ~ group,
    switch = "y"
  ) +
  
  scale_x_continuous(
    name = "Time (years)",
    limits = c(0, 3),
    breaks = seq(0, 3, by = 1),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    name = "Proportion of cases averted",
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.3),
    expand = c(0, 0)
  ) +
  
  scale_fill_viridis_c(
    name = "Vaccine efficacy (%)",
    limits = c(0, 100)
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    strip.placement = "outside",
    axis.title.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    panel.spacing.y = unit(0.05, "cm"),
    panel.spacing.x = unit(0.3, "cm")
  )

# 5. Save & display
ggsave("heatmap_facet_L_axes.png", p, width = 12, height = 9, dpi = 300)
print(p)








######################################################################################


p <- ggplot(df, aes(x = time, y = cases_averted, fill = vaccine_efficacy)) +
  geom_tile() +
  
  # L-shaped axes lines
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  geom_vline(xintercept = 0, color = "black", size = 0.5) +
  
  facet_grid(
    deployment ~ group,
    switch = "y"
  ) +
  
  scale_x_continuous(
    name = "Time (years)",
    limits = c(0, 3),
    breaks = seq(0, 3, by = 1),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    name = "Proportion of cases averted",
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.3),
    expand = c(0, 0)
  ) +
  
  scale_fill_viridis_c(
    name = "Vaccine efficacy (%)",
    limits = c(0, 100)
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    strip.placement = "outside",
    axis.title.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", size = 0.5),      # <-- Tick lines
    axis.ticks.length = unit(0.2, "cm"),                         # <-- Tick length
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    panel.spacing.y = unit(0.05, "cm"),
    panel.spacing.x = unit(0.3, "cm")
  )

#ggsave("heatmap_with_ticks.png", p, width = 12, height = 9, dpi = 300)
print(p)


#########################################################################


###############################
# Faceted Heatmap with L-shaped Axes & Ticks
###############################

library(tidyverse)
library(ggplot2)
library(viridis)
library(grid)   # for unit()

# 1. Sample data grid
df <- expand.grid(
  time = seq(0, 3, by = 0.05),           # Time in years
  cases_averted = seq(0, 1, by = 0.05),  # Proportion of cases averted
  group = c("IG1", "IG2", "P3", "RG"),   
  deployment = c("New", "12M", "24M", "36M")
)

# 2. Vaccine efficacy calculation
df$vaccine_efficacy <- with(df,
                            pmin(
                              100,
                              100 *
                                (1 - 0.7 * (time / 3)) *           # linear decay
                                (1 - cases_averted^1.2) *
                                c(1.0, 1.1, 0.95, 0.9)[match(group, c("IG1","IG2","P3","RG"))] *  # group effect
                                c(1.0, 0.95, 0.9, 0.85)[match(deployment, c("New","12M","24M","36M"))]  # deployment effect
                            )
)

# 3. Factor ordering for facets
df <- df %>%
  mutate(
    group = factor(group, levels = c("IG1", "IG2", "P3", "RG")),
    deployment = factor(deployment, levels = c("36M", "24M", "12M", "New"))
  )

# 4. Plot
p <- ggplot(df, aes(x = time, y = cases_averted, fill = vaccine_efficacy)) +
  geom_tile() +
  
  # L-shaped axes at 0
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  geom_vline(xintercept = 0, color = "black", size = 0.5) +
  
  # Faceting
  facet_grid(
    deployment ~ group,
    switch = "y"
  ) +
  
  # X-axis scale
  scale_x_continuous(
    name = "Time (years)",
    limits = c(0, 3),
    breaks = seq(0, 3, by = 1),
    expand = c(0, 0)
  ) +
  
  # Y-axis scale
  scale_y_continuous(
    name = "Proportion of cases averted",
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.3),
    expand = c(0, 0)
  ) +
  
  # Fill (vaccine efficacy)
  scale_fill_viridis_c(
    name = "Cases Averted (%)",
    limits = c(0, 100)
  ) +
  
  # Theme settings
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    strip.placement = "outside",
    axis.title.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", size = 0.5),      # Tick lines visible
    axis.ticks.length = unit(0.2, "cm"),                         # Tick size
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    panel.spacing.y = unit(0.05, "cm"),
    panel.spacing.x = unit(0.3, "cm")
  )

# 5. Save and display
#ggsave("heatmap_with_ticks.png", p, width = 12, height = 9, dpi = 300)
print(p)








#############################################two plots





library(tidyverse)
library(ggplot2)
library(viridis)
library(grid)   # for unit()

# 1. Sample data grid
df <- expand.grid(
  time = seq(0, 3, by = 0.05),           
  cases_averted = seq(0, 1, by = 0.05),  
  group = c("IG1", "IG2", "P3", "RG"),   
  deployment = c("New", "12M", "24M", "36M")
)

# 2. Vaccine efficacy calculation
df$vaccine_efficacy <- with(df,
                            pmin(
                              100,
                              100 *
                                (1 - 0.7 * (time / 3)) *           
                                (1 - cases_averted^1.2) *
                                c(1.0, 1.1, 0.95, 0.9)[match(group, c("IG1","IG2","P3","RG"))] *
                                c(1.0, 0.95, 0.9, 0.85)[match(deployment, c("New","12M","24M","36M"))]
                            )
)

# 3. Factor ordering for facets
df <- df %>%
  mutate(
    group = factor(group, levels = c("IG1", "IG2", "P3", "RG")),
    deployment = factor(deployment, levels = c("36M", "24M", "12M", "New"))
  )

# 4. Define max infectious mosquitoes (example scaling)
# Here we map cases_averted (0–1) to infectious mosquitoes (0–1000)
df$infectious_mosq <- df$cases_averted * 1000  # linear mapping

# 5. Plot with secondary y-axis
p <- ggplot(df, aes(x = time, y = cases_averted, fill = vaccine_efficacy)) +
  geom_tile() +
  
  # L-shaped axes at 0
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  geom_vline(xintercept = 0, color = "black", size = 0.5) +
  
  # Faceting
  facet_grid(deployment ~ group, switch = "y") +
  
  # X-axis
  scale_x_continuous(
    name = "Time (years)",
    limits = c(0, 3),
    breaks = seq(0, 3, by = 1),
    expand = c(0, 0)
  ) +
  
  # Y-axis with secondary axis
  scale_y_continuous(
    name = "Proportion of cases averted",
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = c(0, 0),
    sec.axis = sec_axis(~ . * 1000, name = "Infectious mosquitoes")  # linear scaling
  ) +
  
  # Fill
  scale_fill_viridis_c(
    name = "Vaccine efficacy (%)",
    limits = c(0, 100)
  ) +
  
  # Theme
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    strip.placement = "outside",
    axis.title.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    panel.spacing.y = unit(0.05, "cm"),
    panel.spacing.x = unit(0.3, "cm")
  )

print(p)





###########33333333


library(tidyverse)
library(ggplot2)
library(viridis)
library(grid)   # for unit()

# 1. Sample data grid
df <- expand.grid(
  time = seq(0, 3, by = 0.05),           
  cases_averted = seq(0, 1, by = 0.05),  
  group = c("IG1", "IG2", "P3", "RG"),   
  deployment = c("New", "12M", "24M", "36M")
)

# 2. Vaccine efficacy calculation
df$vaccine_efficacy <- with(df,
                            pmin(
                              100,
                              100 *
                                (1 - 0.7 * (time / 3)) *           
                                (1 - cases_averted^1.2) *
                                c(1.0, 1.1, 0.95, 0.9)[match(group, c("IG1","IG2","P3","RG"))] *
                                c(1.0, 0.95, 0.9, 0.85)[match(deployment, c("New","12M","24M","36M"))]
                            )
)

# 3. Factor ordering for facets
df <- df %>%
  mutate(
    group = factor(group, levels = c("IG1", "IG2", "P3", "RG")),
    deployment = factor(deployment, levels = c("36M", "24M", "12M", "New"))
  )

# 4. Single line data for proportion infectious mosquitoes
line_df <- data.frame(
  time = seq(0, 3, by = 0.05),
  prop_infectious = seq(1, 0.2, length.out = length(seq(0, 3, by = 0.05)))  # single line
)

# 5. Plot
p <- ggplot(df, aes(x = time, y = cases_averted, fill = vaccine_efficacy)) +
  geom_tile() +
  
  # Single line plot (override default aes)
  geom_line(data = line_df, aes(x = time, y = prop_infectious),
            inherit.aes = FALSE, color = "red", size = 1) +
  
  # L-shaped axes
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  geom_vline(xintercept = 0, color = "black", size = 0.5) +
  
  # Faceting
  facet_grid(deployment ~ group, switch = "y") +
  
  # X-axis
  scale_x_continuous(
    name = "Time (years)",
    limits = c(0, 3),
    breaks = seq(0, 3, by = 1),
    expand = c(0, 0)
  ) +
  
  # Primary Y-axis with secondary axis
  scale_y_continuous(
    name = "Proportion of cases averted",
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.3),
    expand = c(0, 0),
    sec.axis = sec_axis(~ ., name = "Proportion of infectious mosquitoes", breaks = seq(0, 1, by = 0.3))
  ) +
  
  # Fill
  scale_fill_viridis_c(
    name = "Vaccine efficacy (%)",
    limits = c(0, 100)
  ) +
  
  # Theme
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    strip.placement = "outside",
    axis.title.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    panel.spacing.y = unit(0.05, "cm"),
    panel.spacing.x = unit(0.3, "cm")
  )

print(p)
