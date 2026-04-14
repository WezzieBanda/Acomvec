library(dplyr)
library(rstan)
library(readxl)
library(ggplot2)


mydata <- read_excel("benin trial data results/Benin_NNP_hut_trial_raw data (2).xlsx", 
                 sheet = "Combined2")

mydata <- mydata1 %>% 
  filter(Treatment == "IG2")

stan_data <- list(
  S = nrow(mydata),
  X_survive = mydata$total - mydata$tot_dead,
  N_caught = mydata$total,
  X_succfed = mydata$bf_live,
  site = as.numeric(as.factor(mydata$hut)), 
  S_sim = 100,
  theta_sim = seq(0.01, 0.99, length.out = 100)
)

model <- stan_model("feeding_survival_model_simple.stan")
model <- stan_model("benin trial data results/stan/feeding_survival_model_simple.stan.stan")

fit <- sampling(model, data = stan_data, iter = 2000, chains = 1, seed = 123)
posterior <- rstan::extract(fit)




# posterior is the object extracted from Stan
theta_post <- posterior$theta        # matrix: iterations x S
P_fed_post <- posterior$P_fed_sim   # matrix: iterations x S_sim



library(dplyr)

# Mean and 95% credible interval for theta
theta_summary <- apply(theta_post, 2, function(x){
  c(mean = mean(x), 
    lower = quantile(x, 0.025), 
    upper = quantile(x, 0.975))
}) %>% t() %>% as.data.frame()

theta_summary$obs <- 1:nrow(theta_summary)
theta_summary$hut <- mydata$hut







# Assuming theta_sim is 0.01 to 0.99
theta_sim <- seq(0.01, 0.99, length.out = 100)
a_mean <- mean(posterior$a)
b_mean <- mean(posterior$b)

P_fed_curve <- 1 - exp(b_mean * (1 - exp(a_mean * theta_sim)) / a_mean)
plot_data <- data.frame(theta = theta_sim, P_fed = P_fed_curve)













library(ggplot2)

#

P_fed_ci <- apply(posterior$P_fed_sim, 2, function(x) {
  c(lower = quantile(x, 0.025),
    upper = quantile(x, 0.975))
}) %>% t() %>% as.data.frame()

plot_data$lower <- P_fed_ci$lower
plot_data$upper <- P_fed_ci$upper






ggplot() +
  geom_point(data = mydata, 
             aes(x = (total - tot_dead)/total * 100, 
                 y = ((total - tot_dead)/total) * (bf_live/total) * 100),
             color = "gray", size = 3, alpha = 0.6) +
  geom_line(data = plot_data, aes(x = theta*100, y = P_fed*100), color = "gray", size = 1.2) +
  geom_ribbon(data = plot_data, aes(x = theta*100, ymin = lower*100, ymax = upper*100), 
              alpha = 0.2, fill = "gray") +
  labs(x = "Percent mosquito survival ",
       y = "Percent successfully fed)",
       title = "Relationship between mosquito survival and feeding (IG2)") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))










ggplot() +
  geom_point(data = mydata, 
             aes(x = (total - tot_dead)/total * 100, 
                 y = ((total - tot_dead)/total) * (bf_live/total) * 100,
                 size = total),  # point size based on total mosquitoes
             color = "gray", alpha = 0.6) +
  geom_line(data = plot_data, aes(x = theta*100, y = P_fed*100), color = "gray", size = 1.2) +
  geom_ribbon(data = plot_data, aes(x = theta*100, ymin = lower*100, ymax = upper*100), 
              alpha = 0.2, fill = "gray") +
  scale_size_continuous(name = "Total mosquitoes", range = c(2,6)) +  # control point size
  labs(x = "Percent mosquito survival",
       y = "Percent successfully fed",
       title = "Relationship between mosquito survival and feeding (IG2)") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))





ggplot() +
  geom_point(data = mydata, 
             aes(x = (total - tot_dead)/total * 100, 
                 y = ((total - tot_dead)/total) * (bf_live/total) * 100,
                 size = total),  # point size based on total mosquitoes
             color = "gray", alpha = 0.6) +
  geom_line(data = plot_data, 
            aes(x = theta*100, y = P_fed*100), 
            color = "gray", size = 1.2) +
  geom_ribbon(data = plot_data, 
              aes(x = theta*100, ymin = lower*100, ymax = upper*100), 
              alpha = 0.2, fill = "gray") +
  scale_size_continuous(name = "Total mosquitoes", range = c(2,6)) +
  labs(x = "Percent mosquito survival",
       y = "Percent successfully fed",
       title = "Relationship between mosquito survival and feeding (IG2)") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))













ggplot() +
  geom_point(data = mydata, 
             aes(x = (total - tot_dead)/total * 100, 
                 y = ((total - tot_dead)/total) * (bf_live/total) * 100,
                 size = total),  # point size based on total mosquitoes
             color = "gray", alpha = 0.6) +
  geom_line(data = plot_data, 
            aes(x = theta*100, y = P_fed*100), 
            color = "gray", size = 1.2) +
  geom_ribbon(data = plot_data, 
              aes(x = theta*100, ymin = lower*100, ymax = upper*100), 
              alpha = 0.2, fill = "gray") +
  scale_size_continuous(
    name = "Total mosquitoes",
    range = c(2, 6),
    breaks = scales::pretty_breaks()(range(mydata$total))  # auto breaks without 0
  ) +
  labs(
    x = "Percent mosquito survival",
    y = "Percent successfully fed",
    title = "Relationship between mosquito survival and feeding (IG2)"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )
































##############Combinealltreatments 


library(dplyr)
library(rstan)
library(readxl)
library(ggplot2)

#---------------------------
# 1️⃣ Load Data
#---------------------------
mydata <- read_excel("benin trial data results/Benin_NNP_hut_trial_raw data (2).xlsx", 
                     sheet = "Combined2")

#---------------------------
# 2️⃣ Define treatments to analyze
#---------------------------
treatments <- c("IG1", "IG2", "P3", "RG", "UT")

#---------------------------
# 3️⃣ Compile Stan model once
#---------------------------
model <- stan_model("benin trial data results/stan/feeding_survival_model_simple.stan.stan")

#---------------------------
# 4️⃣ Run Stan model for each treatment
#---------------------------
plot_list <- list()

for (trt in treatments) {
  message("Running model for treatment: ", trt)
  
  dat <- mydata %>% filter(Treatment == trt)
  
  stan_data <- list(
    S = nrow(dat),
    X_survive = dat$total - dat$tot_dead,
    N_caught = dat$total,
    X_succfed = dat$bf_live,
    site = as.numeric(as.factor(dat$hut)), 
    S_sim = 100,
    theta_sim = seq(0.01, 0.99, length.out = 100)
  )
  
  # Fit model
  fit <- sampling(model, data = stan_data, iter = 2000, chains = 1, seed = 123)
  posterior <- rstan::extract(fit)
  
  # Extract posterior predictions
  P_fed_ci <- apply(posterior$P_fed_sim, 2, function(x) {
    c(lower = quantile(x, 0.025), upper = quantile(x, 0.975))
  }) %>% t() %>% as.data.frame()
  
  a_mean <- mean(posterior$a)
  b_mean <- mean(posterior$b)
  theta_sim <- seq(0.01, 0.99, length.out = 100)
  P_fed_curve <- 1 - exp(b_mean * (1 - exp(a_mean * theta_sim)) / a_mean)
  
  plot_data <- data.frame(
    theta = theta_sim,
    P_fed = P_fed_curve,
    lower = P_fed_ci$lower,
    upper = P_fed_ci$upper,
    Treatment = trt
  )
  
  dat$Treatment <- trt
  plot_list[[trt]] <- list(plot_data = plot_data, obs_data = dat)
}

#---------------------------
# 5️⃣ Combine results from all treatments
#---------------------------
plot_all <- do.call(rbind, lapply(plot_list, function(x) x$plot_data))
obs_all  <- do.call(rbind, lapply(plot_list, function(x) x$obs_data))

#---------------------------
# 6️⃣ Plot all treatments together (5 panels)
#---------------------------
ggplot() +
  geom_point(data = obs_all, 
             aes(x = (total - tot_dead)/total * 100,
                 y = ((total - tot_dead)/total) * (bf_live/total) * 100,
                 size = total),
             color = "gray", alpha = 0.6) +
  geom_line(data = plot_all, 
            aes(x = theta*100, y = P_fed*100), 
            color = "gray", size = 1.2) +
  geom_ribbon(data = plot_all, 
              aes(x = theta*100, ymin = lower*100, ymax = upper*100),
              alpha = 0.2, fill = "gray") +
  facet_wrap(~Treatment, ncol = 3) +  # adjust ncol to layout style
  scale_size_continuous(name = "Total mosquitoes", range = c(2,6)) +
  labs(
    x = "Percent of mosquito survival in EHT",
    y = "Percent of successfully fed mosquitoes",
    title = "Relationship between mosquito survival and feeding across treatments"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.text = element_text(face = "bold")
  )






















ggplot() +
  geom_point(data = obs_all, 
             aes(x = (total - tot_dead)/total * 100,
                 y = ((total - tot_dead)/total) * (bf_live/total) * 100,
                 size = total),
             color = "gray", alpha = 0.6) +
  geom_line(data = plot_all, 
            aes(x = theta*100, y = P_fed*100), 
            color = "gray", size = 1.2) +
  geom_ribbon(data = plot_all, 
              aes(x = theta*100, ymin = lower*100, ymax = upper*100),
              alpha = 0.2, fill = "gray") +
  facet_wrap(~Treatment, ncol = 3) +
  scale_size_continuous(name = "Total mosquitoes", range = c(2,6)) +
  labs(
    x = "Mosquito survival in EHT (%)",
    y = "successfully fed mosquitoes (%)",
    title = " mosquito survival and feeding"
  ) +
  theme_minimal(base_size = 12) +  # base text size for uniformity
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.text = element_text(face = "bold", size = 10),   # facet labels
    axis.title = element_text(size = 12, face = "bold"),   # axis labels
    axis.text = element_text(size = 12),                   # axis tick labels
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5) # centered title
  )










































# Create size bins
obs_all <- obs_all %>%
  mutate(size_bin = cut(
    total,
    breaks = c(0, 30, 60, 100, 150),
    labels = c("1-30", "31-60", "61-100", "101-150"),
    include.lowest = TRUE
  ))

# Plot with custom size bins
ggplot() +
  geom_point(data = obs_all, 
             aes(x = (total - tot_dead)/total * 100,
                 y = ((total - tot_dead)/total) * (bf_live/total) * 100,
                 size = size_bin),  # use binned sizes
             color = "gray", alpha = 0.5) +
  geom_line(data = plot_all, 
            aes(x = theta*100, y = P_fed*100), 
            color = "gray", size = 1.2) +
  geom_ribbon(data = plot_all, 
              aes(x = theta*100, ymin = lower*100, ymax = upper*100),
              alpha = 0.2, fill = "gray") +
  facet_wrap(~Treatment, ncol = 3) +
  scale_size_manual(
    name = "Total mosquitoes",
    values = c("1-30" = 2, "31-60" = 4, "61-100" = 6, "101-150" = 8) # point sizes
  ) +
  labs(
    x = "Mosquito survival in EHT (%)",
    y = "successfully fed mosquitoes (%)",
    title = " mosquito survival and feeding"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.text = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )


























# Create size bins
obs_all <- obs_all %>%
  mutate(size_bin = cut(
    total,
    breaks = c(0, 30, 60, 100, 150),
    labels = c("1-30", "31-60", "61-100", "101-150"),
    include.lowest = TRUE
  ))

# Ensure plot_all has a Treatment column matching obs_all
if(!"Treatment" %in% names(plot_all)){
  plot_all$Treatment <- obs_all$Treatment[1]  # or assign correctly for all treatments
}

# Plot with 95% CI ribbon
ggplot() +
  # 1️⃣ Ribbon for 95% credible interval (drawn first)
  geom_ribbon(data = plot_all, 
              aes(x = theta*100, ymin = lower*100, ymax = upper*100),
              fill = "gray", alpha = 0.2) +
  
  # 2️⃣ Median predicted line
  geom_line(data = plot_all, 
            aes(x = theta*100, y = P_fed*100), 
            color = "gray20", size = 1.2) +
  
  # 3️⃣ Observed points
  geom_point(data = obs_all, 
             aes(x = (total - tot_dead)/total * 100,
                 y = ((total - tot_dead)/total) * (bf_live/total) * 100,
                 size = size_bin),
             shape = 21, fill = NA, color = "gray20", stroke = 1, alpha = 0.7) +
  
  # 4️⃣ Facet by Treatment
  facet_wrap(~Treatment, ncol = 3) +
  
  # 5️⃣ Point sizes
  scale_size_manual(
    name = "Total mosquitoes",
    values = c("1-30" = 2, "31-60" = 4, "61-100" = 6, "101-150" = 8)
  ) +
  
  # 6️⃣ Labels and theme
  labs(
    x = "Mosquito survival in EHT (%)",
    y = "Successfully fed mosquitoes (%)",
    title = "Mosquito survival and feeding"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.text = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )



































