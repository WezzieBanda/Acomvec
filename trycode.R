#usingoldmodel


library(rstan)

#packageVersion("rstan")
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)



df <- read_excel("Wezzie-Tom/EHT DATA/Benin Trial/benin trial data results/Benin_NNP_hut_trial_raw data (2).xlsx", 
                 sheet = "Combined2") #benindata




df$Treatment <- factor(df$Treatment, levels = c("UT", "IG1", "IG2", "P3","RG"))
df$treatment_id <- as.numeric(df$Treatment)
# Convert categorical variables to index values
df$treatment_id <- as.numeric(as.factor(df$Treatment))
df$age_id <- as.numeric(as.factor(df$Age))
df$time_id <- as.numeric(as.factor(df$Time))




###########model stan


# Prepare data for Stan
stan_data <- list(
  N = nrow(df),
  total_dead = df$unf_live,
  total = df$total,
  A = length(unique(df$age_id)),
  T = length(unique(df$treatment_id)),
  age = df$age_id,
  treatment = df$treatment_id,
  time = df$Time
)



# Stan model code as string
stan_model_code <- "
data {
  int<lower=1> N;
  int<lower=0> total_dead[N];
  int<lower=0> total[N];
  int<lower=1> A;
  int<lower=1> T;
  int<lower=1> age[N];
  int<lower=1> treatment[N];
  int<lower=1> time[N];
}
parameters {
  real alpha;
  vector[T] beta_treatment;
  vector[A] beta_age;
  real beta_time;
}
model {
  alpha ~ normal(0, 5);
  beta_treatment ~ normal(0, 2);
  beta_age ~ normal(0, 2);
  beta_time ~ normal(0, 2);

  for (i in 1:N) {
    real eta = alpha + beta_treatment[treatment[i]] + beta_age[age[i]] + beta_time * time[i];
    total_dead[i] ~ binomial_logit(total[i], eta);
  }
}
"

# Compile and run the model
fit <- stan(
  model_code = stan_model_code,
  data = stan_data,
  chains = 1,
  iter = 2000,
  control = list(max_treedepth = 15)
  #seed = 123,
  #cores = 2
)

# Print a summary
print(fit, pars = c("alpha", "beta_treatment", "beta_age", "beta_time"))










# Extract posterior samples
posterior <- rstan::extract(fit)

# Number of posterior draws
n_draws <- length(posterior$alpha)

# Unique treatment and age IDs
treatment_ids <- 1:stan_data$T
age_ids <- 1:stan_data$A

# Make a grid of all combinations
combo_grid <- expand.grid(treatment = treatment_ids, age = age_ids)

# Compute posterior predicted probability for each combo
results <- data.frame()

for (i in 1:nrow(combo_grid)) {
  t <- combo_grid$treatment[i]
  a <- combo_grid$age[i]
  
  # Linear predictor for all posterior draws
  eta <- posterior$alpha + posterior$beta_treatment[, t] + posterior$beta_age[, a] + posterior$beta_time * mean(stan_data$time)
  
  # Convert to probability using logistic function
  p <- plogis(eta)
  
  # Summarize posterior
  summary_row <- data.frame(
    treatment = t,
    age = a,
    mean = mean(p),
    lower = quantile(p, 0.025),
    upper = quantile(p, 0.975)
  )
  
  results <- rbind(results, summary_row)
}

# Optional: map treatment index back to names
treatment_labels <- levels(as.factor(df$Treatment))
results$treatment <- treatment_labels[results$treatment]

# View the estimates
print(results)

library(writexl)
write_xlsx(results, "beninunflive.xlsx")




























#letmedper timenow

#usingoldmodel


library(rstan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(writexl)

# ---- LOAD DATA ----
df <- read_excel("Wezzie-Tom/EHT DATA/Benin Trial/benin trial data results/Benin_NNP_hut_trial_raw data (2).xlsx",
                 sheet = "Combined2")



df <- read_excel("Wezzie-Tom/EHT DATA/Benin Trial/benin trial data results/NewNets.xlsx", 
                       sheet = "Aged")

df$Treatment <- factor(df$Treatment, levels = c( "IG1", "IG2", "P3","RG"))

df$treatment_id <- as.numeric(df$Treatment)
df$age_id <- as.numeric(as.factor(df$Age))
df$time_id <- as.numeric(as.factor(df$Time))

# ---- STAN DATA ----
stan_data <- list(
  N = nrow(df),
  total_dead = df$tot_72h_dead,
  total = df$total,
  A = length(unique(df$age_id)),
  T = length(unique(df$treatment_id)),
  age = df$age_id,
  treatment = df$treatment_id,
  time = df$Time
)

# ---- STAN MODEL ----
stan_model_code <- "
data {
  int<lower=1> N;
  int<lower=0> total_dead[N];
  int<lower=0> total[N];
  int<lower=1> A;
  int<lower=1> T;
  int<lower=1> age[N];
  int<lower=1> treatment[N];
  int<lower=1> time[N];
}
parameters {
  real alpha;
  vector[T] beta_treatment;
  vector[A] beta_age;
  real beta_time;
}
model {
  alpha ~ normal(0, 5);
  beta_treatment ~ normal(0, 2);
  beta_age ~ normal(0, 2);
  beta_time ~ normal(0, 2);

  for (i in 1:N) {
    real eta = alpha + beta_treatment[treatment[i]] + beta_age[age[i]] + beta_time * time[i];
    total_dead[i] ~ binomial_logit(total[i], eta);
  }
}
"

fit <- stan(
  model_code = stan_model_code,
  data = stan_data,
  chains = 1,
  iter = 2000,
  control = list(max_treedepth = 15)
)

print(fit, pars = c("alpha", "beta_treatment", "beta_age", "beta_time"))

# ---- EXTRACT POSTERIOR ----
posterior <- rstan::extract(fit)

# ---- PREPARE GRID FOR Time × Treatment ----
treatment_ids <- 1:stan_data$T
time_ids <- sort(unique(stan_data$time))

combo_grid <- expand.grid(treatment = treatment_ids, time = time_ids)

# ---- AVERAGE AGE EFFECT (RECOMMENDED) ----
avg_age_effect <- apply(posterior$beta_age, 1, mean)

# ---- POSTERIOR PREDICTION ----
results_time <- data.frame()

for (i in 1:nrow(combo_grid)) {
  t <- combo_grid$treatment[i]
  tm <- combo_grid$time[i]
  
  eta <- posterior$alpha +
    posterior$beta_treatment[, t] +
    avg_age_effect +
    posterior$beta_time * tm
  
  p <- plogis(eta)
  
  summary_row <- data.frame(
    treatment = t,
    time = tm,
    mean = mean(p),
    lower = quantile(p, 0.025),
    upper = quantile(p, 0.975)
  )
  
  results_time <- rbind(results_time, summary_row)
}

# ---- ADD TREATMENT LABELS ----
treatment_labels <- levels(df$Treatment)
results_time$treatment <- treatment_labels[results_time$treatment]

# ---- REORDER COLUMNS ----
results_time <- results_time %>% select(treatment, time, mean, lower, upper)

# ---- SAVE TO EXCEL ----
write_xlsx(results_time, "benin_deadliveeAGED_by_time.xlsx")

# ---- VIEW OUTPUT ----
print(results_time)
