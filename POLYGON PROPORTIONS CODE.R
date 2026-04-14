library(rstan)

#packageVersion("rstan")
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
options(mc.cores = parallel::detectCores())
# To avoid recompilation of unchanged Stan programs, we recommend calling
rstan_options(auto_write = TRUE)



df<- read_excel("Wezzie-Tom/EHT DATA/Benin Trial/benin trial data results/NewNets.xlsx", 
                sheet = "Aged")




#View(df)

df$Treatment <- factor(df$Treatment, levels = c( "IG1", "IG2", "P3","RG"))
df$treatment_id <- as.numeric(df$Treatment)
# Convert categorical variables to index values
df$treatment_id <- as.numeric(as.factor(df$Treatment))
df$age_id <- as.numeric(as.factor(df$Age))
df$time_id <- as.numeric(as.factor(df$Time))




###########model stan


# Prepare data for Stan
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
  real alpha;               // intercept
  vector[T] beta_treatment; // treatment effect
  vector[A] beta_age;       // age effect
  real beta_time;           // time slope
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





# Unique treatment and time IDs
treatment_ids <- 1:stan_data$T
time_values <- sort(unique(df$time_id))

# Make a grid of all combinations
combo_grid <- expand.grid(treatment = treatment_ids, time = time_values)




results_time <- data.frame()
fixed_age <- 1   # baseline age group

for (i in 1:nrow(combo_grid)) {
  t <- combo_grid$treatment[i]
  tm <- combo_grid$time[i]
  
  # Linear predictor (include fixed age)
  eta <- posterior$alpha + 
    posterior$beta_treatment[, t] + 
    posterior$beta_age[, fixed_age] +
    posterior$beta_time * tm
  
  # Convert to probability
  p <- plogis(eta)
  
  # Summarise using MEDIAN
  summary_row <- data.frame(
    treatment = t,
    time = tm,
    median = median(p),
    lower = quantile(p, 0.025),
    upper = quantile(p, 0.975)
  )
  
  results_time <- rbind(results_time, summary_row)
}

# Map treatment IDs to names
results_time$treatment <- levels(df$Treatment)[results_time$treatment]

# Optional: map time IDs to original labels
results_time$time <- factor(results_time$time,
                            labels = sort(unique(df$Time)))

# View results
print(results_time)


library(writexl)

# Save posterior summary results per treatment per time to Excel
write_xlsx(results_time, "medianmortalityAGED.xlsx")




#######New mortality






df<- read_excel("Wezzie-Tom/EHT DATA/Benin Trial/benin trial data results/NewNets.xlsx", 
                sheet = "New")




#View(df)

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
  total_dead = df$tot_72h_dead,
  total = df$total,
  A = length(unique(df$age_id)),
  T = length(unique(df$treatment_id)),
  age = df$age_id,
  treatment = df$treatment_id,
  time = df$Time
)


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
  real alpha;               // intercept
  vector[T] beta_treatment; // treatment effect
  vector[A] beta_age;       // age effect
  real beta_time;           // time slope
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





# Unique treatment and time IDs
treatment_ids <- 1:stan_data$T
time_values <- sort(unique(df$time_id))

# Make a grid of all combinations
combo_grid <- expand.grid(treatment = treatment_ids, time = time_values)




results_time <- data.frame()
fixed_age <- 1   # baseline age group

for (i in 1:nrow(combo_grid)) {
  t <- combo_grid$treatment[i]
  tm <- combo_grid$time[i]
  
  # Linear predictor (include fixed age)
  eta <- posterior$alpha + 
    posterior$beta_treatment[, t] + 
    posterior$beta_age[, fixed_age] +
    posterior$beta_time * tm
  
  # Convert to probability
  p <- plogis(eta)
  
  # Summarise using MEDIAN
  summary_row <- data.frame(
    treatment = t,
    time = tm,
    median = median(p),
    lower = quantile(p, 0.025),
    upper = quantile(p, 0.975)
  )
  
  results_time <- rbind(results_time, summary_row)
}

# Map treatment IDs to names
results_time$treatment <- levels(df$Treatment)[results_time$treatment]

# Optional: map time IDs to original labels
results_time$time <- factor(results_time$time,
                            labels = sort(unique(df$Time)))

# View results
print(results_time)


library(writexl)

# Save posterior summary results per treatment per time to Excel
write_xlsx(results_time, "medianmortalityNEW.xlsx")










#########################################################bflive







df<- read_excel("Wezzie-Tom/EHT DATA/Benin Trial/benin trial data results/NewNets.xlsx", 
                sheet = "New")




#View(df)

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
  total_dead = df$bf_live,
  total = df$total,
  A = length(unique(df$age_id)),
  T = length(unique(df$treatment_id)),
  age = df$age_id,
  treatment = df$treatment_id,
  time = df$Time
)


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
  real alpha;               // intercept
  vector[T] beta_treatment; // treatment effect
  vector[A] beta_age;       // age effect
  real beta_time;           // time slope
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





# Unique treatment and time IDs
treatment_ids <- 1:stan_data$T
time_values <- sort(unique(df$time_id))

# Make a grid of all combinations
combo_grid <- expand.grid(treatment = treatment_ids, time = time_values)




results_time <- data.frame()
fixed_age <- 1   # baseline age group

for (i in 1:nrow(combo_grid)) {
  t <- combo_grid$treatment[i]
  tm <- combo_grid$time[i]
  
  # Linear predictor (include fixed age)
  eta <- posterior$alpha + 
    posterior$beta_treatment[, t] + 
    posterior$beta_age[, fixed_age] +
    posterior$beta_time * tm
  
  # Convert to probability
  p <- plogis(eta)
  
  # Summarise using MEDIAN
  summary_row <- data.frame(
    treatment = t,
    time = tm,
    median = median(p),
    lower = quantile(p, 0.025),
    upper = quantile(p, 0.975)
  )
  
  results_time <- rbind(results_time, summary_row)
}

# Map treatment IDs to names
results_time$treatment <- levels(df$Treatment)[results_time$treatment]

# Optional: map time IDs to original labels
results_time$time <- factor(results_time$time,
                            labels = sort(unique(df$Time)))

# View results
print(results_time)


library(writexl)

# Save posterior summary results per treatment per time to Excel
write_xlsx(results_time, "medianbloodfedNEW.xlsx")













df<- read_excel("Wezzie-Tom/EHT DATA/Benin Trial/benin trial data results/NewNets.xlsx", 
                sheet = "New")




#View(df)

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
  total_dead = df$tot_72h_dead,
  total = df$total,
  A = length(unique(df$age_id)),
  T = length(unique(df$treatment_id)),
  age = df$age_id,
  treatment = df$treatment_id,
  time = df$Time
)


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
  real alpha;               // intercept
  vector[T] beta_treatment; // treatment effect
  vector[A] beta_age;       // age effect
  real beta_time;           // time slope
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





# Unique treatment and time IDs
treatment_ids <- 1:stan_data$T
time_values <- sort(unique(df$time_id))

# Make a grid of all combinations
combo_grid <- expand.grid(treatment = treatment_ids, time = time_values)




results_time <- data.frame()
fixed_age <- 1   # baseline age group

for (i in 1:nrow(combo_grid)) {
  t <- combo_grid$treatment[i]
  tm <- combo_grid$time[i]
  
  # Linear predictor (include fixed age)
  eta <- posterior$alpha + 
    posterior$beta_treatment[, t] + 
    posterior$beta_age[, fixed_age] +
    posterior$beta_time * tm
  
  # Convert to probability
  p <- plogis(eta)
  
  # Summarise using MEDIAN
  summary_row <- data.frame(
    treatment = t,
    time = tm,
    median = median(p),
    lower = quantile(p, 0.025),
    upper = quantile(p, 0.975)
  )
  
  results_time <- rbind(results_time, summary_row)
}

# Map treatment IDs to names
results_time$treatment <- levels(df$Treatment)[results_time$treatment]

# Optional: map time IDs to original labels
results_time$time <- factor(results_time$time,
                            labels = sort(unique(df$Time)))

# View results
print(results_time)


library(writexl)

# Save posterior summary results per treatment per time to Excel
write_xlsx(results_time, "medianmortalityNEW.xlsx")










#########################################################bflive







df<- read_excel("Wezzie-Tom/EHT DATA/Benin Trial/benin trial data results/NewNets.xlsx", 
                sheet = "Aged")




#View(df)

df$Treatment <- factor(df$Treatment, levels = c("IG1", "IG2", "P3","RG"))
df$treatment_id <- as.numeric(df$Treatment)
# Convert categorical variables to index values
df$treatment_id <- as.numeric(as.factor(df$Treatment))
df$age_id <- as.numeric(as.factor(df$Age))
df$time_id <- as.numeric(as.factor(df$Time))




###########model stan


# Prepare data for Stan
stan_data <- list(
  N = nrow(df),
  total_dead = df$bf_live,
  total = df$total,
  A = length(unique(df$age_id)),
  T = length(unique(df$treatment_id)),
  age = df$age_id,
  treatment = df$treatment_id,
  time = df$Time
)


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
  real alpha;               // intercept
  vector[T] beta_treatment; // treatment effect
  vector[A] beta_age;       // age effect
  real beta_time;           // time slope
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





# Unique treatment and time IDs
treatment_ids <- 1:stan_data$T
time_values <- sort(unique(df$time_id))

# Make a grid of all combinations
combo_grid <- expand.grid(treatment = treatment_ids, time = time_values)




results_time <- data.frame()
fixed_age <- 1   # baseline age group

for (i in 1:nrow(combo_grid)) {
  t <- combo_grid$treatment[i]
  tm <- combo_grid$time[i]
  
  # Linear predictor (include fixed age)
  eta <- posterior$alpha + 
    posterior$beta_treatment[, t] + 
    posterior$beta_age[, fixed_age] +
    posterior$beta_time * tm
  
  # Convert to probability
  p <- plogis(eta)
  
  # Summarise using MEDIAN
  summary_row <- data.frame(
    treatment = t,
    time = tm,
    median = median(p),
    lower = quantile(p, 0.025),
    upper = quantile(p, 0.975)
  )
  
  results_time <- rbind(results_time, summary_row)
}

# Map treatment IDs to names
results_time$treatment <- levels(df$Treatment)[results_time$treatment]

# Optional: map time IDs to original labels
results_time$time <- factor(results_time$time,
                            labels = sort(unique(df$Time)))

# View results
print(results_time)


library(writexl)

# Save posterior summary results per treatment per time to Excel
write_xlsx(results_time, "medianbloodfedAGED.xlsx")










#########################################################unfedlive







df<- read_excel("Wezzie-Tom/EHT DATA/Benin Trial/benin trial data results/NewNets.xlsx", 
                sheet = "Aged")




#View(df)

df$Treatment <- factor(df$Treatment, levels = c("IG1", "IG2", "P3","RG"))
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
  real alpha;               // intercept
  vector[T] beta_treatment; // treatment effect
  vector[A] beta_age;       // age effect
  real beta_time;           // time slope
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





# Unique treatment and time IDs
treatment_ids <- 1:stan_data$T
time_values <- sort(unique(df$time_id))

# Make a grid of all combinations
combo_grid <- expand.grid(treatment = treatment_ids, time = time_values)




results_time <- data.frame()
fixed_age <- 1   # baseline age group

for (i in 1:nrow(combo_grid)) {
  t <- combo_grid$treatment[i]
  tm <- combo_grid$time[i]
  
  # Linear predictor (include fixed age)
  eta <- posterior$alpha + 
    posterior$beta_treatment[, t] + 
    posterior$beta_age[, fixed_age] +
    posterior$beta_time * tm
  
  # Convert to probability
  p <- plogis(eta)
  
  # Summarise using MEDIAN
  summary_row <- data.frame(
    treatment = t,
    time = tm,
    median = median(p),
    lower = quantile(p, 0.025),
    upper = quantile(p, 0.975)
  )
  
  results_time <- rbind(results_time, summary_row)
}

# Map treatment IDs to names
results_time$treatment <- levels(df$Treatment)[results_time$treatment]

# Optional: map time IDs to original labels
results_time$time <- factor(results_time$time,
                            labels = sort(unique(df$Time)))

# View results
print(results_time)


library(writexl)

# Save posterior summary results per treatment per time to Excel
write_xlsx(results_time, "medianunfedAGED.xlsx")
















df<- read_excel("Wezzie-Tom/EHT DATA/Benin Trial/benin trial data results/NewNets.xlsx", 
                sheet = "New")




#View(df)

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
  real alpha;               // intercept
  vector[T] beta_treatment; // treatment effect
  vector[A] beta_age;       // age effect
  real beta_time;           // time slope
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





# Unique treatment and time IDs
treatment_ids <- 1:stan_data$T
time_values <- sort(unique(df$time_id))

# Make a grid of all combinations
combo_grid <- expand.grid(treatment = treatment_ids, time = time_values)




results_time <- data.frame()
fixed_age <- 1   # baseline age group

for (i in 1:nrow(combo_grid)) {
  t <- combo_grid$treatment[i]
  tm <- combo_grid$time[i]
  
  # Linear predictor (include fixed age)
  eta <- posterior$alpha + 
    posterior$beta_treatment[, t] + 
    posterior$beta_age[, fixed_age] +
    posterior$beta_time * tm
  
  # Convert to probability
  p <- plogis(eta)
  
  # Summarise using MEDIAN
  summary_row <- data.frame(
    treatment = t,
    time = tm,
    median = median(p),
    lower = quantile(p, 0.025),
    upper = quantile(p, 0.975)
  )
  
  results_time <- rbind(results_time, summary_row)
}

# Map treatment IDs to names
results_time$treatment <- levels(df$Treatment)[results_time$treatment]

# Optional: map time IDs to original labels
results_time$time <- factor(results_time$time,
                            labels = sort(unique(df$Time)))

# View results
print(results_time)


library(writexl)

# Save posterior summary results per treatment per time to Excel
write_xlsx(results_time, "medianunfedNew.xlsx")
























##########################deterencemedian



df<- read_excel("Wezzie-Tom/EHT DATA/Benin Trial/benin trial data results/NewNets.xlsx", 
                sheet = "New")


###########################
# 📌 Prepare Categorical IDs
###########################
df <- df %>%
  mutate(
    treatment_id = as.numeric(factor(Treatment, levels = c("UT","P3","IG1","IG2","RG"))),
    hut_id       = as.numeric(factor(hut)),
    sleeper_id   = as.numeric(factor(sleeper)),
    age_id       = as.numeric(factor(Age)),
    time_id      = as.numeric(factor(Time))
  )

###########################
# 📌 Stan Data List
###########################
stan_data <- list(
  N = nrow(df),
  total = df$total,
  T = length(unique(df$treatment_id)),
  H = length(unique(df$hut_id)),
  S = length(unique(df$sleeper_id)),
  A = length(unique(df$age_id)),
  Time_levels = length(unique(df$time_id)),
  treatment = df$treatment_id,
  hut       = df$hut_id,
  sleeper   = df$sleeper_id,
  age       = df$age_id,
  time      = df$time_id
)

###########################
# 📌 Stan Model
###########################
stan_model_code <- "
data {
  int<lower=1> N;
  int<lower=0> total[N];
  int<lower=1> T;
  int<lower=1> H;
  int<lower=1> S;
  int<lower=1> A;
  int<lower=1> Time_levels;
  int<lower=1> treatment[N];
  int<lower=1> hut[N];
  int<lower=1> sleeper[N];
  int<lower=1> age[N];
  int<lower=1, upper=Time_levels> time[N];
}
parameters {
  real alpha;
  vector[T] beta_treatment;
  vector[A] beta_age;
  matrix[A, T] beta_age_treat;
  vector[H] beta_hut;
  vector[S] beta_sleeper;
  vector[Time_levels] beta_time;
  real<lower=0> phi;
}
model {
  alpha ~ normal(0,5);
  beta_treatment ~ normal(0,2);
  beta_age ~ normal(0,2);
  to_vector(beta_age_treat) ~ normal(0,2);
  beta_hut ~ normal(0,2);
  beta_sleeper ~ normal(0,2);
  beta_time ~ normal(0,1);
  phi ~ gamma(2,0.1);

  for (i in 1:N) {
    real eta = alpha +
               beta_treatment[treatment[i]] +
               beta_age[age[i]] +
               beta_age_treat[age[i], treatment[i]] +
               beta_hut[hut[i]] +
               beta_sleeper[sleeper[i]] +
               beta_time[time[i]];

    total[i] ~ neg_binomial_2_log(eta, phi);
  }
}
"

###########################
# 📌 Fit Stan Model
###########################
fit <- rstan::stan(
  model_code = stan_model_code,
  data       = stan_data,
  iter       = 2000,
  chains     = 1,
  control    = list(max_treedepth = 12, adapt_delta = 0.95)
)

###########################
# 📌 Posterior Predictions
###########################
posterior <- rstan::extract(fit)

predicted_totals <- sapply(1:nrow(df), function(i) {
  eta_i <- posterior$alpha +
    posterior$beta_treatment[, df$treatment_id[i]] +
    posterior$beta_age[, df$age_id[i]] +
    posterior$beta_age_treat[, df$age_id[i], df$treatment_id[i]] +
    posterior$beta_hut[, df$hut_id[i]] +
    posterior$beta_sleeper[, df$sleeper_id[i]] +
    posterior$beta_time[, df$time_id[i]]
  exp(eta_i)
})

df$predicted_total <- colMeans(predicted_totals)

###########################
# 📌 Deterrence per Treatment × Time (NO AGE)
###########################
control_totals <- df %>% filter(Treatment == "UT")

deterrence_df <- df %>%
  filter(Treatment != "UT") %>%
  group_by(Treatment, Time) %>%
  summarise(
    treated_mean = mean(predicted_total),
    control_mean = mean(
      control_totals$predicted_total[
        control_totals$Time == Time
      ]
    ),
    deterrence = 1 - treated_mean / control_mean,
    .groups = "drop"
  )

print(deterrence_df, n = 50)





library(writexl)

# Save posterior summary results per treatment per time to Excel
write_xlsx(deterrence_df, "deterNew.xlsx")














##########################deterencemedian



df<- read_excel("Wezzie-Tom/EHT DATA/Benin Trial/benin trial data results/NewNets.xlsx", 
                sheet = "Aged")


###########################
# 📌 Prepare Categorical IDs
###########################
df <- df %>%
  mutate(
    treatment_id = as.numeric(factor(Treatment, levels = c("P3","IG1","IG2","RG"))),
    hut_id       = as.numeric(factor(hut)),
    sleeper_id   = as.numeric(factor(sleeper)),
    age_id       = as.numeric(factor(Age)),
    time_id      = as.numeric(factor(Time))
  )

###########################
# 📌 Stan Data List
###########################
stan_data <- list(
  N = nrow(df),
  total = df$total,
  T = length(unique(df$treatment_id)),
  H = length(unique(df$hut_id)),
  S = length(unique(df$sleeper_id)),
  A = length(unique(df$age_id)),
  Time_levels = length(unique(df$time_id)),
  treatment = df$treatment_id,
  hut       = df$hut_id,
  sleeper   = df$sleeper_id,
  age       = df$age_id,
  time      = df$time_id
)

###########################
# 📌 Stan Model
###########################
stan_model_code <- "
data {
  int<lower=1> N;
  int<lower=0> total[N];
  int<lower=1> T;
  int<lower=1> H;
  int<lower=1> S;
  int<lower=1> A;
  int<lower=1> Time_levels;
  int<lower=1> treatment[N];
  int<lower=1> hut[N];
  int<lower=1> sleeper[N];
  int<lower=1> age[N];
  int<lower=1, upper=Time_levels> time[N];
}
parameters {
  real alpha;
  vector[T] beta_treatment;
  vector[A] beta_age;
  matrix[A, T] beta_age_treat;
  vector[H] beta_hut;
  vector[S] beta_sleeper;
  vector[Time_levels] beta_time;
  real<lower=0> phi;
}
model {
  alpha ~ normal(0,5);
  beta_treatment ~ normal(0,2);
  beta_age ~ normal(0,2);
  to_vector(beta_age_treat) ~ normal(0,2);
  beta_hut ~ normal(0,2);
  beta_sleeper ~ normal(0,2);
  beta_time ~ normal(0,1);
  phi ~ gamma(2,0.1);

  for (i in 1:N) {
    real eta = alpha +
               beta_treatment[treatment[i]] +
               beta_age[age[i]] +
               beta_age_treat[age[i], treatment[i]] +
               beta_hut[hut[i]] +
               beta_sleeper[sleeper[i]] +
               beta_time[time[i]];

    total[i] ~ neg_binomial_2_log(eta, phi);
  }
}
"

###########################
# 📌 Fit Stan Model
###########################
fit <- rstan::stan(
  model_code = stan_model_code,
  data       = stan_data,
  iter       = 2000,
  chains     = 1,
  control    = list(max_treedepth = 12, adapt_delta = 0.95)
)

###########################
# 📌 Posterior Predictions
###########################
posterior <- rstan::extract(fit)

predicted_totals <- sapply(1:nrow(df), function(i) {
  eta_i <- posterior$alpha +
    posterior$beta_treatment[, df$treatment_id[i]] +
    posterior$beta_age[, df$age_id[i]] +
    posterior$beta_age_treat[, df$age_id[i], df$treatment_id[i]] +
    posterior$beta_hut[, df$hut_id[i]] +
    posterior$beta_sleeper[, df$sleeper_id[i]] +
    posterior$beta_time[, df$time_id[i]]
  exp(eta_i)
})

df$predicted_total <- colMeans(predicted_totals)

###########################
# 📌 Deterrence per Treatment × Time (NO AGE)
###########################
control_totals <- df %>% filter(Treatment == "UT")

deterrence_df <- df %>%
  filter(Treatment != "UT") %>%
  group_by(Treatment, Time) %>%
  summarise(
    treated_mean = mean(predicted_total),
    control_mean = mean(
      control_totals$predicted_total[
        control_totals$Time == Time
      ]
    ),
    deterrence = 1 - treated_mean / control_mean,
    .groups = "drop"
  )

print(deterrence_df, n = 50)





library(writexl)

# Save posterior summary results per treatment per time to Excel
write_xlsx(deterrence_df, "deterAGED.xlsx")









##########################################################################################################################################################
##########################################################################################################################################################
################################################################################################POLYGON PLOTS##########################################################
##########################################################################################################################################################
##########################################################################################################################################################

















###############OLD NETS  IG1

# --- Libraries ---
library(RColorBrewer)

# --- Time points (x-axis) ---
time <- c(0, 1,2,3)   # e.g. 0, 10, 20 washes
n_smooth <- 300
#x_smooth <- seq(min(time), max(time), length.out = n_smooth)




# --- Example probabilities (replace with your real data if needed) ---
Killed    <- c(0.113, 0.0810, 0.0580,  0.0416) #, 0.025874814)
Deterred  <- c(0.357,0.333, 0.288, 0.086) #, 0.288)
Repelled  <- c(0.582, 0.485, 0.415, 0.366) #, 0.393831024)
Fed       <- c(0.305, 0.434, 0.526, 0.593) #, 0.292294162)

# --- Smooth interpolation using natural spline (works with 3 points) ---
#smooth_fun <- function(x, y) spline(x, y, xout = x_smooth, method = "natural")$y
Killed_s   <- smooth_fun(time, Killed)
Deterred_s <- smooth_fun(time, Deterred)
Repelled_s <- smooth_fun(time, Repelled)
Fed_s      <- smooth_fun(time, Fed)

# --- Normalize so probabilities sum to 1 at each x ---
total <- Killed_s + Deterred_s + Repelled_s + Fed_s
Killed_s   <- Killed_s / total
Deterred_s <- Deterred_s / total
Repelled_s <- Repelled_s / total
Fed_s      <- Fed_s / total

# --- Cumulative layers for stacking (bottom-up) ---
y0 <- rep(0, n_smooth)
y1 <- Fed_s                        # bottom band: Successfully blood fed (red)
y2 <- Fed_s + Repelled_s           # next: Exited without feeding (yellow)
y3 <- Fed_s + Repelled_s + Deterred_s  # next: Deterred (green)
y4 <- Fed_s + Repelled_s + Deterred_s + Killed_s  # top: Killed (blue)

# --- Pastel colors with transparency (order matches y1..y4) ---
my.cols <- c(
  adjustcolor("#E15759", alpha.f = 0.8),  # Fed (red) - bottom
  adjustcolor("#FAD02E", alpha.f = 0.6),  # Exited without feeding (yellow)
  adjustcolor("#7AD151", alpha.f = 0.6),  # Deterred (green)
  adjustcolor("#6BAED6", alpha.f = 0.5)   # Killed (blue) - top
)

# --- Plot setup ---
par(bg = "white", mar = c(5,6,3,2), las = 1)
plot(x_smooth, x_smooth, type = "n", ylim = c(0,1), xlim = range(time),
     xlab = "Net Age", ylab = "Probable outcome of \n a single feeding attempt (%)",
     main = "Interceptor", cex.lab = 1.25, cex.main = 1.35, cex.axis = 1.05,
     frame.plot = FALSE)

# --- Draw stacked smooth polygons (bottom -> top) ---
polygon(c(x_smooth, rev(x_smooth)), c(y0, rev(y1)), col = my.cols[1], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y1, rev(y2)), col = my.cols[2], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y2, rev(y3)), col = my.cols[3], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y3, rev(y4)), col = my.cols[4], border = NA)

# --- Compute label positions robustly from the actual band positions ---
# pick an index near the right-middle of the plot so labels sit in the wider part of each band
idx <- round(n_smooth * 0.65)         # ~65% along x range (adjust if you prefer left/right)
x_lbl <- x_smooth[idx]

# For each band compute a y position inside the band (midpoint between band's lower & upper boundaries)
# Bands: (y0 -> y1) = Fed (red), (y1 -> y2) = Exited without feeding (yellow),
#        (y2 -> y3) = Deterred (green), (y3 -> y4) = Killed (blue)
fed_y    <- (y0[idx] + y1[idx]) / 2
exit_y   <- (y1[idx] + y2[idx]) / 2
deter_y  <- (y2[idx] + y3[idx]) / 2
killed_y <- (y3[idx] + y4[idx]) / 2

# Slight horizontal nudges so overlapping labels don't collide
x_nudge <- 0.12 * (max(time) - min(time))  # 12% of x-range
x_lbl_fed    <- x_lbl + x_nudge * 0.6
x_lbl_exit   <- x_lbl - x_nudge * 0.1
x_lbl_deter  <- x_lbl - x_nudge * 0.6
x_lbl_killed <- x_lbl + x_nudge * 0.2

# --- Add the labels with contrasting text colors chosen for readability inside each band ---
text(x = x_lbl_killed, y = killed_y, labels = "Killed", col = "#03396C", cex = 1.05, font = 2)
text(x = x_lbl_deter,  y = deter_y,  labels = "Deterred", col = "#1B5E20", cex = 1.05, font = 2)
text(x = x_lbl_exit,   y = exit_y,   labels = "Exited without feeding", col = "#8A4B00", cex = 1.00)
text(x = x_lbl_fed,    y = fed_y,    labels = "Successfully blood fed", col = "#7A0314", cex = 1.00)

# --- Y-axis as percentage ticks and soft frame outline ---
#axis(2, at = seq(0,1,0.25), labels = seq(0,100,25), las = 1)
#box(bty = "l", col = "gray70")


# --- Y-axis as percentage ---
#axis(2, at = seq(0,1,0.25), labels = seq(0,100,25))

# --- Frame outline ---
box(bty="l", col="gray70")



###ends here






###############OLD NETS  IG2

# --- Libraries ---
library(RColorBrewer)

# --- Time points (x-axis) ---
time <- c(0, 1,2,3)   # e.g. 0, 10, 20 washes
n_smooth <- 300
#x_smooth <- seq(min(time), max(time), length.out = n_smooth)


# --- Example probabilities (replace with your real data if needed) ---
Killed    <- c(0.307, 0.259, 0.218,  0.184) #, 0.025874814)
Deterred  <- c(0.372,0.33, 0.28, 0.109) #, 0.288)
Repelled  <- c(0.474, 0.438, 0.406, 0.380) #, 0.393831024)
Fed       <- c(0.219, 0.304, 0.375, 0.436) #, 0.292294162)

# --- Smooth interpolation using natural spline (works with 3 points) ---
#smooth_fun <- function(x, y) spline(x, y, xout = x_smooth, method = "natural")$y
Killed_s   <- smooth_fun(time, Killed)
Deterred_s <- smooth_fun(time, Deterred)
Repelled_s <- smooth_fun(time, Repelled)
Fed_s      <- smooth_fun(time, Fed)

# --- Normalize so probabilities sum to 1 at each x ---
total <- Killed_s + Deterred_s + Repelled_s + Fed_s
Killed_s   <- Killed_s / total
Deterred_s <- Deterred_s / total
Repelled_s <- Repelled_s / total
Fed_s      <- Fed_s / total

# --- Cumulative layers for stacking (bottom-up) ---
y0 <- rep(0, n_smooth)
y1 <- Fed_s                        # bottom band: Successfully blood fed (red)
y2 <- Fed_s + Repelled_s           # next: Exited without feeding (yellow)
y3 <- Fed_s + Repelled_s + Deterred_s  # next: Deterred (green)
y4 <- Fed_s + Repelled_s + Deterred_s + Killed_s  # top: Killed (blue)

# --- Pastel colors with transparency (order matches y1..y4) ---
my.cols <- c(
  adjustcolor("#E15759", alpha.f = 0.8),  # Fed (red) - bottom
  adjustcolor("#FAD02E", alpha.f = 0.6),  # Exited without feeding (yellow)
  adjustcolor("#7AD151", alpha.f = 0.6),  # Deterred (green)
  adjustcolor("#6BAED6", alpha.f = 0.5)   # Killed (blue) - top
)

# --- Plot setup ---
par(bg = "white", mar = c(5,6,3,2), las = 1)
plot(x_smooth, x_smooth, type = "n", ylim = c(0,1), xlim = range(time),
     xlab = "Net Age", ylab = "Probable outcome of \n a single feeding attempt (%)",
     main = "Interceptor G2", cex.lab = 1.25, cex.main = 1.35, cex.axis = 1.05,
     frame.plot = FALSE)

# --- Draw stacked smooth polygons (bottom -> top) ---
polygon(c(x_smooth, rev(x_smooth)), c(y0, rev(y1)), col = my.cols[1], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y1, rev(y2)), col = my.cols[2], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y2, rev(y3)), col = my.cols[3], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y3, rev(y4)), col = my.cols[4], border = NA)

# --- Compute label positions robustly from the actual band positions ---
# pick an index near the right-middle of the plot so labels sit in the wider part of each band
idx <- round(n_smooth * 0.65)         # ~65% along x range (adjust if you prefer left/right)
x_lbl <- x_smooth[idx]

# For each band compute a y position inside the band (midpoint between band's lower & upper boundaries)
# Bands: (y0 -> y1) = Fed (red), (y1 -> y2) = Exited without feeding (yellow),
#        (y2 -> y3) = Deterred (green), (y3 -> y4) = Killed (blue)
fed_y    <- (y0[idx] + y1[idx]) / 2
exit_y   <- (y1[idx] + y2[idx]) / 2
deter_y  <- (y2[idx] + y3[idx]) / 2
killed_y <- (y3[idx] + y4[idx]) / 2

# Slight horizontal nudges so overlapping labels don't collide
x_nudge <- 0.12 * (max(time) - min(time))  # 12% of x-range
x_lbl_fed    <- x_lbl + x_nudge * 0.6
x_lbl_exit   <- x_lbl - x_nudge * 0.1
x_lbl_deter  <- x_lbl - x_nudge * 0.6
x_lbl_killed <- x_lbl + x_nudge * 0.2

# --- Add the labels with contrasting text colors chosen for readability inside each band ---
text(x = x_lbl_killed, y = killed_y, labels = "Killed", col = "#03396C", cex = 1.05, font = 2)
text(x = x_lbl_deter,  y = deter_y,  labels = "Deterred", col = "#1B5E20", cex = 1.05, font = 2)
text(x = x_lbl_exit,   y = exit_y,   labels = "Exited without feeding", col = "#8A4B00", cex = 1.00)
text(x = x_lbl_fed,    y = fed_y,    labels = "Successfully blood fed", col = "#7A0314", cex = 1.00)

# --- Y-axis as percentage ticks and soft frame outline ---
#axis(2, at = seq(0,1,0.25), labels = seq(0,100,25), las = 1)
#box(bty = "l", col = "gray70")


# --- Y-axis as percentage ---
#axis(2, at = seq(0,1,0.25), labels = seq(0,100,25))

# --- Frame outline ---
box(bty="l", col="gray70")








###############OLD NETS  P3

# --- Libraries ---
library(RColorBrewer)

# --- Time points (x-axis) ---
time <- c(0, 1,2,3)   # e.g. 0, 10, 20 washes
n_smooth <- 300
x_smooth <- seq(min(time), max(time), length.out = n_smooth)



# --- Example probabilities (replace with your real data if needed) ---
Killed    <- c(0.189, 0.156, 0.128,  0.106) #, 0.025874814)
Deterred  <- c(0.42,0.216, 0.18, 0.086) #, 0.288)
Repelled  <- c(0.640, 0.574, 0.516, 0.468) #, 0.393831024)
Fed       <- c(0.163, 0.269, 0.355, 0.427) #, 0.292294162)

# --- Smooth interpolation using natural spline (works with 3 points) ---
#smooth_fun <- function(x, y) spline(x, y, xout = x_smooth, method = "natural")$y
Killed_s   <- smooth_fun(time, Killed)
Deterred_s <- smooth_fun(time, Deterred)
Repelled_s <- smooth_fun(time, Repelled)
Fed_s      <- smooth_fun(time, Fed)

# --- Normalize so probabilities sum to 1 at each x ---
total <- Killed_s + Deterred_s + Repelled_s + Fed_s
Killed_s   <- Killed_s / total
Deterred_s <- Deterred_s / total
Repelled_s <- Repelled_s / total
Fed_s      <- Fed_s / total

# --- Cumulative layers for stacking (bottom-up) ---
y0 <- rep(0, n_smooth)
y1 <- Fed_s                        # bottom band: Successfully blood fed (red)
y2 <- Fed_s + Repelled_s           # next: Exited without feeding (yellow)
y3 <- Fed_s + Repelled_s + Deterred_s  # next: Deterred (green)
y4 <- Fed_s + Repelled_s + Deterred_s + Killed_s  # top: Killed (blue)

# --- Pastel colors with transparency (order matches y1..y4) ---
my.cols <- c(
  adjustcolor("#E15759", alpha.f = 0.8),  # Fed (red) - bottom
  adjustcolor("#FAD02E", alpha.f = 0.6),  # Exited without feeding (yellow)
  adjustcolor("#7AD151", alpha.f = 0.6),  # Deterred (green)
  adjustcolor("#6BAED6", alpha.f = 0.5)   # Killed (blue) - top
)

# --- Plot setup ---
par(bg = "white", mar = c(5,6,3,2), las = 1)
plot(x_smooth, x_smooth, type = "n", ylim = c(0,1), xlim = range(time),
     xlab = "Net Age", ylab = "Probable outcome of\n a single feeding attempt (%)",
     main = "Olyset Plus", cex.lab = 1.25, cex.main = 1.35, cex.axis = 1.05,
     frame.plot = FALSE)

# --- Draw stacked smooth polygons (bottom -> top) ---
polygon(c(x_smooth, rev(x_smooth)), c(y0, rev(y1)), col = my.cols[1], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y1, rev(y2)), col = my.cols[2], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y2, rev(y3)), col = my.cols[3], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y3, rev(y4)), col = my.cols[4], border = NA)

# --- Compute label positions robustly from the actual band positions ---
# pick an index near the right-middle of the plot so labels sit in the wider part of each band
idx <- round(n_smooth * 0.65)         # ~65% along x range (adjust if you prefer left/right)
x_lbl <- x_smooth[idx]

# For each band compute a y position inside the band (midpoint between band's lower & upper boundaries)
# Bands: (y0 -> y1) = Fed (red), (y1 -> y2) = Exited without feeding (yellow),
#        (y2 -> y3) = Deterred (green), (y3 -> y4) = Killed (blue)
fed_y    <- (y0[idx] + y1[idx]) / 2
exit_y   <- (y1[idx] + y2[idx]) / 2
deter_y  <- (y2[idx] + y3[idx]) / 2
killed_y <- (y3[idx] + y4[idx]) / 2

# Slight horizontal nudges so overlapping labels don't collide
x_nudge <- 0.12 * (max(time) - min(time))  # 12% of x-range
x_lbl_fed    <- x_lbl + x_nudge * 0.6
x_lbl_exit   <- x_lbl - x_nudge * 0.1
x_lbl_deter  <- x_lbl - x_nudge * 0.6
x_lbl_killed <- x_lbl + x_nudge * 0.2

# --- Add the labels with contrasting text colors chosen for readability inside each band ---
text(x = x_lbl_killed, y = killed_y, labels = "Killed", col = "#03396C", cex = 1.05, font = 2)
text(x = x_lbl_deter,  y = deter_y,  labels = "Deterred", col = "#1B5E20", cex = 1.05, font = 2)
text(x = x_lbl_exit,   y = exit_y,   labels = "Exited without feeding", col = "#8A4B00", cex = 1.00)
text(x = x_lbl_fed,    y = fed_y,    labels = "Successfully blood fed", col = "#7A0314", cex = 1.00)

# --- Y-axis as percentage ticks and soft frame outline ---
#axis(2, at = seq(0,1,0.25), labels = seq(0,100,25), las = 1)
#box(bty = "l", col = "gray70")


# --- Y-axis as percentage ---
#axis(2, at = seq(0,1,0.25), labels = seq(0,100,25))

# --- Frame outline ---
box(bty="l", col="gray70")





###############OLD NETS  RG

# --- Libraries ---
library(RColorBrewer)

# --- Time points (x-axis) ---
time <- c(0, 1,2,3)   # e.g. 0, 10, 20 washes
n_smooth <- 300
#x_smooth <- seq(min(time), max(time), length.out = n_smooth)



# --- Example probabilities (replace with your real data if needed) ---
Killed    <- c(0.299, 0.247, 0.204,  0.168) #, 0.025874814)
Deterred  <- c(0,0,0,0) #,0.216, 0.18, 0.086) #, 0.288)
Repelled  <- c(0.493, 0.448, 0.412, 0.382) #, 0.393831024)
Fed       <- c(0.208, 0.305, 0.385, 0.451) #, 0.292294162)

# --- Smooth interpolation using natural spline (works with 3 points) ---
smooth_fun <- function(x, y) spline(x, y, xout = x_smooth, method = "natural")$y
Killed_s   <- smooth_fun(time, Killed)
Deterred_s <- smooth_fun(time, Deterred)
Repelled_s <- smooth_fun(time, Repelled)
Fed_s      <- smooth_fun(time, Fed)

# --- Normalize so probabilities sum to 1 at each x ---
total <- Killed_s + Deterred_s + Repelled_s + Fed_s
Killed_s   <- Killed_s / total
Deterred_s <- Deterred_s / total
Repelled_s <- Repelled_s / total
Fed_s      <- Fed_s / total

# --- Cumulative layers for stacking (bottom-up) ---
y0 <- rep(0, n_smooth)
y1 <- Fed_s                        # bottom band: Successfully blood fed (red)
y2 <- Fed_s + Repelled_s           # next: Exited without feeding (yellow)
y3 <- Fed_s + Repelled_s + Deterred_s  # next: Deterred (green)
y4 <- Fed_s + Repelled_s + Deterred_s + Killed_s  # top: Killed (blue)

# --- Pastel colors with transparency (order matches y1..y4) ---
my.cols <- c(
  adjustcolor("#E15759", alpha.f = 0.8),  # Fed (red) - bottom
  adjustcolor("#FAD02E", alpha.f = 0.6),  # Exited without feeding (yellow)
  adjustcolor("#7AD151", alpha.f = 0.6),  # Deterred (green)
  adjustcolor("#6BAED6", alpha.f = 0.5)   # Killed (blue) - top
)

# --- Plot setup ---
par(bg = "white", mar = c(5,6,3,2), las = 1)
plot(x_smooth, x_smooth, type = "n", ylim = c(0,1), xlim = range(time),
     xlab = "Net Age", ylab = "Probable outcome of \n a single feeding attempt (%)",
     main = "Royal Guard", cex.lab = 1.25, cex.main = 1.35, cex.axis = 1.05,
     frame.plot = FALSE)

# --- Draw stacked smooth polygons (bottom -> top) ---
polygon(c(x_smooth, rev(x_smooth)), c(y0, rev(y1)), col = my.cols[1], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y1, rev(y2)), col = my.cols[2], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y2, rev(y3)), col = my.cols[3], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y3, rev(y4)), col = my.cols[4], border = NA)

# --- Compute label positions robustly from the actual band positions ---
# pick an index near the right-middle of the plot so labels sit in the wider part of each band
idx <- round(n_smooth * 0.65)         # ~65% along x range (adjust if you prefer left/right)
x_lbl <- x_smooth[idx]

# For each band compute a y position inside the band (midpoint between band's lower & upper boundaries)
# Bands: (y0 -> y1) = Fed (red), (y1 -> y2) = Exited without feeding (yellow),
#        (y2 -> y3) = Deterred (green), (y3 -> y4) = Killed (blue)
fed_y    <- (y0[idx] + y1[idx]) / 2
exit_y   <- (y1[idx] + y2[idx]) / 2
deter_y  <- (y2[idx] + y3[idx]) / 2
killed_y <- (y3[idx] + y4[idx]) / 2

# Slight horizontal nudges so overlapping labels don't collide
x_nudge <- 0.12 * (max(time) - min(time))  # 12% of x-range
x_lbl_fed    <- x_lbl + x_nudge * 0.6
x_lbl_exit   <- x_lbl - x_nudge * 0.1
x_lbl_deter  <- x_lbl - x_nudge * 0.6
x_lbl_killed <- x_lbl + x_nudge * 0.2

# --- Add the labels with contrasting text colors chosen for readability inside each band ---
text(x = x_lbl_killed, y = killed_y, labels = "Killed", col = "#03396C", cex = 1.05, font = 2)
#text(x = x_lbl_deter,  y = deter_y,  labels = "Deterred", col = "#1B5E20", cex = 1.05, font = 2)
text(x = x_lbl_exit,   y = exit_y,   labels = "Exited without feeding", col = "#8A4B00", cex = 1.00)
text(x = x_lbl_fed,    y = fed_y,    labels = "Successfully blood fed", col = "#7A0314", cex = 1.00)

# --- Y-axis as percentage ticks and soft frame outline ---
#axis(2, at = seq(0,1,0.25), labels = seq(0,100,25), las = 1)
#box(bty = "l", col = "gray70")


# --- Y-axis as percentage ---
#axis(2, at = seq(0,1,0.25), labels = seq(0,100,25))

# --- Frame outline ---
box(bty="l", col="gray70")




























##########################################################################################################################################################
##########################################################################################################################################################
################################################################################################TANZANIA##########################################################
##########################################################################################################################################################
##########################################################################################################################################################

















###############OLD NETS  IG1

# --- Libraries ---
library(RColorBrewer)

# --- Time points (x-axis) ---
time <- c(0, 1,2,3)   # e.g. 0, 10, 20 washes
n_smooth <- 300
#x_smooth <- seq(min(time), max(time), length.out = n_smooth)




# --- Example probabilities (replace with your real data if needed) ---
Killed    <- c(0.113, 0.0946, 0.0688,  0.0501) #, 0.025874814)
Deterred  <- c(0.15,0.12, 0.095, 0.02) #, 0.288)
Repelled  <- c(0.72, 0.589, 0.494, 0.425) #, 0.393831024)
Fed       <- c(0.305, 0.434, 0.526, 0.593) #, 0.292294162)

# --- Smooth interpolation using natural spline (works with 3 points) ---
#smooth_fun <- function(x, y) spline(x, y, xout = x_smooth, method = "natural")$y
Killed_s   <- smooth_fun(time, Killed)
Deterred_s <- smooth_fun(time, Deterred)
Repelled_s <- smooth_fun(time, Repelled)
Fed_s      <- smooth_fun(time, Fed)

# --- Normalize so probabilities sum to 1 at each x ---
total <- Killed_s + Deterred_s + Repelled_s + Fed_s
Killed_s   <- Killed_s / total
Deterred_s <- Deterred_s / total
Repelled_s <- Repelled_s / total
Fed_s      <- Fed_s / total

# --- Cumulative layers for stacking (bottom-up) ---
y0 <- rep(0, n_smooth)
y1 <- Fed_s                        # bottom band: Successfully blood fed (red)
y2 <- Fed_s + Repelled_s           # next: Exited without feeding (yellow)
y3 <- Fed_s + Repelled_s + Deterred_s  # next: Deterred (green)
y4 <- Fed_s + Repelled_s + Deterred_s + Killed_s  # top: Killed (blue)

# --- Pastel colors with transparency (order matches y1..y4) ---
my.cols <- c(
  adjustcolor("#E15759", alpha.f = 0.8),  # Fed (red) - bottom
  adjustcolor("#FAD02E", alpha.f = 0.6),  # Exited without feeding (yellow)
  adjustcolor("#7AD151", alpha.f = 0.6),  # Deterred (green)
  adjustcolor("#6BAED6", alpha.f = 0.5)   # Killed (blue) - top
)

# --- Plot setup ---
par(bg = "white", mar = c(5,6,3,2), las = 1)
plot(x_smooth, x_smooth, type = "n", ylim = c(0,1), xlim = range(time),
     xlab = "Net Age", ylab = "Probable outcome of \n a single feeding attempt (%)",
     main = "Interceptor", cex.lab = 1.25, cex.main = 1.35, cex.axis = 1.05,
     frame.plot = FALSE)

# --- Draw stacked smooth polygons (bottom -> top) ---
polygon(c(x_smooth, rev(x_smooth)), c(y0, rev(y1)), col = my.cols[1], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y1, rev(y2)), col = my.cols[2], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y2, rev(y3)), col = my.cols[3], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y3, rev(y4)), col = my.cols[4], border = NA)

# --- Compute label positions robustly from the actual band positions ---
# pick an index near the right-middle of the plot so labels sit in the wider part of each band
idx <- round(n_smooth * 0.65)         # ~65% along x range (adjust if you prefer left/right)
x_lbl <- x_smooth[idx]

# For each band compute a y position inside the band (midpoint between band's lower & upper boundaries)
# Bands: (y0 -> y1) = Fed (red), (y1 -> y2) = Exited without feeding (yellow),
#        (y2 -> y3) = Deterred (green), (y3 -> y4) = Killed (blue)
fed_y    <- (y0[idx] + y1[idx]) / 2
exit_y   <- (y1[idx] + y2[idx]) / 2
deter_y  <- (y2[idx] + y3[idx]) / 2
killed_y <- (y3[idx] + y4[idx]) / 2

# Slight horizontal nudges so overlapping labels don't collide
x_nudge <- 0.12 * (max(time) - min(time))  # 12% of x-range
x_lbl_fed    <- x_lbl + x_nudge * 0.6
x_lbl_exit   <- x_lbl - x_nudge * 0.1
x_lbl_deter  <- x_lbl - x_nudge * 0.6
x_lbl_killed <- x_lbl + x_nudge * 0.2

# --- Add the labels with contrasting text colors chosen for readability inside each band ---
text(x = x_lbl_killed, y = killed_y, labels = "Killed", col = "#03396C", cex = 1.05, font = 2)
text(x = x_lbl_deter,  y = deter_y,  labels = "Deterred", col = "#1B5E20", cex = 1.05, font = 2)
text(x = x_lbl_exit,   y = exit_y,   labels = "Exited without feeding", col = "#8A4B00", cex = 1.00)
text(x = x_lbl_fed,    y = fed_y,    labels = "Successfully blood fed", col = "#7A0314", cex = 1.00)

# --- Y-axis as percentage ticks and soft frame outline ---
#axis(2, at = seq(0,1,0.25), labels = seq(0,100,25), las = 1)
#box(bty = "l", col = "gray70")


# --- Y-axis as percentage ---
#axis(2, at = seq(0,1,0.25), labels = seq(0,100,25))

# --- Frame outline ---
box(bty="l", col="gray70")



###ends here






###############OLD NETS  IG2

# --- Libraries ---
library(RColorBrewer)

# --- Time points (x-axis) ---
time <- c(0, 1,2,3)   # e.g. 0, 10, 20 washes
n_smooth <- 300
#x_smooth <- seq(min(time), max(time), length.out = n_smooth)



# --- Example probabilities (replace with your real data if needed) ---
Killed    <- c(0.28, 0.222, 0.176,  0.140) #, 0.025874814)
Deterred  <- c(0.15, 0.134, 0.061, 0.047) #, 0.288)
Repelled  <- c(0.59, 0.518, 0.461, 0.415) #, 0.393831024)
Fed       <- c(0.13, 0.260, 0.363, 0.445) #, 0.292294162)

# --- Smooth interpolation using natural spline (works with 3 points) ---
#smooth_fun <- function(x, y) spline(x, y, xout = x_smooth, method = "natural")$y
Killed_s   <- smooth_fun(time, Killed)
Deterred_s <- smooth_fun(time, Deterred)
Repelled_s <- smooth_fun(time, Repelled)
Fed_s      <- smooth_fun(time, Fed)

# --- Normalize so probabilities sum to 1 at each x ---
total <- Killed_s + Deterred_s + Repelled_s + Fed_s
Killed_s   <- Killed_s / total
Deterred_s <- Deterred_s / total
Repelled_s <- Repelled_s / total
Fed_s      <- Fed_s / total

# --- Cumulative layers for stacking (bottom-up) ---
y0 <- rep(0, n_smooth)
y1 <- Fed_s                        # bottom band: Successfully blood fed (red)
y2 <- Fed_s + Repelled_s           # next: Exited without feeding (yellow)
y3 <- Fed_s + Repelled_s + Deterred_s  # next: Deterred (green)
y4 <- Fed_s + Repelled_s + Deterred_s + Killed_s  # top: Killed (blue)

# --- Pastel colors with transparency (order matches y1..y4) ---
my.cols <- c(
  adjustcolor("#E15759", alpha.f = 0.8),  # Fed (red) - bottom
  adjustcolor("#FAD02E", alpha.f = 0.6),  # Exited without feeding (yellow)
  adjustcolor("#7AD151", alpha.f = 0.6),  # Deterred (green)
  adjustcolor("#6BAED6", alpha.f = 0.5)   # Killed (blue) - top
)

# --- Plot setup ---
par(bg = "white", mar = c(5,6,3,2), las = 1)
plot(x_smooth, x_smooth, type = "n", ylim = c(0,1), xlim = range(time),
     xlab = "Net Age", ylab = "Probable outcome of \n a single feeding attempt (%)",
     main = "Interceptor G2", cex.lab = 1.25, cex.main = 1.35, cex.axis = 1.05,
     frame.plot = FALSE)

# --- Draw stacked smooth polygons (bottom -> top) ---
polygon(c(x_smooth, rev(x_smooth)), c(y0, rev(y1)), col = my.cols[1], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y1, rev(y2)), col = my.cols[2], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y2, rev(y3)), col = my.cols[3], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y3, rev(y4)), col = my.cols[4], border = NA)

# --- Compute label positions robustly from the actual band positions ---
# pick an index near the right-middle of the plot so labels sit in the wider part of each band
idx <- round(n_smooth * 0.65)         # ~65% along x range (adjust if you prefer left/right)
x_lbl <- x_smooth[idx]

# For each band compute a y position inside the band (midpoint between band's lower & upper boundaries)
# Bands: (y0 -> y1) = Fed (red), (y1 -> y2) = Exited without feeding (yellow),
#        (y2 -> y3) = Deterred (green), (y3 -> y4) = Killed (blue)
fed_y    <- (y0[idx] + y1[idx]) / 2
exit_y   <- (y1[idx] + y2[idx]) / 2
deter_y  <- (y2[idx] + y3[idx]) / 2
killed_y <- (y3[idx] + y4[idx]) / 2

# Slight horizontal nudges so overlapping labels don't collide
x_nudge <- 0.12 * (max(time) - min(time))  # 12% of x-range
x_lbl_fed    <- x_lbl + x_nudge * 0.6
x_lbl_exit   <- x_lbl - x_nudge * 0.1
x_lbl_deter  <- x_lbl - x_nudge * 0.6
x_lbl_killed <- x_lbl + x_nudge * 0.2

# --- Add the labels with contrasting text colors chosen for readability inside each band ---
text(x = x_lbl_killed, y = killed_y, labels = "Killed", col = "#03396C", cex = 1.05, font = 2)
text(x = x_lbl_deter,  y = deter_y,  labels = "Deterred", col = "#1B5E20", cex = 1.05, font = 2)
text(x = x_lbl_exit,   y = exit_y,   labels = "Exited without feeding", col = "#8A4B00", cex = 1.00)
text(x = x_lbl_fed,    y = fed_y,    labels = "Successfully blood fed", col = "#7A0314", cex = 1.00)

# --- Y-axis as percentage ticks and soft frame outline ---
#axis(2, at = seq(0,1,0.25), labels = seq(0,100,25), las = 1)
#box(bty = "l", col = "gray70")


# --- Y-axis as percentage ---
#axis(2, at = seq(0,1,0.25), labels = seq(0,100,25))

# --- Frame outline ---
box(bty="l", col="gray70")








###############OLD NETS  P3

# --- Libraries ---
library(RColorBrewer)

# --- Time points (x-axis) ---
time <- c(0, 1,2,3)   # e.g. 0, 10, 20 washes
n_smooth <- 300
#x_smooth <- seq(min(time), max(time), length.out = n_smooth)



# --- Example probabilities (replace with your real data if needed) ---
Killed    <- c(0.27, 0.200, 0.148,  0.109) #, 0.025874814)
Deterred  <- c(0,0 ,0, 0 ) #,0.216, 0.18, 0.086) #, 0.288)
Repelled  <- c(0.58, 0.492, 0.426, 0.378) #, 0.393831024)
Fed       <- c(0.15, 0.309, 0.426, 0.513) #, 0.292294162)

# --- Smooth interpolation using natural spline (works with 3 points) ---
#smooth_fun <- function(x, y) spline(x, y, xout = x_smooth, method = "natural")$y
Killed_s   <- smooth_fun(time, Killed)
Deterred_s <- smooth_fun(time, Deterred)
Repelled_s <- smooth_fun(time, Repelled)
Fed_s      <- smooth_fun(time, Fed)

# --- Normalize so probabilities sum to 1 at each x ---
total <- Killed_s + Deterred_s + Repelled_s + Fed_s
Killed_s   <- Killed_s / total
Deterred_s <- Deterred_s / total
Repelled_s <- Repelled_s / total
Fed_s      <- Fed_s / total

# --- Cumulative layers for stacking (bottom-up) ---
y0 <- rep(0, n_smooth)
y1 <- Fed_s                        # bottom band: Successfully blood fed (red)
y2 <- Fed_s + Repelled_s           # next: Exited without feeding (yellow)
y3 <- Fed_s + Repelled_s + Deterred_s  # next: Deterred (green)
y4 <- Fed_s + Repelled_s + Deterred_s + Killed_s  # top: Killed (blue)

# --- Pastel colors with transparency (order matches y1..y4) ---
my.cols <- c(
  adjustcolor("#E15759", alpha.f = 0.8),  # Fed (red) - bottom
  adjustcolor("#FAD02E", alpha.f = 0.6),  # Exited without feeding (yellow)
  adjustcolor("#7AD151", alpha.f = 0.6),  # Deterred (green)
  adjustcolor("#6BAED6", alpha.f = 0.5)   # Killed (blue) - top
)

# --- Plot setup ---
par(bg = "white", mar = c(5,6,3,2), las = 1)
plot(x_smooth, x_smooth, type = "n", ylim = c(0,1), xlim = range(time),
     xlab = "Net Age", ylab = "Probable outcome of\n a single feeding attempt (%)",
     main = "Olyset Plus", cex.lab = 1.25, cex.main = 1.35, cex.axis = 1.05,
     frame.plot = FALSE)

# --- Draw stacked smooth polygons (bottom -> top) ---
polygon(c(x_smooth, rev(x_smooth)), c(y0, rev(y1)), col = my.cols[1], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y1, rev(y2)), col = my.cols[2], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y2, rev(y3)), col = my.cols[3], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y3, rev(y4)), col = my.cols[4], border = NA)

# --- Compute label positions robustly from the actual band positions ---
# pick an index near the right-middle of the plot so labels sit in the wider part of each band
idx <- round(n_smooth * 0.65)         # ~65% along x range (adjust if you prefer left/right)
x_lbl <- x_smooth[idx]

# For each band compute a y position inside the band (midpoint between band's lower & upper boundaries)
# Bands: (y0 -> y1) = Fed (red), (y1 -> y2) = Exited without feeding (yellow),
#        (y2 -> y3) = Deterred (green), (y3 -> y4) = Killed (blue)
fed_y    <- (y0[idx] + y1[idx]) / 2
exit_y   <- (y1[idx] + y2[idx]) / 2
deter_y  <- (y2[idx] + y3[idx]) / 2
killed_y <- (y3[idx] + y4[idx]) / 2

# Slight horizontal nudges so overlapping labels don't collide
x_nudge <- 0.12 * (max(time) - min(time))  # 12% of x-range
x_lbl_fed    <- x_lbl + x_nudge * 0.6
x_lbl_exit   <- x_lbl - x_nudge * 0.1
x_lbl_deter  <- x_lbl - x_nudge * 0.6
x_lbl_killed <- x_lbl + x_nudge * 0.2

# --- Add the labels with contrasting text colors chosen for readability inside each band ---
text(x = x_lbl_killed, y = killed_y, labels = "Killed", col = "#03396C", cex = 1.05, font = 2)
#text(x = x_lbl_deter,  y = deter_y,  labels = "Deterred", col = "#1B5E20", cex = 1.05, font = 2)
text(x = x_lbl_exit,   y = exit_y,   labels = "Exited without feeding", col = "#8A4B00", cex = 1.00)
text(x = x_lbl_fed,    y = fed_y,    labels = "Successfully blood fed", col = "#7A0314", cex = 1.00)

# --- Y-axis as percentage ticks and soft frame outline ---
#axis(2, at = seq(0,1,0.25), labels = seq(0,100,25), las = 1)
#box(bty = "l", col = "gray70")


# --- Y-axis as percentage ---
#axis(2, at = seq(0,1,0.25), labels = seq(0,100,25))

# --- Frame outline ---
box(bty="l", col="gray70")





###############OLD NETS  RG

# --- Libraries ---
library(RColorBrewer)

# --- Time points (x-axis) ---
time <- c(0, 1,2,3)   # e.g. 0, 10, 20 washes
n_smooth <- 300
#x_smooth <- seq(min(time), max(time), length.out = n_smooth)



# --- Example probabilities (replace with your real data if needed) ---
Killed    <- c(0.7,0.45 , 0.32,  0.09) #, 0.025874814)
Deterred  <- c(0,0,0,0) #,0.216, 0.18, 0.086) #, 0.288)
Fed  <- c(0.06, 0.098, 0.11, 0.123) #, 0.393831024)
Repelled       <- 1-(killed+Fed) #, 0.292294162)

# --- Smooth interpolation using natural spline (works with 3 points) ---
smooth_fun <- function(x, y) spline(x, y, xout = x_smooth, method = "natural")$y
Killed_s   <- smooth_fun(time, Killed)
Deterred_s <- smooth_fun(time, Deterred)
Repelled_s <- smooth_fun(time, Repelled)
Fed_s      <- smooth_fun(time, Fed)

# --- Normalize so probabilities sum to 1 at each x ---
total <- Killed_s + Deterred_s + Repelled_s + Fed_s
Killed_s   <- Killed_s / total
Deterred_s <- Deterred_s / total
Repelled_s <- Repelled_s / total
Fed_s      <- Fed_s / total

# --- Cumulative layers for stacking (bottom-up) ---
y0 <- rep(0, n_smooth)
y1 <- Fed_s                        # bottom band: Successfully blood fed (red)
y2 <- Fed_s + Repelled_s           # next: Exited without feeding (yellow)
y3 <- Fed_s + Repelled_s + Deterred_s  # next: Deterred (green)
y4 <- Fed_s + Repelled_s + Deterred_s + Killed_s  # top: Killed (blue)

# --- Pastel colors with transparency (order matches y1..y4) ---
my.cols <- c(
  adjustcolor("#E15759", alpha.f = 0.8),  # Fed (red) - bottom
  adjustcolor("#FAD02E", alpha.f = 0.6),  # Exited without feeding (yellow)
  adjustcolor("#7AD151", alpha.f = 0.6),  # Deterred (green)
  adjustcolor("#6BAED6", alpha.f = 0.5)   # Killed (blue) - top
)

# --- Plot setup ---
par(bg = "white", mar = c(5,6,3,2), las = 1)
plot(x_smooth, x_smooth, type = "n", ylim = c(0,1), xlim = range(time),
     xlab = "Net Age", ylab = "Probable outcome of \n a single feeding attempt (%)",
     main = "Royal Guard", cex.lab = 1.25, cex.main = 1.35, cex.axis = 1.05,
     frame.plot = FALSE)

# --- Draw stacked smooth polygons (bottom -> top) ---
polygon(c(x_smooth, rev(x_smooth)), c(y0, rev(y1)), col = my.cols[1], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y1, rev(y2)), col = my.cols[2], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y2, rev(y3)), col = my.cols[3], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y3, rev(y4)), col = my.cols[4], border = NA)

# --- Compute label positions robustly from the actual band positions ---
# pick an index near the right-middle of the plot so labels sit in the wider part of each band
idx <- round(n_smooth * 0.65)         # ~65% along x range (adjust if you prefer left/right)
x_lbl <- x_smooth[idx]

# For each band compute a y position inside the band (midpoint between band's lower & upper boundaries)
# Bands: (y0 -> y1) = Fed (red), (y1 -> y2) = Exited without feeding (yellow),
#        (y2 -> y3) = Deterred (green), (y3 -> y4) = Killed (blue)
fed_y    <- (y0[idx] + y1[idx]) / 2
exit_y   <- (y1[idx] + y2[idx]) / 2
deter_y  <- (y2[idx] + y3[idx]) / 2
killed_y <- (y3[idx] + y4[idx]) / 2

# Slight horizontal nudges so overlapping labels don't collide
x_nudge <- 0.12 * (max(time) - min(time))  # 12% of x-range
x_lbl_fed    <- x_lbl + x_nudge * 0.6
x_lbl_exit   <- x_lbl - x_nudge * 0.1
x_lbl_deter  <- x_lbl - x_nudge * 0.6
x_lbl_killed <- x_lbl + x_nudge * 0.2

# --- Add the labels with contrasting text colors chosen for readability inside each band ---
text(x = x_lbl_killed, y = killed_y, labels = "Killed", col = "#03396C", cex = 1.05, font = 2)
#text(x = x_lbl_deter,  y = deter_y,  labels = "Deterred", col = "#1B5E20", cex = 1.05, font = 2)
text(x = x_lbl_exit,   y = exit_y,   labels = "Exited without feeding", col = "#8A4B00", cex = 1.00)
text(x = x_lbl_fed,    y = fed_y,    labels = "Successfully blood fed", col = "#7A0314", cex = 1.00)

# --- Y-axis as percentage ticks and soft frame outline ---
#axis(2, at = seq(0,1,0.25), labels = seq(0,100,25), las = 1)
#box(bty = "l", col = "gray70")


# --- Y-axis as percentage ---
#axis(2, at = seq(0,1,0.25), labels = seq(0,100,25))

# --- Frame outline ---
box(bty="l", col="gray70")



































####################################TANZANIA TRIAL



# --- Libraries ---
library(RColorBrewer)

# --- Time points (x-axis) ---
time <- c(0, 1,2,3)   # e.g. 0, 10, 20 washes
n_smooth <- 300


# =========================
# RESET WORKSPACE (optional but helpful)
# =========================
#rm(list = ls())

# =========================
# LIBRARIES
# =========================
library(dplyr)

# =========================
# TIME POINTS
# =========================
time <- c(0, 1, 2, 3)

# =========================
# SMOOTH GRID (THIS WAS MISSING)
# =========================
x_smooth <- seq(min(time), max(time), length.out = 300)

# =========================
# DATA (USE CONSISTENT NAMES)
0.332775239
0.241200643
0.236239984
0.219162581


Deterred <- c(0.372, 0.33, 0.28, 0.109)
Repelled <- c(0.298, 0.330, 0.311, 0.279)
Killed  <- c(0.332775239, 0.241200643, 0.236239984, 0.219162581)
Fed  <- c(0.137, 0.257, 0.388, 0.506)
# =========================
# SMOOTH FUNCTION
# =========================
smooth_fun <- function(x, y) {
  spline(x, y, xout = x_smooth, method = "natural")$y
}

# =========================
# APPLY SMOOTHING
# =========================
Killed_s   <- smooth_fun(time, Killed)
Deterred_s <- smooth_fun(time, Deterred)
Repelled_s <- smooth_fun(time, Repelled)
Fed_s      <- smooth_fun(time, Fed)

# =========================
# CHECK (VERY IMPORTANT)
# =========================
round(Killed + Deterred + Repelled + Fed, 3)








# --- Cumulative layers for stacking (bottom-up) ---
y0 <- rep(0, n_smooth)
y1 <- Fed_s                        # bottom band: Successfully blood fed (red)
y2 <- Fed_s + Repelled_s           # next: Exited without feeding (yellow)
y3 <- Fed_s + Repelled_s + Deterred_s  # next: Deterred (green)
y4 <- Fed_s + Repelled_s + Deterred_s + Killed_s  # top: Killed (blue)

# --- Pastel colors with transparency (order matches y1..y4) ---
my.cols <- c(
  adjustcolor("#E15759", alpha.f = 0.8),  # Fed (red) - bottom
  adjustcolor("#FAD02E", alpha.f = 0.6),  # Exited without feeding (yellow)
  adjustcolor("#7AD151", alpha.f = 0.6),  # Deterred (green)
  adjustcolor("#6BAED6", alpha.f = 0.5)   # Killed (blue) - top
)

# --- Plot setup ---
par(bg = "white", mar = c(5,6,3,2), las = 1)
plot(x_smooth, x_smooth, type = "n", ylim = c(0,1), xlim = range(time),
     xlab = "Net Age (years)", ylab = "Probable outcome of \n a single feeding attempt",
     main = "Royal Guard", cex.lab = 1.25, cex.main = 1.35, cex.axis = 1.05,
     frame.plot = FALSE)

# --- Draw stacked smooth polygons (bottom -> top) ---
polygon(c(x_smooth, rev(x_smooth)), c(y0, rev(y1)), col = my.cols[1], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y1, rev(y2)), col = my.cols[2], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y2, rev(y3)), col = my.cols[3], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y3, rev(y4)), col = my.cols[4], border = NA)

# --- Compute label positions robustly from the actual band positions ---
# pick an index near the right-middle of the plot so labels sit in the wider part of each band
idx <- round(n_smooth * 0.65)         # ~65% along x range (adjust if you prefer left/right)
x_lbl <- x_smooth[idx]

# For each band compute a y position inside the band (midpoint between band's lower & upper boundaries)
# Bands: (y0 -> y1) = Fed (red), (y1 -> y2) = Exited without feeding (yellow),
#        (y2 -> y3) = Deterred (green), (y3 -> y4) = Killed (blue)
fed_y    <- (y0[idx] + y1[idx]) / 2
exit_y   <- (y1[idx] + y2[idx]) / 2
deter_y  <- (y2[idx] + y3[idx]) / 2
killed_y <- (y3[idx] + y4[idx]) / 2

# Slight horizontal nudges so overlapping labels don't collide
x_nudge <- 0.12 * (max(time) - min(time))  # 12% of x-range
x_lbl_fed    <- x_lbl + x_nudge * 0.6
x_lbl_exit   <- x_lbl - x_nudge * 0.1
x_lbl_deter  <- x_lbl - x_nudge * 0.6
x_lbl_killed <- x_lbl + x_nudge * 0.2

# --- Add the labels with contrasting text colors chosen for readability inside each band ---
text(x = x_lbl_killed, y = killed_y, labels = "Killed", col = "#03396C", cex = 1.05, font = 2)
#text(x = x_lbl_deter,  y = deter_y,  labels = "Deterred", col = "#1B5E20", cex = 1.05, font = 2)
text(x = x_lbl_exit,   y = exit_y,   labels = "Exited without feeding", col = "#8A4B00", cex = 1.00)
text(x = x_lbl_fed,    y = fed_y,    labels = "Successfully blood fed", col = "#7A0314", cex = 1.00)

# --- Y-axis as percentage ticks and soft frame outline ---
#axis(2, at = seq(0,1,0.25), labels = seq(0,100,25), las = 1)
#box(bty = "l", col = "gray70")


# --- Y-axis as percentage ---
#axis(2, at = seq(0,1,0.25), labels = seq(0,100,25))

# --- Frame outline ---
box(bty="l", col="gray70")







###############usingthe functionlike


# --- Parameters ---
gamma <- 0.363
gamma2<- 0.0472
n_smooth <- 300
x_smooth <- seq(0, 3, length.out = n_smooth)   # time 0–3

# --- Exponential decay formulas ---
Killed_s   <- 0.19* exp(-gamma * x_smooth)
Deterred_s <- 0.352* exp(-gamma * x_smooth)
Repelled_s <- (0.46 - 0.24) * exp(-gamma2 * x_smooth) + 0.24
Fed_s      <- 1 - Killed_s - Repelled_s   # remaining probability

# --- (Optional but safe) normalize ---
total <- Killed_s + Deterred_s + Repelled_s + Fed_s
Killed_s   <- Killed_s / total
Deterred_s <- Deterred_s / total
Repelled_s <- Repelled_s / total
Fed_s      <- Fed_s / total

# --- Cumulative stacking ---
y0 <- rep(0, n_smooth)
y1 <- Fed_s
y2 <- Fed_s + Repelled_s
y3 <- Fed_s + Repelled_s + Deterred_s
y4 <- Fed_s + Repelled_s + Deterred_s + Killed_s

# --- Colors ---
my.cols <- c(
  adjustcolor("#E15759", alpha.f = 0.8),  # Fed
  adjustcolor("#FAD02E", alpha.f = 0.6),  # Repelled
  adjustcolor("#7AD151", alpha.f = 0.6),  # Deterred
  adjustcolor("#6BAED6", alpha.f = 0.5)   # Killed
)

# --- Plot ---
par(bg = "white", mar = c(5,6,3,2), las = 1)
plot(x_smooth, x_smooth, type = "n",
     ylim = c(0,1), xlim = c(0,3),
     xlab = "Net Age (Years)",
     ylab = "Probable outcome of\n a single feeding attempt",
     main = "Pyrethroid-only", #pyriproxyfen net",
     frame.plot = FALSE)

# --- Smooth stacked polygons ---
polygon(c(x_smooth, rev(x_smooth)), c(y0, rev(y1)), col = my.cols[1], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y1, rev(y2)), col = my.cols[2], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y2, rev(y3)), col = my.cols[3], border = NA)
polygon(c(x_smooth, rev(x_smooth)), c(y3, rev(y4)), col = my.cols[4], border = NA)

# --- Label positions (automatic) ---
idx <- round(n_smooth * 0.65)
x_lbl <- x_smooth[idx]

fed_y    <- Fed_s[idx] / 2
repel_y  <- Fed_s[idx] + Repelled_s[idx] / 2
deter_y  <- Fed_s[idx] + Repelled_s[idx] + Deterred_s[idx] / 2
killed_y <- Fed_s[idx] + Repelled_s[idx] + Deterred_s[idx] + Killed_s[idx] / 2

# --- Labels ---
text(x_lbl, killed_y, "Killed", col = "#03396C", font = 2)
text(x_lbl, deter_y, "Deterred", col = "#1B5E20", font = 2)
text(x_lbl, repel_y, "Repelled", col = "#8A4B00")
text(x_lbl, fed_y, "Successfully blood fed", col = "#7A0314")

# --- Box ---
box(bty="l", col="gray70")











#####################################################################TANZANIA

