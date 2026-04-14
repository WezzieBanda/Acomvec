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

# Compute posterior predicted probabilities for each combo
results_time <- data.frame()

for (i in 1:nrow(combo_grid)) {
  t <- combo_grid$treatment[i]
  tm <- combo_grid$time[i]
  
  # Linear predictor for all posterior draws
  eta <- posterior$alpha + 
    posterior$beta_treatment[, t] + 
    posterior$beta_time * tm
  
  # Convert to probability using logistic function
  p <- plogis(eta)
  
  # Summarize posterior distribution
  summary_row <- data.frame(
    treatment = t,
    time = tm,
    mean = mean(p),
    lower = quantile(p, 0.025),
    upper = quantile(p, 0.975)
  )
  
  results_time <- rbind(results_time, summary_row)
}

# Map treatment index back to names
treatment_labels <- levels(df$Treatment)
results_time$treatment <- treatment_labels[results_time$treatment]

# Optional: label times
results_time$time <- factor(results_time$time,
                            labels = sort(unique(df$Time)))

# View the summarized results
print(results_time)

library(writexl)

# Save posterior summary results per treatment per time to Excel
write_xlsx(results_time, "motold.xlsx")















#polygon plot



# --- Libraries ---
library(RColorBrewer)

# --- Time points (x-axis) ---
time <- c(1, 2, 3,4)   # e.g. 0, 10, 20 washes
n_smooth <- 300
x_smooth <- seq(min(time), max(time), length.out = n_smooth)




#another set


# --- Example probabilities (replace with your real data if needed) ---
Killed    <- c(0.118365097, 0.087817133, 0.056130924, 0.03641781) #, 0.105202576)
Deterred  <- c(0.339, 0.339, 0.339, 0.339) #, 0.0193)
Repelled  <- c(0.567303213, 0.583074198, 0.699270331, 0.554160257) #, 0.617404221)
Fed       <- c(0.31433169, 0.329108669, 0.244598745, 0.409421933) 

#realIG1 set
# --- Example probabilities (replace with your real data if needed) ---
Killed    <- c(0.078239329, 0.067355741, 0.036092184, 0.025929481) #, 0.105202576)
Deterred  <- c(0.339, 0.339, 0.339, 0.339) #, 0.0193)
Repelled  <- c(0.374987424, 0.44721791, 0.449630823, 0.394562103) #, 0.617404221)
Fed       <- c(0.207773247, 0.252426349, 0.157276993, 0.291508416) 

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
par(bg = "white", mar = c(5,5,3,2), las = 1)
plot(x_smooth, x_smooth, type = "n", ylim = c(0,1), xlim = range(time),
     xlab = "Time ", ylab = "Probability (%)",
     main = "Pyrethroid-only Net", cex.lab = 1.25, cex.main = 1.35, cex.axis = 1.05,
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


###trying others



















####################IG2################################################
############################################################################

# --- Libraries ---
library(RColorBrewer)

# --- Time points (x-axis) ---
time <- c(1, 2, 3,4)   # e.g. 0, 10, 20 washes
n_smooth <- 300
x_smooth <- seq(min(time), max(time), length.out = n_smooth)


#another set IG2


# --- Example probabilities (replace with your real data if needed) ---
Killed    <- c(0.315539497, 0.240175734, 0.225942114, 0.194298385) #, 0.105202576)
Deterred  <- c(0.336, 0.372, 0.28, 0.109) #, 0.0193)
Repelled  <- c(0.453319954, 0.56507877, 0.578755715, 0.41980714) #, 0.617404221)
Fed       <- c(0.231140548, 0.194745496, 0.19530217, 0.385894475) 



#REALIG
# --- Example probabilities (replace with your real data if needed) ---
Killed    <- c(0.209518226, 0.150830361, 0.162678322, 0.173119861) #, 0.105202576)
Deterred  <- c(0.336, 0.372, 0.28, 0.109) #, 0.0193)
Repelled  <- c(0.30100445, 0.354869468, 0.416704115, 0.374048162) #, 0.617404221)
Fed       <- c(0.153477324, 0.122300171, 0.140617563, 0.343831977) 

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
par(bg = "white", mar = c(5,5,3,2), las = 1)
plot(x_smooth, x_smooth, type = "n", ylim = c(0,1), xlim = range(time),
     xlab = "Time ", ylab = "Probability (%)",
     main = "Pyrethroid-Pyrrole Net", cex.lab = 1.25, cex.main = 1.35, cex.axis = 1.05,
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


###trying others
























####################P3################################################
############################################################################

# --- Libraries ---
library(RColorBrewer)

# --- Time points (x-axis) ---
time <- c(1, 2, 3,4)   # e.g. 0, 10, 20 washes
n_smooth <- 300
x_smooth <- seq(min(time), max(time), length.out = n_smooth)


#another set P3
R,S,D,
0.664184314	0.154951607	0.180864079	0.449
0.62705913	0.195090727	0.177850142	0.216
0.64465229	0.197585226	0.157762484	0.086
0.612830536	0.291327124	0.09584234	0.18

##Another ser

# --- Example probabilities (replace with your real data if needed) ---
Killed    <- c(0.180864079, 0.177850142, 0.157762484, 0.09584234) #, 0.105202576)
Deterred  <- c(0.449, 0.216, 0.086, 0.18) #, 0.0193)
Repelled  <- c(0.664184314, 0.62705913, 0.64465229, 0.612830536) #, 0.617404221)
Fed       <- c(0.154951607, 0.195090727, 0.197585226, 0.291327124) 

#pr real set

Killed    <- c(0.099656107, 0.139434512, 0.144194911, 0.078590719) #, 0.105202576)
Deterred  <- c(0.449, 0.216, 0.086, 0.18) #, 0.0193)
Repelled  <- c(0.365965557, 0.491614358, 0.589212193, 0.502521039) #, 0.617404221)
Fed       <- c(0.085378335, 0.15295113, 0.180592896, 0.238888242) 



#h interpolation using natural spline (works with 3 points) ---
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
par(bg = "white", mar = c(5,5,3,2), las = 1)
plot(x_smooth, x_smooth, type = "n", ylim = c(0,1), xlim = range(time),
     xlab = "Time ", ylab = "Probability (%)",
     main = "Pyrethroid-PBO Net", cex.lab = 1.25, cex.main = 1.35, cex.axis = 1.05,
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














####################RG################################################
############################################################################

# --- Libraries ---
library(RColorBrewer)

# --- Time points (x-axis) ---
time <- c(1, 2, 3,4)   # e.g. 0, 10, 20 washes
n_smooth <- 300
x_smooth <- seq(min(time), max(time), length.out = n_smooth)




##Another ser


Killed    <- c(0.301035734, 0.17961269, 0.14834662, 0.107335176) #, 0.105202576)
Deterred  <- c(0, 0.1554, 0.081, 0.0193) #, 0.0193)
Repelled  <- c(0.490784752, 0.665895573, 0.706693018, 0.630287511) #, 0.617404221)
Fed       <- c(0.208179514, 0.154491736, 0.144960362, 0.262377313) 

#RG real set


Killed    <- c(0.301035734, 0.151700878, 0.136330544, 0.105263607) #, 0.105202576)
Deterred  <- c(0, 0.1554, 0.081, 0.0193) #, 0.0193)
Repelled  <- c(0.490784752, 0.562415401, 0.649450884, 0.618122962) #, 0.617404221)
Fed       <- c(0.208179514, 0.130483721, 0.133218572, 0.25731343) 



#h interpolation using natural spline (works with 3 points) ---
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
par(bg = "white", mar = c(5,5,3,2), las = 1)
plot(x_smooth, x_smooth, type = "n", ylim = c(0,1), xlim = range(time),
     xlab = "Time ", ylab = "Probability (%)",
     main = "Pyrethroid-pyriproxyfen Net", cex.lab = 1.25, cex.main = 1.35, cex.axis = 1.05,
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











































#####################now lwtmedo pertime



