# Load necessary libraries
library(dplyr)  # for data manipulation (not needed for basic operations)
#library(Metrics)  # for handling logit functions, if needed

# Given parameters
lp <-  0.24495673  #l_p value
rho <- -3.01       # ρ_p value
mu <- -2.43        # μ_p value
tau <- 0.5         # τ value

# Calculate logit(γ_p)
logit_gamma_p <- mu + rho * (lp - tau)

# Convert from logit to probability (γ_p)
gamma_p <- 1 / (1 + exp(-logit_gamma_p))

# Calculate the half-life HW
HW <- log(2) / gamma_p

# Print the results
cat("logit(γ_p):", logit_gamma_p, "\n")
cat("γ_p:", gamma_p, "\n")
cat("Half-life (HW):", HW, "\n")

