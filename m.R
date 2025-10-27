library(readr)
library(rstan)
library(loo)
library(dplyr)  
library(openxlsx)

# Define the DETA_ijv function to calculate order intervals and DETA array
# T: Matrix of response times (rows = students, columns = items)
# S: Vector of interval breakpoints (quantiles of response times)
DETA_ijv <- function(T, S){
  n_student <- nrow(T)  # Number of students
  n_item <- ncol(T)     # Number of items
  V <- length(S)        # Number of intervals (V intervals correspond to V+1 breakpoints)
  
  # Initialize DETA array (dimensions: students × items × intervals) and order interval matrix
  DETA <- array(0, dim = c(n_student, n_item, V))
  order_interval <- matrix(, n_student, n_item)
  
  # Iterate over each student and item to determine interval membership
  for(i in 1:n_student){
    for(j in 1:n_item){
      # Find the first interval where response time T[i,j] is ≤ the breakpoint S
      v <- which(T[i, j] - S <= 0)[1]
      order_interval[i, j] <- v  # Record the interval index for (student i, item j)
      DETA[i, j, v] <- 1         # One-hot encode interval membership in DETA array
    }}
  return(list(order = order_interval, DETA = DETA))  # Return results as a list
}

# Load response data (resp.csv) and response time data (rt.csv)
res <- read_csv("D:\\res.csv")
rt <- read_csv("D:\\rt.csv")

# Convert data to matrices, excluding the first column (assumed to be non-data columns like IDs)
rt <- as.matrix(rt[2:length(rt)])
res <- as.matrix(res[2:length(res)])

# Check dimensions of rt and res matrices to confirm data loading
dim(rt)
dim(res)

# Extract sample size (N = number of students) and item count (J = number of items)
size <- dim(rt)
N <- size[1]
J <- size[2]

# Set number of intervals (V) and random seed for reproducibility
V <- 5
set.seed(V)

# Print current working directory and interval count for verification
print(getwd())
print(paste0('V=', V))

# Calculate interval breakpoints: quantiles of response times (exclude 0th and 100th percentiles)
S <- quantile(rt, probs = seq(0, 1, length = V + 1))[-c(1, V + 1)]

# Verify length of breakpoints and print values
length(S)
S

# Extend breakpoints with max response time + 0.1 to cover upper bound of longest response times
s <- c(S, max(rt) + 0.1)
s

# Save breakpoint vectors to CSV files for later reference
write.csv(s, "s.csv")
write.csv(S, "SS.csv")

# Use DETA_ijv function to get interval assignments for each (student, item) pair
OR <- DETA_ijv(rt, s)$order
write.csv(OR, "OR_Pieces.csv")

# Prepare input data list for the main Stan model
data <- list(
  N = N,          # Number of students
  J = J,          # Number of items
  V = V,          # Number of intervals
  rt = rt,        # Response time matrix
  resp = res,     # Response data matrix
  s = s           # Extended breakpoint vector
)
data  # Print data list to confirm structure

# Define initial parameter values for the main Stan model (single chain)
init <- list(list(
  theta = rnorm(N),          # Student ability parameters (normal prior)
  a = runif(J, 0.5, 1.5),    # Item discrimination parameters (uniform prior)
  b = rnorm(J),              # Item difficulty parameters (normal prior)
  mu_b = 0,                  # Mean of item difficulty distribution (fixed prior mean)
  sigma2_b = 1,              # Variance of item difficulty distribution (fixed prior variance)
  phi = runif(J, 0.5, 1.5),  # Item speed parameters (uniform prior)
  zeta = rnorm(J),           # Item time threshold parameters (normal prior)
  tau = rnorm(N),            # Student speed parameters (normal prior)
  varphi = 0.5,              # Scale parameter for response time model (fixed initial value)
  gamma = rep(1, J),         # Item interval weight parameters (fixed initial value = 1)
  lambda = rep(1, V)         # Interval weight parameters (fixed initial value = 1)
))

init  # Print initial values to confirm

# Record start time to calculate model runtime
fir <- Sys.time()
fir

# Run main Stan model: 1 chain, 6000 total iterations (2000 warmup + 4000 sampling)
mcmc.out <- stan(
  file = "m.stan", 
  data = data, 
  chains = 1, 
  iter = 6000, 
  init = init, 
  warmup = 2000
)

# Record end time and calculate runtime
end <- Sys.time()
end - fir

# Extract posterior samples from the main model
samples <- extract(mcmc.out) 

# Save posterior samples to CSV (filename includes sample size and item count)
write.csv(samples, paste0('samples','-', N, '-', J, '.csv'))

# Generate summary statistics for model parameters
par.summary <- summary(mcmc.out)
all.params <- as.data.frame(par.summary)

# Save parameter summaries to CSV
write.csv(all.params, paste0('allparams','-', N, '-', J, '.csv'))

# Re-load posterior samples (for reprocessing if needed)
samples <- read_csv(paste0('samples','-', N, '-', J, '.csv'))

# Extract specific parameter subsets from posterior samples
mcmc_theta <- samples[, 1:N + 1]                          # Student ability (theta)
mcmc_tau <- samples[, (N + 2):(2*N + 1)]                  # Student speed (tau)
mcmc_a <- samples[, (2*N + 2):(2*N + J + 1)]              # Item discrimination (a)
mcmc_b <- samples[, (2*N + J + 2):(2*N + 2*J + 1)]        # Item difficulty (b)
mcmc_gamma <- samples[, (2*N + 2*J + 2):(2*N + 3*J + 1)]  # Item interval weights (gamma)
mcmc_phi <- samples[, (3*J + 2*N + 2):(2*N + 4*J + 1)]    # Item speed (phi)
mcmc_zeta <- samples[, (2*N + 4*J + 2):(5*J + 2*N + 1)]   # Item time thresholds (zeta)
mcmc_varphi <- samples[, 2*N + 5*J + 2]                   # Response time scale (varphi)
mcmc_lambda <- samples[, (ncol(samples) - V):(ncol(samples)-1)]  # Interval weights (lambda)

# ------------------------------
# Conditional sampling 1: Estimate response parameters GIVEN response times
# ------------------------------
# Calculate EAP (Expected A Posteriori) estimates for fixed parameters
eap_a <- colMeans(mcmc_a)        # EAP of item discrimination (a)
eap_b <- colMeans(mcmc_b)        # EAP of item difficulty (b)
eap_th <- colMeans(mcmc_theta)   # EAP of student ability (theta)
eap_gamma <- colMeans(mcmc_gamma)# EAP of item interval weights (gamma)
eap_lambda <- colMeans(mcmc_lambda)# EAP of interval weights (lambda)
eap_phi <- colMeans(mcmc_phi)    # EAP of item speed (phi)
eap_tau <- colMeans(mcmc_tau)    # EAP of student speed (tau)
eap_zeta <- colMeans(mcmc_zeta)  # EAP of item time thresholds (zeta)
eap_varphi <- colMeans(mcmc_varphi)# EAP of response time scale (varphi)

# Define parameter names for EAP results (for readability)
param_names <- c(
  paste0("a_", 1:J),            # Item discrimination (a_1 to a_J)
  paste0("b_", 1:J),            # Item difficulty (b_1 to b_J)
  paste0("theta_", 1:N),        # Student ability (theta_1 to theta_N)
  paste0("gamma_", 1:J),        # Item interval weights (gamma_1 to gamma_J)
  paste0("lambda_", 1:V),       # Interval weights (lambda_1 to lambda_V)
  paste0("phi_", 1:J),          # Item speed (phi_1 to phi_J)
  paste0("tau_", 1:N),          # Student speed (tau_1 to tau_N)
  paste0("zeta_", 1:J),         # Item time thresholds (zeta_1 to zeta_J)
  "varphi"                      # Response time scale (varphi)
)

# Combine EAP values into a single vector
param_values <- c(
  eap_a,
  eap_b,
  eap_th,
  eap_gamma,
  eap_lambda,
  eap_phi,
  eap_tau,
  eap_zeta,
  eap_varphi
)

# Create data frame of EAP estimates (parameter name + value)
eap_df <- data.frame(
  parameter = param_names,
  value = param_values
)

# Save EAP results to CSV (filename includes sample size and item count)
write.csv(
  eap_df, 
  paste0('all_eap_results-', N, '-', J, '.csv'), 
  row.names = FALSE
)

# Prepare input data for conditional model (fixed parameters = EAP estimates)
data1 <- list(
  N = N,
  J = J,
  V = V,
  rt = rt,
  resp = res,
  s = s,
  constant_gamma = eap_gamma,  # Fixed: EAP of gamma
  constant_phi = eap_phi,      # Fixed: EAP of phi
  constant_zeta = eap_zeta,    # Fixed: EAP of zeta
  constant_varphi = eap_varphi,# Fixed: EAP of varphi
  constant_lambda = eap_lambda # Fixed: EAP of lambda
)

# Define initial values for conditional model 1 (estimating theta, a, b, mu_b, sigma2_b, tau)
init1 <- list(list(
  theta = rnorm(N),          # Student ability (theta)
  a = runif(J, 0.5, 1.5),    # Item discrimination (a)
  b = rnorm(J),              # Item difficulty (b)
  mu_b = 0,                  # Mean of item difficulty distribution
  sigma2_b = 1,              # Variance of item difficulty distribution
  tau = rnorm(N)             # Student speed (tau)
))

# Run conditional Stan model 1: 1 chain, 4000 total iterations (2000 warmup + 2000 sampling)
mcmc.out1 <- stan(
  file = "m_resp_give_rt.stan", 
  data = data1, 
  chains = 1, 
  iter = 4000, 
  init = init1, 
  warmup = 2000
)

# Extract posterior samples from conditional model 1
samples1 <- extract(mcmc.out1)

# Plot traceplots for item parameters (a, b) to check convergence (exclude warmup)
traceplot(mcmc.out1, pars = c("b", "a"), inc_warmup = FALSE)

# Save samples and parameter summaries from conditional model 1
write.csv(samples1, "samples_res_give_rt.csv")
par.summary1 <- summary(mcmc.out1)
all.params2 <- as.data.frame(par.summary1)
write.csv(all.params2, "allparams_res_give_rt.csv")

# ------------------------------
# Conditional sampling 2: Estimate response time parameters GIVEN responses
# ------------------------------
# Prepare input data for conditional model 2 (fixed parameters = EAP estimates of a, b)
data2 <- list(
  N = N,
  J = J,
  V = V,
  rt = rt,
  resp = res,
  s = s,
  constant_a = eap_a,   # Fixed: EAP of item discrimination (a)
  constant_b = eap_b    # Fixed: EAP of item difficulty (b)
)

# Define initial values for conditional model 2 (estimating theta, mu_b, sigma2_b, phi, zeta, tau, varphi, gamma, lambda)
init2 <- list(list(
  theta = rnorm(N),          # Student ability (theta)
  mu_b = 0,                  # Mean of item difficulty distribution
  sigma2_b = 1,              # Variance of item difficulty distribution
  phi = runif(J, 0.5, 1.5),  # Item speed (phi)
  zeta = rnorm(J),           # Item time thresholds (zeta)
  tau = rnorm(N),            # Student speed (tau)
  varphi = 0.5,              # Response time scale (varphi)
  gamma = rep(1, J),         # Item interval weights (gamma)
  lambda = rep(1, V)         # Interval weights (lambda)
))

# Run conditional Stan model 2: 1 chain, 8000 total iterations (2000 warmup + 6000 sampling)
mcmc.out2 <- stan(
  file = "m_rt_give_resp.stan", 
  data = data2, 
  chains = 1, 
  iter = 8000, 
  init = init2, 
  warmup = 2000
)

# Extract posterior samples from conditional model 2
samples2 <- extract(mcmc.out2)

# Plot traceplot for item interval weights (gamma) to check convergence (exclude warmup)
traceplot(mcmc.out2, pars = c("gamma"), inc_warmup = FALSE)

# Save samples and parameter summaries from conditional model 2
write.csv(samples2, "samples_rt_give_resp.csv")
par.summary2 <- summary(mcmc.out2)
all.params3 <- as.data.frame(par.summary2)
write.csv(all.params3, "allparams_rt_give_resp.csv")