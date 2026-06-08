library(rstan)

# read the data file
dat = read.table("simul_data_hw6_model1.txt", header=T, sep="\t")

allSubjs = unique(dat$subjID)  # all subject IDs
N = length(allSubjs)      # number of subjects
T = table(dat$subjID)[1]  # number of trials per subject 

choice  <- array(-1, c(N, T))
outcome <- array(0, c(N, T))

for (i in 1:N) {
  curSubj = allSubjs[i]
  tmp     = subset(dat, subjID == curSubj)
  choice[i, 1:T] <- tmp$choice
  outcome[i, 1:T] <- tmp$outcome
}

dataList <- list(
  N       = N,
  T       = T,
  Tsubj   = rep(T, N),
  choice  = choice,
  outcome = outcome
)

# run!
output_nonhier <- stan(
  file = "hw6_model1_nonhier.stan",
  data = dataList,
  pars = c("alpha", "beta"),
  iter = 2000,
  warmup = 1000,
  chains = 2,
  cores = 2
)

# traceplot
traceplot(output_nonhier)

# print summary
print(output_nonhier)

# extract Stan fit object (parameters)
parameters_nonhier <- rstan::extract(output_nonhier)

alpha_mean_nonhier <- apply(parameters_nonhier$alpha, 2, mean)
alpha_sd_nonhier <- apply(parameters_nonhier$alpha, 2, sd)
beta_mean_nonhier <- apply(parameters_nonhier$beta, 2, mean)
beta_sd_nonhier <- apply(parameters_nonhier$beta, 2, sd)

plot(simul_pars$alpha, alpha_mean_nonhier, xlim = c(0, 0.5), ylim = c(0, 0.5))
abline(0, 1)
arrows(x0=simul_pars$alpha,y0=alpha_mean_nonhier - alpha_sd_nonhier,
       y1=alpha_mean_nonhier + alpha_sd_nonhier,
       length=0.02, angle=90, code=3)

plot(simul_pars$beta, beta_mean_nonhier, xlim = c(0, 4), ylim = c(0, 4))
abline(0, 1)
arrows(x0=simul_pars$beta, y0=beta_mean_nonhier - beta_sd_nonhier,
       y1=beta_mean_nonhier + beta_sd_nonhier,
       length=0.05, angle=90, code=3)