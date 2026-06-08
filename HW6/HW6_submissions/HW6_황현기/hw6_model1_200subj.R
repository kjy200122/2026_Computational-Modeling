# hw6_model1_exec.R
# Programmed by Woo-Young Ahn (wahn55@snu.ac.kr)
 
# rm(list=ls())  # remove all variables

library(rstan)

# read the data file
dat = read.table("simul_data_hw6_model1_200subj.txt", header=TRUE, sep="\t")

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
output_200subj = stan("hw6_model1_200subj.stan", data = dataList, pars = c("mu_p", "sigma", "alpha", "beta"),
              iter = 2000, warmup=1000, chains=2, cores=2)

# variational inference
# m = stan_model("hw6_model1.stan")
# output_vb = vb(m, data = dataList, pars = c("alpha", "beta"), #c("mu_alpha", "mu_beta", "alpha", "beta"),)

# traceplot
traceplot(output_200subj, pars = c("alpha[1]", "alpha[50]", "alpha[100]", "alpha[150]", "alpha[200]"))
traceplot(output_200subj, pars = c("beta[1]", "beta[50]", "beta[100]", "beta[150]", "beta[200]"))

# print summary
print(output_200subj)

# extract Stan fit object (parameters)
parameters_200subj <- rstan::extract(output_200subj)

alpha_mean_200subj = apply(parameters_200subj$alpha, 2, mean)
alpha_sd_200subj = apply(parameters_200subj$alpha, 2, sd)
beta_mean_200subj = apply(parameters_200subj$beta, 2, mean)
beta_sd_200subj = apply(parameters_200subj$beta, 2, sd)

# plot posterior distributions
hist(parameters_200subj$alpha[, 1])
hist(parameters_200subj$beta[, 1])
hist(parameters_200subj$mu_p[ ,1])
hist(parameters_200subj$sigma[ ,1])

plot(simul_pars$alpha, alpha_mean_200subj, xlim=c(0, 0.5), ylim=c(0, 0.5)); abline(0,1)
arrows(x0=simul_pars$alpha, y0= alpha_mean_200subj - alpha_sd_200subj, y1= alpha_mean_200subj + alpha_sd_200subj, length=0.02, angle=90, code=3)

plot(simul_pars$beta, beta_mean_200subj, xlim=c(0, 4), ylim=c(0, 4)); abline(0,1)
arrows(x0=simul_pars$beta, y0=beta_mean_200subj - beta_sd_200subj, y1=beta_mean_200subj + beta_sd_200subj, length=0.05, angle=90, code=3)

