# hw6_model1_exec.R
# Programmed by Woo-Young Ahn (wahn55@snu.ac.kr)
 
rm(list=ls())  # remove all variables

library(rstan)
set.seed(2026)
# generate data using R files
source("simulate_hw6_model1.R")

# read the data file
dat = read.table("simul_data_hw6_model1.txt", header=T, sep="\t")
allSubjs = unique(dat$subjID)  # all subject IDs
N = length(allSubjs)      # number of subjects
T = table(dat$subjID)[1]  # number of trials per subject 

choice  <- array(-1, c(N, T))
outcome <- array(0, c(N, T))

print(N)
print(T)
print(mean(dat$outcome))

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
output = stan("hw6_model_indiv.stan", data = dataList, pars = c("alpha", "beta"),
              iter = 2000, warmup=1000, chains=2, cores=2)

# save fitted model
# saveRDS(output, "output_indiv.rds")

# load saved model
#output = readRDS("output_indiv.rds")

# variational inference
#m = stan_model("hw6_model1.stan")
#output_vb = vb(m, data = dataList, pars = c("alpha", "beta"), #c("mu_alpha", "mu_beta", "alpha", "beta"),
#)

# traceplot
traceplot(output)

# print summary
print(output)

rstan::check_divergences(output)

# check parameter with Rhat > 1.01
summary_output = summary(output)$summary
print(rownames(summary_output)[summary_output[, "Rhat"] > 1.01]) #NO parameters with Rhat>1.01


# extract Stan fit object (parameters)
parameters <- rstan::extract(output)

#1.(2a) posterior distribution of individual and group parameters
hist(parameters$mu_alpha, main = "mu_alpha", xlab="values")
hist(parameters$mu_beta, main = "mu_beta", xlab="values")

alpha_mean = apply(parameters$alpha, 2, mean)
alpha_sd = apply(parameters$alpha, 2, sd)

beta_mean = apply(parameters$beta, 2, mean)
beta_sd = apply(parameters$beta, 2, sd)

hist(x= alpha_mean, breaks=10, main="mean of individual alpha")
hist(x= beta_mean, breaks=10, main="mean of individual beta")

plot(simul_pars$alpha, 
     alpha_mean, 
     xlim=c(0, 0.5), 
     ylim=c(0, 0.5),
     xlab="true alpha",
     ylab="estimated alpha",
     main = "alpha"
)

abline(0,1)
arrows(x0=simul_pars$alpha, y0= alpha_mean - alpha_sd, y1= alpha_mean + alpha_sd, length=0.02, angle=90, code=3)

plot(simul_pars$beta, 
     beta_mean, 
     xlim=c(0, 4), 
     ylim=c(0, 4),
     xlab="true beta",
     ylab="estimated beta",
     main="beta")
abline(0,1)
arrows(x0=simul_pars$beta, y0=beta_mean-beta_sd, y1=beta_mean+beta_sd, length=0.05, angle=90, code=3)

recovered_alpha = sum((simul_pars$alpha < alpha_mean+alpha_sd)&(alpha_mean-alpha_sd < simul_pars$alpha))
recovered_beta = sum((simul_pars$beta < beta_mean+beta_sd)&(beta_mean-beta_sd < simul_pars$beta))

print(recovered_alpha)
print(recovered_beta)
print(mean((alpha_mean - simul_pars$alpha)^2))
print(mean((beta_mean - simul_pars$beta)^2))
