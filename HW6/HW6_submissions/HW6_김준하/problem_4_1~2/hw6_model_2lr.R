# hw6_model1_exec.R
# Programmed by Woo-Young Ahn (wahn55@snu.ac.kr)
 
rm(list=ls())  # remove all variables

source("simulate_hw6_model_2lr.R")
library(rstan)
# read the data file
dat = read.table("simul_data_hw6_model_2lr.txt", header=T, sep="\t")

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
output = stan("hw6_model_2lr.stan", data = dataList, pars = c("mu_alpha_pos", "mu_alpha_neg", "mu_beta", "sigma", "alpha_pos", "alpha_neg", "beta"),
              iter = 2000, warmup=1000, chains=2, cores=2)

# # save fitted model
# saveRDS(output, "output_2lr.rds")

# # load saved output
# output = readRDS("output_2lr.rds")

# variational inference
# m = stan_model("hw6_model1.stan")
# output_vb = vb(m, data = dataList, pars = c("alpha", "beta"), #c("mu_alpha", "mu_beta", "alpha", "beta"),
# )
# traceplot
traceplot(output)

# print summary
print(output)

rstan::check_divergences(output)

# check parameter with Rhat > 1.01
summary_output = summary(output)$summary
print(rownames(summary_output)[summary_output[, "Rhat"] > 1.01]) #NO parameters with Rhat>1.01
print(summary_output[summary_output[, "Rhat"]>1.01, "Rhat"])

# extract Stan fit object (parameters)
parameters <- rstan::extract(output)

# posterior distribution of individual and group parameters
hist(parameters$mu_alpha_pos, main = "mu_alpha_pos", xlab="values")
hist(parameters$mu_alpha_neg, main = "mu_alpha_neg", xlab="values")
hist(parameters$mu_beta, main = "mu_beta", xlab="values")

alpha_pos_mean = apply(parameters$alpha_pos, 2, mean)
alpha_pos_sd = apply(parameters$alpha_pos, 2, sd)

alpha_neg_mean = apply(parameters$alpha_neg, 2, mean)
alpha_neg_sd = apply(parameters$alpha_neg, 2, sd)

beta_mean = apply(parameters$beta, 2, mean)
beta_sd = apply(parameters$beta, 2, sd)

hist(x= alpha_pos_mean, breaks=5, main="mean of individual alpha_pos")
hist(x= alpha_neg_mean, breaks=10, main="mean of individual alpha_neg")
hist(x= beta_mean, breaks=5, main="mean of individual beta")

plot(simul_pars$alpha_pos, 
     alpha_pos_mean, 
     xlim=c(0, 0.5), 
     ylim=c(0, 0.5),
     xlab="true alpha_pos",
     ylab="estimated alpha_pos",
     main = "alpha_pos"
)

## alpha_pos
abline(0,1)
arrows(x0=simul_pars$alpha_pos, 
       y0= alpha_pos_mean - alpha_pos_sd, 
       y1= alpha_pos_mean + alpha_pos_sd, 
       length=0.02, angle=90, code=3)

## alpha_neg
plot(simul_pars$alpha_neg, 
     alpha_neg_mean, 
     xlim=c(0, 0.5), 
     ylim=c(0, 0.5),
     xlab="true alpha_neg",
     ylab="estimated alpha_neg",
     main = "alpha_neg"
)

abline(0,1)
arrows(x0=simul_pars$alpha_neg, 
       y0= alpha_neg_mean - alpha_neg_sd, 
       y1= alpha_neg_mean + alpha_neg_sd, 
       length=0.02, angle=90, code=3)

## beta
plot(simul_pars$beta, 
     beta_mean, 
     xlim=c(0, 4), 
     ylim=c(0, 4),
     xlab="true beta",
     ylab="estimated beta",
     main="beta")
abline(0,1)
arrows(x0=simul_pars$beta, 
       y0=beta_mean-beta_sd, 
       y1=beta_mean+beta_sd, 
       length=0.05, angle=90, code=3)

recovered_alpha_pos = sum((simul_pars$alpha_pos < alpha_pos_mean+alpha_pos_sd)
                          &(alpha_pos_mean-alpha_pos_sd < simul_pars$alpha_pos))

recovered_alpha_neg = sum((simul_pars$alpha_neg < alpha_neg_mean+alpha_neg_sd)
                          &(alpha_neg_mean-alpha_neg_sd < simul_pars$alpha_neg))

recovered_beta = sum((simul_pars$beta < beta_mean+beta_sd)&(beta_mean-beta_sd < simul_pars$beta))

print(recovered_alpha_pos)
print(recovered_alpha_neg)
print(recovered_beta)
print(mean((alpha_pos_mean - simul_pars$alpha_pos)^2))
print(mean((alpha_neg_mean - simul_pars$alpha_neg)^2))
print(mean((beta_mean - simul_pars$beta)^2))

