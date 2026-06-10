# hw6_model1_exec.R
# Programmed by Woo-Young Ahn (wahn55@snu.ac.kr)
 
#rm(list=ls())  # remove all variables

library(rstan)

# read the data file
dat = read.table("simul_data_hw6_Q4_2.txt", header=T, sep="\t")

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
output_Q4_2 = stan("hw6_Q4.stan", data = dataList, pars = c("mu_alpha_pos", "mu_alpha_neg", "mu_beta", "sigma", 
                                                          "alpha_pos", "alpha_neg", "beta"),
              iter = 2000, warmup=1000, chains=2, cores=2)

# # variational inference
# m = stan_model("hw6_model1.stan")
# output_vb = vb(m, data = dataList, pars = c("alpha", "beta"), #c("mu_alpha", "mu_beta", "alpha", "beta"),
# )

save(output_Q4_2, file = "output_Q4_2.RData")

# traceplot
traceplot(output_Q4_2)
traceplot(output_Q4_2, pars = c("alpha"))
traceplot(output_Q4_2, pars = c("beta"))


# print summary
print(output_Q4_2)

# extract Stan fit object (parameters)
parameters <- rstan::extract(output_Q4_2)

alpha_pos_mean = apply(parameters$alpha_pos, 2, mean)
alpha_pos_sd = apply(parameters$alpha_pos, 2, sd)
alpha_neg_mean = apply(parameters$alpha_neg, 2, mean)
alpha_neg_sd = apply(parameters$alpha_neg, 2, sd)

beta_mean = apply(parameters$beta, 2, mean)
beta_sd = apply(parameters$beta, 2, sd)

plot(simul_pars$alpha_pos, alpha_pos_mean, xlim=c(0, 0.5), ylim=c(0, 0.5)); abline(0,1)
arrows(x0=simul_pars$alpha_pos, y0= alpha_pos_mean - alpha_pos_sd, y1= alpha_pos_mean + alpha_pos_sd, length=0.02, angle=90, code=3)
plot(simul_pars$alpha_neg, alpha_neg_mean, xlim=c(0, 0.5), ylim=c(0, 0.5)); abline(0,1)
arrows(x0=simul_pars$alpha_neg, y0= alpha_neg_mean - alpha_neg_sd, y1= alpha_neg_mean + alpha_neg_sd, length=0.02, angle=90, code=3)

plot(simul_pars$beta, beta_mean, xlim=c(0, 4), ylim=c(0, 4)); abline(0,1)
arrows(x0=simul_pars$beta, y0=beta_mean-beta_sd, y1=beta_mean+beta_sd, length=0.05, angle=90, code=3)

