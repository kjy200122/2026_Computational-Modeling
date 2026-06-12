# hw6_model1_exec.R
# Programmed by Woo-Young Ahn (wahn55@snu.ac.kr)
 
#rm(list=ls())  # remove all variables

library(rstan)

# read the data file
dat = read.table("simul_data_hw6_model1_rein.txt", header=T, sep="\t")

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
output = stan("hw6_model1_rein.stan", data = dataList, pars = c("alpha_pos", "alpha_neg", "beta"),
              iter = 2000, warmup=1000, chains=2, cores=2)

# variational inference
m = stan_model("hw6_model1_rein.stan")
output_vb = vb(m, data = dataList, pars = c("alpha", "beta"),
)
# traceplot
traceplot(output)

# print summary
print(output)

# extract Stan fit object (parameters)
parameters <- rstan::extract(output)

alpha_pos_mean <- apply(parameters$alpha_pos, 2, mean)
alpha_pos_sd   <- apply(parameters$alpha_pos, 2, sd)
alpha_neg_mean <- apply(parameters$alpha_neg, 2, mean)
alpha_neg_sd   <- apply(parameters$alpha_neg, 2, sd)
beta_mean      <- apply(parameters$beta, 2, mean)
beta_sd        <- apply(parameters$beta, 2, sd)

# scatter plot - alpha_pos
plot(simul_pars$alpha_pos, alpha_pos_mean, xlim=c(0,0.5), ylim=c(0,0.5)); abline(0,1)
arrows(x0=simul_pars$alpha_pos, y0=alpha_pos_mean-alpha_pos_sd, y1=alpha_pos_mean+alpha_pos_sd, length=0.02, angle=90, code=3)

# scatter plot - alpha_neg
plot(simul_pars$alpha_neg, alpha_neg_mean, xlim=c(0,0.5), ylim=c(0,0.5)); abline(0,1)
arrows(x0=simul_pars$alpha_neg, y0=alpha_neg_mean-alpha_neg_sd, y1=alpha_neg_mean+alpha_neg_sd, length=0.02, angle=90, code=3)

plot(simul_pars$beta, beta_mean, xlim=c(0, 4), ylim=c(0, 4)); abline(0,1)
arrows(x0=simul_pars$beta, y0=beta_mean-beta_sd, y1=beta_mean+beta_sd, length=0.05, angle=90, code=3)

