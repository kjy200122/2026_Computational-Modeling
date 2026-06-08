# hw6_model1_exec.R
# Programmed by Woo-Young Ahn (wahn55@snu.ac.kr)
 
#rm(list=ls())  # remove all variables

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
output = stan("hw6_model1.stan", data = dataList, pars = c("mu_p", "sigma", "alpha_pos", "alpha_neg", "beta"), #c("mu_alpha", "mu_beta", "alpha", "beta"),
              iter = 2000, warmup=1000, chains=2, cores=2)

# variational inference
m = stan_model("hw6_model1.stan")

# traceplot
traceplot(output)

# print summary
print(output)

# extract Stan fit object (parameters)
parameters <- rstan::extract(output)

mu_alpha_pos <- mean(pnorm(parameters$mu_p[, 1]))
mu_alpha_neg <- mean(pnorm(parameters$mu_p[, 2]))
mu_beta      <- mean(pnorm(parameters$mu_p[, 3])) * 5

alpha_pos_mean <- apply(parameters$alpha_pos, 2, mean)
alpha_pos_sd   <- apply(parameters$alpha_pos, 2, sd)
alpha_neg_mean <- apply(parameters$alpha_neg, 2, mean)
alpha_neg_sd   <- apply(parameters$alpha_neg, 2, sd)
beta_mean      <- apply(parameters$beta,      2, mean)
beta_sd        <- apply(parameters$beta,      2, sd)


plot(simul_pars$alpha_pos, alpha_pos_mean,
     xlim=c(0,0.5), ylim=c(0,0.5),
     xlab="True Alpha_pos",
     ylab="Estimated Alpha_pos (posterior mean)",
     main="Parameter Recovery: Alpha_pos (num 3)"); abline(0,1)
arrows(x0=simul_pars$alpha_pos,
       y0=alpha_pos_mean - alpha_pos_sd,
       y1=alpha_pos_mean + alpha_pos_sd,
       length=0.02, angle=90, code=3)

plot(simul_pars$alpha_neg, alpha_neg_mean,
     xlim=c(0,0.6), ylim=c(0,0.6),
     xlab="True Alpha_neg",
     ylab="Estimated Alpha_neg (posterior mean)",
     main="Parameter Recovery: Alpha_neg (num 3)"); abline(0,1)
arrows(x0=simul_pars$alpha_neg,
       y0=alpha_neg_mean - alpha_neg_sd,
       y1=alpha_neg_mean + alpha_neg_sd,
       length=0.02, angle=90, code=3)

# beta plot
plot(simul_pars$beta, beta_mean,
     xlim=c(0,4), ylim=c(0,4),
     xlab="True Beta",
     ylab="Estimated Beta (posterior mean)",
     main="Parameter Recovery: Beta (num 3)"); abline(0,1)
arrows(x0=simul_pars$beta,
       y0=beta_mean - beta_sd,
       y1=beta_mean + beta_sd,
       length=0.05, angle=90, code=3)



### 
print(beta_sd)
cat("mu_alpha_pos:", round(mu_alpha_pos, 3), "(true: 0.20)\n")
cat("mu_alpha_neg:", round(mu_alpha_neg, 3), "(true: 0.30)\n")
cat("mu_beta:     ", round(mu_beta,      3), "(true: 2.00)\n")
