library(rstan)
source("simulate_hw6_model2.R")

dat = read.table("simul_data_hw6_model2.txt", header=T, sep="\t")

allSubjs = unique(dat$subjID)
N = length(allSubjs)
T = table(dat$subjID)[1]

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

# Run Stan
output_m2 = stan("hw6_model2.stan", data = dataList,
                 pars = c("alpha_pos", "alpha_neg", "beta"),
                 iter = 2000, warmup=1000, chains=2, cores=2)

print(output_m2)

# Extract posterior
parameters_m2 <- rstan::extract(output_m2)

alpha_pos_mean = apply(parameters_m2$alpha_pos, 2, mean)
alpha_pos_sd   = apply(parameters_m2$alpha_pos, 2, sd)
alpha_neg_mean = apply(parameters_m2$alpha_neg, 2, mean)
alpha_neg_sd   = apply(parameters_m2$alpha_neg, 2, sd)
beta_mean_m2   = apply(parameters_m2$beta, 2, mean)
beta_sd_m2     = apply(parameters_m2$beta, 2, sd)

#Scatter plot
par(mfrow=c(1,3))

plot(simul_pars$alpha_pos, alpha_pos_mean,
     xlim=c(0, 0.5), ylim=c(0, 0.5),
     xlab="True alpha_pos", ylab="Estimated alpha_pos",
     main="Q4: alpha_pos recovery")
abline(0, 1, col="red")
arrows(x0=simul_pars$alpha_pos, y0=alpha_pos_mean - alpha_pos_sd,
       y1=alpha_pos_mean + alpha_pos_sd,
       length=0.02, angle=90, code=3)

plot(simul_pars$alpha_neg, alpha_neg_mean,
     xlim=c(0, 0.6), ylim=c(0, 0.6),
     xlab="True alpha_neg", ylab="Estimated alpha_neg",
     main="Q4: alpha_neg recovery")
abline(0, 1, col="red")
arrows(x0=simul_pars$alpha_neg, y0=alpha_neg_mean - alpha_neg_sd,
       y1=alpha_neg_mean + alpha_neg_sd,
       length=0.02, angle=90, code=3)

plot(simul_pars$beta, beta_mean_m2,
     xlim=c(0, 4), ylim=c(0, 4),
     xlab="True beta", ylab="Estimated beta",
     main="Q4: beta recovery")
abline(0, 1, col="red")
arrows(x0=simul_pars$beta, y0=beta_mean_m2 - beta_sd_m2,
       y1=beta_mean_m2 + beta_sd_m2,
       length=0.05, angle=90, code=3)

par(mfrow=c(1,1))
