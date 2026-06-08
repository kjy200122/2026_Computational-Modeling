library(rstan)
library(hBayesDM)

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

data <- read.table("ra_exampleData.txt", header = TRUE, sep = "\t")

# data of 5 subjects
data_subs <- subset(data, subjID %in% c(2,3,4,6,7))

# hBayesDM fit
fit <- ra_prospect(
  data = data_subs,
  niter = 4000,
  nwarmup = 1000,
  nchain = 4,
  ncore = 4 
)

print(fit$fit)

# posterior samples
post <- rstan::extract(fit$fit)

# posterior means
posterior_summary <- data.frame(
  subjID = c(2,3,4,6,7),
  rho_post = apply(post$rho, 2, mean),
  lambda_post = apply(post$lambda, 2, mean),
  tau_post = apply(post$tau, 2, mean)
)

mle_summary <- read.csv("HW2_MLE_results.csv")
comparison <- merge(mle_summary, posterior_summary, by = "subjID")

# posterior distributions
par(mfrow = c(5,3))

for(i in 1:5){
  hist(post$rho[,i],
       main = paste("Subject", comparison$subjID[i], "rho"),
       xlab = "rho")
  
  hist(post$lambda[,i],
       main = paste("Subject", comparison$subjID[i], "lambda"),
       xlab = "lambda")
  
  hist(post$tau[,i],
       main = paste("Subject", comparison$subjID[i], "tau"),
       xlab = "tau")
}

# plot
par(mfrow = c(1,3))

rho_lim <- range(c(comparison$rho,
                   comparison$rho_post))
plot(comparison$rho,
     comparison$rho_post,
     xlab = "MLE rho",
     ylab = "Posterior mean rho",
     main = "rho",
     xlim = rho_lim,
     ylim = rho_lim)
abline(0, 1, lty=2)


lambda_lim <- range(c(comparison$lambda,
                      comparison$lambda_post))
plot(comparison$lambda,
     comparison$lambda_post,
     xlab = "MLE lambda",
     ylab = "Posterior mean lambda",
     main = "lambda",
     xlim = lambda_lim,
     ylim = lambda_lim)
abline(0, 1, lty=2)

tau_lim <- range(c(comparison$tau,
                   comparison$tau_post))
plot(comparison$tau,
     comparison$tau_post,
     xlab = "MLE tau",
     ylab = "Posterior mean tau",
     main = "tau",
     xlim = tau_lim,
     ylim = tau_lim)
abline(0, 1, lty=2)