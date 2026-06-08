rm(list=ls())
library(hBayesDM)

set.seed(202604)  # set a seed number for replication

### hBayes fitting ###
df <- read.delim("ra_exampleData.txt")

hba_output = ra_prospect(data=df,
                     niter = 4000,
                     nwarmup = 1000,
                     nchain = 4,
                     ncore = 4)

### check hBayes results ###
print(hba_output$allIndPars)
print(hba_output$fit)

plot(hba_output, type="trace", fontSize=11)

# convergence check for group level params
rstan::traceplot(hba_output$fit, pars = c("mu_rho", "mu_lambda", "mu_tau",
                                          "sigma[1]", "sigma[2]", "sigma[3]"))
# convergence check for individual params
rstan::traceplot(hba_output$fit, pars = c("rho", "lambda", "tau"))

plotInd(hba_output, "rho", show_density = TRUE)
plotInd(hba_output, "lambda", show_density = TRUE)
plotInd(hba_output, "tau", show_density = TRUE)

### MLE fitting (from HW2) ###

set.seed(202604)  # set a seed number for replication

N = 5  # number of subjects
T = 140 # number of trials per subject

ra_prospect_mle <- function(param, n_trials, data){
  # param[1]: rho (risk aversion) 
  # param[2]: lambda (loss aversion)
  # param[3]: tau (inverse temperatue)
  
  sum_minusLL = 0
  rho = param[1]
  lambda = param[2]
  tau = param[3]
  
  gain = data$gain
  loss = data$loss
  cert = data$cert
  gamble = data$gamble
  
  for (t in 1:n_trials) {
    evSafe   = cert[t]^rho
    evGamble = 0.5*(gain[t]^rho - lambda*abs(loss[t])^rho) 
    pGamble  = 1 / (1 + exp(tau*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    tmp_minusLL = -log(pGamble)*gamble[t] - log(1-pGamble)*(1-gamble[t])  # -LL of trial t
    sum_minusLL = sum_minusLL + tmp_minusLL
  }
  return(sum_minusLL)
}

mle_iter <- function(model, n_params, param_low, param_high, n_iter, n_trials, data){
  
  mle_model = list(value=Inf)
  
  for (iter in 1:n_iter) {
    
    param_init <- runif(n_params, param_low, param_high)
    
    temp_model <- optim(param_init, 
                        model, method="L-BFGS-B", 
                        lower = param_low, upper = param_high, n_trials=n_trials, data=data)
    
    if (temp_model$value < mle_model$value){
      mle_model <- temp_model
    }
  }
  return(mle_model)
}

pars_sub_lst = c()
LL_sub_lst = c()

for (i in unique(df$subjID)) {
  df_sub = df[df$subjID == i, ]
  
  param_low <- c(0, 0, 0); param_high <- c(2, 10, 5)
  
  mle_ra_prospect = mle_iter(model=ra_prospect_mle, n_params=3, param_low=param_low, param_high=param_high, 
                             n_iter=100, n_trials=T, data=df_sub)
  pars_sub_lst[[length(pars_sub_lst)+1]] <- mle_ra_prospect$par
  LL_sub_lst[[length(LL_sub_lst)+1]] <- -mle_ra_prospect$value
}

mle_result_df <- data.frame(
  subjID = unique(df$subjID),
  rho    = sapply(pars_sub_lst, function(x) round(x[1],3)),
  lambda = sapply(pars_sub_lst, function(x) round(x[2],3)),
  tau    = sapply(pars_sub_lst, function(x) round(x[3],3)),
  LL = round(unlist(LL_sub_lst), 3)
)

### compare HBA and MLE estimates ###

# rho
plot(x=mle_result_df$rho, y=hba_output$allIndPars$rho,
             xlim=c(0.65, 1.05), ylim=c(0.65, 1.05), main = "rho",
             xlab = "MLE estimates", ylab="Posterior mean(HBA)")
abline(a=0, b=1, col="red",lty=2)

# lambda
plot(x=mle_result_df$lambda, y=hba_output$allIndPars$lambda, 
     xlim=c(0.5, 2.6), ylim=c(0.5, 2.6), 
     main = "lambda",
     xlab = "MLE estimates", ylab="Posterior mean(HBA)")
abline(a=0, b=1, col="red",lty=2)

# tau
plot(x=mle_result_df$tau, y=hba_output$allIndPars$tau, 
     xlim=c(1, 5), ylim=c(1, 5),
     main = "tau",
     xlab = "MLE estimates", ylab="Posterior mean(HBA)")
abline(a=0, b=1, col="red",lty=2)

