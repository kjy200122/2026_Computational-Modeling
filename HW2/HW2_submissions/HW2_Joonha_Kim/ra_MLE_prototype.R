rm(list=ls())  # clear workspace
graphics.off() # close all figures

set.seed(202604)  # set a seed number for replication

N = 5  # number of subjects
T = 140 # number of trials per subject

ra_prospect <- function(param, n_trials, gain, loss, cert, gamble){
  # param[1]: rho (risk aversion) 
  # param[2]: lambda (loss aversion)
  # param[3]: tau (inverse temperatue)
  
  sum_minusLL = 0
  rho = param[1]
  lambda = param[2]
  tau = param[3]
  for (t in 1:n_trials) {
    evSafe   = cert[t]^rho
    evGamble = 0.5*(gain[t]^rho - lambda*abs(loss[t])^rho) 
    pGamble  = 1 / (1 + exp(tau*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    tmp_minusLL = -log(pGamble)*gamble[t] - log(1-pGamble)*(1-gamble[t])  # -LL of trial t
    sum_minusLL = sum_minusLL + tmp_minusLL
  }
  sum_minusLL
}

mle_iter <- function(model, n_params, param_low, param_high, n_iter, n_trials, df_sub){
  
  mle_model = list(value=Inf)
  
  for (iter in 1:n_iter) {
    
    param_init <- runif(n_params, param_low, param_high)
    
    temp_model <- optim(param_init, 
                       model, method="L-BFGS-B", 
                       lower = param_low, upper = param_high, n_trials=n_trials,
                       gain=df_sub$gain, loss=df_sub$loss, cert = df_sub$cert, gamble = df_sub$gamble)
    
    if (temp_model$value < mle_model$value){
      mle_model <- temp_model
    }
  }
  mle_model
}

df <- read.delim("ra_exampleData.txt")

### subject 2 ###
df_sub = df[df$subjID == 2, ]
param_low <- c(0, 0, 0); param_high <- c(2, 10, 5)
mle_ra_prospect <- mle_iter(model=ra_prospect, n_params=3, param_low=param_low, param_high=param_high, 
                              n_iter=1000, n_trials=T, df_sub=df_sub)
print(mle_ra_prospect$par)

pars_sub_lst = c()
LL_sub_lst = c()

### all subjects ###
for (i in unique(df$subjID)) {
  df_sub = df[df$subjID == i, ]
  
  param_low <- c(0, 0, 0); param_high <- c(2, 10, 5)
  
  mle_ra_prospect = mle_iter(model=ra_prospect, n_params=3, param_low=param_low, param_high=param_high, 
                               n_iter=1000, n_trials=T, df_sub=df_sub)
  pars_sub_lst[[length(pars_sub_lst)+1]] <- mle_ra_prospect$par
  LL_sub_lst[[length(LL_sub_lst)+1]] <- -mle_ra_prospect$value
  }

ra_prospect_df <- data.frame(
  subjID = unique(df$subjID),
  rho    = sapply(pars_sub_lst, function(x) round(x[1],3)),
  lambda = sapply(pars_sub_lst, function(x) round(x[2],3)),
  tau    = sapply(pars_sub_lst, function(x) round(x[3],3)),
  LL = round(unlist(LL_sub_lst), 3)
)

ra_prospect_df$AIC <- round(-2*ra_prospect_df$LL + 2*3, 3)
ra_prospect_df$BIC <- round(-2*ra_prospect_df$LL + log(T)*3, 3)
print(ra_prospect_df)

### without lambda (ra_noLA) ###
ra_noLA <- function(param, n_trials, gain, loss, cert, gamble){
  # param[1]: rho (risk aversion) 
  # param[2]: tau (inverse temperatue)
  rho = param[1]
  tau = param[2]
    
  sum_minusLL = 0
  for (t in 1:n_trials) {
    evSafe   = cert[t]^rho
    evGamble = 0.5*(gain[t]^rho - abs(loss[t])^rho)
    pGamble  = 1 / (1 + exp(tau*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    tmp_minusLL = -log(pGamble)*gamble[t] - log(1-pGamble)*(1-gamble[t])  # -LL of trial t
    sum_minusLL = sum_minusLL + tmp_minusLL 
  }
  sum_minusLL
}


pars_sub_lst = c()
LL_sub_lst = c()

for (i in unique(df$subjID)) {
  df_sub = df[df$subjID == i, ]
  
  param_low <- c(0, 0); param_high <- c(2, 5)
  
  mle_ra_noLA = mle_iter(model=ra_noLA, n_params=2, param_low=param_low, param_high=param_high, 
                             n_iter=1000, n_trials=T, df_sub=df_sub)
  pars_sub_lst[[length(pars_sub_lst)+1]] <- mle_ra_noLA$par
  LL_sub_lst[[length(LL_sub_lst)+1]] <- -mle_ra_noLA$value
}

ra_noLA_df <- data.frame(
  subjID = unique(df$subjID),
  rho    = sapply(pars_sub_lst, function(x) round(x[1],3)),
  tau    = sapply(pars_sub_lst, function(x) round(x[2],3)),
  LL = round(unlist(LL_sub_lst), 3)
)

ra_noLA_df$AIC <- round(-2*ra_noLA_df$LL + 2*2, 3)
ra_noLA_df$BIC <- round(-2*ra_noLA_df$LL + log(T)*2, 3)
print(ra_noLA_df)


### without rho (ra_noRA) ###

ra_noRA <- function(param, n_trials, gain, loss, cert, gamble){
  # param[1]: lambda (loss aversion) 
  # param[2]: tau (inverse temperatue)
  lambda = param[1]
  tau = param[2]
  
  sum_minusLL = 0
  for (t in 1:n_trials) {
    evSafe   = cert[t]
    evGamble = 0.5*(gain[t] - lambda*abs(loss[t]))
    pGamble  = 1 / (1 + exp(tau*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    tmp_minusLL = -log(pGamble)*gamble[t] - log(1-pGamble)*(1-gamble[t])  # -LL of trial t
    sum_minusLL = sum_minusLL + tmp_minusLL 
  }
  sum_minusLL
}

pars_sub_lst = c()
LL_sub_lst = c()

for (i in unique(df$subjID)) {
  df_sub = df[df$subjID == i, ]
  
  param_low <- c(0, 0); param_high <- c(2, 5)
  
  mle_ra_noRA = mle_iter(model=ra_noRA, n_params=2, param_low=param_low, param_high=param_high, 
                         n_iter=1000, n_trials=T, df_sub=df_sub)
  pars_sub_lst[[length(pars_sub_lst)+1]] <- mle_ra_noRA$par
  LL_sub_lst[[length(LL_sub_lst)+1]] <- -mle_ra_noRA$value
}

ra_noRA_df <- data.frame(
  subjID = unique(df$subjID),
  lambda    = sapply(pars_sub_lst,function(x) round(x[1],3)),
  tau    = sapply(pars_sub_lst, function(x) round(x[2],3)),
  LL = round(unlist(LL_sub_lst), 3)
)

ra_noRA_df$AIC <- round(-2*ra_noRA_df$LL + 2*2, 3)
ra_noRA_df$BIC <- round(-2*ra_noRA_df$LL + log(T)*2, 3)

print(ra_prospect_df)
print(ra_noLA_df)
print(ra_noRA_df)

ra_prospect_AIC = sum(ra_prospect_df$AIC)
ra_prospect_BIC = sum(ra_prospect_df$BIC)
ra_noLA_AIC = sum(ra_noLA_df$AIC)
ra_noLA_BIC = sum(ra_noLA_df$BIC)
ra_noRA_AIC = sum(ra_noRA_df$AIC)
ra_noRA_BIC = sum(ra_noRA_df$BIC)

model_comp <- data.frame(
  model = c("ra_prospect", "ra_noLA", "ra_noRA"),
  AIC   = c(ra_prospect_AIC, ra_noLA_AIC, ra_noRA_AIC),
  BIC   = c(ra_prospect_BIC, ra_noLA_BIC, ra_noRA_BIC)
)

print(model_comp)
write.csv(ra_prospect_df, "./results/ra_prospect.csv")
write.csv(ra_noLA_df, "./results/ra_noLA.csv")
write.csv(ra_noRA_df, "./results/ra_noRA.csv")
write.csv(model_comp, "./results/ra_model_comp.csv")