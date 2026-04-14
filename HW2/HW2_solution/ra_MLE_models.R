rm(list=ls())  # clear workspace
graphics.off() # close all figures

set.seed(20250403)  # set a seed number for replication
source('ra_MLE_prototype.R')

# Loading data
data = read.table("ra_exampleData.txt", header = TRUE)

N = 5  # number of subjects
trial = 140 # number of trials per subject
subjList = unique(data$subjID)

k_ra_prospect <- 3 # number of parameters
k_ra_noLA <- 2 # number of parameters
k_ra_noRA <- 2 # number of parameters

# Generate random uniform numbers between 0 and 1 to use as initials for the optim procedure
param_init_prospect <- runif(3)
param_init <- runif(2)

param_low_prospect <- c(0, 0, 0); param_up_prospect <- c(2, 10, 5);  # lower and upper bounds of prospect model (0<rho<2, 0<lambda<10, 0<tau<5)
param_low_noLA <- c(0, 0); param_up_noLA <- c(2, 5);  # lower and upper bounds of noLA model (0<rho<2, 0<tau<5)
param_low_noRA <- c(0, 0); param_up_noRA <- c(10, 5);  # lower and upper bounds of noRA model (0<lambda<10, 0<tau<5)

####### Functions for AIC & BIC #######
# Compute AIC = -2*log(lik) + 2*K
compute_AIC <- function(log_lik, k){
  AIC <- 2*log_lik + 2*k
  return(AIC)
}

# Compute BIC = -2*log(lik) + K*log(N)
compute_BIC <- function(log_lik, k, N){
  BIC <- 2*log_lik + k*log(N)
  return(BIC)
}

##########################
## MLE                  ##
##########################

# dataframe to save parameter estimates 
subjlabelList <- paste('sub', subjList, sep = '')
# ra_prospect model
df_estimates_prospect <- data.frame(matrix(nrow = N, ncol = 3), row.names = subjlabelList)
colnames(df_estimates_prospect) <- c('rho', 'lambda','tau')
# ra_noLA model
df_estimates_noLA <- data.frame(matrix(nrow = N, ncol = 2), row.names = subjlabelList)
colnames(df_estimates_noLA) <- c('rho', 'tau')
# ra_noRA model
df_estimates_noRA <- data.frame(matrix(nrow = N, ncol = 2), row.names = subjlabelList)
colnames(df_estimates_noRA) <- c('lambda', 'tau')

# dataframe to save AIC & BIC
# ra_prospect model
modelcomp_prospect <- data.frame(matrix(nrow = N, ncol = 3), row.names = subjlabelList)
colnames(modelcomp_prospect) <- c('AIC', 'BIC')
# ra_noLA model
modelcomp_noLA <- data.frame(matrix(nrow = N, ncol = 2), row.names = subjlabelList)
colnames(modelcomp_noLA) <- c('AIC', 'BIC')
# ra_noRA model
modelcomp_noRA <- data.frame(matrix(nrow = N, ncol = 2), row.names = subjlabelList)
colnames(modelcomp_noRA) <- c('AIC', 'BIC')

# run MLE for all subjects
for (s in 1:N){ 
  subject = subjList[s] 
  subjlabel = subjlabelList[s] # ex. sub1, sub2,... 
  subdata = data[data$subjID == subject,] # data of each subjects

  # Call general purpose optimization rountine
  mle_model_prospect <- optim(param_init_prospect, mle_ra_prospect, method = "L-BFGS-B", lower = param_low_prospect, upper = param_up_prospect, T = trial, cert = subdata$cert, gain = subdata$gain, loss = subdata$loss, gamble = subdata$gamble)
  mle_model_noLA <- optim(param_init, mle_ra_noLA, method = "L-BFGS-B", lower = param_low_noLA, upper = param_up_noLA, T = trial, cert = subdata$cert, gain = subdata$gain, loss = subdata$loss, gamble = subdata$gamble)
  mle_model_noRA <- optim(param_init, mle_ra_noRA, method = "L-BFGS-B", lower = param_low_noRA, upper = param_up_noRA, T = trial, cert = subdata$cert, gain = subdata$gain, loss = subdata$loss, gamble = subdata$gamble)
  
  # Try many different inits to escape from the local maxima
  for (i in 1:100) {
    # 초기값을 각 모델의 범위에 맞춰 생성
    param_init_prospect <- runif(3, min = param_low_prospect, max = param_up_prospect)
    param_init_noLA     <- runif(2, min = param_low_noLA,     max = param_up_noLA)
    param_init_noRA     <- runif(2, min = param_low_noRA,     max = param_up_noRA)
  
    # 최적화 실행
    temp_ra_prospect <- optim(param_init_prospect, mle_ra_prospect, method = 'L-BFGS-B',
                              lower = param_low_prospect, upper = param_up_prospect,
                              T = trial, cert = subdata$cert, gain = subdata$gain,
                              loss = subdata$loss, gamble = subdata$gamble)
  
    temp_ra_noLA <- optim(param_init_noLA, mle_ra_noLA, method = 'L-BFGS-B',
                          lower = param_low_noLA, upper = param_up_noLA,
                          T = trial, cert = subdata$cert, gain = subdata$gain,
                          loss = subdata$loss, gamble = subdata$gamble)
  
    temp_ra_noRA <- optim(param_init_noRA, mle_ra_noRA, method = 'L-BFGS-B',
                          lower = param_low_noRA, upper = param_up_noRA,
                          T = trial, cert = subdata$cert, gain = subdata$gain,
                          loss = subdata$loss, gamble = subdata$gamble)
  
    # 최적 결과로 갱신
    if (temp_ra_prospect$value < mle_model_prospect$value) mle_model_prospect <- temp_ra_prospect
    if (temp_ra_noLA$value     < mle_model_noLA$value)     mle_model_noLA     <- temp_ra_noLA
    if (temp_ra_noRA$value     < mle_model_noRA$value)     mle_model_noRA     <- temp_ra_noRA
  }

  parm_ra_prospect <- mle_model_prospect$par
  parm_ra_noLA <- mle_model_noLA$par
  parm_ra_noRA <- mle_model_noRA$par
  
  df_estimates_prospect[subjlabel,] <- parm_ra_prospect
  df_estimates_noLA[subjlabel,] <- parm_ra_noLA
  df_estimates_noRA[subjlabel,] <- parm_ra_noRA
  
  modelcomp_prospect[subjlabel,'AIC'] <- compute_AIC(mle_model_prospect$value, k_ra_prospect)
  modelcomp_noLA[subjlabel,'AIC'] <- compute_AIC(mle_model_noLA$value, k_ra_noLA)
  modelcomp_noRA[subjlabel,'AIC'] <- compute_AIC(mle_model_noRA$value, k_ra_noRA)
  
  modelcomp_prospect[subjlabel,'BIC'] <- compute_BIC(mle_model_prospect$value, k_ra_prospect, trial)
  modelcomp_noLA[subjlabel,'BIC'] <- compute_BIC(mle_model_noLA$value, k_ra_noLA, trial)
  modelcomp_noRA[subjlabel,'BIC'] <- compute_BIC(mle_model_noRA$value, k_ra_noRA, trial)
}

modelcomp_summary <- data.frame(row.names = c('ra prospect','ra noLA', 'ra noRA'),
                                'AIC' = c(mean(modelcomp_prospect$AIC), mean(modelcomp_noLA$AIC), mean(modelcomp_noRA$AIC)),
                                'sd (AIC)' = c(sd(modelcomp_prospect$AIC), sd(modelcomp_noLA$AIC),  sd(modelcomp_noRA$AIC)),
                                'BIC' = c(mean(modelcomp_prospect$BIC), mean(modelcomp_noLA$BIC), mean(modelcomp_noRA$BIC)),
                                'sd (BIC)' = c(sd(modelcomp_prospect$BIC), sd(modelcomp_noLA$BIC),  sd(modelcomp_noRA$BIC)))


# save for HW2 Q2
df_prospect <- df_estimates_prospect
df_noLA <- df_estimates_noLA
df_noRA <- df_estimates_noRA
save(df_prospect, file = "HW2_ra_prospect.RData")
save(df_noLA, file = "HW2_ra_noLA.RData")
save(df_noRA, file = "HW2_ra_noRA.RData")
