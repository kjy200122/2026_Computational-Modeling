## relying on HW1 codes

N = 5  # number of subjects
T = 140 # number of trials per subject  

# read_txt
data <- read.table("ra_exampleData.txt", header=TRUE, sep="\t")
data
subj_2 <- subset(data, subjID == 2)
subj_3 <- subset(data, subjID == 3)
subj_4 <- subset(data, subjID == 4)
subj_6 <- subset(data, subjID == 6)
subj_7 <- subset(data, subjID == 7)

#### ra_noLA ####
ra_noLA <- function(param, gain, loss, cert, gamble){
  rho    <- param[1]
  lambda <- 1           
  tau    <- param[2]    
  
  sum_minusLL = 0
  for (t in 1:T) {
    evSafe   = cert[t]^rho
    evGamble = 0.5*(gain[t]^rho - lambda*abs(loss[t])^rho) 
    pGamble  = 1 / (1 + exp(tau*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001
    tmp_minusLL = -log(pGamble)*gamble[t] - log(1-pGamble)*(1-gamble[t])
    sum_minusLL = sum_minusLL + tmp_minusLL
  }
  return(sum_minusLL)
}

# parameter bound 
param_low_noLA <- c(0, 0)
param_up_noLA  <- c(2, 5)

# subject 2
param_init <- c(runif(1, 0.1, 2), runif(1, 0.1, 5))
result_noLA_subj2 <- optim(param_init, ra_noLA, method="L-BFGS-B",
                           lower=param_low_noLA, upper=param_up_noLA,
                           gain=subj_2$gain, loss=subj_2$loss,
                           cert=subj_2$cert, gamble=subj_2$gamble)
for (i in 1:100) {
  param_init <- c(runif(1, 0.1, 2), runif(1, 0.1, 5))
  temp <- optim(param_init, ra_noLA, method="L-BFGS-B",
                lower=param_low_noLA, upper=param_up_noLA,
                gain=subj_2$gain, loss=subj_2$loss,
                cert=subj_2$cert, gamble=subj_2$gamble)
  if (temp$value < result_noLA_subj2$value) result_noLA_subj2 <- temp
}

# subject 3
param_init <- c(runif(1, 0.1, 2), runif(1, 0.1, 5))
result_noLA_subj3 <- optim(param_init, ra_noLA, method="L-BFGS-B",
                           lower=param_low_noLA, upper=param_up_noLA,
                           gain=subj_3$gain, loss=subj_3$loss,
                           cert=subj_3$cert, gamble=subj_3$gamble)
for (i in 1:100) {
  param_init <- c(runif(1, 0.1, 2), runif(1, 0.1, 5))
  temp <- optim(param_init, ra_noLA, method="L-BFGS-B",
                lower=param_low_noLA, upper=param_up_noLA,
                gain=subj_3$gain, loss=subj_3$loss,
                cert=subj_3$cert, gamble=subj_3$gamble)
  if (temp$value < result_noLA_subj3$value) result_noLA_subj3 <- temp
}

# subject 4
param_init <- c(runif(1, 0.1, 2), runif(1, 0.1, 5))
result_noLA_subj4 <- optim(param_init, ra_noLA, method="L-BFGS-B",
                           lower=param_low_noLA, upper=param_up_noLA,
                           gain=subj_4$gain, loss=subj_4$loss,
                           cert=subj_4$cert, gamble=subj_4$gamble)
for (i in 1:100) {
  param_init <- c(runif(1, 0.1, 2), runif(1, 0.1, 5))
  temp <- optim(param_init, ra_noLA, method="L-BFGS-B",
                lower=param_low_noLA, upper=param_up_noLA,
                gain=subj_4$gain, loss=subj_4$loss,
                cert=subj_4$cert, gamble=subj_4$gamble)
  if (temp$value < result_noLA_subj4$value) result_noLA_subj4 <- temp
}

# subject 6
param_init <- c(runif(1, 0.1, 2), runif(1, 0.1, 5))
result_noLA_subj6 <- optim(param_init, ra_noLA, method="L-BFGS-B",
                           lower=param_low_noLA, upper=param_up_noLA,
                           gain=subj_6$gain, loss=subj_6$loss,
                           cert=subj_6$cert, gamble=subj_6$gamble)
for (i in 1:100) {
  param_init <- c(runif(1, 0.1, 2), runif(1, 0.1, 5))
  temp <- optim(param_init, ra_noLA, method="L-BFGS-B",
                lower=param_low_noLA, upper=param_up_noLA,
                gain=subj_6$gain, loss=subj_6$loss,
                cert=subj_6$cert, gamble=subj_6$gamble)
  if (temp$value < result_noLA_subj6$value) result_noLA_subj6 <- temp
}

# subject 7
param_init <- c(runif(1, 0.1, 2), runif(1, 0.1, 5))
result_noLA_subj7 <- optim(param_init, ra_noLA, method="L-BFGS-B",
                           lower=param_low_noLA, upper=param_up_noLA,
                           gain=subj_7$gain, loss=subj_7$loss,
                           cert=subj_7$cert, gamble=subj_7$gamble)
for (i in 1:100) {
  param_init <- c(runif(1, 0.1, 2), runif(1, 0.1, 5))
  temp <- optim(param_init, ra_noLA, method="L-BFGS-B",
                lower=param_low_noLA, upper=param_up_noLA,
                gain=subj_7$gain, loss=subj_7$loss,
                cert=subj_7$cert, gamble=subj_7$gamble)
  if (temp$value < result_noLA_subj7$value) result_noLA_subj7 <- temp
}

# save result: ra_noLA
all_params_noLA <- matrix(NA, nrow=5, ncol=2)
colnames(all_params_noLA) <- c("rho", "tau")
rownames(all_params_noLA) <- c("subj_2", "subj_3", "subj_4", "subj_6", "subj_7")

all_params_noLA[1, ] <- result_noLA_subj2$par
all_params_noLA[2, ] <- result_noLA_subj3$par
all_params_noLA[3, ] <- result_noLA_subj4$par
all_params_noLA[4, ] <- result_noLA_subj6$par
all_params_noLA[5, ] <- result_noLA_subj7$par

print("--- ra_noLA ---")
print(all_params_noLA)

# AIC, BIC
N_trials <- 140 
k_noLA <- 2
AIC_noLA <- c(
  2*result_noLA_subj2$value + 2*k_noLA,
  2*result_noLA_subj3$value + 2*k_noLA,
  2*result_noLA_subj4$value + 2*k_noLA,
  2*result_noLA_subj6$value + 2*k_noLA,
  2*result_noLA_subj7$value + 2*k_noLA
)
BIC_noLA <- c(
  2*result_noLA_subj2$value + k_noLA*log(N_trials),
  2*result_noLA_subj3$value + k_noLA*log(N_trials),
  2*result_noLA_subj4$value + k_noLA*log(N_trials),
  2*result_noLA_subj6$value + k_noLA*log(N_trials),
  2*result_noLA_subj7$value + k_noLA*log(N_trials)
)

# ra_noLA: per-subject table
noLA_table <- data.frame(
  Subject = c("subj_2", "subj_3", "subj_4", "subj_6", "subj_7"),
  AIC = round(AIC_noLA, 3),
  BIC = round(BIC_noLA, 3)
)
print("--- ra_noLA (per subject) ---")
print(noLA_table)

# ra_noLA: all-subject
sum_AIC <- c(sum(AIC_noLA))
sum_BIC <- c(sum(BIC_noLA))

model_comparison <- data.frame(
  Model = ("ra_noLA"),
  sum_AIC = round(sum_AIC, 3),
  sum_BIC = round(sum_BIC, 3)
)

print("--- Model Comparison (summed across subjects) ---")
print(model_comparison)
