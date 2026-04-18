rm(list=ls())  # clear workspace
graphics.off() # close all figures

N = 5  # number of subjects
T = 140 # number of trials per subject

#Open file and save the data
dat <- read.delim("ra_exampleData.txt", header = TRUE)
sub_dat <- list() #Save the data according to people number
for (i in dat$subjID){
  sub_dat[[as.character(i)]] <- subset(dat, subjID == i)
}

#parameter initialization
param3_init <- runif(3)
param_gam_low <- c(0,0,0);param_gam_up <- c(2,10,5);   # free parameters: rho, tau, lambda




#Mle equation for gamble

mle_gam <- function(param, gain, loss, cert, gamble){
  
  rho <- param[1]
  lambda <- param[2]
  tau <- param[3]
  
  sum_minusLL = 0  # sum of minus log likelihood. Initialize
  
  for (t in 1:T) {
  
    # evSafe: expected value of a certain (safe) option
    # evGamble: expected value of a risky option (gamble)
    # pGamble   # probability of choosing a gamble on each trial
    # free parameters: rho, tau, lambda
  evSafe   = cert[t]^rho #확실한 보상
  evGamble = 0.5*(gain[t]^rho - lambda*abs(loss[t])^rho) #도박 보상
  pGamble  = 1 / (1 + exp(tau*(evSafe - evGamble))) #도박을 선택할 확률. evGamble이 evSafe보다 클수록 1에 가까워짐
  pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
  
  # 확률 무한대 발산 방지
  if (pGamble <= 0) pGamble <- 10e-6
  if (pGamble >= 1) pGamble <- 1 - 10e-6
  
  tmp_minusLL = -log(pGamble)*gamble[t] - log(1-pGamble)*(1-gamble[t])  # gamble을 했으면 -log(pGamble), 안 했으면 그 반대 확률
  sum_minusLL <- sum_minusLL + tmp_minusLL
  }
  sum_minusLL
}


#Data for one subjects(subject 2)

#MLE calculation
gain_2 <- sub_dat[["2"]]$gain# gain[t] & loss[t]: gain and loss outcome of a risky option on trial t
loss_2 <- sub_dat[["2"]]$loss
cert_2 <- sub_dat[["2"]]$cert # cert[t]: outcome of a certain option on trial t
gamble_2 <- sub_dat[["2"]]$gamble# gamble[t]: choice on trial t. gamble[t]=1 --> chose a risky option. 0 --> chose a certain option.

param3_init <- runif(3);

mle_model_gam <- optim(param3_init, mle_gam, method="L-BFGS-B", lower=param_gam_low, upper=param_gam_up, gain = gain_2, loss = loss_2,cert = cert_2, gamble = gamble_2)

for (i in 1:100) {
  param3_init <- runif(3)

  temp_gam <- optim(param3_init, mle_gam, method="L-BFGS-B", lower=param_gam_low, upper=param_gam_up, gain = gain_2, loss = loss_2,cert = cert_2, gamble = gamble_2)
  if(temp_gam$value < mle_model_gam$value) mle_model_gam <- temp_gam
}

parm_gam <- mle_model_gam$par #parameters

p_prd_gam <- rep(NA, T) #predictability

for (t in 1:T) {
  evSafe   <- cert_2[t]^parm_gam[1]
  evGamble <- 0.5 * (gain_2[t]^parm_gam[1] - parm_gam[2] * abs(loss_2[t])^parm_gam[1])
  p_prd_gam[t] <- 1 / (1 + exp(parm_gam[3] * (evSafe - evGamble)))
  p_prd_gam[t] <- p_prd_gam[t] * 0.9998 + 0.0001
}


r2_gam <- 1 - sum((gamble_2 - p_prd_gam)^2) / sum((gamble_2 - mean(gamble_2))^2)

#question 1 answer
result_sub2 <- data.frame(
  Subject = 2,
  rho     = mle_model_gam$par[1],
  lambda  = mle_model_gam$par[2],
  tau     = mle_model_gam$par[3],
  minusLL = mle_model_gam$value,
  R2 = r2_gam
)
print(result_sub2)



# for multiple subjects

subj_names <- names(sub_dat) #list which saves the name

#MLE calculation

for (number in 1:N){
  
  n <- subj_names[number] #each subject name
  
  gain_n <- sub_dat[[n]]$gain# gain[t] & loss[t]: gain and loss outcome of a risky option on trial t
  loss_n <- sub_dat[[n]]$loss
  cert_n <- sub_dat[[n]]$cert# cert[t]: outcome of a certain option on trial t
  gamble_n <- sub_dat[[n]]$gamble# gamble[t]: choice on trial t. gamble[t]=1 --> chose a risky option. 0 --> chose a certain option.

param3_init <- runif(3);

mle_model_gam <- optim(param3_init, mle_gam, method="L-BFGS-B", lower=param_gam_low, upper=param_gam_up, gain = gain_n, loss = loss_n,cert = cert_n, gamble = gamble_n)

for (i in 1:100) {
  param3_init <- runif(3)
  
  temp_gam <- optim(param3_init, mle_gam, method="L-BFGS-B", lower=param_gam_low, upper=param_gam_up, gain = gain_n, loss = loss_n,cert = cert_n, gamble = gamble_n)
  if(temp_gam$value < mle_model_gam$value) mle_model_gam <- temp_gam
}

parm_gam <- mle_model_gam$par #parameter

p_prd_gam <- rep(NA, T) #predictability

for (t in 1:T) {
  evSafe   <- cert_n[t]^parm_gam[1]
  evGamble <- 0.5 * (gain_n[t]^parm_gam[1] - parm_gam[2] * abs(loss_n[t])^parm_gam[1])
  p_prd_gam[t] <- 1 / (1 + exp(parm_gam[3] * (evSafe - evGamble)))
  p_prd_gam[t] <- p_prd_gam[t] * 0.9998 + 0.0001
}


r2_gam <- 1 - sum((gamble_n - p_prd_gam)^2) / sum((gamble_n - mean(gamble_n))^2)

#Part 1. If you want know the result for Each 
"result_all <- rbind(result_all, data.frame(
  Subject = n,
  rho     = mle_model_gam$par[1],
  lambda  = mle_model_gam$par[2],
  tau     = mle_model_gam$par[3],
  minusLL = mle_model_gam$value,
  R2      = r2_gam
))"

#Part 2. If you want to know the all
result_all <- data.frame(
  rho     = mle_model_gam$par[1],
  lambda  = mle_model_gam$par[2],
  tau     = mle_model_gam$par[3],
  minusLL = mle_model_gam$value,
  R2 = r2_gam
)

}
print(result_all)

