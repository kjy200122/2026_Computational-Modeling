#if (file.exists(".RData")) file.remove(".RData")
#install.packages("hBayesDM")
#Sys.setenv(TAR = "C:/rtools42/usr/bin/tar.exe")
#Sys.setenv(TAR_OPTIONS = "--force-local")
#remotes::install_github("CCS-Lab/hBayesDM", subdir="R", force=TRUE)


library("rstan")
library(hBayesDM)

library(ggplot2) #For graph
library(gridExtra)

####################### problem 1 ########################
setwd("C:/HW5") # Set your direction
ra_data <- read.table("ra_exampleData.txt", header = TRUE)
colnames(ra_data)  # 컬럼명 확인 후 결과 붙여넣어줘
head(ra_data)


# There are few problem, so we have to reconstruct ra_data (name : ra_data_fixed)
ra_data_fixed <- data.frame(
  subjID = ra_data$subjID,
  gain   = ra_data$gain,
  loss   = abs(ra_data$loss),   # 음수 → 양수로 변환 (hBayesDM 요구사항)
  cert   = ra_data$cert,
  gamble = ra_data$gamble
)
write.table(ra_data_fixed, "ra_exampleData_fixed.txt", sep = "\t", row.names = FALSE, quote = FALSE) # Save file in my local folder


output1 <- ra_prospect(data = ra_data_fixed, niter = 4000, nwarmup = 1000, nchain  = 4, ncore  = 4)

output1$allIndPars
output1$fit
plot(output1, type="trace", fontSize=11)
plot(output1, type="trace", inc_warmup=T)
plot(output1)
plotInd(output1, "ep")
print(output1)

##### (1)
plot(output1, type="trace", fontSize=11)
plot(output1)
plotInd(output1, "rho") 
plotInd(output1, "lambda")
plotInd(output1, "tau")


bayes_result <- output1$allIndPars

mle_result <- data.frame(
  subjID = c(2, 3, 4, 6, 7),
  rho    = c(1.0199, 0.6927, 0.7944, 0.8656, 0.9616),
  lambda = c(0.7828, 2.4747, 1.0765, 0.9823, 1.4085),
  tau    = c(1.0071, 4.7421, 1.0355, 3.1299, 2.3468)
) # we already calculated in Hw 2

compare_model <- data.frame(
  subjID       = c(2, 3, 4, 6, 7),
  mle_rho      = mle_result$rho,
  mle_lambda   = mle_result$lambda,
  mle_tau      = mle_result$tau,
  
  bayes_rho    = bayes_result$rho,
  bayes_lambda = bayes_result$lambda,
  bayes_tau    = bayes_result$tau
)
print(compare_model)

# rho
plot(compare_model$mle_rho, compare_model$bayes_rho,
     main = "rho",
     xlab = "MLE", ylab = "Posterior Mean",
     pch = 19, col = "blue",
     xlim = c(0, 2), ylim = c(0, 2))
abline(0, 1, lty = 2, col = "gray")


# lambda
plot(compare_model$mle_lambda, compare_model$bayes_lambda,
     main = "lambda",
     xlab = "MLE", ylab = "Posterior Mean",
     pch = 19, col = "red",
     xlim = c(0, 5), ylim = c(0, 5))
abline(0, 1, lty = 2, col = "gray")

# tau
plot(compare_model$mle_tau, compare_model$bayes_tau,
     main = "tau",
     xlab = "MLE", ylab = "Posterior Mean",
     pch = 19, col = "green",
     xlim = c(0, 6), ylim = c(0, 6))
abline(0, 1, lty = 2, col = "gray")



####################### problem 2 ########################
##### (a)##############
BIC_prospect <- 580.434 # We already calculated BIC in HW 2
BIC_noLA     <- 651.866

BF_12 <- exp(- (BIC_prospect - BIC_noLA) / 2)
print(BF_12) # 3.245351e+15

##### (b)##############
#BF_12 > 100, So we conclude ra_prospect (model 1) is more fit than ra_noLA (model 2) 

##### (c)##############
Prior_odds <-  0.6/0.4 # P(M1)/P(M2) = P(ra_prospect) / P(ra_noLA) = 0.6/0.4 = 1.5
Posterior_odds <- BF_12 * Prior_odds
print(Posterior_odds) #  = BF_12 * 1.5 = about 4.86 * 10^15


#### (d) #############3
#P(M1 / D) = Posterior odds / (1 + Posterior odds) 
pro_M1_D = ((Posterior_odds) / ( 1 + Posterior_odds))
print(pro_M1_D) # = 1


###3 (E) #####
BIC_POW1 = 2.093 # Also, We already calculated BIC in HW 2
BIC_POW2 = 4.173
BIC_EXP1 = 2.393
BIC_EXP2 = 4.266
BIC_EXPOW = 6.252
BIC_HYP1 = 2.152
BIC_HYP2 = 4.205

BIC_models <- c(BIC_POW1, BIC_POW2, BIC_EXP1, BIC_EXP2, BIC_EXPOW, BIC_HYP1, BIC_HYP2)
model_names <- c("POW1", "POW2", "EXP1", "EXP2", "EXPOW", "HYP1", "HYP2")


for (i in 2:7){
  BF_memory <- exp(- (BIC_POW1 - BIC_models[i]) / 2) #BF12 -> 2 : Models, 1 : POW1
  cat(model_names[i], ':', round(BF_memory, 3), '\n')
}



################ problem 3 ############3
### (1)
# Change 'ra_prospect_singleSubj.stan'

parameters <- rstan::extract(output)

# apply(..., 2, mean): 각 피험자(열)에 대해 mean 계산
rho_means    <- apply(parameters$rho,    2, mean)
lambda_means <- apply(parameters$lambda, 2, mean)
tau_means    <- apply(parameters$tau,    2, mean)

# 결과 데이터프레임으로 정리
results <- data.frame(
  subjID = allSubjs,
  rho    = round(rho_means,    4),
  lambda = round(lambda_means, 4),
  tau    = round(tau_means,    4)
)

print(results)