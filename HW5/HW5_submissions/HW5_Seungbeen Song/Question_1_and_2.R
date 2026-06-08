Sys.setenv(TMPDIR = "C:/Rtmp")
Sys.setenv(TEMP = "C:/Rtmp")
Sys.setenv(TMP = "C:/Rtmp")

install.packages("remotes")

remotes::install_github(
  "CCS-Lab/hBayesDM",
  subdir = "R",
  dependencies = TRUE,
  upgrade = "never",
  build_vignettes = FALSE
)

remotes::install_github(
  "CCS-Lab/hBayesDM",
  subdir = "R",
  dependencies = TRUE,
  upgrade = "never",
  build = FALSE
)
library(hBayesDM)

list.files(system.file("stan_files", package = "hBayesDM"))

pkg_path <- system.file(package = "hBayesDM")
pkg_path


# `devtools` is required to install hBayesDM from GitHub
if (!require(devtools)) install.packages("devtools")

devtools::install_github("CCS-Lab/hBayesDM", subdir="R")

# ======================================================================= #

########### Question 1 ########### 

data <- read.table("ra_exampleData.txt", header = TRUE)

### 1. hBayesDM ###
output = ra_prospect(data, niter=4000, nwarmup=1000, nchain=4, ncore=4)

# results
output$allIndPars

# visualizations
plotInd(output, "rho")
plotInd(output, "lambda")
plotInd(output, "tau")
plot(output, type="trace", fontSize=11) 

### 2. MLE (directly from HW2: Quesion 2) ###
N = 5  # number of subjects
T = 140 # number of trials per subject

param_init <-runif(3) # initial values
param_low <- c(0,0,0); param_up <- c(2,10,5) # parameter bounds

data <- read.table("ra_exampleData.txt", header = TRUE)

ra_prospect <- function (param, cert, gain, loss, gamble, T){
  sum_minusLL = 0  # sum of minus log likelihood. Initialize.
  
  for (t in 1:T) {
    # cert[t]: outcome of a certain option on trial t
    # gain[t] & loss[t]: gain and loss outcome of a risky option on trial t
    # gamble[t]: choice on trial t. gamble[t]=1 --> chose a risky option. 0 --> chose a certain option.
    # 
    # evSafe: expected value of a certain (safe) option
    # evGamble: expected value of a risky option (gamble)
    # pGamble   # probability of choosing a gamble on each trial
    # free parameters: rho, tau, lambda
    
    evSafe   = cert[t]^param[1]
    evGamble = 0.5*(gain[t]^param[1] - param[2]*abs(loss[t])^param[1]) 
    pGamble  = 1 / (1 + exp(param[3]*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    tmp_minusLL = -log(pGamble)*gamble[t] - log(1-pGamble)*(1-gamble[t])  # -LL of trial t
    sum_minusLL = sum_minusLL + tmp_minusLL
  }
  
  sum_minusLL
}

# Run MLE for all subjects 

all_parameters <- data.frame(
  subjID = unique(data$subjID),
  rho = numeric(N),
  lambda = numeric(N),
  tau = numeric(N)
)

idx = 1

for (i in unique(data$subjID)){
  curr_data <- data[data$subjID == i, ]
  param_init <- runif(3)
  
  # Run MLE
  rp_mle <- optim(param_init, ra_prospect, method="L-BFGS-B", lower=param_low, upper=param_up, cert=curr_data$cert, gain=curr_data$gain, loss=curr_data$loss, gamble=curr_data$gamble, T=140)
  
  # Try many different initial values
  for (j in 1:100){
    param_init <- runif(3)
    
    rp_temp <- optim(param_init, ra_prospect, method="L-BFGS-B", lower=param_low, upper=param_up, cert=curr_data$cert, gain=curr_data$gain, loss=curr_data$loss, gamble=curr_data$gamble, T=140)
    
    if (rp_temp$value < rp_mle$value) rp_mle <- rp_temp
  }
  
  all_parameters$rho[idx] <- rp_mle$par[1]
  all_parameters$lambda[idx] <- rp_mle$par[2]
  all_parameters$tau[idx] <- rp_mle$par[3]
  
  idx = idx + 1
  
}

print(all_parameters)

### Compare posterior means from hBayesDM and MLE estimates
mle <- all_parameters
post <- output$allIndPars

names(mle) <- c("subjID", "rho_mle", "lambda_mle", "tau_mle")
names(post) <- c("subjID", "rho_post", "lambda_post", "tau_post")

compare <- merge(mle, post, by = "subjID")

# plot for rho
plot(
  compare$rho_mle,
  compare$rho_post,
  xlab = "MLE estimate",
  ylab = "Posterior mean",
  main = "MLE vs Posterior Mean: rho"
)

abline(0, 1, lty = 2)

text(
  compare$rho_mle,
  compare$rho_post,
  labels = compare$subjID,
  pos = 3
)

# plot for lambda
plot(
  compare$lambda_mle,
  compare$lambda_post,
  xlab = "MLE estimate",
  ylab = "Posterior mean",
  main = "MLE vs Posterior Mean: lambda"
)

abline(0, 1, lty = 2)

text(
  compare$lambda_mle,
  compare$lambda_post,
  labels = compare$subjID,
  pos = 3
)

# plot for tau
plot(
  compare$tau_mle,
  compare$tau_post,
  xlab = "MLE estimate",
  ylab = "Posterior mean",
  main = "MLE vs Posterior Mean: tau"
)

abline(0, 1, lty = 2)

text(
  compare$tau_mle,
  compare$tau_post,
  labels = compare$subjID,
  pos = 3
)


########### Question 2 ########### 

### (1)
# BIC values straight from HW2, Question 3
BIC1 <- 580.4341 # ra_prospect
BIC2 <- 651.8656 #ra_noLA

BF12 <- exp((BIC2 - BIC1) / 2)
BF12


### (3)
prior_odds = 0.6 / 0.4
prior_odds

post_odds = BF12 * prior_odds
post_odds

### (4)
post_M1 <- post_odds / (1 + post_odds)
post_M1

### (5)
options(digits = 6)

BIC_df <- data.frame(
  model = c("POW1", "POW2", "EXP1", "EXP2", "EXPOW", "HYP1", "HYP2"),
  BIC = c(475.292, 477.371, 524.532, 489.335, 479.451, 482.954, 481.634)
) # BIC values straight from HW2, Question 1
BIC_df

n <- nrow(BIC_df)
all_pairs_BF <- data.frame()

for (i in 1:(n-1)) {
  for (j in (i+1):n){
    
    BICi <- BIC_df$BIC[i]
    BICj <- BIC_df$BIC[j]
    
    BFij <- exp((BIC2 - BIC1) / 2)
    
    all_pairs_BF <- rbind(
      all_pairs_BF, 
      data.frame(model_i = BIC_df$model[i],
                 model_j = BIC_df$model[j],
                 BF_i_over_j = BFij
                 )
    )
  }
}

all_pairs_BF

