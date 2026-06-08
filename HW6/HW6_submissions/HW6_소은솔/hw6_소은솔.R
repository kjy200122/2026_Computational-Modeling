#################################################################
###HW6 소은솔 ##############
################################################################
#Q1-1 data simulation###########################################

# Simulation parameters
seed <- 2026    # do not change the seed number!
num_subjs  <- 30 # number of subjects
num_trials <- 300 # number of trials per subject
pr_correct_option1 <- 0.5  # reward probability in option 1
pr_correct_option2 <- 0.8  # reward probability in option 2

# Set seed
set.seed(seed)   # always set a seed number for this homework!

# True parameters 
simul_pars <- data.frame(alpha = rnorm(num_subjs, 0.20, 0.08),
                         beta = rnorm(num_subjs, 2.00, 0.70),
                         subjID  = 1:num_subjs)

# For storing simulated choice data for all subjects
all_data <- NULL

for (i in 1:num_subjs) {
  # Individual-level (i.e. per subject) parameter values
  alpha <- simul_pars$alpha[i]
  beta <- simul_pars$beta[i]
  
  # geneate payoff structure for each subject
  # Defaults for the two options 
  # option 1: 50% win (+1), 50% loss (-1)
  # option 2: 80% win (+1), 20% loss (-1)
  payoff_option1 = rbinom(size=1, n = num_trials, prob = pr_correct_option1)
  payoff_option2 = rbinom(size=1, n = num_trials, prob = pr_correct_option2)
  
  # Replace 0 with -1
  payoff_option1[payoff_option1 == 0] = -1   # if 0 --> replace it with -1
  payoff_option2[payoff_option2 == 0] = -1   # if 0 --> replace it with -1
  payoff_both = data.frame(payoff_option1, payoff_option2)
  
  # For storing simulated data for current subject
  # subjID = subject ID
  # trial = trial number
  # choice = choice made on each trial (1 or 2)
  # outcome = outcome reveived on each trial (1 or -1)
  tmp_data = data.frame( subjID=NULL, trial=NULL, choice=NULL, outcome=NULL)
  
  # initialize some variables
  sv = c(0, 0)  # stimulus value of two options
  
  for (t in 1:num_trials)  {
    # Prob of choosing option 2
    prob_choose2 = 1 / (1 + exp(beta * (sv[1] - sv[2])))  # exploration/exploitation parameter is set to 1
    
    # choice
    choice = rbinom(size=1, n = 1, prob = prob_choose2 )
    choice = choice + 1  # 0 or 1 --> 1 (option 1) or 2 (option 2)
    
    # outcome
    outcome = payoff_both[t, choice]
    
    # after receiving outcome (feedback), update sv[t+1]
    # prediction error (PE)
    PE = outcome - sv[choice]
    
    # update stimulus value (sv) of the chosen option
    sv[choice] = sv[choice] + alpha * (outcome - sv[choice] )
    
    # append simulated task/response to subject data
    tmp_data[t, "subjID"] = i
    tmp_data[t, "trial"] = t
    tmp_data[t, "choice"] = choice
    tmp_data[t, "outcome"]    = outcome
  } # end of t loop
  # Append current subject with all subjects' data
  all_data = rbind(all_data, tmp_data)
}

# Write out data
write.table(all_data, file = "~/Downloads/simul_data_hw6_model1.txt", row.names = F, col.names = T, sep = "\t")

################################################################
#Q1-2 parameter recovery###########################################

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
output = stan("hw6_model1.stan", data = dataList, pars = c("mu_alpha", "mu_beta", "alpha", "beta"), 
              iter = 2000, warmup=1000, chains=2, cores=2)

# variational inference
m = stan_model("hw6_model1.stan")
output_vb = vb(m, data = dataList, pars = c("mu_alpha", "mu_beta", "alpha", "beta")
)
# traceplot
traceplot(output)

# print summary of parameters 
print(output)

# extract Stan fit object (parameters)
parameters <- rstan::extract(output)

########################################################################
#######2a visualisation of posterior distributions 

#group-level parameters 
stan_hist(output, pars = c("mu_alpha", "mu_beta"))

#individual level parameters 
stan_hist(output, pars = "alpha")
stan_hist(output, pars = "beta")

######2b scatterplot of true and estimated parameters 

alpha_mean = apply(parameters$alpha, 2, mean)
alpha_sd   = apply(parameters$alpha, 2, sd)
beta_mean  = apply(parameters$beta,  2, mean)
beta_sd    = apply(parameters$beta,  2, sd)

plot(simul_pars$alpha, alpha_mean, xlim=c(0, 0.5), ylim=c(0, 0.5),
     xlab="True alpha", ylab="Estimated alpha (posterior mean ± 1 SD)",
     main="Alpha recovery", pch=19)
abline(0, 1, lty=2)
arrows(x0=simul_pars$alpha, y0=alpha_mean-alpha_sd, y1=alpha_mean+alpha_sd,
       length=0.02, angle=90, code=3)

plot(simul_pars$beta, beta_mean, xlim=c(0, 4), ylim=c(0, 4),
     xlab="True beta", ylab="Estimated beta (posterior mean ± 1 SD)",
     main="Beta recovery", pch=19)
abline(0, 1, lty=2)
arrows(x0=simul_pars$beta, y0=beta_mean-beta_sd, y1=beta_mean+beta_sd,
       length=0.05, angle=90, code=3)

########################################################################
#Q2 comparing hierarchical/non-hierarchical model outputs###############

dat = read.table("simul_data_hw6_model1.txt", header=T, sep="\t")
allSubjs = unique(dat$subjID)
N = length(allSubjs)
T = table(dat$subjID)[1]
choice  <- array(-1, c(N, T))
outcome <- array(0, c(N, T))
for (i in 1:N) {
  curSubj = allSubjs[i]
  tmp     = subset(dat, subjID == curSubj)
  choice[i, 1:T]  <- tmp$choice
  outcome[i, 1:T] <- tmp$outcome
}
dataList <- list(
  N = N, T = T, Tsubj = rep(T, N),
  choice = choice, outcome = outcome
)

# Run the NON-hierarchical model
output2 = stan("hw6_model1_nh.stan", data = dataList,
               pars = c("alpha", "beta"),
               iter = 2000, warmup = 1000, chains = 2, cores = 2)

#arbitrary alpha/beta parameters for convergence check 
traceplot(output2, pars = c("alpha[1]", "alpha[15]", "beta[1]", "beta[15]"))

#summary
print(output2)

#extraction for recovery plots 
parameters2 <- rstan::extract(output2)
alpha_mean2 = apply(parameters2$alpha, 2, mean)
alpha_sd2   = apply(parameters2$alpha, 2, sd)
beta_mean2  = apply(parameters2$beta,  2, mean)
beta_sd2    = apply(parameters2$beta,  2, sd)

# Same scatterplots as 2b
plot(simul_pars$alpha, alpha_mean2, xlim=c(0, 0.5), ylim=c(0, 1),
     xlab="True alpha", ylab="Estimated alpha (non-hierarchical)",
     main="Alpha recovery (non-hierarchical)", pch=19)
abline(0, 1, lty=2)
arrows(x0=simul_pars$alpha, y0=alpha_mean2-alpha_sd2, y1=alpha_mean2+alpha_sd2,
       length=0.02, angle=90, code=3)

plot(simul_pars$beta, beta_mean2, xlim=c(0, 4), ylim=c(-1, 4),
     xlab="True beta", ylab="Estimated beta (non-hierarchical)",
     main="Beta recovery (non-hierarchical)", pch=19)
abline(0, 1, lty=2)
arrows(x0=simul_pars$beta, y0=beta_mean2-beta_sd2, y1=beta_mean2+beta_sd2,
       length=0.05, angle=90, code=3)

########################################################################
#Q3 computing shrinkage in HBA##########################################
rm(list=ls())
# Simulation parameters
rm(list = ls()) 
seed <- 2026    # do not change the seed number!
num_subjs  <- 200 # number of subjects
num_trials <- 100 # number of trials per subject
pr_correct_option1 <- 0.5  # reward probability in option 1
pr_correct_option2 <- 0.8  # reward probability in option 2

# Set seed
set.seed(seed)   # always set a seed number for this homework!

# True parameters 
simul_pars <- data.frame(alpha = rnorm(num_subjs, 0.20, 0.08),
                         beta = rnorm(num_subjs, 2.00, 0.70),
                         subjID  = 1:num_subjs)

# For storing simulated choice data for all subjects
all_data <- NULL

for (i in 1:num_subjs) {
  # Individual-level (i.e. per subject) parameter values
  alpha <- simul_pars$alpha[i]
  beta <- simul_pars$beta[i]
  
  # geneate payoff structure for each subject
  # Defaults for the two options 
  # option 1: 50% win (+1), 50% loss (-1)
  # option 2: 80% win (+1), 20% loss (-1)
  payoff_option1 = rbinom(size=1, n = num_trials, prob = pr_correct_option1)
  payoff_option2 = rbinom(size=1, n = num_trials, prob = pr_correct_option2)
  
  # Replace 0 with -1
  payoff_option1[payoff_option1 == 0] = -1   # if 0 --> replace it with -1
  payoff_option2[payoff_option2 == 0] = -1   # if 0 --> replace it with -1
  payoff_both = data.frame(payoff_option1, payoff_option2)
  
  # For storing simulated data for current subject
  # subjID = subject ID
  # trial = trial number
  # choice = choice made on each trial (1 or 2)
  # outcome = outcome reveived on each trial (1 or -1)
  tmp_data = data.frame( subjID=NULL, trial=NULL, choice=NULL, outcome=NULL)
  
  # initialize some variables
  sv = c(0, 0)  # stimulus value of two options
  
  for (t in 1:num_trials)  {
    # Prob of choosing option 2
    prob_choose2 = 1 / (1 + exp(beta * (sv[1] - sv[2])))  # exploration/exploitation parameter is set to 1
    
    # choice
    choice = rbinom(size=1, n = 1, prob = prob_choose2 )
    choice = choice + 1  # 0 or 1 --> 1 (option 1) or 2 (option 2)
    
    # outcome
    outcome = payoff_both[t, choice]
    
    # after receiving outcome (feedback), update sv[t+1]
    # prediction error (PE)
    PE = outcome - sv[choice]
    
    # update stimulus value (sv) of the chosen option
    sv[choice] = sv[choice] + alpha * (outcome - sv[choice] )
    
    # append simulated task/response to subject data
    tmp_data[t, "subjID"] = i
    tmp_data[t, "trial"] = t
    tmp_data[t, "choice"] = choice
    tmp_data[t, "outcome"]    = outcome
  } # end of t loop
  # Append current subject with all subjects' data
  all_data = rbind(all_data, tmp_data)
}

# Write out data
write.table(all_data, file = "simul_data_hw6_model2.txt", row.names = F, col.names = T, sep = "\t")

########################################################################
#Q3-(2) posterior distribution of individual and group parameters#######


dat = read.table("simul_data_hw6_model2.txt", header=T, sep="\t")

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
output = stan("hw6_model1.stan", data = dataList, pars = c("mu_alpha", "mu_beta", "alpha", "beta"), 
              iter = 2000, warmup=1000, chains=2, cores=2)

# variational inference
m = stan_model("hw6_model1.stan")
output_vb = vb(m, data = dataList, pars = c("mu_alpha", "mu_beta", "alpha", "beta")
)
# traceplot
traceplot(output)

# print summary of parameters 
print(output)

# extract Stan fit object (parameters)
parameters <- rstan::extract(output)

########################################################################
#visualisation#######

#group-level paramaeters 
stan_hist(output, pars = c("mu_alpha", "mu_beta"))

#individual level parameters 
stan_hist(output, pars = paste0("beta[", 1:30, "]"))
stan_hist(output, pars = paste0("alpha[", 1:30, "]"))

########################################################################
#Q3-(3) scatter plot####################################################

alpha_mean = apply(parameters$alpha, 2, mean)
alpha_sd   = apply(parameters$alpha, 2, sd)
beta_mean  = apply(parameters$beta,  2, mean)
beta_sd    = apply(parameters$beta,  2, sd)

plot(simul_pars$alpha, alpha_mean, xlim=c(0, 0.5), ylim=c(0, 0.5),
     xlab="True alpha", ylab="Estimated alpha (posterior mean ± 1 SD)",
     main="Alpha recovery", pch=19)
abline(0, 1, lty=2)
arrows(x0=simul_pars$alpha, y0=alpha_mean-alpha_sd, y1=alpha_mean+alpha_sd,
       length=0.02, angle=90, code=3)

plot(simul_pars$beta, beta_mean, xlim=c(0, 4), ylim=c(0, 4),
     xlab="True beta", ylab="Estimated beta (posterior mean ± 1 SD)",
     main="Beta recovery", pch=19)
abline(0, 1, lty=2)
arrows(x0=simul_pars$beta, y0=beta_mean-beta_sd, y1=beta_mean+beta_sd,
       length=0.05, angle=90, code=3)

########################################################################
#Q4-(1) posterior distribution using alpha_pos and alpha_neg ###########
rm(list=ls())

#data simulation 
seed <- 2026
num_subjs  <- 30
num_trials <- 100
pr_correct_option1 <- 0.5
pr_correct_option2 <- 0.8
set.seed(seed)

simul_pars = data.frame(alpha_pos = rnorm(num_subjs, 0.20, 0.08),
                        alpha_neg = rnorm(num_subjs, 0.30, 0.10),
                        beta      = rnorm(num_subjs, 2.00, 0.70),
                        subjID    = 1:num_subjs)

all_data <- NULL
for (i in 1:num_subjs) {
  alpha_pos <- simul_pars$alpha_pos[i]
  alpha_neg <- simul_pars$alpha_neg[i]
  beta      <- simul_pars$beta[i]
  
  payoff_option1 = rbinom(size = 1, n = num_trials, prob = pr_correct_option1)
  payoff_option2 = rbinom(size = 1, n = num_trials, prob = pr_correct_option2)
  payoff_option1[payoff_option1 == 0] = -1
  payoff_option2[payoff_option2 == 0] = -1
  payoff_both = data.frame(payoff_option1, payoff_option2)
  
  tmp_data = data.frame(subjID = NULL, trial = NULL, choice = NULL, outcome = NULL)
  sv = c(0, 0)
  
  for (t in 1:num_trials) {
    prob_choose2 = 1 / (1 + exp(beta * (sv[1] - sv[2])))
    choice = rbinom(size = 1, n = 1, prob = prob_choose2) + 1
    outcome = payoff_both[t, choice]
    
    PE = outcome - sv[choice]
    
    # KEY CHANGE: use different learning rate based on PE sign
    if (PE >= 0) {
      sv[choice] = sv[choice] + alpha_pos * PE
    } else {
      sv[choice] = sv[choice] + alpha_neg * PE
    }
    
    tmp_data[t, "subjID"]  = i
    tmp_data[t, "trial"]   = t
    tmp_data[t, "choice"]  = choice
    tmp_data[t, "outcome"] = outcome
  }
  all_data = rbind(all_data, tmp_data)
}

write.table(all_data, file = "simul_data_hw6_model3.txt",
            row.names = FALSE, col.names = TRUE, sep = "\t")

########################################################################
#########posterior distribution######################################### 

dat = read.table("simul_data_hw6_model3.txt", header = T, sep = "\t")
allSubjs = unique(dat$subjID)
N = length(allSubjs)
T = table(dat$subjID)[1]
choice  <- array(-1, c(N, T))
outcome <- array(0,  c(N, T))
for (i in 1:N) {
  curSubj = allSubjs[i]
  tmp     = subset(dat, subjID == curSubj)
  choice[i, 1:T]  <- tmp$choice
  outcome[i, 1:T] <- tmp$outcome
}
dataList <- list(
  N = N, T = T, Tsubj = rep(T, N),
  choice = choice, outcome = outcome
)

# Run
output3 = stan("hw6_model2.stan", data = dataList,
               pars = c("mu_alpha_pos", "mu_alpha_neg", "mu_beta",
                        "alpha_pos", "alpha_neg", "beta"),
               iter = 2000, warmup = 1000, chains = 2, cores = 2)

# Diagnostics and summary (this is part 4a / posterior summary)
print(output3, pars = c("mu_alpha_pos", "mu_alpha_neg", "mu_beta",
                        "alpha_pos", "alpha_neg", "beta"))
traceplot(output3, pars = c("mu_alpha_pos", "mu_alpha_neg", "mu_beta"))

# Group-level (hyper) parameters
stan_hist(output3, pars = c("mu_alpha_pos", "mu_alpha_neg", "mu_beta"))

# Individual-level parameters
stan_hist(output3, pars = "alpha_pos")
stan_hist(output3, pars = "alpha_neg")
stan_hist(output3, pars = "beta")

########################################################################
#Q4-(2) scatter plots #################################################

# Extract for recovery plots
parameters3 <- rstan::extract(output3)

alpha_pos_mean = apply(parameters3$alpha_pos, 2, mean)
alpha_pos_sd   = apply(parameters3$alpha_pos, 2, sd)
alpha_neg_mean = apply(parameters3$alpha_neg, 2, mean)
alpha_neg_sd   = apply(parameters3$alpha_neg, 2, sd)
beta_mean      = apply(parameters3$beta, 2, mean)
beta_sd        = apply(parameters3$beta, 2, sd)

# Three recovery scatterplots

plot(simul_pars$alpha_pos, alpha_pos_mean, xlim = c(0, 0.5), ylim = c(0, 0.5),
     xlab = "True alpha_pos", ylab = "Estimated alpha_pos",
     main = "Alpha_pos recovery", pch = 19)
abline(0, 1, lty = 2)
arrows(x0 = simul_pars$alpha_pos,
       y0 = alpha_pos_mean - alpha_pos_sd, y1 = alpha_pos_mean + alpha_pos_sd,
       length = 0.02, angle = 90, code = 3)

plot(simul_pars$alpha_neg, alpha_neg_mean, xlim = c(0, 0.6), ylim = c(0, 0.6),
     xlab = "True alpha_neg", ylab = "Estimated alpha_neg",
     main = "Alpha_neg recovery", pch = 19)
abline(0, 1, lty = 2)
arrows(x0 = simul_pars$alpha_neg,
       y0 = alpha_neg_mean - alpha_neg_sd, y1 = alpha_neg_mean + alpha_neg_sd,
       length = 0.02, angle = 90, code = 3)

plot(simul_pars$beta, beta_mean, xlim = c(0, 4), ylim = c(0, 4),
     xlab = "True beta", ylab = "Estimated beta",
     main = "Beta recovery", pch = 19)
abline(0, 1, lty = 2)
arrows(x0 = simul_pars$beta,
       y0 = beta_mean - beta_sd, y1 = beta_mean + beta_sd,
       length = 0.05, angle = 90, code = 3)

########################################################################
#####Q4-(3) changing # of trials and prob of options ###################

seed <- 2026
num_subjs  <- 30
num_trials <- 300
pr_correct_option1 <- 0.3
pr_correct_option2 <- 0.6
set.seed(seed)

simul_pars = data.frame(alpha_pos = rnorm(num_subjs, 0.20, 0.08),
                        alpha_neg = rnorm(num_subjs, 0.30, 0.10),
                        beta      = rnorm(num_subjs, 2.00, 0.70),
                        subjID    = 1:num_subjs)

all_data <- NULL
for (i in 1:num_subjs) {
  alpha_pos <- simul_pars$alpha_pos[i]
  alpha_neg <- simul_pars$alpha_neg[i]
  beta      <- simul_pars$beta[i]
  
  payoff_option1 = rbinom(size = 1, n = num_trials, prob = pr_correct_option1)
  payoff_option2 = rbinom(size = 1, n = num_trials, prob = pr_correct_option2)
  payoff_option1[payoff_option1 == 0] = -1
  payoff_option2[payoff_option2 == 0] = -1
  payoff_both = data.frame(payoff_option1, payoff_option2)
  
  tmp_data = data.frame(subjID = NULL, trial = NULL, choice = NULL, outcome = NULL)
  sv = c(0, 0)
  
  for (t in 1:num_trials) {
    prob_choose2 = 1 / (1 + exp(beta * (sv[1] - sv[2])))
    choice = rbinom(size = 1, n = 1, prob = prob_choose2) + 1
    outcome = payoff_both[t, choice]
    
    PE = outcome - sv[choice]
    
    # KEY CHANGE: use different learning rate based on PE sign
    if (PE >= 0) {
      sv[choice] = sv[choice] + alpha_pos * PE
    } else {
      sv[choice] = sv[choice] + alpha_neg * PE
    }
    
    tmp_data[t, "subjID"]  = i
    tmp_data[t, "trial"]   = t
    tmp_data[t, "choice"]  = choice
    tmp_data[t, "outcome"] = outcome
  }
  all_data = rbind(all_data, tmp_data)
}

write.table(all_data, file = "simul_data_hw6_model3_1.txt",
            row.names = FALSE, col.names = TRUE, sep = "\t")

########################################################################
#########posterior distribution######################################### 

dat = read.table("simul_data_hw6_model3_1.txt", header = T, sep = "\t")
allSubjs = unique(dat$subjID)
N = length(allSubjs)
T = table(dat$subjID)[1]
choice  <- array(-1, c(N, T))
outcome <- array(0,  c(N, T))
for (i in 1:N) {
  curSubj = allSubjs[i]
  tmp     = subset(dat, subjID == curSubj)
  choice[i, 1:T]  <- tmp$choice
  outcome[i, 1:T] <- tmp$outcome
}
dataList <- list(
  N = N, T = T, Tsubj = rep(T, N),
  choice = choice, outcome = outcome
)

# Run
output4 = stan("hw6_model2.stan", data = dataList,
               pars = c("mu_alpha_pos", "mu_alpha_neg", "mu_beta",
                        "alpha_pos", "alpha_neg", "beta"),
               iter = 2000, warmup = 1000, chains = 2, cores = 2)

# Diagnostics and summary (this is part 4a / posterior summary)
print(output4, pars = c("mu_alpha_pos", "mu_alpha_neg", "mu_beta",
                        "alpha_pos", "alpha_neg", "beta"))
traceplot(output4, pars = c("mu_alpha_pos", "mu_alpha_neg", "mu_beta"))

# Group-level (hyper) parameters
stan_hist(output4, pars = c("mu_alpha_pos", "mu_alpha_neg", "mu_beta"))

# Individual-level parameters
stan_hist(output4, pars = "alpha_pos")
stan_hist(output4, pars = "alpha_neg")
stan_hist(output4, pars = "beta")

########################################################################
#scatter plots #########################################################

# Extract for recovery plots
parameters3 <- rstan::extract(output4)

alpha_pos_mean = apply(parameters3$alpha_pos, 2, mean)
alpha_pos_sd   = apply(parameters3$alpha_pos, 2, sd)
alpha_neg_mean = apply(parameters3$alpha_neg, 2, mean)
alpha_neg_sd   = apply(parameters3$alpha_neg, 2, sd)
beta_mean      = apply(parameters3$beta, 2, mean)
beta_sd        = apply(parameters3$beta, 2, sd)

# Three recovery scatterplots

plot(simul_pars$alpha_pos, alpha_pos_mean, xlim = c(0, 0.5), ylim = c(0, 0.7),
     xlab = "True alpha_pos", ylab = "Estimated alpha_pos",
     main = "Alpha_pos recovery", pch = 19)
abline(0, 1, lty = 2)
arrows(x0 = simul_pars$alpha_pos,
       y0 = alpha_pos_mean - alpha_pos_sd, y1 = alpha_pos_mean + alpha_pos_sd,
       length = 0.02, angle = 90, code = 3)

plot(simul_pars$alpha_neg, alpha_neg_mean, xlim = c(0, 0.6), ylim = c(0, 0.6),
     xlab = "True alpha_neg", ylab = "Estimated alpha_neg",
     main = "Alpha_neg recovery", pch = 19)
abline(0, 1, lty = 2)
arrows(x0 = simul_pars$alpha_neg,
       y0 = alpha_neg_mean - alpha_neg_sd, y1 = alpha_neg_mean + alpha_neg_sd,
       length = 0.02, angle = 90, code = 3)

plot(simul_pars$beta, beta_mean, xlim = c(0, 4), ylim = c(-1, 4),
     xlab = "True beta", ylab = "Estimated beta",
     main = "Beta recovery", pch = 19)
abline(0, 1, lty = 2)
arrows(x0 = simul_pars$beta,
       y0 = beta_mean - beta_sd, y1 = beta_mean + beta_sd,
       length = 0.05, angle = 90, code = 3)


########################################################################
#Q5 (1) bonus question##################################################

rm(list = ls())
graphics.off()

setwd("~/Downloads/HW2")
data <- read.table("ra_exampleData.txt", header = TRUE)
all_subjects <- unique(data$subjID)

#likelihood function from HW2
ra_prospect <- function(params, cert, gain, loss, gamble) {
  rho    <- params[1]
  lambda <- params[2]
  tau    <- params[3]
  T <- length(cert)
  sum_minusLL <- 0
  for (t in 1:T) {
    evSafe   <- cert[t]^rho
    evGamble <- 0.5 * (gain[t]^rho - lambda * abs(loss[t])^rho)
    pGamble  <- 1 / (1 + exp(tau * (evSafe - evGamble)))
    pGamble  <- pGamble * 0.9998 + 0.0001
    tmp_minusLL <- -log(pGamble) * gamble[t] - log(1 - pGamble) * (1 - gamble[t])
    sum_minusLL <- sum_minusLL + tmp_minusLL
  }
  return(sum_minusLL)
}

#grids and priors
n_grid <- 10
rho_grid    <- seq(0, 2,  length.out = n_grid)
lambda_grid <- seq(0, 10, length.out = n_grid)
tau_grid    <- seq(0, 5,  length.out = n_grid)

#to store the marginal posterior afterwards
results <- list()

# Loop over subjects
for (s in all_subjects) {
  subj_data <- subset(data, subjID == s)
  
  log_posterior <- array(NA, dim = c(n_grid, n_grid, n_grid))
  
  for (i in 1:n_grid) {
    for (j in 1:n_grid) {
      for (k in 1:n_grid) {
        params <- c(rho_grid[i], lambda_grid[j], tau_grid[k])
        log_lik <- -ra_prospect(params,
                                cert   = subj_data$cert,
                                gain   = subj_data$gain,
                                loss   = subj_data$loss,
                                gamble = subj_data$gamble)
        # Uniform priors → log prior is 0 (constant, drops out in normalization)
        log_posterior[i, j, k] <- log_lik
      }
    }
  }
  
  # Normalize in log space
  log_posterior_shifted <- log_posterior - max(log_posterior)
  posterior_unnorm <- exp(log_posterior_shifted)
  posterior <- posterior_unnorm / sum(posterior_unnorm)
  
  # Marginalize
  results[[as.character(s)]] <- list(
    rho    = apply(posterior, 1, sum),
    lambda = apply(posterior, 2, sum),
    tau    = apply(posterior, 3, sum)
  )
  
}

#histogram for posterior distributions 
par(mfrow = c(5, 3), mar = c(4, 4, 2, 1))

for (s in all_subjects) {
  r <- results[[as.character(s)]]
  
  plot(rho_grid, r$rho, type = "h", lwd = 5, lend = 1,
       xlab = expression(rho), ylab = "P",
       main = paste("Subject", s, "- rho"))
  points(rho_grid, r$rho, pch = 19)
  
  plot(lambda_grid, r$lambda, type = "h", lwd = 5, lend = 1,
       xlab = expression(lambda), ylab = "P",
       main = paste("Subject", s, "- lambda"))
  points(lambda_grid, r$lambda, pch = 19)
  
  plot(tau_grid, r$tau, type = "h", lwd = 5, lend = 1,
       xlab = expression(tau), ylab = "P",
       main = paste("Subject", s, "- tau"))
  points(tau_grid, r$tau, pch = 19)
}

par(mfrow = c(1, 1))


########################################################################
#Q5 (2) bonus question##################################################

#same as Q5-1 but using 30 grids per parameter 
rm(list = ls())
graphics.off()

setwd("~/Downloads/HW2")
data <- read.table("ra_exampleData.txt", header = TRUE)
all_subjects <- unique(data$subjID)

ra_prospect <- function(params, cert, gain, loss, gamble) {
  rho    <- params[1]
  lambda <- params[2]
  tau    <- params[3]
  T <- length(cert)
  sum_minusLL <- 0
  for (t in 1:T) {
    evSafe   <- cert[t]^rho
    evGamble <- 0.5 * (gain[t]^rho - lambda * abs(loss[t])^rho)
    pGamble  <- 1 / (1 + exp(tau * (evSafe - evGamble)))
    pGamble  <- pGamble * 0.9998 + 0.0001
    tmp_minusLL <- -log(pGamble) * gamble[t] - log(1 - pGamble) * (1 - gamble[t])
    sum_minusLL <- sum_minusLL + tmp_minusLL
  }
  return(sum_minusLL)
}

#grids and priors
n_grid <- 30
rho_grid    <- seq(0, 2,  length.out = n_grid)
lambda_grid <- seq(0, 10, length.out = n_grid)
tau_grid    <- seq(0, 5,  length.out = n_grid)

#to store the marginal posterior afterwards
results <- list()

# Loop over subjects
for (s in all_subjects) {
  subj_data <- subset(data, subjID == s)
  
  log_posterior <- array(NA, dim = c(n_grid, n_grid, n_grid))
  
  for (i in 1:n_grid) {
    for (j in 1:n_grid) {
      for (k in 1:n_grid) {
        params <- c(rho_grid[i], lambda_grid[j], tau_grid[k])
        log_lik <- -ra_prospect(params,
                                cert   = subj_data$cert,
                                gain   = subj_data$gain,
                                loss   = subj_data$loss,
                                gamble = subj_data$gamble)
        # Uniform priors → log prior is 0 (constant, drops out in normalization)
        log_posterior[i, j, k] <- log_lik
      }
    }
  }
  
  # Normalize in log space
  log_posterior_shifted <- log_posterior - max(log_posterior)
  posterior_unnorm <- exp(log_posterior_shifted)
  posterior <- posterior_unnorm / sum(posterior_unnorm)
  
  # Marginalize
  results[[as.character(s)]] <- list(
    rho    = apply(posterior, 1, sum),
    lambda = apply(posterior, 2, sum),
    tau    = apply(posterior, 3, sum)
  )
  
}

#histogram for posterior distributions 
par(mfrow = c(5, 3), mar = c(4, 4, 2, 1))

for (s in all_subjects) {
  r <- results[[as.character(s)]]
  
  plot(rho_grid, r$rho, type = "h", lwd = 5, lend = 1,
       xlab = expression(rho), ylab = "P",
       main = paste("Subject", s, "- rho"))
  points(rho_grid, r$rho, pch = 19)
  
  plot(lambda_grid, r$lambda, type = "h", lwd = 5, lend = 1,
       xlab = expression(lambda), ylab = "P",
       main = paste("Subject", s, "- lambda"))
  points(lambda_grid, r$lambda, pch = 19)
  
  plot(tau_grid, r$tau, type = "h", lwd = 5, lend = 1,
       xlab = expression(tau), ylab = "P",
       main = paste("Subject", s, "- tau"))
  points(tau_grid, r$tau, pch = 19)
}

par(mfrow = c(1, 1))


