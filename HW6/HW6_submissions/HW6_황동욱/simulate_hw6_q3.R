# Same as simulate_hw6_model1.R, but with num_subjs=200, num_trials=100

# Simulation parameters
seed <- 2026
num_subjs  <- 200   # changed
num_trials <- 100   # changed
pr_correct_option1 <- 0.5
pr_correct_option2 <- 0.8

set.seed(seed)

simul_pars <- data.frame(alpha = rnorm(num_subjs, 0.20, 0.08),
                         beta  = rnorm(num_subjs, 2.00, 0.70),
                         subjID = 1:num_subjs)

all_data <- NULL

for (i in 1:num_subjs) {
  alpha <- simul_pars$alpha[i]
  beta  <- simul_pars$beta[i]

  payoff_option1 = rbinom(size=1, n = num_trials, prob = pr_correct_option1)
  payoff_option2 = rbinom(size=1, n = num_trials, prob = pr_correct_option2)

  payoff_option1[payoff_option1 == 0] = -1
  payoff_option2[payoff_option2 == 0] = -1
  payoff_both = data.frame(payoff_option1, payoff_option2)

  tmp_data = data.frame(subjID=NULL, trial=NULL, choice=NULL, outcome=NULL)

  sv = c(0, 0)

  for (t in 1:num_trials) {
    prob_choose2 = 1 / (1 + exp(beta * (sv[1] - sv[2])))

    choice = rbinom(size=1, n = 1, prob = prob_choose2)
    choice = choice + 1

    outcome = payoff_both[t, choice]

    PE = outcome - sv[choice]
    sv[choice] = sv[choice] + alpha * (outcome - sv[choice])

    tmp_data[t, "subjID"]  = i
    tmp_data[t, "trial"]   = t
    tmp_data[t, "choice"]  = choice
    tmp_data[t, "outcome"] = outcome
  }
  all_data = rbind(all_data, tmp_data)
}

write.table(all_data, file = "simul_data_hw6_q3.txt",
            row.names = F, col.names = T, sep = "\t")
