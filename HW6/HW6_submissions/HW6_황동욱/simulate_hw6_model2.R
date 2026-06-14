# Simulation for Q4: model with two learning rates (alpha_pos, alpha_neg)

# Simulation parameters
seed <- 2026
num_subjs  <- 30
num_trials <- 300             # Q4-(1),(2) setting
pr_correct_option1 <- 0.3
pr_correct_option2 <- 0.6

# For Q4-(3), change the three lines above to:
#   num_trials <- 300
#   pr_correct_option1 <- 0.3
#   pr_correct_option2 <- 0.6

set.seed(seed)

# True parameters (alpha_pos, alpha_neg, beta)
simul_pars <- data.frame(alpha_pos = rnorm(num_subjs, 0.20, 0.08),
                         alpha_neg = rnorm(num_subjs, 0.30, 0.10),
                         beta      = rnorm(num_subjs, 2.00, 0.70),
                         subjID    = 1:num_subjs)

all_data <- NULL

for (i in 1:num_subjs) {
  alpha_pos <- simul_pars$alpha_pos[i]
  alpha_neg <- simul_pars$alpha_neg[i]
  beta      <- simul_pars$beta[i]

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

    # prediction error
    PE = outcome - sv[choice]

    # update with different learning rates depending on PE sign
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

write.table(all_data, file = "simul_data_hw6_model2.txt",
            row.names = F, col.names = T, sep = "\t")
