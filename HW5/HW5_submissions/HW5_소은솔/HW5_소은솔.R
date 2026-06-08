########################################################################
# HW5 - Eunsol So
########################################################################
rm(list = ls()); graphics.off()
library(rstan)
library(hBayesDM)
set.seed(08826)
setwd("~/Downloads/HW5_소은솔")
source("HDIofMCMC.R")

dat      <- read.table("ra_exampleData.txt", header = TRUE, sep = "\t")
allSubjs <- sort(unique(dat$subjID))         # subjects 2, 3, 4, 6, 7
N        <- length(allSubjs)
T        <- as.numeric(table(dat$subjID)[1]) # trials per subject (=140)
numPars  <- 3

########################################################################
# QUESTION 1
########################################################################

# ----MLE estimates for the 5 subjects ----
ra_prospect <- function(params, cert, gain, loss, gamble){
  rho <- params[1]; lambda <- params[2]; tau <- params[3]
  T = 140
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
param_prospect_low <- c(0, 0, 0); param_prospect_up <- c(2, 10, 5)

mle_summary <- data.frame()
for (sid in allSubjs) {
  subj_data <- subset(dat, subjID == sid)
  best <- optim(runif(3), ra_prospect, method = "L-BFGS-B",
                lower = param_prospect_low, upper = param_prospect_up,
                cert = subj_data$cert, gain = subj_data$gain,
                loss = subj_data$loss, gamble = subj_data$gamble)
  for (i in 1:100) {
    temp <- optim(runif(3), ra_prospect, method = "L-BFGS-B",
                  lower = param_prospect_low, upper = param_prospect_up,
                  cert = subj_data$cert, gain = subj_data$gain,
                  loss = subj_data$loss, gamble = subj_data$gamble)
    if (temp$value < best$value) best <- temp
  }
  parm <- round(best$par, 3)
  mle_summary <- rbind(mle_summary, data.frame(
    subjID = sid, rho_mle = parm[1], lambda_mle = parm[2], tau_mle = parm[3]))
}
print(mle_summary)

# ----hBayesDM hierarchical fit of the same 5 subjects ----
subData <- subset(dat, subjID %in% allSubjs)
write.table(subData, "ra_fiveSubj.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)

output_hbayes <- hBayesDM::ra_prospect(
  data = "ra_fiveSubj.txt",
  niter = 4000, nwarmup = 1000, nchain = 4, ncore = 4)

# --- convergence check for the hBayesDM model ---
plot(output_hbayes, type = "trace")    
print(output_hbayes$fit)

# posterior distributions of the 5 subjects
plotInd(output_hbayes, "rho")
plotInd(output_hbayes, "lambda")
plotInd(output_hbayes, "tau")

hb <- output_hbayes$allIndPars
print(hb)

# ----comparison graphs: MLE vs hBayesDM, per parameter ----
cmp <- merge(mle_summary, hb, by = "subjID")

plot_param <- function(mle_vals, hb_vals, ids, param_name) {
  yr <- range(c(mle_vals, hb_vals)); x <- seq_along(ids)
  plot(x, mle_vals, pch = 19, col = "tomato",
       xlim = c(0.5, length(ids) + 0.5), ylim = yr,
       xaxt = "n", xlab = "Subject", ylab = "Parameter value", main = param_name)
  axis(1, at = x, labels = ids)
  points(x, hb_vals, pch = 17, col = "steelblue")
  segments(x, mle_vals, x, hb_vals, col = "gray60")
  legend("topright", legend = c("MLE", "hBayesDM"),
         pch = c(19, 17), col = c("tomato", "steelblue"), bty = "n")
}
par(mfrow = c(1, 3))
plot_param(cmp$rho_mle,    cmp$rho,    cmp$subjID, "rho (risk aversion)")
plot_param(cmp$lambda_mle, cmp$lambda, cmp$subjID, "lambda (loss aversion)")
plot_param(cmp$tau_mle,    cmp$tau,    cmp$subjID, "tau (inverse temp.)")
par(mfrow = c(1, 1))
print(cmp)


########################################################################
# QUESTION 2 -- Bayes factors from BIC values
########################################################################
BIC1 <- 580.433     # ra_prospect (model 1), summed BIC from HW2
BIC2 <- 651.865     # ra_noLA     (model 2), summed BIC from HW2

BF12 <- exp((BIC2 - BIC1) / 2)
print(BF12)

prior_odds     <- 0.6 / 0.4
posterior_odds <- prior_odds * BF12
print(posterior_odds)

pM1_D <- posterior_odds / (1 + posterior_odds)
format(pM1_D, digits = 17)

# (2e) seven memory-retention models
model_name <- c("POW2", "EXP2", "POW1", "EXP1", "EXPOW", "HYP1", "HYP2")
BIC_mem    <- c(477.371, 489.335, 475.292, 524.532, 479.451, 482.954, 481.634)
best       <- which.min(BIC_mem)
BF_vs_best <- exp((BIC_mem - BIC_mem[best]) / 2)
postProb   <- exp(-(BIC_mem - min(BIC_mem)) / 2)
postProb   <- postProb / sum(postProb)
cat("(2e) best model =", model_name[best], "\n")
print(data.frame(model = model_name, BIC = BIC_mem,
                 BF_best_vs_model = round(BF_vs_best, 3),
                 postProb = round(postProb, 4)))


########################################################################
# QUESTION 3 -- single-subject Stan model fit to each of the 5 subjects
########################################################################
# build [N, T] matrices -- one row per subject
gain_m   <- matrix(0, N, T); loss_m   <- matrix(0, N, T)
cert_m   <- matrix(0, N, T); gamble_m <- matrix(0, N, T)
for (i in 1:N) {
  tmp <- subset(dat, subjID == allSubjs[i])
  gain_m[i, ]   <- tmp$gain
  loss_m[i, ]   <- abs(tmp$loss)        # absolute value
  cert_m[i, ]   <- tmp$cert
  gamble_m[i, ] <- tmp$gamble
}
dataList <- list(N = N, T = T,
                 gain = gain_m, loss = loss_m,
                 cert = cert_m, gamble = gamble_m)

# run the modified model for all 5 subjects
output <- stan("ra_prospect_multiSubj.stan", data = dataList,
               pars = c("rho", "lambda", "tau"),
               iter = 4000, warmup = 1000, chains = 4, cores = 4)

rstan::traceplot(output)
print(output)


# extract -- 
parameters <- rstan::extract(output)

# ---- posterior distributions: 5 histograms per parameter ----
cols <- c("skyblue", "salmon", "lightgreen", "orange", "pink")
for (p in c("rho", "lambda", "tau")) {
  par(mfrow = c(1, 5))
  for (s in 1:N) {
    hist(parameters[[p]][, s],
         main = paste0(p, "_sub", allSubjs[s]),
         xlab = p, ylab = "frequency", col = cols[s])
  }
}
par(mfrow = c(1, 1))

# ---- Q3(2): posterior means + 95% HDI ----
result_stan <- data.frame()
for (s in 1:N) {
  hdi_rho    <- HDIofMCMC(parameters$rho[, s],    credMass = 0.95)
  hdi_lambda <- HDIofMCMC(parameters$lambda[, s], credMass = 0.95)
  hdi_tau    <- HDIofMCMC(parameters$tau[, s],    credMass = 0.95)
  result_stan <- rbind(result_stan, data.frame(
    RHO       = mean(parameters$rho[, s]),
    RHO_HDIlo = hdi_rho[1],    RHO_HDIhi = hdi_rho[2],
    LAMBDA    = mean(parameters$lambda[, s]),
    LAM_HDIlo = hdi_lambda[1], LAM_HDIhi = hdi_lambda[2],
    TAU       = mean(parameters$tau[, s]),
    TAU_HDIlo = hdi_tau[1],    TAU_HDIhi = hdi_tau[2]))
}
rownames(result_stan) <- paste0("subject", allSubjs)
print(round(result_stan, 3))
