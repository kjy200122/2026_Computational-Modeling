rm(list=ls())  # remove all variables 

library(rstan)
library(hBayesDM)

# read_txt
data <- read.table("ra_exampleData.txt", header=TRUE, sep="\t")
data
subj_2 <- subset(data, subjID == 2)
subj_3 <- subset(data, subjID == 3)
subj_4 <- subset(data, subjID == 4)
subj_6 <- subset(data, subjID == 6)
subj_7 <- subset(data, subjID == 7)

#HW2 MLE estimation
mle_results <- data.frame(
  subjID = c(2, 3, 4, 6, 7),
  rho    = c(1.0199, 0.6927, 0.7944, 0.8656, 0.9616),
  lambda = c(0.7828, 2.4747, 1.0765, 0.9823, 1.4085),
  tau    = c(1.0071, 4.7421, 1.0355, 3.1299, 2.3468)
)

output1 = ra_prospect(data="ra_exampleData.txt", niter=4000, nwarmup=1000, nchain=4, ncore=4)

plot(output1)

plotInd(output1, "rho")
plotInd(output1, "lambda")
plotInd(output1, "tau")

print(output1$allIndPars) #bayes result
print(mle_results)

#plot MLE vs Bayesian
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))  # 1 row, 3 cols layout

params <- c("rho", "lambda", "tau")
subj_ids <- c(2, 3, 4, 6, 7)

for (p in params) {
  xvals <- mle_results[[p]]       # MLE values
  yvals <- output1$allIndPars[[p]]     # Bayesian posterior means
  
  # Find a symmetric axis range (so y=x line is centered)
  lim <- range(c(xvals, yvals)) + c(-0.1, 0.1)
  
  plot(xvals, yvals,
       xlim = lim, ylim = lim,
       pch = 19, col = "darkblue", cex = 1.5,
       xlab = paste("MLE", p),
       ylab = paste("hBayesDM posterior mean", p),
       main = paste(p, ": MLE vs Bayesian"))
  
  # y = x identity line (gray dashed)
  abline(0, 1, col = "gray60", lty = 2, lwd = 1.5)
  
  # Label each point with subject ID
  text(xvals, yvals, labels = subj_ids, pos = 4, cex = 0.8)
}

