# code for devtools installation due to the error

# remove.packages(c("StanHeaders", "rstan"))
remove.packages("rstan")
if (file.exists(".RData")) file.remove(".RData")
# Install and load rstan
install.packages("rstan", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
library("rstan") #Load the rstan package

#Install and load hBayesDM
remove.packages("hBayesDM")
install.packages("D:/대학원/2026-봄/수업/임상심리세미나/HW5/hBayesDM_1.2.1.tar.gz", repos = NULL, type = "source")
library("hBayesDM")

# Programmed by Woo-Young Ahn (wahn55@snu.ac.kr)
# read the data file
dat = read.table("ra_exampleData.txt", header=T, sep="\t")

#Change the format for ra_prospect
dat_hb <- dat[, c("subjID", "gain", "loss", "cert", "gamble")]

#make a file
write.table(dat_hb, "ra_exampleData_hbayesdm.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# run!
output = ra_prospect(data = "ra_exampleData_hbayesdm.txt", niter = 4000, nwarmup = 1000, nchain = 4, ncore = 4)

#Plotting
#plot(output, type="trace", fontSize=11)   # traceplot of hyper parameters. Set font size 11.
#plot(output)
plotInd(output,"rho") #individual plot

# individual posterior mean of rho
rho_result <- data.frame(
  subjID = output$allIndPars$subjID,
  rho_postmean = output$allIndPars$rho
)
rho_result


#For plotting MLE and Bayesian
# 직접 입력
rho_mle <- c(1.019878, 0.692663, 0.794402, 0.865562, 0.961615)
rho_postmean <- c(0.9488431, 0.6962745, 0.7812755, 0.8816324, 0.9279043)
subj <- c(2, 3, 4, 6, 7)

# 산점도
plot(rho_mle, rho_postmean,
     xlab = "MLE rho",
     ylab = "Bayesian posterior mean rho",
     pch = 19,
     col = "blue",
     xlim = range(c(rho_mle, rho_postmean)),
     ylim = range(c(rho_mle, rho_postmean)),
     main = "MLE vs Bayesian posterior Mean")

# 기준선
abline(0, 1, lty = 2, col = "red")

# subject 번호 표시
text(rho_mle, rho_postmean, labels = subj, pos = 3, cex = 0.9)
