setwd("~/Documents/GitHub/2026_Computational-Modeling/HW6/HW6_solutions/Q5")

###############################################################################
## Grid-based posterior + summary table  (seed = 08826)
###############################################################################
rm(list = ls())
set.seed(08826)

library(ggplot2)
library(gridExtra)

## ① Load data ---------------------------------------------------------------
df           <- read.table("ra_exampleData.txt", header = TRUE, sep = "\t")
allSubjs     <- sort(unique(df$subjID))
N            <- length(allSubjs)
T            <- 140                             # 140 trials per subject

gain_mat     <- loss_mat <- cert_mat <- gamble_mat <- matrix(NA, N, T)
for (i in seq_len(N)) {
  idx                 <- df$subjID == allSubjs[i]
  gain_mat  [i, ]     <- df$gain   [idx]
  loss_mat  [i, ]     <- abs(df$loss[idx])
  cert_mat  [i, ]     <- df$cert   [idx]
  gamble_mat[i, ]     <- df$gamble [idx]
}

## ② Minus-log-likelihood function -------------------------------------------
minusLL <- function(par, subj) {
  rho <- par[1]; tau <- par[2]; lambda <- par[3]
  g   <- gain_mat   [subj, ];  l <- loss_mat[subj, ]
  c   <- cert_mat   [subj, ];  y <- gamble_mat[subj, ]
  
  evSafe <- c^rho
  evGamb <- 0.5 * (g^rho - lambda * l^rho)
  pGamb  <- 1 / (1 + exp(tau * (evSafe - evGamb)))
  pGamb  <- pGamb * 0.9998 + 0.0001          # numerical stabilisation
  
  -sum( log(pGamb)*y + log(1 - pGamb)*(1 - y) )
}

## ③ Grid-posterior function --------------------------------------------------
grid_post <- function(subj, K = 10) {
  rho_g <- seq(0, 2,  length.out = K)
  tau_g <- seq(0, 5,  length.out = K)
  lam_g <- seq(0, 10, length.out = K)
  
  grd   <- expand.grid(rho = rho_g, tau = tau_g, lambda = lam_g,
                       KEEP.OUT.ATTRS = FALSE)
  lik   <- exp(-apply(grd, 1, minusLL, subj = subj))      # likelihood
  post  <- lik / sum(lik)                                 # normalise
  cbind(grd, post)
}

## ④ Compute posterior means for every subject -------------------------------
get_means <- function(K = 10) {
  res <- data.frame(
    subjID  = allSubjs,
    rho     = NA_real_,
    lambda  = NA_real_,
    tau     = NA_real_
  )
  
  for (s in seq_len(N)) {
    pst <- grid_post(s, K)
    
    # marginals & expectations
    marg_rho    <- rowsum(pst$post, pst$rho)    [,1]   # Σ_{τ,λ}
    marg_lambda <- rowsum(pst$post, pst$lambda) [,1]   # Σ_{ρ,τ}
    marg_tau    <- rowsum(pst$post, pst$tau)    [,1]   # Σ_{ρ,λ}
    
    rho_g <- sort(unique(pst$rho))
    lam_g <- sort(unique(pst$lambda))
    tau_g <- sort(unique(pst$tau))
    
    res$rho   [s] <- sum(rho_g * marg_rho)
    res$lambda[s] <- sum(lam_g * marg_lambda)
    res$tau   [s] <- sum(tau_g * marg_tau)
  }
  res
}

## ⑤ Tables for 10-grid and 30-grid ------------------------------------------
tab10 <- get_means(10)
tab30 <- get_means(30)

cat("### Posterior means (10-grid)\n");  print(tab10, row.names = FALSE); cat("\n")
cat("### Posterior means (30-grid)\n");  print(tab30, row.names = FALSE); cat("\n")

## ⑥ Save as CSV --------------------------------------------------------------
write.csv(tab10, "posterior_means_grid10.csv", row.names = FALSE)
write.csv(tab30, "posterior_means_grid30.csv", row.names = FALSE)
cat(">> Two CSV files have been saved in the working directory.\n")

## ⑦ (Optional) save / show bar plots ----------------------------------------
cols_param <- c(rho = "#C3A6FF", lambda = "#8ECFFF", tau = "#FFB3C6")

make_bar <- function(df_grid, par, subj, K) {
  step   <- diff(sort(unique(df_grid[[par]])))[1]
  ggplot(df_grid, aes_string(par, "post")) +
    geom_bar(stat = "identity", width = step * 0.8, fill = cols_param[[par]]) +
    labs(title = sprintf("%s – Subj%d", par, allSubjs[subj]),
         x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(hjust = .5, size = 9),
          axis.text.x  = element_text(angle = 0, vjust = 0.5))
}

plot_overview <- function(K) {
  plist <- vector("list", N * 3);  id <- 1
  for (s in seq_len(N)) {
    pst <- grid_post(s, K)
    plist[[id]] <- make_bar(pst, "rho",    s, K); id <- id + 1
    plist[[id]] <- make_bar(pst, "lambda", s, K); id <- id + 1
    plist[[id]] <- make_bar(pst, "tau",    s, K); id <- id + 1
  }
  grid.arrange(grobs = plist,
               layout_matrix = matrix(seq_len(N * 3), nrow = N, byrow = TRUE),
               top = sprintf("Grid approximation – %d grid / parameter", K))
}

# Uncomment to view in the console
plot_overview(10)
plot_overview(30)
