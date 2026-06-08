# Install pacakges
if (!require(rstan)) install.packages('rstan')
if (!require(hBayesDM)) install.packages('hBayesDM')

library(rstan)
library(hBayesDM)


# Fit models using hBayesDM
output_ra_prospect = ra_prospect('./ra_exampleData.txt', 4000, 1000, 4, 4)

# Check effective sample size (n_eff) and convergence (Rhat)
output_ra_prospect$fit

# Check convergence with trace plots
plot(output_ra_prospect, type="trace")
rstan::traceplot(output_ra_prospect$fit, pars = 'rho')
rstan::traceplot(output_ra_prospect$fit, pars = 'lambda')
rstan::traceplot(output_ra_prospect$fit, pars = 'tau')

# Save output
save(output_ra_prospect, file = './hBayesDM_output/hBayesDM_ra_prospect.RData')
