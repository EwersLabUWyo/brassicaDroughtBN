# simBNprior.R
# R version 3.3.1 (2016-06-21)
# June 7, 2017. Mallory B. Lai.
# Creating a small simulated BN for learning how BN can be used
# as to inform other Bayesian models as priors. 

#-----------------------------------------------------------------------
library(bnlearn)
library(rjags)
#-----------------------------------------------------------------------



modelStr="

model {
  csup ~ dcat(sp);
  cdiam ~ dnorm(mu[csup], 1/sigma^2);
}

"

modelSpec <- textConnection("modelStr")

sp <- c(0.5, 0.5)
mu <- c(6.1, 6.25)
sigma <- 0.05
jags.data <- list(sp = sp, mu = mu, sigma = sigma, cdiam = 6.20)
model1 <- jags.model(modelSpec, data = jags.data, n.chains = 2)



