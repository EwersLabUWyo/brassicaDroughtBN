# tinySimBN.R
# R version 3.3.1 (2016-06-21)
# June 4, 2017. Mallory B. Lai.
# Creating a tiny simulated BN for learning purposes. 

#-----------------------------------------------------------------------
library(bnlearn)
#-----------------------------------------------------------------------

### Discrete network:

# Create parent node x.
x <- as.factor(rep(c(1, 2, 3), each = 50))

# Create child node y. 
y <- as.factor(c(sample(1:3, 50, replace = T, prob = c(.7, .2, .1)), 
                 sample(2:4, 50, replace = T, prob = c(.4, .5, .1)), 
                 sample(7:9, 50, replace = T, prob = c(.3, .3, .4))))

# Combine into data frame. 
net <- data.frame(x, y)

# Use a tabu search to learn the network. 
bn <- tabu(net, score = "bde", iss = 10, tabu = 50)

# Plot the network. 
plot(bn)

# Fit the parameters conditional on its structure. 
fit <- bn.fit(bn, net, method = "bayes")

# Plot barchart of P(y|x).
bn.fit.barchart(fit$y)

# List probabilities associated with barchart. 
fit$y





