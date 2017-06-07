# tinySimBN.R
# R version 3.3.1 (2016-06-21)
# June 4, 2017. Mallory B. Lai.
# Creating a tiny simulated BN for learning purposes. 

#-----------------------------------------------------------------------
library(bnlearn)
#-----------------------------------------------------------------------

### Discrete network x and y:
  
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

### Discrete network x, y, and z:

  # Create parent node x.
  x <- as.factor(rep(c(1, 2, 3), each = 50))
  
  # Create child node y. 
  y <- as.factor(c(sample(1:3, 50, replace = T, prob = c(.7, .2, .1)), 
                   sample(2:4, 50, replace = T, prob = c(.4, .5, .1)), 
                   sample(7:9, 50, replace = T, prob = c(.3, .3, .4))))
  
  # Create child node z. 
  z <- as.factor(c(sample(3:5, 100, replace = T, prob = c(.2, .6, .2)), 
                   sample(6:8, 50, replace = T, prob = c(.5, .1, .4))))
  
  
  # Combine into data frame. 
  net <- data.frame(x, y, z)
  
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
  
  # Plot barchart of P(z|x).
  bn.fit.barchart(fit$z)
  
  # List probabilities associated with barchart. 
  fit$z
  
# Re-run with iss > 10. 
  # Use a tabu search to learn the network. 
  bn <- tabu(net, score = "bde", iss = 15, tabu = 50)
  
  # Plot the network. 
  plot(bn)
  
  # Fit the parameters conditional on its structure. 
  fit <- bn.fit(bn, net, method = "bayes")
  
  # List probabilities associated with barchart. 
  fit$z # Notice the probabilities get more flat.

### Discrete network w, x, and y:
  
  # Create parent node x.
  w <- as.factor(rep(c(0, 1), times = c(100, 50)))
  
  # Create parent and child node x.
  x <- as.factor(rep(c(1, 2, 3), each = 50))
  
  # Create child node y. 
  y <- as.factor(c(sample(1:3, 50, replace = T, prob = c(.7, .2, .1)), 
                   sample(2:4, 50, replace = T, prob = c(.4, .5, .1)), 
                   sample(7:9, 50, replace = T, prob = c(.3, .3, .4))))
  
  # Combine into data frame. 
  net <- data.frame(w, x, y)
  
  # Use a tabu search to learn the network. 
  bn <- tabu(net, score = "bde", iss = 5, tabu = 50)
  
  # Plot the network. 
  plot(bn)
  
  # Fit the parameters conditional on its structure. 
  fit <- bn.fit(bn, net, method = "bayes")
  
  # Plot barchart of P(y|x).
  bn.fit.barchart(fit$y)
  
  # List probabilities associated with barchart. 
  fit$y
  
  # Plot barchart of P(x|w).
  bn.fit.barchart(fit$x)
  
  # List probabilities associated with barchart. 
  fit$x
  
### Continuous network x, y, and z:
  
  # Create parent node x.
  x <- rnorm(150, mean = 5, sd = 3)
  
  # Create child node y. 
  y <- rnorm(150, mean = x + 5, sd = 2)
  
  # Create child node y. 
  z <- rnorm(150, mean = x - 3, sd = 1)
  
  # Combine into data frame. 
  net <- data.frame(x, y, z)
  
  # Use a tabu search to learn the network. 
  bn <- tabu(net, tabu = 50)
  
  # Plot the network. 
  plot(bn)
  
  # Fit the parameters conditional on its structure. 
  fit <- bn.fit(bn, net, method = "mle")
  
  # List probabilities. 
  fit$y
  
### Continuous network x, y, and z:
  
  # Create parent node x.
  x <- rnorm(150, mean = 5, sd = 3)
  
  # Create child node y. 
  y <- rnorm(150, mean = x + 5, sd = 2)
  
  # Create child node y. 
  z <- rnorm(150, mean = y + 1, sd = 1)
  
  # Combine into data frame. 
  net <- data.frame(x, y, z)
  
  # Use a tabu search to learn the network. 
  bn <- tabu(net, tabu = 50)
  
  # Plot the network. 
  plot(bn)
  
  # Fit the parameters conditional on its structure. 
  fit <- bn.fit(bn, net, method = "mle")
  
  # List probabilities. 
  fit$y
  
