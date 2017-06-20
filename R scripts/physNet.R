# physNet.R
# R version 3.3.1 (2016-06-21)
# June 15, 2017. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Creating Bayesian network from Lina's droughted phys experiments. 

#-----------------------------------------------------------------------
library(bnlearn)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)
#-----------------------------------------------------------------------

# Read in well-watered data. 
ww <- fread(file.choose())
 
# Remove blank cells. 
ww <- ww[1:33, 1:29]

# Keep only Photosynthesis, stomatal conductance, fluorescence, 
# soil moisture, total NSC, and total Starch. 
ww <- ww[, c(7, 11, 15, 16, 25, 28)]

# Read in droughted data. 
d <- fread(file.choose())

# Remove blank cells. 
d <- d[1:48, 1:29]

# Keep only Photosynthesis, stomatal conductance, fluorescence, 
# soil moisture, total NSC, and total Starch. 
d <- d[, c(7, 12, 16, 17, 25, 28)]

phys <- rbind(ww, d)
physVar <- names(phys)
phys <- lapply(phys, function(x){str_replace(x, "N/A", "0")})
phys <- lapply(phys, function(x){as.numeric(as.character(x))})
phys <- data.frame(matrix(unlist(phys), nrow = 81))

colnames(phys) <- physVar



discPhys <- discretize(phys, method = 'interval', breaks = 5)

wh <- data.frame(from = c("CONDUCTANCE", "MIDDAY SOIL MOISTURE (%)"), 
                 to = c("PHOTOSYNTHESIS", "CONDUCTANCE"))
bl <- tiers2blacklist(list(colnames(discPhys)[1], 
                           colnames(discPhys)[-1]))
bn <- rsmax2(discPhys, restrict = "aracne", maximize = "tabu", 
             score = "bde", maximize.args = list(iss = 15), 
             whitelist = wh, blacklist = bl)
plot(bn)



nodes <- names(discPhys)
start <- random.graph(nodes = nodes, method = "ic-dag", num = 100, 
                      every = 3)
netlist <- suppressWarnings(lapply(start, function(net){
  tabu(discPhys, score = "bde", tabu = 50, iss = 10)
}))


rnd <- custom.strength(netlist, nodes = nodes)
modelAvg <- rnd[(rnd$strength > .85) & (rnd$direction >= .5), ]
avg.start <- averaged.network(rnd, threshold = .85)

plot(avg.start)


