# simBNprior.R
# R version 3.3.1 (2016-06-21)
# June 7, 2017. Mallory B. Lai.
# Creating a small simulated BN for learning how BN can be used
# as to inform other Bayesian models as priors. 

#-----------------------------------------------------------------------
library(bnlearn)
library(rjags)
#-----------------------------------------------------------------------

pest.dag <- model2network("[PR][CL][G1|PR:CL][G2|G1][TR|G1][LO|G2:TR]")


###################################################
### code chunk number 12: chapter3.rnw:416-422
###################################################
dat0 <- list(p.PR = c(0.7, 0.2, 0.1),
             a.CL = 3, b.CL = 1,
             g.G1 = c(1, 3, 10),
             k.G2 = 10,
             m.TR = 5, s.TR = 2.5,
             r.LO = 1/3, d.LO = 1)


###################################################
### code chunk number 13: chapter3.rnw:474-475
###################################################
set.seed(123)


###################################################
### code chunk number 14: pest.PR123
###################################################
exp.loss  <- rep(NA, 3)
names(exp.loss) <- paste("PR=", 1:3, sep = "")
qua.loss <- exp.loss
for (PR in 1:3) {
  dat1 <- dat0
  dat1$PR <- PR
  mopest <- jags.model(file = "inclu.pest.jam", data = dat1,
                       quiet = TRUE)
  update(mopest, 3000)
  sipest <- 
    coda.samples(model = mopest, variable.names = "LO",
                 n.iter  =  50000)
  summa <- summary(sipest)
  exp.loss[PR] <- summa$statistics["Mean"]
  qua.loss[PR] <- summa$quantiles["75%"]
}#FOR
mean3 <- mean(sipest[[1]][, "LO"])
round(c(exp.loss, MEAN = mean(exp.loss)), 1)


###################################################
### code chunk number 15: pest.PR1
###################################################
set.seed(567)
nbs <- 50000
PR <- 1
g <- function(pr) c(1, 3, 10)[pr] 
CL <- rbeta(nbs, 3, 1)
G1 <- rpois(nbs, CL * g(PR))
G2 <- rpois(nbs, G1 * 10)
il <- function(x) { 
  exp((x - 5) / 2.5)/(1 + exp((x - 5) / 2.5))
}#IL
TR <- rbinom(nbs, 1, il(G1))
x.lo <- G2 * (1 - (1-1/3)*TR)
LO <- rchisq(nbs, 1, ncp = x.lo)
round(mean(LO), 1)


###################################################
### code chunk number 16: pest.PRquantile
###################################################
round(qua.loss)


#### Simulate data from p.76
set.seed(567)
nbs <- 50000
PR <- sample(1:3, nbs, prob = c(.7, .2, .1), replace = T)
g <- function(pr)c(1, 3, 10)[pr]
CL <- rbeta(nbs, 3,1)
G1 <- rpois(nbs, CL * g(PR))
G2 <- rpois(nbs, G1 * 10)
il <- function(x){
  exp((x - 5) / 2.5)/(1 + exp((x-5) / 2.5))
}
TR <- rbinom(nbs, 1, il(G1))
x.lo <- G2 * (1 - (1 - 1/3) * TR)
LO <- rchisq(nbs, 1, ncp = x.lo)

crop <- data.frame(CL = CL, G1 = G1, G2 = G2, TR = TR, LO = LO, PR = PR)

# Learn structure from simulated data. 
crop$G1 <- as.numeric(crop$G1)
crop$G2 <- as.numeric(crop$G2)
crop$TR <- as.numeric(crop$TR)
crop$PR <- as.numeric(crop$PR)

cropDisc <- discretize(crop, method = "interval", 
                       breaks = c(5, 5, 5, 2, 5, 3))

cropBN <- tabu(cropDisc, score = "bde", iss = 150, tabu = 50)
plot(cropBN)

nodes <- names(cropDisc)
start <- random.graph(nodes = nodes, method = "ic-dag", num = 100, 
                      every = 3)
netlist <- suppressWarnings(lapply(start, function(net){
  tabu(cropDisc, score = "bde", tabu = 50, iss = 50)
}))


rnd <- custom.strength(netlist, nodes = nodes)
modelAvg <- rnd[(rnd$strength > .85) & (rnd$direction >= .5), ]
avg.start <- averaged.network(rnd, threshold = .85)

plot(avg.start)

cropBN <- rsmax2(cropDisc, restrict = "si.hiton.pc", test = "x2", maximize = "tabu",
       score = "bde", maximize.args = list(iss = 15))
plot(cropBN)



#### MOST PROMISING:

nodes <- names(cropDisc)
start <- random.graph(nodes = nodes, method = "ic-dag", num = 100, 
                      every = 3)
netlist <- suppressWarnings(lapply(start, function(net){
  rsmax2(cropDisc, restrict = "aracne", maximize = "tabu",
         score = "bde", maximize.args = list(iss = 5))
}))


rnd <- custom.strength(netlist, nodes = nodes)
modelAvg <- rnd[(rnd$strength > .85) & (rnd$direction >= .5), ]
avg.start <- averaged.network(rnd, threshold = .85)

plot(avg.start)


avg.start$arcs <- rbind(avg.start$arcs, c("TR", "LO"))
plot(avg.start)

bnparam <- bn.fit(avg.start, cropDisc, method = "bayes")
bn.fit.barchart(bnparam$G1, xlab = "P()")

bnparam$G1
bnparam$CL

g.one <- cpdist(bnparam, nodes = ("G1"), evidence = (PR == "[0.998,1.67]"))
x <- as.data.frame(table(g.one))
x

g.one <- cpdist(bnparam, nodes = ("LO"), evidence = (PR == "[0.998,1.67]"))
x <- as.data.frame(table(g.one))
x

g.one <- cpdist(bnparam, nodes = ("TR"), evidence = (PR == "[0.998,1.67]"))
x <- as.data.frame(table(g.one))
x

g.one <- cpdist(bnparam, nodes = ("G2"), evidence = (PR == "[0.998,1.67]"))
x <- as.data.frame(table(g.one))
x

g.one <- cpdist(bnparam, nodes = ("G1"), evidence = (PR == "[0.998,1.67]") | 
                  (PR == "(1.67,2.33]") | (PR == "(2.33,3]") & ((CL == "[0.02,0.217]") |
                  (CL == "(0.217,0.413]") | (CL == "(0.413,0.608]") | 
                  (CL == "(0.608,0.804]") | (CL == "(0.804,1]")))
x <- as.data.frame(table(g.one))
x
