## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE,
                      eval = TRUE)

## ------------------------------------------------------------------------
library(iAdapt)
std.nTTP = 0.15 # standard deviation of nTTP value

coh.size = 3 # number pts per dose
ntox <- 3 # Number of unique toxicities
d <- 6 # Number of dose levels
N <- 25 # maximum number of patients

# Variance of the efficacy endpoints used for stage 2 randomization, assumed known and constant across doses.
v <- rep(0.01, 6) 

K = 2 # for LRT

# Stopping rule: if dose 1 is the only safe dose, allocate up to 9 pts before ending the trial to collect more information
stop.rule <- 9 

# Dose-efficacy curve
m = c(10, 20, 30, 40, 70, 90)

#### Define the weight matrix, from Du et al.(2019)
W <- matrix(c(0, 0.5, 0.75, 1.0, 1.5, # Burden weight for grades 1-4 for toxicity 1
              0, 0.5, 0.75, 1.0, 1.5, # Burden weight for grades 1-4 for toxicity 2
              0, 0.0, 0.00, 0.5, 1), ## Burden weight for grades 1-4 for toxicity 3
            nrow = ntox, byrow = T)

#### Define an array to hold toxicitiy probabilities 
data("TOX"); TOX # laod sample TOX array

## Grade at which an AE is defined as DLT for each toxicity type
grade.thresh = c(3, 3, 4)

## ------------------------------------------------------------------------
# Obtain the expected mean toxicity score at each dose level
tru_mnTTP = get.thresh(dose = d, ntox = ntox, W = W, TOX = TOX)

# Get the expected probability of a DLT at each dose level
pDLT = dlt.prob(dose = d, TOX = TOX, ntox = ntox, grade.thresh = grade.thresh) 

# Specify hypotheses
h1 = 0.35
h2 = 0.10

## ---- echo=FALSE, fig.width=6, fig.height=5------------------------------
# plot mnTTP
plot(x = 1:d,
     y = tru_mnTTP,
     ylim = c(0, 0.5),
     type = 'b',
     ylab = NA, 
     xlab = NA)
# plot corresponsing DLT rate
lines(x = 1:d,
      y = pDLT,
      col = "red",
      type = 'b')
# plot benchmark value of 0.4 from Chiuzan
abline(h = 0.4,
       lty = 2)
# plot vertical line at intersection of benchmark with dose-DLT
abline(v = 5.6,
       lty = 4)
# plot vertical line at intersection of vertical line and dose-nTTP
abline(h = 0.35,
       lty = 3)
# labels
title(main = "Dose-level toxicity",
       ylab = "Toxicity",
       xlab = "Dose")
# curve labels
text(x = 1.8,
     y = 0.1,
     labels = "mean nTTP",
     cex = 0.7)
text(x = 4,
     y = 0.1,
     labels = "DLT probability",
     cex = 0.7,
     col = "red")

## ------------------------------------------------------------------------
set.seed(3)
tox.profile.nTTP(dose = d, p1 = h1, p2 = h2, K = K, coh.size = coh.size, ntox = ntox, W = W, TOX = TOX, std.nTTP = std.nTTP)

## ------------------------------------------------------------------------
set.seed(3)
safe.dose.nTTP(dose = d, p1 = h1, p2 = h2, K = K, coh.size = coh.size, W = W, TOX = TOX, ntox = ntox, std.nTTP = std.nTTP) 

## ------------------------------------------------------------------------
set.seed(3)
eff.stg1.nTTP(dose = d, p1 = h1, p2 = h2, K = K, m = m, v = v, coh.size = coh.size, nbb = 100, W = W, TOX = TOX, ntox = ntox, std.nTTP = std.nTTP) 

## ------------------------------------------------------------------------
set.seed(3)
rand.stg2.nTTP(dose = d, p1 = h1, p2 = h2, K = K, coh.size = coh.size, m = m, v = v, N = N, stop.rule = stop.rule, cohort = 1, samedose = TRUE, nbb = 100, W = W, TOX = TOX, ntox = ntox, std.nTTP = std.nTTP) 

## ---- results='hide'-----------------------------------------------------
set.seed(3)
sims = 1e2 # number of trials to simulate

simulations = sim.trials.nTTP(numsims = sims,
                              dose = d,
                              p1 = h1,
                              p2 = h2,
                              K = K,
                              coh.size = coh.size,
                              m = m,
                              v = v,
                              N = N,
                              W = W,
                              TOX = TOX,
                              ntox = ntox,
                              std.nTTP = std.nTTP)

## ------------------------------------------------------------------------
head(simulations$safe.d)
head(simulations$sim.Y)
head(simulations$sim.d)

## ------------------------------------------------------------------------
colSums(simulations$safe.d) / sims

## ------------------------------------------------------------------------
sim.tables = sim.summary(simulations)

