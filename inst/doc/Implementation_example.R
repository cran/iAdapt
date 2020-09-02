## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
library(iAdapt)

# Acceptable (p_yes) and unacceptable (p_no) DLT rates used for establishing safety
p_no <- 0.40                                     
p_yes <- 0.15    

# Likelihood-ratio (LR) threshold
K <- 2                                          

# Cohort size used in stage 1
coh.size <- 3 

# number of observed DLTs
x <- 1


## ------------------------------------------------------------------------
LRtox(coh.size, x, p_no, p_yes, K)

## ------------------------------------------------------------------------
LRtox(coh.size, ndlt = 2, p_no, p_yes, K)

## ------------------------------------------------------------------------
ntox = 3 # three different types of toxicity 
coh.size = 3 # number of patients enrolled per dose

# Observed AE grades for each patient on tested dose
obs = data.frame(tox1 = c(0, 1, 1),
                 tox2 = c(1, 0, 0),
                 tox3 = c(2, 0, 1))

# Toxicity burden weight matrix
W = matrix(c(0, 0.5, 0.75, 1.0, 1.5, # Burden weight for grades 0-4 for toxicity 1
             0, 0.5, 0.75, 1.0, 1.5, # Burden weight for grades 0-4 for toxicity 2
             0, 0.00, 0.00, 0.5, 1), # Burden weight for grades 0-4 for toxicity 3
           nrow = ntox, byrow = TRUE) 

# Acceptable (p2) and unacceptable nTTP values
p1 <- 0.35                                     
p2 <- 0.10       

LRtox.nTTP(obs, ntox, coh.size, W, p1, p2, K = 2, std.nTTP = 0.15) 

## ------------------------------------------------------------------------
y.eff <- c(9, 1, 0, 34, 10, 27, 38, 42, 60, 75, 48, 62)
d.safe <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4)
rand.prob(y.eff, d.safe)

