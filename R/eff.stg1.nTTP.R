#' @title Generates efficacy outcomes for stage 1 when using nTTP to measure toxicity
#' 
#' @description Function \code{eff.stg1.nTTP()} uses a beta-binomial distribution to generate 
#' outcomes (Ys) corresponding to acceptable dose assignments from stage 1. 
#' 
#' @return List of efficacy outcomes for subjects enrolled during stage 1 (dose-escalation)
#' \itemize{
#' \item Y.safe - vector of efficacy outcomes for each subject assigned to an acceptable safe dose
#' \item d.safe - vector of dose allocation for each subject assigned to an acceptable safe dose
#' \item tox.safe - number of dose-limiting toxicities for each safe dose level
#' \item Y.alloc - vector of efficacy outcomes for all subjects from stage 1 
#' (acceptable and unsafe doses)
#' \item d.alloc - vector of dose allocation for all subjects from stage 1 
#' (acceptable and unsafe doses)
#' \item all_nttp - all observed nTTP values
#' }
#' 
#' @param dose  number of doses to be tested (scalar)
#' @param p1  toxicity under null (unsafe nTTP). Values range from 0 - 1.
#' @param p2  toxicity under alternative (safe nTTP). Values range from 0 - 1; p1 > p2
#' @param K  threshold for LR. Takes integer values: 1,2,...(recommended K=2)
#' @param coh.size  cohort size (number of patients) per dose (Stage 1) 
#' @param m  vector of mean efficacies per dose. Values range from 0 - 100. 
#' (e.g, T cell persistence - values b/w 5 and 80 per cent) 
#' @param v  vector of efficacy variances per dose. Values range from 0 - 1. (e.g., 0.01)
#' @param nbb  binomial parameter (default = 100 cells per patient)
#' @param ntox  number (integer) of different toxicity types
#' @param W  matrix defining burden weight of each grade level for all toxicity types. 
#' The dimensions are ntox rows by 4 columns (for grades 0-4). 
#' See Ezzalfani et al. (2013) for details.
#' @param TOX  matrix array of toxicity probabilities. There should be ntox matrices. 
#' Each matrix represents one toxicity type, where probabilities of each toxicity grade are 
#' specified across each dose. Each matrix has the same dimensions: 
#' n rows, representing number of doses, and 5 columns (for grades 0-4). 
#' Probabilities across each dose (rows) must sum to 1. 
#' See Ezzalfani et al. (2013) for details.
#' @param std.nTTP the standard deviation of nTTP scores at each dose level 
#' (assumed constant across doses) 
#' 
#' @examples 
#' # Number of pre-specified dose levels
#' dose <- 6      
#' 
#' # Acceptable (p2) and unacceptable nTTP values
#' p1 <- 0.35                                     
#' p2 <- 0.10    
#' 
#' # Likelihood-ratio (LR) threshold
#' K <- 2                                          
#' 
#' # Cohort size used in stage 1
#' coh.size <- 3 
#' 
#' # Efficacy (equal) variance per dose
#' v <- rep(0.01, 6)
#' 
#' # Dose-efficacy curve
#' m = c(10, 20, 30, 40, 70, 90)
#' 
#' # Number of toxicity types
#' ntox <- 3
#' 
#' # Toxicity burden weight matrix
#' W = matrix(c(0, 0.5, 0.75, 1.0, 1.5, # Burden weight for grades 0-4 for toxicity 1
#'              0, 0.5, 0.75, 1.0, 1.5, # Burden weight for grades 0-4 for toxicity 2
#'              0, 0.00, 0.00, 0.5, 1), # Burden weight for grades 0-4 for toxicity 3
#'            nrow = ntox, byrow = TRUE)
#'            
#' 
#' # standard deviation of nTTP values
#' std.nTTP = 0.15
#' 
#' # Array of toxicity event probabilities
#' TOX <- array(NA, c(dose, 5, ntox)) 
#' 
#' TOX[, , 1] = matrix(c(0.823, 0.152, 0.022, 0.002, 0.001,
#'                       0.791, 0.172, 0.032, 0.004, 0.001,
#'                       0.758, 0.180, 0.043, 0.010, 0.009,
#'                       0.685, 0.190, 0.068, 0.044, 0.013,
#'                       0.662, 0.200, 0.078, 0.046, 0.014,
#'                       0.605, 0.223, 0.082, 0.070, 0.020),
#'                     nrow = 6, byrow = TRUE)
#' TOX[, , 2] = matrix(c(0.970, 0.027, 0.002, 0.001, 0.000,
#'                       0.968, 0.029, 0.002, 0.001, 0.000,
#'                       0.813, 0.172, 0.006, 0.009, 0.000,
#'                       0.762, 0.183, 0.041, 0.010, 0.004,
#'                       0.671, 0.205, 0.108, 0.011, 0.005,
#'                       0.397, 0.258, 0.277, 0.060, 0.008),
#'                     nrow = 6, byrow = TRUE)
#' TOX[, , 3] = matrix(c(0.930, 0.060, 0.005, 0.001, 0.004,
#'                       0.917, 0.070, 0.007, 0.001, 0.005,
#'                       0.652, 0.280, 0.010, 0.021, 0.037,
#'                       0.536, 0.209, 0.031, 0.090, 0.134,
#'                       0.015, 0.134, 0.240, 0.335, 0.276,
#'                       0.005, 0.052, 0.224, 0.372, 0.347),
#'                     nrow = 6, byrow = TRUE)
#' 
#' eff.stg1.nTTP(dose = dose, p1 = p1, p2 = p2, K = K, coh.size = coh.size, 
#' m = m, v = v, nbb = 100, W = W, TOX = TOX, ntox = ntox, std.nTTP = std.nTTP) 
#' 
#' 
#' @export


eff.stg1.nTTP <- function(dose, p1, p2, K, coh.size, m, v, nbb = 100, W, TOX, ntox, std.nTTP) {
  
  res <- safe.dose.nTTP(dose, p1, p2, K, coh.size, W, TOX, ntox, std.nTTP)
  d.alloc <- res$alloc.total
  val.safe <- res$alloc.safe
  Y.safe <- d.safe <- tox.safe <- Y.alloc <- NULL
  n1 <- res$n1
  all_nttp <- res$all_nttp
  
  # Simulate efficacy outcomes for all doses that enrolled patients
  for (i in 1:length(d.alloc)) {
    ab <- beta.ab(m[d.alloc[i]]/100, v[d.alloc[i]])
    p <- stats::rbeta(1, ab$a, ab$b)
    Y.alloc[i] <- 100 * stats::rbinom(1, nbb, p)/nbb
  }
  
  if (length(val.safe) > 2) { # if more than one dose found safe
    d.safe <- sort(rep(val.safe[, 1], coh.size))
    tox.safe <- res$alloc.safe[, 2]
    Y.safe <- Y.alloc[1:length(d.safe)]
  } else if (length(val.safe) == 2) { # if exactly one dose (the first) found safe
    d.safe <- sort(rep(val.safe[1], coh.size))
    tox.safe <- res$alloc.safe[2]
    Y.safe <- Y.alloc[1:length(d.safe)]
  } else { # if no doses found safe (first dose is too toxic)
    Y.safe <- d.safe <- NULL
    tox.safe <- res$alloc.safe[, 2]
  }
  
  return(list(Y.safe = Y.safe, 
              d.safe = d.safe, 
              tox.safe = tox.safe, 
              n1 = n1, 
              Y.alloc = Y.alloc, 
              d.alloc = d.alloc, 
              all_nttp = all_nttp))
}
