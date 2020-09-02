#' @title Identify safe/acceptable doses from stage 1 based on nTTP scores. 
#' 
#' @description Function \code{safe.dose.nTTP()} distinguishes acceptable from unacceptable doses
#' 
#' @return List of the following objects:
#' \itemize{
#' \item alloc.safe - matrix of assignments only for acceptable doses (to be used in stage 2) 
#' and their corresponding toxicities
#' \item alloc.total - vector of all dose assignments from stage 1 
#' \item n1 - total number of subjects allocated in stage 1
#' \item all_nttp - all observed nTTP values
#' }
#'          
#' @param dose  number of doses to be tested (scalar)
#' @param p1  toxicity under null (unsafe nTTP). Values range from 0 - 1.
#' @param p2  toxicity under alternative (safe nTTP). Values range from 0 - 1; p1 > p2
#' @param K  threshold for LR. Takes integer values: 1,2,...(recommended K=2)
#' @param coh.size  cohort size (number of patients) per dose (Stage 1) 
#' @param ntox  number (integer) of different toxicity types
#' @param W  matrix defining burden weight of each grade level for all toxicity types. 
#' The dimensions are ntox rows by 4 columns (for grades 0-4). 
#' See Ezzalfani et al. (2013) for details.
#' @param TOX  matrix array of toxicity probabilities. There should be ntox matrices. 
#' Each matrix represents one toxicity type, where probabilities of each toxicity grade 
#' are specified across each dose. Each matrix has the same dimensions: n rows, representing 
#' number of doses, and 5 columns (for grades 0-4). 
#' Probabilities across each dose (rows) must sum to 1. 
#' See Ezzalfani et al. (2013) for details.
#' @param std.nTTP the standard deviation of nTTP scores at each dose level (constant across doses) 
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
#' # Number of toxicity types
#' ntox <- 3
#' 
#' # Standard deviation of nTTP values
#' std.nTTP = 0.15
#' 
#' # Toxicity burden weight matrix
#' W = matrix(c(0, 0.5, 0.75, 1.0, 1.5, # Burden weight for grades 0-4 for toxicity 1
#'              0, 0.5, 0.75, 1.0, 1.5, # Burden weight for grades 0-4 for toxicity 2
#'              0, 0.00, 0.00, 0.5, 1), # Burden weight for grades 0-4 for toxicity 3
#'            nrow = ntox, byrow = TRUE)
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
#'                       nrow = 6, byrow = TRUE)
#' TOX[, , 2] = matrix(c(0.970, 0.027, 0.002, 0.001, 0.000,
#'                       0.968, 0.029, 0.002, 0.001, 0.000,
#'                       0.813, 0.172, 0.006, 0.009, 0.000,
#'                       0.762, 0.183, 0.041, 0.010, 0.004,
#'                       0.671, 0.205, 0.108, 0.011, 0.005,
#'                       0.397, 0.258, 0.277, 0.060, 0.008),
#'                       nrow = 6, byrow = TRUE)
#' TOX[, , 3] = matrix(c(0.930, 0.060, 0.005, 0.001, 0.004,
#'                       0.917, 0.070, 0.007, 0.001, 0.005,
#'                       0.652, 0.280, 0.010, 0.021, 0.037,
#'                       0.536, 0.209, 0.031, 0.090, 0.134,
#'                       0.015, 0.134, 0.240, 0.335, 0.276,
#'                       0.005, 0.052, 0.224, 0.372, 0.347),
#'                       nrow = 6, byrow = TRUE)
#' 
#' safe.dose.nTTP(dose = dose, p1 = p1, p2 = p2, K = K, coh.size = coh.size, 
#' W = W, TOX = TOX, ntox = ntox, std.nTTP = std.nTTP) 
#' 
#' 
#' @export


safe.dose.nTTP <- function(dose, p1, p2, K, coh.size, W, TOX, ntox, std.nTTP = 0.15) {
  
  # simulate 
  tox <- tox.profile.nTTP(dose = dose, p1 = p1, p2 = p2, K = K, 
                          coh.size = coh.size, ntox = ntox, W = W, TOX = TOX, std.nTTP = std.nTTP)
  res <- tox$mnTTP
  all_nttp <- tox$all_nTTP # all observed nTTP in stage 1
  
  alloc.total <- sort(rep(res[, 1], coh.size)) # dose allocation for all patients (stage 1)
  n1 <- nrow(res) * coh.size # sample size in stage 1
  
  unsafe.dose <- which(res[, 4] <= (1/K))
  if (length(unsafe.dose) == 0) {
    alloc.safe <- res[, 1:2]
  } else {
    alloc.safe <- res[res[, 1] != unsafe.dose, 1:2]
  }
  return(list(alloc.safe = alloc.safe, 
              alloc.total = alloc.total, 
              n1 = n1, 
              all_nttp = all_nttp))
}
