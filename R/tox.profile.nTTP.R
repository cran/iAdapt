#' @title Generate nTTPs toxicity scores and the likelihhood-ratio (LR) per dose 
#' 
#' @description The normalized total toxicity profiles (nTTP) are calculated by 
#' combining multiple toxicity grades and their weights. The nTTPs are considered 
#' a quasi-continuous toxicity measure that follows a normal distribution truncated to [0, 1].
#' The likelihood ratio per dose are based on nTTP toxicity.
#' 
#' @return 
#' \itemize{
#' \item mnTTP - 4-column matrix containing dose assignment, mean nTTP at each dose, 
#' cohort number, and likelihood ratio.
#' \item all_nttp - all observed nTTP values
#' }
#' 
#' @param dose  number of doses to be tested (scalar)
#' @param p1  toxicity under null (unsafe nTTP). Values range from 0 - 1.
#' @param p2  toxicity under alternative (safe nTTP). Values range from 0 - 1; p1 > p2
#' @param K  threshold for LR. Takes integer values: 1,2,...(recommended K=2)
#' @param coh.size  cohort size (number of patients) per dose (Stage 1) 
#' @param ntox  number (integer) of different toxicity types (e.g, hematological, neurological, GI)
#' @param W  matrix defines burden weight of each grade level for all toxicity types. 
#' The dimensions are ntox rows by 5 columns (for grades 0-4). See Ezzalfani et al. (2013) for details.
#' @param TOX  matrix array of toxicity probabilities. There should be ntox matrices. 
#' Each matrix represents one toxicity type, where probabilities of each toxicity grade 
#' are specified across each dose. Each matrix has the same dimensions: n rows, representing 
#' number of doses, and 5 columns (for grades 0-4). Probabilities across each dose (rows) must sum to 1. 
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
#'              nrow = ntox, byrow = TRUE)
#'            
#' 
#' # Array of toxicity event probabilities
#' TOX = array(NA, c(dose, 5, ntox)) 
#' 
#' TOX[, , 1] = matrix(c(0.823, 0.152, 0.022, 0.002, 0.001,  #prob of tox for dose 1 and tox type 1
#'                       0.791, 0.172, 0.032, 0.004, 0.001,  #prob of tox for dose 2 and tox type 1
#'                       0.758, 0.180, 0.043, 0.010, 0.009,  #prob of tox for dose 3 and tox type 1
#'                       0.685, 0.190, 0.068, 0.044, 0.013,  #prob of tox for dose 4 and tox type 1
#'                       0.662, 0.200, 0.078, 0.046, 0.014,  #prob of tox for dose 5 and tox type 1
#'                       0.605, 0.223, 0.082, 0.070, 0.020), #prob of tox for dose 6 and tox type 1
#'                       nrow = 6, byrow = TRUE)
#' TOX[, , 2] = matrix(c(0.970, 0.027, 0.002, 0.001, 0.000,  #prob of tox for dose 1 and tox type 2
#'                       0.968, 0.029, 0.002, 0.001, 0.000,  #prob of tox for dose 2 and tox type 2
#'                       0.813, 0.172, 0.006, 0.009, 0.000,  #prob of tox for dose 3 and tox type 2
#'                       0.762, 0.183, 0.041, 0.010, 0.004,  #prob of tox for dose 4 and tox type 2
#'                       0.671, 0.205, 0.108, 0.011, 0.005,  #prob of tox for dose 5 and tox type 2
#'                       0.397, 0.258, 0.277, 0.060, 0.008), #prob of tox for dose 6 and tox type 2
#'                       nrow = 6, byrow = TRUE)
#' TOX[, , 3] = matrix(c(0.930, 0.060, 0.005, 0.001, 0.004,  #prob of tox for dose 1 and tox type 3
#'                       0.917, 0.070, 0.007, 0.001, 0.005,  #prob of tox for dose 2 and tox type 3
#'                       0.652, 0.280, 0.010, 0.021, 0.037,  #prob of tox for dose 3 and tox type 3
#'                       0.536, 0.209, 0.031, 0.090, 0.134,  #prob of tox for dose 4 and tox type 3
#'                       0.015, 0.134, 0.240, 0.335, 0.276,  #prob of tox for dose 5 and tox type 3
#'                       0.005, 0.052, 0.224, 0.372, 0.347), #prob of tox for dose 6 and tox type 3
#'                       nrow = 6, byrow = TRUE)
#' 
#' tox.profile.nTTP(dose = dose, 
#' p1 = p1, 
#' p2 = p2, 
#' K = K, 
#' coh.size = coh.size, 
#' ntox = ntox, 
#' W = W, 
#' TOX = TOX, 
#' std.nTTP = std.nTTP)
#' 
#' 
#' @export

tox.profile.nTTP <- function(dose, p1, p2, K, coh.size, ntox, W, TOX, std.nTTP = 0.15){ 
  
  dose   <- c(1:dose) # vector of counts up to number of doses given
  stop   <- 0
  cohort <- 0
  i      <- 1
  x      <- c()
  nttp  <- NULL
  
  # bounds for nTTP (truncated normal distribution)
  a = 0
  b = 1
  
  while ((stop == 0) & (i <= length(dose))) {
    cohort <- cohort + 1                                        # current cohort corresponding to dose
    
    # nTTPs for all patients (size coh.size) on dose i
    coh.nttp  <- replicate(coh.size, nTTP.indiv.sim(W = W, 
                                                    TOX = TOX, 
                                                    ntox = ntox, 
                                                    dose = dose[i]))	# nTTPs for that dose based on tox prob
    
    # Calculate LR
    l.p2   <- prod(sapply(coh.nttp, FUN = function(i){ dnorm((i - p2)/std.nTTP) })) / 
      (std.nTTP*(pnorm((b - p2)/std.nTTP) - pnorm((a - p2)/std.nTTP)))^coh.size # likelihood of acceptable/alternative hypothesis 
    l.p1   <- prod(sapply(coh.nttp, FUN = function(i){ dnorm((i - p1)/std.nTTP) })) / 
      (std.nTTP*(pnorm((b - p1)/std.nTTP) - pnorm((a - p1)/std.nTTP)))^coh.size # likelihood of unacceptable/null hypothesis
    LR     <- round(l.p2/l.p1, 2)    
    
    x <- c(x, dose[i], mean(coh.nttp), cohort, LR)                              			
    
    # list of observed nTTP
    nttp = append(nttp, coh.nttp) 
    
    if (LR <= (1/K)) {       # stop escalation
      stop <- 1
    } else if (LR > (1/K)) { # escalate to next dose i + 1
      i <- i + 1	 
    }
    
  }       
  return(list(mnTTP = matrix(x, ncol = 4, byrow = TRUE), 
              all_nTTP = nttp))
} 
