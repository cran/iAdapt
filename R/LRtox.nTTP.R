#' @title Calculates likelihood of safety for single dose, using nTTP
#' 
#' @description (nTTP) Function \code{LRtox.nTTP()} calculates the likelihood of safety for a single dose 
#' and designates whether to escalate to the next dose (safe) or stop dose escalation and move onto stage 2 (unsafe).
#' 
#' @return List object that gives the likelihood ratio of safety and indicates whether to escalate to the 
#' next highest dose level, or stop dose escalation and move onto stage 2.
#'         
#' @param tox_grades  data frame of observed AE grades for each patient (rows) across all toxicity types (columns).
#' e.g. for one patient, grades for 3 toxicity types might be c(3, 2, 4), where they experienced a grade 3 AE for tox type 1,
#' grade 2 AE for tox type 2, etc.
#' @param ntox number (integer) of different toxicity types
#' @param coh.size  cohort size (number of patients) per dose (Stage 1) 
#' @param W  matrix defining burden weight of each grade level for all toxicity types. 
#' The dimensions are ntox rows by 5 columns (for grades 0-4).
#' See Ezzalfani et al. (2013) for details.
#' @param p1  toxicity under null (unsafe DLT rate). Values range from 0 - 1.
#' @param p2  toxicity under alternative (safe DLT rate). Values range from 0 - 1; p1 > p2
#' @param K  threshold for LR. Takes integer values: 1,2,...(recommended K=2)
#' @param std.nTTP the standard deviation of nTTP scores at each dose level (constant across doses) 
#' 
#' @importFrom Matrix sparseMatrix
#' 
#' @examples 
#' ntox = 3 # three different types of toxicity 
#' coh.size = 3 # number of patients enrolled per dose
#' 
#' # Observed AE grades for each patient on tested dose
#' obs = data.frame(tox1 = c(3, 2, 4),
#'                  tox2 = c(1, 1, 2),
#'                  tox3 = c(2, 3, 3))
#'                 
#' # Toxicity burden weight matrix
#' W = matrix(c(0, 0.5, 0.75, 1.0, 1.5, # Burden weight for grades 0-4 for toxicity 1
#'              0, 0.5, 0.75, 1.0, 1.5, # Burden weight for grades 0-4 for toxicity 2
#'              0, 0.00, 0.00, 0.5, 1), # Burden weight for grades 0-4 for toxicity 3
#'              nrow = ntox, byrow = TRUE) 
#'              
#' # Acceptable (p2) and unacceptable nTTP values
#' p1 <- 0.35                                     
#' p2 <- 0.10       
#'              
#' LRtox.nTTP(obs, ntox, coh.size, W, p1, p2, K = 2, std.nTTP = 0.15)                                
#' 
#' @export


LRtox.nTTP <- function(tox_grades, ntox, coh.size, W, p1, p2, K = 2, std.nTTP = 0.15) {
  
  # Calculate nTTPs for cohort
  # function to calculate nTTP for one observed patient
  cc = function(obs, ntox, W) {
    x = sparseMatrix(i = obs + 1, j = 1:3, dims = c(5, 3))
    scores = as.matrix(W %*% x)
    scores = diag(scores)
    thetamax = sum(W[, 5]^2)
    nttp = sqrt(sum(scores^2) / thetamax)
    
    return(nttp)
  }
  
  coh.nttp = apply(tox_grades, 1, cc, ntox = ntox, W = W)
  
  # Bounds for truncated normal
  a = 0; b = 1
  
  # Calculate LR
  l.p2   <- prod(sapply(coh.nttp, FUN = function(i){ dnorm((i - p2)/std.nTTP) })) / 
    (std.nTTP*(pnorm((b - p2)/std.nTTP) - pnorm((a - p2)/std.nTTP)))^coh.size # likelihood of acceptable/alternative hypothesis 
  l.p1   <- prod(sapply(coh.nttp, FUN = function(i){ dnorm((i - p1)/std.nTTP) })) / 
    (std.nTTP*(pnorm((b - p1)/std.nTTP) - pnorm((a - p1)/std.nTTP)))^coh.size # likelihood of unacceptable/null hypothesis
  LR     <- round(l.p2/l.p1, 2)  
  
  if (LR > 1/K) {
    print("Safe/Escalate") 
  } else {
    print("Unsafe/Stop")
  }
  
  return(list(LR = LR))
}
