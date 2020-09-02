#' @title Calculate DLT probability corresponding to average nTTP for each dose
#' 
#' @description Calculate DLT probability corresponding to average nTTP for each dose
#' 
#' @return ptox - Vector of DLT probabilities per dose.
#'      
#' @param dose  number of doses to be tested (scalar)
#' @param ntox number (integer) of different toxicity types (e.g, hematological, neurological, GI)              
#' @param TOX  matrix array of toxicity probabilities. There should be ntox matrices. 
#' Each matrix represents one toxicity type, where probabilities of each toxicity grade 
#' are specified across each dose. Each matrix has the same dimensions: n rows, representing 
#' number of doses, and 5 columns (for grades 0-4, since the probability of a grade 0 event 
#' may not be 0). Probabilities across each dose (rows) must sum to 1. 
#' See Ezzalfani et al. (2013) for details.
#' @param grade.thresh grade (0-4) at which each toxicity type qualifies as a DLT
#' 
#' @examples
#' # Number of test doses
#' dose = 6
#' 
#' # Number of toxicity types
#' ntox <- 3
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
#' # Grades at which each tox type qualifies as DLT
#' grade.thresh = c(3, 3, 4)
#' 
#' dlt.prob(dose = dose, ntox = ntox, TOX = TOX, grade.thresh = grade.thresh)
#' 
#' @export

dlt.prob <- function(dose, ntox, TOX, grade.thresh){ 
  
  cp <- function(p) 
  {
    ev <- do.call(expand.grid, replicate(length(p), 0:1, simplify = FALSE)) # all permutations (rows) of DLT for each tox type
    pe <- apply(ev, 1, function(x) prod(p*(x == 1) + (1 - p)*(x == 0))) # Law of total prob ; DLT occurs with probability p 
    tapply(pe,rowSums(ev), sum) # sums probabilities for each group (total # DLTs across tox types)
  }
  
  ## probability of DLT for each dose at each tox type
  dlt.matrix <- matrix(0, nrow = nrow(TOX), ncol = dose)
  for (i in 1:ntox) {
    xxx <- TOX[, -(1:grade.thresh[i]), i] # for the dose and tox type, give only the cols corresponding to the grades that qualify as DLT
    xxx <- cbind(xxx, 0) # ensures that xxx is a matrix (instead of vector)
    dlt.matrix[i, ] <- apply(xxx, 1, sum) # why sum?
  }
  ## probability of DLT for each dose across tox types
  ptox <- NULL
  for (i in 1:dose) {
    ptox[i] <- 1 - cp(dlt.matrix[1:ntox, i])[1] #Ptox is probability of DLT at each dose level
  } 
  
  return(ptox)
}
