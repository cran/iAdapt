#' @title Obtain average nTTP at each dose level
#' 
#' @description Obtain average nTTP at each dose level
#' 
#' @return Vector of average nTTP for each dose level.
#'          
#' @param dose  number of doses to be tested (scalar)
#' @param ntox  number (integer) of different toxicity types (e.g, hematological, neurological, GI)
#' @param W  matrix defines burden weight of each grade level for all toxicity types. 
#' The dimensions are ntox rows by 4 columns (for grades 0-4). See Ezzalfani et al. (2013) for details.
#' @param TOX  matrix array of toxicity probabilities. There should be ntox matrices. 
#' Each matrix represents one toxicity type, where probabilities of each toxicity grade 
#' are specified across each dose. Each matrix has the same dimensions: n rows, representing 
#' number of doses, and 5 columns (for grades 0-4). 
#' Probabilities across each dose (rows) must sum to 1. 
#' See Ezzalfani et al. (2013) for details.
#' 
#' @examples
#' # Number of test doses
#' dose = 6
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
#' get.thresh(dose = dose, ntox = ntox, W = W, TOX = TOX)
#'    
#' @export
#' 

get.thresh <- function(dose, ntox, W, TOX){
  
  W = W[, 2:5] # quick fix for grade 0 column input 
  thetamax = sum(W[, 4]^2)
  
  grades <- seq(from = 0, to = 4, by = 1)
  tox.type <- matrix(rep(grades, each = ntox), ncol = ntox, byrow = TRUE)
  tox.type <- as.data.frame(tox.type)
  possible_outcomes <- expand.grid(tox.type[, 1:ntox]) #Permutes all possible AE grades into df
  
  # weights corresponding to AE grades
  mapped.weight <- NA
  vec <- NULL
  for (i in 1:ntox) {
    for (k in 1:nrow(possible_outcomes)) {
      mapped.weight[k] <- ifelse(possible_outcomes[k, i] > 0, 
                                 W[i, possible_outcomes[k, i]], 
                                 0)
    }
    vec <- cbind(vec, mapped.weight)
  }
  
  
  mapped_data <- matrix(vec, ncol = ntox, byrow = FALSE)
  scores <- apply(mapped_data^2, 1, sum) # Calculate toxicity scores for each profile in data set
  normalized_scores <- sqrt(scores/thetamax) # Calculate nTTP
  
  prob.data <- array(NA, c(nrow(possible_outcomes), ntox, dose))
  
  # probabilities corresponding to AE grades
  for (j in 1:dose) {
    for (i in 1:nrow(possible_outcomes)) {
      for (k in 1:ntox) {
        tox.prob <- TOX[j, possible_outcomes[i, k] + 1, k] 
        prob.data[i, k, j] <- tox.prob
      }
    }
  }
  
  prob.tox <- apply(prob.data[,,], c(1, 3), prod)
  prob.scores <- normalized_scores*prob.tox
  thresh <- apply(prob.scores, 2, sum)
  return(thresh)
}
