#' @title Simulate full trial (both stages) x times when using nTTP to measure toxicity
#' 
#' @description Results are displayed in a matrix format, where each row represents one trial simulation  
#' 
#' @return List of the following objects:
#' \itemize{
#' \item sim.Y - estimated efficacy per each dose assignment 
#' \item sim.d - dose assignment for each patient in the trial 
#' }
#'         
#'          
#' @param dose  number of doses to be tested (scalar)
#' @param ntox  number (integer) of different toxicity types
#' @param W  matrix defining burden weight of each grade level for all toxicity types. 
#' The dimensions are ntox rows by 5 columns (for grades 0-4).
#' See Ezzalfani et al. (2013) for details.
#' @param TOX  matrix array of toxicity probabilities. There should be ntox matrices. 
#' Each matrix represents one toxicity type, where probabilities of each toxicity 
#' grade are specified across each dose. Each matrix has the same dimensions: 
#' n rows, representing number of doses, and 5 columns (for grades 0-4). 
#' Probabilities across each dose (rows) must sum to 1. See Ezzalfani et al. (2013) for details.
#' 
#' 


nTTP.indiv.sim <- function(W, TOX, ntox, dose){
  
  W = W[, 2:5] # quick fix for grade 0 column input 
  
  # Simulate grade observed for each toxicity type, based on TOX probabilities
  Tox <- NA
  for (k in 1:ntox) {
    # Tox vector of observed toxicity's grades for each toxicity given the current dose
    Tox[k] = sample(0:4, size = 1, replace = TRUE, prob = TOX[dose, , k])
  }
  
  # Corresponding weights
  toxscores <- NA
  for (k in 1:ntox) {
    # ifelse statement bc W does not have zero weights for grade 0
    toxscores[k] <- ifelse(Tox[k] == 0, 
                           0, # weight for grade 0
                           W[k, Tox[k]]) # weight of AE given grade, tox type
  }
  
  # nTTP calculation
  thetamax = sum(W[, 4]^2)
  nTTP <- sqrt(sum(toxscores^2) / thetamax) # The observed nTTP for the patient 
  
  return(nTTP)
}
