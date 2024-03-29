#' @title Simulate full trial (both stages) x times
#' 
#' @description Results are displayed in a matrix format, where each row 
#' represents one trial simulation.
#' 
#' @return List of the following objects:
#' \itemize{
#' \item sim.Y - estimated efficacy per each dose assignment 
#' \item sim.d - dose assignment for each patient in the trial 
#' \item safe.d - indicator of whether dose was declared safe
#' }
#'         
#'          
#' @param numsims  number of simulated trials
#' @param dose  number of doses to be tested (scalar)
#' @param dose.tox  vector of true toxicities for each dose. Values range from 0 - 1.
#' @param p1  toxicity under null (unsafe DLT rate). Values range from 0 - 1.
#' @param p2  toxicity under alternative (safe DLT rate). Values range from 0 - 1; p1 > p2
#' @param K  threshold for LR. Takes integer values: 1,2,...(recommended K=2)
#' @param coh.size  cohort size (number of patients) per dose (Stage 1) 
#' @param m  vector of mean efficacies per dose. Values range from 0 - 100. 
#' (e.g, T cell persistence - values b/w 5 and 80 per cent) 
#' @param v  vector of efficacy variances per dose. Values range from 0 - 1. (e.g., 0.01)
#' @param nbb  binomial parameter (default = 100 cells per patient)
#' @param N  total sample size for stages 1&2
#' @param stop.rule  if only dose 1 safe, allocate up to 9 (default) patients at dose 1 
#' to collect more info.
#' @param cohort cohort size (number of patients) per dose (Stage 2). Default is 1.
#' @param samedose designates whether the next patient is allocated to the same dose 
#' as the previous patient. Default is TRUE. Function adjusts accordingly.
#' 
#' @examples
#' # Number of pre-specified dose levels
#' dose <- 5
#' 
#' # Vector of true toxicities associated with each dose
#' dose.tox <- c(0.05, 0.10, 0.20, 0.35, 0.45)       
#' 
#' # Acceptable (p_yes) and unacceptable (p_no) DLT rates used for establishing safety
#' p_no <- 0.40                                     
#' p_yes <- 0.15    
#' 
#' # Likelihood-ratio (LR) threshold
#' K <- 2                                          
#' 
#' # Cohort size used in stage 1
#' coh.size <- 3 
#' 
#' # Vector of true mean efficacies per dose (here mean percent persistence per dose)
#' m <- c(5, 15, 40, 65, 80)   # MUST BE THE SAME LENGTH AS dose.tox                  
#' 
#' # Efficacy(equal) variance per dose
#' v <- rep(0.01, 5) 
#' 
#' # Total sample size (stages 1&2)                            
#' N <- 25                                        
#' 
#' # Stopping rule: if dose 1 is the only safe dose, allocate up to 9 pts.
#' stop.rule <- 9 
#' 
#' sim.trials(numsims = 10, dose, dose.tox, p1 = p_no, p2 = p_yes, K, 
#' coh.size, m, v, N, stop.rule = stop.rule, cohort = 1, samedose = TRUE, nbb = 100)
#' 
#' @export


sim.trials <- function(numsims, dose, dose.tox, p1, p2, K, coh.size, m, v, N, stop.rule = 9, 
                       cohort = 1, samedose = TRUE, nbb = 100){
  
  sim.yk <- sim.dk <- matrix(NA, nrow = numsims, ncol = N)
  sim.doses <- matrix(NA, nrow = numsims, ncol = dose)
  
  ## add progress bar to output
  pb <- utils::txtProgressBar(min = 0,      # Minimum value of the progress bar
                              max = numsims, # Maximum value of the progress bar
                              style = 3,    # Progress bar style (also available style = 1 and style = 2)
                              width = 50,   # Progress bar width. Defaults to getOption("width")
                              char = "=")   # Character used to create the bar
  
  for (i in 1:numsims) {  
    
    fstudy.out <- rand.stg2(dose, dose.tox, p1, p2, K, coh.size, m, v, N, stop.rule, cohort, samedose, nbb)              
    
    # number of doses declared safe, based on stage 1 sample size
    n.safe <- max(fstudy.out$d.final[(fstudy.out$n1 + 1):length(fstudy.out$d.final)], na.rm = TRUE) 
    sim.doses[i,] <- c(rep(1, n.safe), rep(0, dose - n.safe)) # binary vector for dose safety
    
    
    if (length(fstudy.out$Y.final) < N) { # if max sample size was not reached, fill-in with NAs        
      
      sim.yk[i,] <- c(fstudy.out$Y.final, rep(NA, N - length(fstudy.out$Y.final)))
      sim.dk[i,] <- c(fstudy.out$d.final, rep(NA, N - length(fstudy.out$d.final)))
      
    } else {
      
      sim.yk[i,] <- fstudy.out$Y.final
      sim.dk[i,] <- fstudy.out$d.final
    }
    
    # Sets the progress bar to the current state
    utils::setTxtProgressBar(pb, i)
    
    # cat(i,"\n")
  }
  return(list(sim.Y = sim.yk, 
              sim.d = sim.dk, 
              safe.d = sim.doses))
}
