#' @title Generate plots for estimated percent allocation and response per dose.
#' 
#' @description Generate plots for estimated percent allocation and response per dose.
#' 
#' @return Error plots of estimated (1) percent allocation per dose, and 
#' (2) estimated response per dose.
#'          
#' @param sims output from sim.trials
#' 
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
#' # Vector of true mean efficacies per dose (here mean T-cell persistence per dose (%))
#' m <- c(5, 15, 40, 65, 80)   # MUST BE THE SAME LENGTH AS dose.tox                  
#' 
#' # Efficacy (equal) variance per dose
#' v <- rep(0.01, 5) 
#' 
#' # Total sample size (stages 1&2)                            
#' N <- 25                                        
#' 
#' # Stopping rule
#' stop.rule <- 9   
#' 
#' 
#' numsims = 100
#' 
#' set.seed(1)
#' simulations = sim.trials(numsims = numsims, dose, dose.tox, p1 = p_no, p2 = p_yes, 
#'                          K, coh.size, m, v, N, stop.rule = stop.rule, cohort = 1, 
#'                          samedose = TRUE, nbb = 100)
#'                          
#' # sim.plot(simulations)
#' 
#' @export
#' 
#' 

sim.plot <- function(sims) {
  tables = sim.summary(sims)
  
  graphics::par(mfrow = c(1, 2))
  
  ## Percent allocation
  t = tables$pct.treated
  dd = nrow(t)
  graphics::plot(x = t[,1], y = t[,3],
                 ylim = c(0, 100),
                 ylab = "Percent allocation",
                 xlab = "Dose",
                 cex = 1.3,
                 pch = 19,
                 main = "Percent allocation")
  for (i in 1:dd) {
    graphics::lines(x = rep(i, 2), y = t[i, c(2, 4)], lwd = 3)
  }
  
  ## Estimated response
  t = tables$efficacy
  dd = nrow(t)
  graphics::plot(x = t[,1], y = t[,3],
                 ylim = c(0, 100),
                 ylab = "Estimated response",
                 xlab = "Dose",
                 cex = 1.3,
                 pch = 19,
                 main = "Estimated response")
  for (i in 1:dd) {
    graphics::lines(x = rep(i, 2), y = t[i, c(2, 4)], lwd = 3)
  }
  
  # reset graphics parameters
  grDevices::graphics.off()
  
}
