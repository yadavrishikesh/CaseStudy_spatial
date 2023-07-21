
## The posterior predictive distributions F_i=N(m_i, s_i^2) include the observation noise
## variance, so are appropriate to evaluate on "new" observations y.
##
## Mean scores: S({F_i}, {y_i}) := 1/N \sum_{i=1}^N S(F_i, y_i)
##
## All the scores are negatively oriented, i.e. small values are
## better, and all are proper, i.e. the expected score under y~G
## cannot be made smaller than when F=G (the predictive distribution
## is the same as the truth).
## For more scoring details, see Gneiting & Raftery, JASA, 2007.
##
## mspe or MSE: mean squared error, S(F,y) = (y-m)^2
## mae: mean absolute error, S(F,y) = |y-m|
## normlized_PI or IGN: ignorance score, S(F,y) = (y-m)^2 / s^2 / 2 + log(s)
## crps or CRPS: contiunuous ranked probability score, S(F,y) = CRPS(N(m,s^2), y)
##  Interval score: interval score, S(F,y) = IntervalScore([m-q*s, m+q*s], y)
##  mpiw: mean prediction interval width
## picp or Coverage: prediction interval coverage (not a score) = I{y in [m-q*s, m+q*s]}
## q is the 1-alpha/2 quantile of N(0,1), and alpha=0.05
##
## Remark: (y-m)^2 + s^2 is not a proper score, and should not be
## used. It's frighteningly common to see it used in the wild. Resist
## the temptation. (See G&R above for more details.)
#' Compute the summary statistics of the prediction accuracy  
#'
#' @param pred_val a vectors of predicted values
#' @param true_val  a vector of true values, same length as the pred_val
#' @param std_pred a vector of standard deviation of the predicted values, same length as pred_val
#' @param sumry_type  types of summary statistics to be calculates. This includes, (1) mspe: r mean squared prediction scores, (2) v: mean absolute erorr, (3) mpiw:  mean prediction Interval score
#'  (4) picp: ean prediction coverage, (5) normlized_PI: normalized prediction residuals (6) crps: crps-score, (7) norm_pr: normalized prediction residuals
#' @param alpha  by default set to 0.05 for a coverage probability of 1-alpha
#'
#' @return the calculated summary statistics 
eval_scores<- function(pred_val, 
                       true_val, 
                       std_pred,
                       sumry_type="mspe",
                       alpha=0.05
                       ){
  #browser()
  if(sumry_type=="mspe"){ ## Compute the sqerr-score: or mean squared prediction scores
    z <- as.numeric((true_val - pred_val))
    scores <- mean(z^2, na.rm=TRUE)
  } else if(sumry_type=="mae"){   ## Compute the absolute error score:
    z <- as.numeric((true_val - pred_val))
    scores <- mean(abs(z), na.rm=TRUE)
  } else if(sumry_type=="IS"){  ## Compute the mean prediction Interval score:
    hw <- -qnorm(alpha/2) * std_pred
    scores <- mean(2 * hw + (2/alpha) * (((pred_val - hw) - true_val) * (true_val < pred_val - hw) +
                                      (true_val - (pred_val + hw)) * (true_val > pred_val + hw)), na.rm=TRUE)
  } else if(sumry_type=="mpiw"){  ## Compute the mean prediction Interval score:
      hw <- -qnorm(alpha/2) * std_pred
      scores <- mean(2 * hw, na.rm=TRUE)
  } else if(sumry_type=="picp"){ ## Compute the mean prediction coverage:
    hw <- -qnorm(alpha/2) * std_pred
    scores <- mean((pred_val - hw <= true_val) & (true_val <= pred_val + hw), na.rm=TRUE)
  } else if (sumry_type=="normlized_PI") { ## Compute the Compute normalised prediction residuals
    z <- as.numeric((true_val - pred_val) / std_pred)
    scores <- mean(z^2 / 2 + log(std_pred), na.rm=TRUE)
  } else if (sumry_type=="crps") {  ## Compute the crps-score:
    z <- as.numeric((true_val - pred_val) / std_pred)
    scores <-  mean(std_pred * (z *(2 * pnorm(z, 0, 1) - 1) +
                        2 * dnorm(z, 0, 1) - 1/sqrt(pi)), na.rm=TRUE)
  } else { ## Compute normalised prediction residuals
    z <- as.numeric((true_val - pred_val) / std_pred)
    scores <- mean(pnorm(z, 0, 1), na.rm=TRUE)
  }
  return(scores)
}

