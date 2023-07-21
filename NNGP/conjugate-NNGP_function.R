#rm(list=ls())
#' Estimation of Matern parameters and prediction at unknown spatial locations using  NNGP approach
#'
#' @param data_train the observations vectors of length ns_t in the training set, where ns_t is the number if spatial locations in training datasets
#' @param S the location of spatial observations in the training datasets
#' @param S.pred the location of spatial observations in the testing datasets
#' @param init.par initial parameter values for the reparametrized Matern parameters
#' @param discrete.level 
#'
#' @return
#' @export
#'
#' @examples
conjNNGP <- function(data_train, 
                     S,
                     S.pred, 
                     X.0,
                     init.par,
                     n.neighbors=50,
                     k.fold=5,
                     n.omp.threads=1,
                     score.rule="rmspe",
                     n.samples=200,
                     discrete.level = 3
                     ){
  alpha.init<-c(seq(from=0.01, to=1, length.out=discrete.level),
                (init.par[4]/init.par[1])+0.00001)
  phi.init<-c(seq(from=1, to=100, length.out=discrete.level),1/init.par[2])
  nu.init<-c(seq(from=0.1, to=5, length.out=discrete.level),init.par[3])
  theta.alpha <- as.matrix(expand.grid(alpha.init,phi.init, nu.init))
  colnames(theta.alpha) <- c("alpha","phi","nu")
  
  ord <- order(S[,1]+S[,2])
  z<- data_train$z
  x<- data_train$x
  y<- data_train$y
  sim.c <- spConjNNGP(formula=z~1+x+y, 
                      data=data_train,  
                      coords=S,
                      cov.model="matern", 
                      sigma.sq.IG=c(2,0.5*var(z)),
                      n.neighbors=n.neighbors, 
                      ord=ord,
                      theta.alpha=theta.alpha,
                      k.fold = k.fold, 
                      score.rule = "rmspe",
                      fit.rep=TRUE, 
                      n.samples=n.samples,
                      n.omp.threads=n.omp.threads)
  
  matern.pars <- c(sim.c$sigma.sq.hat, sim.c$theta.alpha[2],  
                   sim.c$theta.alpha[3], sim.c$theta.alpha[2]*sim.c$sigma.sq.hat)
  comp.time.estimate <- sim.c$run.time[3]
  
  preds <- predict(sim.c, 
                   X.0=X.0,
                   coords.0=S.pred, 
                   n.omp.threads=n.omp.threads)
  comp.time.predict <- preds$run.time[3]
  
  out <- list(matern.pars = matern.pars,
              comp.time.estimate = comp.time.estimate,
              preds = preds,
              comp.time.predict = comp.time.predict)
  out}
