#' Fits AFT Semiparametric Mixture Cure Model using the KME-KDE method
#'
#' Fits AFT (Accelerated Failure Time) Semiparametric Mixture Cure Model using the KME-KDE (Kaplan-Meier Estimation and Kernel Density Estimator) method.
#'
#' @param time a survival time to event variable.
#' @param event a survival status variable: 1 for event and 0 for censoring.
#' @param survPreds a matrix of survival predictor variable(s).
#' @param curePreds an optional matrix of curing predictor variable(s).
#' @param multiOptim_maxit a number showing the maximum of allowed multiple optimization. The program uses multiple optimization if the convergence of "optim" does not meet.
#' @param multiOptim_reltol a number showing the relative tolerance in continuing multiple optimization procedure.
#' @param optim_reltol a number showing the relative tolerance in continuing of each optimization.
#' @param optim_maxit a number showing the maximum of allowed iterations in each optimization.
#' @param scale a Boolean value that if set to TRUE the program automatically transform predictors using scale function and then back transform the estimated coefficients to the original scale.
#' @param silent a Boolean value which if set to TRUE it prevent from showing output messages.
#' @param conditional a Boolean value which, if set to TRUE it uses an iterative procedure that estimate parameters survival and cure sub-models conditionally on the last estimation of the other one.
#' @param optim_method a string value showing the method of optimization: "Nelder-Mead" and "SANN" are supported.
#' @param optim_init an optional numeric vector of initial values. For example, it could be the estimated coefficients of a previous fit to be used in continuing of optimization.
#'
#' @return is a kmcure object
#'
#' @examples
#'
#' data(hfp)
#'
#' time = hfp$Time
#'
#' event = hfp$Event
#'
#' survPreds = hfp[, c(3:15)]
#' names(survPreds)
#'
#' curePreds = hfp[, c(3:15)]
#' names(curePreds)
#'
#' fit1 = kmcure (time, event, survPreds, curePreds, multiOptim_maxit = 10)
#'
#' names(fit1)
#'
#' if(fit1$exitcode==0){
#' print(fit1$loglik)
#' print(fit1$timeD)
#' print(fit1$coef)
#' }
#'
#' @export
kmcure <- function(time, event, survPreds, curePreds=NULL,
                   multiOptim_maxit = 10, multiOptim_reltol = 0.001,
                   optim_reltol = 1e-8, optim_maxit = 500,
                   scale = FALSE, silent = FALSE, conditional = FALSE,
                   optim_method = "Nelder-Mead", optim_init = NULL){

  settings = list(multiOptim_maxit=multiOptim_maxit, multiOptim_reltol=multiOptim_reltol,
                  optim_reltol=optim_reltol, optim_maxit=optim_maxit,
                  scale=scale, silent=silent, conditional = conditional,
                  optim_method=optim_method, optim_init=optim_init)

  X = scale(survPreds, center = FALSE, scale = TRUE)
  scaleX = attr(X, "scaled:scale")

  if(is.null(curePreds)){
    Z = NULL
    scaleZ = NULL
  }else{
    Z = scale(curePreds, center = FALSE, scale = TRUE)
    scaleZ = attr(Z, "scaled:scale")
  }

  scaleFull = c(1, scaleZ, scaleX)

  if(scale==FALSE){
    X = survPreds
    Z = curePreds
  }

fit = kmekde(time = time, event = event, survPreds = X, curePreds = Z,
             multiOptim_maxit = multiOptim_maxit, multiOptim_reltol = multiOptim_reltol,
             optim_reltol = optim_reltol, optim_maxit = optim_maxit, silent = silent,
             conditional = conditional, optim_method = optim_method, optim_init = optim_init)

if (scale==TRUE){
  fit$stdmat_coef = fit$mat_coef
  stdmat_coef =  fit$mat_coef
  backstdmat_coef = stdmat_coef / scaleFull
  fit$mat_coef = backstdmat_coef
}else{
  mat_coef = fit$mat_coef
  stdmat_coef = mat_coef * scaleFull
  fit$stdmat_coef = stdmat_coef
}

indexMaxLL = which.max(fit$vecLL)
fit$coef = as.matrix(fit$mat_coef [,indexMaxLL])
fit$stdcoef = as.matrix(fit$stdmat_coef [,indexMaxLL])

fit$scale$survPreds = scaleX
fit$scale$curePreds = scaleZ
fit$scale$fullScale = scaleFull

fit$data$time = time
fit$data$event = event
fit$data$survPreds = survPreds
fit$data$curePreds = curePreds

fit$settings = settings

class(fit) = c("kmcure","list")

return(fit)
}
