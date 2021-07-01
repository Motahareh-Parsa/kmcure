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
#' library(kmcure)
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
#' fit = kmcure (time, event, survPreds, curePreds)
#'
#' # bootstrap("fit") # Do bootstrap if you need to estimate SE and p-values
#'
#' summary(fit)
#'
#' @export
kmcure <- function(time, event, survPreds, curePreds=NULL,
                   multiOptim_maxit = 10, multiOptim_reltol = 0.001,
                   optim_reltol = 1e-8, optim_maxit = 500,
                   scale = FALSE, silent = FALSE, conditional = FALSE,
                   optim_method = "Nelder-Mead", optim_init = NULL){

  call <- match.call()

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

fit$call = call

class(fit) = "kmcure"

pUncureModel = predict(fit)

fit$pcure_modelMean = mean(1 - pUncureModel)

return(fit)
}

#' @export
print.kmcure <- function(fit){
  cat("This is a kmcure fit object that is resulted from the following call:\n\n")
  call = fit$call
  print(call)
  if(is.null(fit$boot)){
    cat("\nThe Bootstrap function is Not applied to this kmcure fit object yet. \n")
  }else{
    cat("\nThe Bootstrap function with",fit$boot$repeats ,"replications is applied to this kmcure fit object. \n")
  }
}

#' @export
summary.kmcure <- function(fit){
  result = list()

  ## define p-value function
  pvalue <- function(est,se) {
    names = rownames(est)
    est = as.numeric(est)
    se = as.numeric(se)
    twald = est/se
    pwald = 2*pnorm(abs(twald),lower.tail=F)
    output = cbind(est,se,twald,pwald)
    colnames(output)=c("Estimate","Std. Error","t value","Pr(>|t|)")
    rownames(output) = names
    return(output)
  }

  ## calculate output (coef and p-value)
  repeats = fit$boot$repeats
  if(is.null(repeats)) repeats = 0
  if(repeats < 30){
    output = fit$coef
    colnames(output) = "Estimate"
    result$output = output
  }else{
    est = fit$coef
    covmat = cov(t(fit$boot$coef))
    se = diag(covmat)
    output = pvalue(est, se)
    result$output = output
  }

  ## add also interesting info to the result
  result$cure = round(100*fit$pcure, 2)
  result$cens = round(100*fit$pcens, 2)
  result$cure_modelMean = round(100*fit$pcure_modelMean, 2)
  result$cure95ci = rep(NA,2)
  result$cens95ci = rep(NA,2)
  result$cure_modelMean95ci = rep(NA,2)
  if(repeats>=30){
    result$cure95ci[1] = round(100*quantile(fit$boot$pcure,0.025), 2)
    result$cure95ci[2] = round(100*quantile(fit$boot$pcure,0.975), 2)
    result$cens95ci[1] = round(100*quantile(fit$boot$pcens,0.025), 2)
    result$cens95ci[2] = round(100*quantile(fit$boot$pcens,0.975), 2)
    result$cure_modelMean95ci[1] = round(100*quantile(fit$boot$pcure_modelMean,0.025), 2)
    result$cure_modelMean95ci[2] = round(100*quantile(fit$boot$pcure_modelMean,0.975), 2)
  }
  result$call = fit$call
  result$Rboot = repeats
  result$loglik = fit$loglik
  result$AIC = fit$AIC
  result$BIC = fit$BIC

  ## split coefficients to gamma and beta parts based on the coef rownames
  names = rownames(result$output)
  gammaLength = sum(grepl("Gamma", names))
  gammaPart = as.matrix(result$output[1:gammaLength,])
  rownames(gammaPart) = names[1:gammaLength]
  betaPart = as.matrix(result$output[(gammaLength+1):nrow(result$output),])
  rownames(betaPart) = names[(gammaLength+1):nrow(result$output)]
  if(ncol(gammaPart)==1) colnames(gammaPart) = "Estimate"
  if(ncol(betaPart)==1) colnames(betaPart) = "Estimate"
  result$gammaPart = gammaPart
  result$betaPart = betaPart

  ## assign an appropriate class name to the summary of kmcure object
  class(result) <- 'summary.kmcure'
  return(result)
}

#' @export
print.summary.kmcure <- function(result){
  cat("\nCall:\n")
  print(result$call)

  cat("\nDescription of data:\n")
  if(result$Rboot==0){
    cat("- Censored percentage is ", result$cens, "\n", sep = "")
    cat("- Kaplan-Meier estimation of Cured percentage is ", result$cure, "\n", sep = "")
    cat('- The "kmcure" estimation of Cured percentage is ', result$cure_modelMean, "\n", sep = "")
  }else{
    cat("- Censored percentage is ", result$cens, " and its 95% Bootstrap CI is (", result$cens95ci[1], ", ", result$cens95ci[2], ")\n", sep = "")
    cat("- Kaplan-Meier estimation of Cured percentage is ", result$cure, " and its 95% Bootstrap CI is (", result$cure95ci[1], ", ", result$cure95ci[2], ")\n", sep = "")
    cat('- The "kmcure" estimation of Cured percentage is ', result$cure_modelMean, " and its 95% Bootstrap CI is (", result$cure_modelMean95ci[1], ", ", result$cure_modelMean95ci[2], ")\n", sep = "")
  }

  cat("\nCure probability model:\n")
  print(result$gammaPart)

  cat("\nFailure time distribution model:\n")
  print(result$betaPart)

  cat("\nThe log-likelihood is",result$loglik,"and the AIC is", result$AIC, "and the BIC is", result$BIC, "\n")

  if(result$Rboot==0){
    cat("\nThere are no bootstrap replications to calculate SE and pvalue.\n")
    cat('Please employ the kmcure "bootstrap" function to produce replications.\n')
  }
  if(result$Rboot>0 & result$Rboot<30){
    cat("\nThere is", result$Rboot, "bootstrap replications which is not enough to calculate SE and pvalue.\n")
    cat('Please increase the number of bootstrap replications using the bootstrap function.\n')
  }
  if(result$Rboot>30){
    cat("\nCalculation of SE and pvalues are based on", result$Rboot, "bootstrap replications.\n")
    cat('Nnumber of bootstrap replications can be increased using the bootstrap function.\n')
  }

}


#' Predict probability of being an Un-cure observation
#'
#' @param fit a kmcure fit object
#'
#' @param newdata an optional matrix with same columns of the model curePreds
#'
#' @export
predict.kmcure <- function(fit, newdata = NULL){
  if(is.null(newdata)){
    curePreds = fit$data$curePreds
  }else{
    if(is.null(fit$data$curePreds)){
      curePreds = NULL
      warning("newdata is ignored because the model does not have any curePreds")
    }else{
      newdata = as.matrix(newdata)
      if(ncol(newdata==1)) newdata = t(newdata) # account for a numeric newdata
      curePreds = as.matrix(fit$data$curePreds)
      if(ncol(newdata)==ncol(curePreds)){
        curePreds = newdata
      }else{
        curePreds = fit$data$curePreds
        warning("newdata is ignored because its number of columns is different from the model curePreds")
      }
    }
  }
  names = rownames(fit$coef)
  gammaLength = sum(grepl("Gamma", names))
  gammaCoef = as.matrix(fit$coef[1:gammaLength,])
  rownames(gammaCoef) = names[1:gammaLength]
  if(is.null(curePreds)){
    Z = as.matrix(rep(1, length(event)))
  }else{
    if(rownames(gammaCoef)[1]=="Gamma(Intercept)") Z = as.matrix(cbind(1,curePreds)) else Z = as.matrix(curePreds)
  }
  Zgamma = Z %*% gammaCoef
  pUncureModel = exp(Zgamma)/(1+exp(Zgamma))
  return(pUncureModel)
}

