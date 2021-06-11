#' Fits AFT Semiparametric Mixture Cure Model using the KME-KDE method
#'
#' Fits AFT (Accelerated Failure Time) Semiparametric Mixture Cure Model using the KME-KDE (Kaplan-Meier Estimation and Kernel Density Estimator) method.
#'
#' @param formula is the formula of survival sub-model which define time and event variables using Surv function on the left of a tilde symbol and it define list of survival predictors on the right of the tilde symbol.
#' @param cureform is the formula of cure sub-model which define predictors of cure sub-model on the right of a tilde symbol, if omitted program use just an intercept in the cure sub-model
#' @param multiOptim_maxit is the maximum of allowed multiple optimization. If the convergence of "optim" does not meet, the program does multiple optimization.
#' @param multiOptim_reltol is the relative tolerance in continuing multiple optimization procedure
#' @param optim_reltol is the relative tolerance in continuing of each optimization
#' @param optim_maxit is the maximum of allowed iterations in each optimization
#' @param silent a Boolean value which if set to TRUE it prevent from showing output messages
#' @param data is the input data frame which its columns' names can be called in the formula and the cureform
#'
#' @return is a kmcure object
#'
#' @examples
#' data(hfp)
#' fit = kmcure(formula = Surv(TIME, Event) ~
#' Gender+Smoking+Diabetes+BP+Anaemia+Age+EF+Sodium+Creatinine+Pletelets+CPK,
#' cureform = ~Gender+Smoking+Diabetes+BP+Anaemia+Age+EF+Sodium+Creatinine+Pletelets+CPK,
#' data = hfp)
#' names(fit)
#' View(fit$mat_coef[, c(1, ncol(fit$mat_coef) )])
#' @export
kmcure <- function(formula, cureform, data,
                   multiOptim_maxit = 100, multiOptim_reltol = 0.001,
                   optim_reltol = 1e-8, optim_maxit = 500, silent = FALSE){

  tryCatch({mf <- model.frame(formula,data)},
           error = function(e) {
             stop("The Survival formula is not valid! This maybe duo to misspelled variabe names.")})

  tryCatch({cvars <- all.vars(cureform); Z <- as.matrix(data[,cvars]); colnames(Z) <- c(cvars)},
           error = function(e) {
             stop("The Cure formula is not valid! This maybe duo to misspelled variabe names.")})

  Y <- model.extract(mf,"response")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
  X <- model.matrix(attr(mf,"terms"), mf)

  time <- Y[,1]
  event <- Y[,2]

  betaNames = colnames(X)
  if(  betaNames[1]=="(Intercept)" ){
    X = as.matrix(X[,-1])
    betaNames = betaNames[-1]
    colnames(X) = betaNames
  }

  gammaNames = colnames(Z)
  if( gammaNames[1]=="(Intercept)" ){
    Z = as.matrix(Z[,-1])
    gammaNames = gammaNames[-1]
    colnames(Z) = gammaNames
  }

if(ncol(Z)==0) Z = NULL

# return(list(time=time, event=event, survPreds=X, curePreds=Z))
result = kmekde(time = time, event = event, survPreds = X, curePreds = Z,
             multiOptim_maxit = multiOptim_maxit, multiOptim_reltol = multiOptim_reltol,
             optim_reltol = optim_reltol, optim_maxit = optim_maxit, silent = silent)

result$data$time = time
result$data$event = event
result$data$survPreds = X
result$data$curePreds = Z

class(result) = c("kmcure","list")

return(result)
}
