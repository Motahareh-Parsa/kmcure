#' Title
#'
#' @param formula is the function input
#' @param cureform is the function input
#' @param data is the function input
#' @return is a kmcure object
#'
#' @examples
#' data(hfp)
#' kmcure(11,22,33)
#'
#' @export
kmcure <- function(formula, cureform, data){

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

# return(list(time=time, event=event, survPreds=X, curePreds=Z))
fit = kmekde(time = time, event = event, survPreds = X, curePreds = Z)
return(fit)
}
