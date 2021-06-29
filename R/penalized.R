#' Calculate Penalized AFT cure model which can be used for predictors selection
#'
#' @param time a survival time to event variable.
#' @param event a survival status variable: 1 for event and 0 for censoring.
#' @param survPreds a matrix of survival predictor variable(s).
#' @param curePreds an optional matrix of curing predictor variable(s).
#' @param R the number of Bootstrap iterations to be used in this procedure
#' @param silent a Boolean value which if set to TRUE it prevent from showing output messages.
#' @return a list showing estimated penalized standard coefficients and non-zoro coefficients
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
#' # penalFit = penalized(time, event, survPreds, curePreds, R = 100)
#'
#' # names(penalFit)
#'
#' # penalFit$isNonZeroCoef$lassoBIC
#'
#' @export
penalized <- function(time, event, survPreds, curePreds=NULL, R = 100, silent = FALSE){

  if(!silent){
    cat("Calculation of Penalized AFT cure model is started at", format(Sys.time(), "%H:%M:%S (%Y-%m-%d)."), "\n")
    cat("Please be patient...\n")
  }

  ## standardize covariates
  stdSurvPreds = scale(survPreds, center = FALSE, scale = TRUE)

  if(is.null(curePreds)){
    stdCurePreds = NULL
  }else{
    stdCurePreds = scale(curePreds, center = FALSE, scale = TRUE)
  }

  ## Fit KMEKDE model to estimate coef
  fit = kmekde (time, event, stdSurvPreds, stdCurePreds, multiOptim_maxit = 1,
                conditional = TRUE, silent = TRUE)
  stdCoef = fit$coef

  ## Bootstrap KMEKDE model to estimate covmat
  BootStdCoef = matrix(nrow=length(stdCoef), ncol = R)
  n = length(time)
  r = 1
  while(r<=R){
    rsample = sample(1:n, n, replace=TRUE)
    rtime = time[rsample]
    revent = event[rsample]
    rStdSurvPreds = stdSurvPreds[rsample,]
    if(is.null(curePreds)){
      rStdCurePreds = NULL
    }else{
      rStdCurePreds = stdCurePreds[rsample,]
    }
    rfit = kmcure (rtime, revent, rStdSurvPreds, rStdCurePreds,  multiOptim_maxit = 1,
                   conditional = FALSE, optim_init = stdCoef, silent = TRUE)
    if(rfit$exitcode==0){
      BootStdCoef[, r] = rfit$coef
      r = r + 1
    }
  } # end while loop of r

  stdCovMat = cov(t(BootStdCoef))

  stdCoefPenLassoBIC = suppressWarnings(LSAkmcure(stdCoef, stdCovMat, n,
                                                  type = "lasso" , criteria = "BIC"))
  stdCoefPenLassoAIC = suppressWarnings(LSAkmcure(stdCoef, stdCovMat, n,
                                                  type = "lasso" , criteria = "AIC"))
  stdCoefPenLarBIC = suppressWarnings(LSAkmcure(stdCoef, stdCovMat, n,
                                                  type = "lar" , criteria = "BIC"))
  stdCoefPenLarAIC = suppressWarnings(LSAkmcure(stdCoef, stdCovMat, n,
                                                  type = "lar" , criteria = "AIC"))

  nonZeroStdCoefPenLassoBIC = (stdCoefPenLassoBIC != 0)
  nonZeroStdCoefPenLassoAIC = (stdCoefPenLassoAIC != 0)
  nonZeroStdCoefPenLarBIC = (stdCoefPenLarBIC != 0)
  nonZeroStdCoefPenLarAIC = (stdCoefPenLarAIC != 0)

  outlist = list()

  outlist$stdCoefPenalized$lassoBIC = stdCoefPenLassoBIC
  outlist$stdCoefPenalized$lassoAIC = stdCoefPenLassoAIC
  outlist$stdCoefPenalized$larBIC = stdCoefPenLarBIC
  outlist$stdCoefPenalized$larAIC = stdCoefPenLarAIC

  outlist$isNonZeroCoef$lassoBIC = nonZeroStdCoefPenLassoBIC
  outlist$isNonZeroCoef$lassoAIC = nonZeroStdCoefPenLassoAIC
  outlist$isNonZeroCoef$larBIC = nonZeroStdCoefPenLarBIC
  outlist$isNonZeroCoef$larAIC = nonZeroStdCoefPenLarAIC

  if(!silent) cat("Calculation of Penalized AFT cure model is completed at", format(Sys.time(), "%H:%M:%S (%Y-%m-%d)."), "\n")

  return(outlist)
}




LSAkmcure = function (coef, covmat, n, type , criteria) {

  type = tolower(type[1])
  criteria = toupper(criteria[1])

  ## split coefficients to gamma and beta parts based on the coef rownames
  names = rownames(coef)
  bintercept <- grepl("(Intercept)", names[1])
  betaintercept <- any(grepl("(Intercept)", names[-1]))
  bLength = sum(grepl("Gamma", names))

  ### split coef and covmat to the b and beta parts and then do shrinkage using the lars.lsa function

  # extract b part coef and it's coresponding covariance matrix
  bcoef = coef[1:bLength]
  bcov = covmat[1:bLength, 1:bLength]
  #dim(bcov)

  # extract beta part coef and it's coresponding covariance matrix
  betacoef = coef[(1+bLength):length(coef)]
  betacov = covmat[(1+bLength):length(coef), (1+bLength):length(coef)]

  bSI <- solve(bcov)

  betaSI <- solve(betacov)

  # ---- shrinkage b part using RobMixReg::lars.lsa
  bPart_lsa <- RobMixReg::lars.lsa(bSI, bcoef, bintercept, n, type)

  bPart_ic = bPart_lsa$BIC
  if(criteria=="AIC") bPart_ic = bPart_lsa$AIC

  bPart_BestIndex = which.min(bPart_ic)
  b0Pen = bPart_lsa$beta0[bPart_BestIndex]
  bcoefPen = bPart_lsa$beta[bPart_BestIndex,]

  # ---- shrinkage beta part using RobMixReg::lars.lsa
  betaPart_lsa <- RobMixReg::lars.lsa(betaSI, betacoef, betaintercept, n, type)

  betaPart_ic = betaPart_lsa$BIC
  if(criteria=="AIC") betaPart_ic = betaPart_lsa$AIC

  betaPart_BestIndex = which.min(betaPart_ic)
  betacoefPen = betaPart_lsa$beta[betaPart_BestIndex,]

  thetaPen = c(b0Pen, bcoefPen, betacoefPen)

  thetaPen = as.matrix(thetaPen)
  rownames(thetaPen) = rownames(coef)
  colnames(thetaPen) = "stdCoefPenalized"

  return(thetaPen)
}
