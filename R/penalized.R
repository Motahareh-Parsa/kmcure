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
    cat("Running of the kmcure 'penalized' function is started at", format(Sys.time(), "%H:%M:%S (%Y-%m-%d)."), "\n")
    cat("Please be patient...\n")
  }

  timeS = Sys.time() # start recording time

  ## make sure that time and event are not dataframe to prevent indexing problem in bootstrapping
  time = as.matrix(time)
  event = as.matrix(event)

  ## standardize covariates
  stdSurvPreds = scale(survPreds, center = FALSE, scale = TRUE)
  scaleSurvPreds = attr(stdSurvPreds, "scaled:scale")

  if(is.null(curePreds)){
    stdCurePreds = NULL
    scaleCurePreds = NULL
  }else{
    stdCurePreds = scale(curePreds, center = FALSE, scale = TRUE)
    scaleCurePreds = attr(stdCurePreds, "scaled:scale")
  }

  scaleFull = c(1, scaleCurePreds, scaleSurvPreds)

  ## Fit KMEKDE model to estimate coef
  if(!silent) cat("The primary fitting of the model is started at", format(Sys.time(), "%H:%M:%S ..."), "\n")
  fit = kmekde (time, event, stdSurvPreds, stdCurePreds, multiOptim_maxit = 1,
                conditional = TRUE, silent = TRUE)
  if(!silent) cat("The primary fitting of the model is finished at", format(Sys.time(), "%H:%M:%S ..."), "\n")
  if(fit$exitcode>=2){
    stop("An error occurred in using kmekde in the penalized function! The kmekde exitcode is", fit$exitcode)
  }else{
    stdCoef = fit$coef

    ## Bootstrap KMEKDE model to estimate covmat
    timeS_boot = Sys.time()
    if(!silent) cat("The Bootstrap fittings of the model is started at", format(Sys.time(), "%H:%M:%S ..."), "\n")
    BootStdCoef = matrix(nrow=length(stdCoef), ncol = R)
    n = length(time)
    r = 1
    while(r<=R){
      if(!silent) if(r%%10==0) cat("Bootstrap #", r, "of", R, "is completed...\n")
      rsample = sample(1:n, n, replace=TRUE)
      rtime = time[rsample]
      revent = event[rsample]
      rStdSurvPreds = stdSurvPreds[rsample,]
      if(is.null(curePreds)){
        rStdCurePreds = NULL
      }else{
        rStdCurePreds = stdCurePreds[rsample,]
      }
      rfit = kmekde (rtime, revent, rStdSurvPreds, rStdCurePreds,  multiOptim_maxit = 1,
                     conditional = FALSE, optim_init = stdCoef, silent = TRUE)
      if(rfit$exitcode==0){
        BootStdCoef[, r] = rfit$coef
        r = r + 1
      }
    } # end while loop of r
    timeE_boot = Sys.time()

    stdCovMat = cov(t(BootStdCoef))

    # do LSAkmcure for lasso type shrinkage
    timeS_lsa_lasso = Sys.time()
    stdCoefPenLasso = suppressWarnings(LSAkmcure(stdCoef, stdCovMat, n, type = "lasso"))
    timeE_lsa_lasso = Sys.time()

    # do LSAkmcure for lar type shrinkage
    timeS_lsa_lar = Sys.time()
    stdCoefPenLar = suppressWarnings(LSAkmcure(stdCoef, stdCovMat, n, type = "lar"))
    timeE_lsa_lar = Sys.time()

    # prepare outlist for the penalized function

    outlist = list()

    outlist$isNonZeroCoef$lassoBIC = (stdCoefPenLasso$theta.bic != 0)
    outlist$isNonZeroCoef$lassoAIC = (stdCoefPenLasso$theta.aic != 0)
    outlist$isNonZeroCoef$larBIC = (stdCoefPenLar$theta.bic != 0)
    outlist$isNonZeroCoef$larAIC = (stdCoefPenLar$theta.aic != 0)

    outlist$coefPenalized$lassoBIC = stdCoefPenLasso$theta.bic / scaleFull
    outlist$coefPenalized$lassoAIC = stdCoefPenLasso$theta.aic / scaleFull
    outlist$coefPenalized$larBIC = stdCoefPenLar$theta.bic / scaleFull
    outlist$coefPenalized$larAIC = stdCoefPenLar$theta.aic / scaleFull

    outlist$stdCoefPenalized$lassoBIC = stdCoefPenLasso$theta.bic
    outlist$stdCoefPenalized$lassoAIC = stdCoefPenLasso$theta.aic
    outlist$stdCoefPenalized$larBIC = stdCoefPenLar$theta.bic
    outlist$stdCoefPenalized$larAIC = stdCoefPenLar$theta.aic

    outlist$ordinaryFit$stdCoef = stdCoef
    outlist$ordinaryFit$stdCovMat = stdCovMat

    outlist$scaleFull = scaleFull

    outlist$timeD_fit = fit$timeD
    outlist$timeD_bootstrap = difftime(timeE_boot, timeS_boot, units = "mins")
    outlist$timeD_lsa_lasso = difftime(timeE_lsa_lasso, timeS_lsa_lasso, units = "mins")
    outlist$timeD_lsa_lar = difftime(timeE_lsa_lar, timeS_lsa_lar, units = "mins")

    timeE = Sys.time()
    timeD = difftime(timeE, timeS, units = "mins")
    outlist$timeD = timeD

    outlist$BIC_lasso = stdCoefPenLasso$best.bic
    outlist$BIC_lar = stdCoefPenLar$best.bic
    outlist$AIC_lasso = stdCoefPenLasso$best.aic
    outlist$AIC_lar = stdCoefPenLar$best.aic

    if(!silent) cat("Running of the kmcure 'penalized' function is finished at", format(Sys.time(), "%H:%M:%S (%Y-%m-%d)."), "\n")

    return(outlist)

  } # end else fit$exitcode
} # end penalized function



LSAkmcure = function (coef, covmat, n, type) {

  type = tolower(type[1])

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

  bExtract = extract_larslsa(l.fit = bPart_lsa, beta.ols=bcoef, intercept=bintercept, n=n)

  # ---- shrinkage beta part using RobMixReg::lars.lsa
  betaPart_lsa <- RobMixReg::lars.lsa(betaSI, betacoef, betaintercept, n, type)

  betaExtract = extract_larslsa(l.fit = betaPart_lsa, beta.ols=betacoef, intercept=betaintercept, n=n)

  # prepare objects for output
  theta.ols = as.matrix(c(bExtract$beta.ols, betaExtract$beta.ols))
  theta.bic = as.matrix(c(bExtract$beta.bic, betaExtract$beta.bic))
  theta.aic = as.matrix(c(bExtract$beta.aic, betaExtract$beta.aic))
  best.bic = c(bExtract$bestBIC, betaExtract$bestBIC)
  best.aic = c(bExtract$bestAIC, betaExtract$bestAIC)
  rownames(theta.ols) = rownames(coef); colnames(theta.ols) = "theta.ols"
  rownames(theta.bic) = rownames(coef); colnames(theta.bic) = "theta.bic"
  rownames(theta.aic) = rownames(coef); colnames(theta.aic) = "theta.aic"
  names(best.bic)=c("GammaBestBIC","BetaBestBIC")
  names(best.aic)=c("GammaBestAIC","BetaBestAIC")
  LSAkmcure_outlist = list(theta.ols=theta.ols, theta.bic=theta.bic, theta.aic=theta.aic,
                           best.bic=best.bic, best.aic=best.aic)
  return(LSAkmcure_outlist)
}



extract_larslsa = function(l.fit, beta.ols, intercept, n){
  t1 <- sort(l.fit$BIC, ind=T)

  t2 <- sort(l.fit$AIC, ind=T)

  beta <- l.fit$beta

  if(intercept) {

    beta0 <- l.fit$beta0+beta.ols[1]

    beta.bic <- c(beta0[t1$ix[1]],beta[t1$ix[1],])

    beta.aic <- c(beta0[t2$ix[1]],beta[t2$ix[1],])

  } else {

    beta0 <- l.fit$beta0

    beta.bic <- beta[t1$ix[1],]

    beta.aic <- beta[t2$ix[1],]

  }
  bestBIC = l.fit$BIC[t1$ix[1]]
  bestAIC = l.fit$AIC[t2$ix[1]]
  obj <- list(beta.ols=beta.ols,
              beta.bic=beta.bic, beta.aic = beta.aic,
              bestBIC = bestBIC, bestAIC = bestAIC)
  return(obj)
}

