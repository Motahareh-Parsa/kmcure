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
    stdCoefPenLasso = LSAkmcure(stdCoef, stdCovMat, n=sum(event), type = "lasso")
    timeE_lsa_lasso = Sys.time()

    # do LSAkmcure for lar type shrinkage
    timeS_lsa_lar = Sys.time()
    stdCoefPenLar = LSAkmcure(stdCoef, stdCovMat, n=sum(event), type = "lar")
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
    outlist$AIC_lasso = stdCoefPenLasso$best.aic
    outlist$BIC_lar = stdCoefPenLar$best.bic
    outlist$AIC_lar = stdCoefPenLar$best.aic

    outlist$lambda_BIC_lasso = stdCoefPenLasso$lambda.bic
    outlist$lambda_AIC_lasso = stdCoefPenLasso$lambda.aic
    outlist$lambda_BIC_lar = stdCoefPenLar$lambda.bic
    outlist$lambda_AIC_lar = stdCoefPenLar$lambda.aic

    if(!silent) cat("Running of the kmcure 'penalized' function is finished at", format(Sys.time(), "%H:%M:%S (%Y-%m-%d)."), "\n")

    return(outlist)

  } # end else fit$exitcode
} # end penalized function
