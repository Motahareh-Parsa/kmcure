#' Do bootstrap replications for a "kmcure" fitted object
#'
#' After estimating model parameters using the "kmcure" function we can use the "bootstrap" function to make it possible to estimate standard errors and p-values of the model parameters.
#'
#' @param fitObjName is the name of a kmcure fit object
#' @param R is the number of Bootstrap repeats to be runs
#' @param useFitEstAsBootInit a Boolean value which if set to TRUE it use the Fit Estimated coefficients as Initial values in Bootstrap Estimations
#' @param considerPreviousBoots a Boolean value which if set to FALSE it ignore possible previously estimated bootstrap replications
#' @param silent a Boolean value which if set to TRUE it prevent from showing output messages.
#' @param multiOptim_maxit (NULL means to employ the same setting from the kmcure fit object) a number showing the maximum of allowed multiple optimization. The program uses multiple optimization if the convergence of "optim" does not meet.
#' @param multiOptim_reltol (NULL means to employ the same setting from the kmcure fit object) a number showing the relative tolerance in continuing multiple optimization procedure.
#' @param optim_reltol (NULL means to employ the same setting from the kmcure fit object) a number showing the relative tolerance in continuing of each optimization.
#' @param optim_maxit (NULL means to employ the same setting from the kmcure fit object) a number showing the maximum of allowed iterations in each optimization.
#'
#' @return No direct value is returned instead add a boot section in the inputted kmcure fit object
#'
#' @examples
#' \dontrun{
#'
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
#' fit = kmcure (time, event, survPreds, curePreds) # Fit AFT cure model to data
#'
#' bootstrap("fit") # Do bootstrap by providing the name of the fitted kmcure object
#'
#' summary(fit)
#'
#' }
#' @export
bootstrap <- function(fitObjName,
                      R = 50,
                      useFitEstAsBootInit = TRUE,
                      considerPreviousBoots = TRUE,
                      silent = FALSE,
                      multiOptim_maxit = NULL,
                      multiOptim_reltol = NULL,
                      optim_reltol = NULL,
                      optim_maxit = NULL){

  if(class(fitObjName)!="character") stop('The fitObjName argument is Not the "name" of a kmcure fit object: It must be a string or character object!')
  if(exists(fitObjName)==FALSE) stop('An object with the specified name of fitObjName does Not exist in the R environment! Please, first of all, check its spelling. It must be the name of a valid kmcure object. Please see "help(bootstrap)" for more information.')
  if(class(eval(str2lang(fitObjName)))!="kmcure") stop('The fitObjName argument is Not the name of a valid "kmcure" fit object:  It must refer to a "kmcure" object!')


  fitObj = eval(str2lang(fitObjName))

  if(is.null(multiOptim_maxit)) multiOptim_maxit = fitObj$settings$multiOptim_maxit
  if(is.null(multiOptim_reltol)) multiOptim_reltol = fitObj$settings$multiOptim_reltol
  if(is.null(optim_reltol)) optim_reltol = fitObj$settings$optim_reltol
  if(is.null(optim_maxit)) optim_maxit = fitObj$settings$optim_maxit
  if(useFitEstAsBootInit){
    optim_init = fitObj$coef
  }else{
    optim_init = NULL
  }

  if(considerPreviousBoots==FALSE){
    fitObj$boot = NULL
    if(!silent) cat('The possible previously completed bootstrap results is dropped because "considerPreviousBoots" was set to FALSE!\n')
  }

  r = ncol(fitObj$boot$coef) # number of previously completed bootstrap repeats; start from zero
  if (is.null(r)) r = 0
  if(r>0){ # check whether the bootstrap settings are the same as the previously completed bootstrap
    useFitEstAsBootInit0 = fitObj$boot$settings$useFitEstAsBootInit
    multiOptim_maxit0 = fitObj$boot$settings$multiOptim_maxit
    multiOptim_reltol0 = fitObj$boot$settings$multiOptim_reltol
    optim_reltol0 = fitObj$boot$settings$optim_reltol
    optim_maxit0 = fitObj$boot$settings$optim_maxit
  if(useFitEstAsBootInit!=useFitEstAsBootInit0) stop("The value of useFitEstAsBootInit is ", useFitEstAsBootInit, " so its value is different from the previously considered bootstrap which was ", useFitEstAsBootInit0, " Re-run the bootstrap function by using the same value or by setting the considerPreviousBoots as FALSE!")
  if(multiOptim_maxit!=multiOptim_maxit0) stop("multiOptim_maxit is ", multiOptim_maxit, " so its value is different from the previously considered bootstrap which was ", multiOptim_maxit0, " Re-run the bootstrap function by using the same value or by setting the considerPreviousBoots as FALSE!")
  if(multiOptim_reltol!=multiOptim_reltol0) stop("multiOptim_reltol is ", multiOptim_reltol, " so its value is different from the previously considered bootstrap which was ", multiOptim_reltol0, " Re-run the bootstrap function by using the same value or by setting the considerPreviousBoots as FALSE!")
  if(optim_reltol!=optim_reltol0) stop("optim_reltol is ", optim_reltol, " so its value is different from the previously considered bootstrap which was ", optim_reltol0, " Re-run the bootstrap function by using the same value or by setting the considerPreviousBoots as FALSE!")
  if(optim_maxit!=optim_maxit0) stop("optim_maxit is ", optim_maxit, " so its value is different from the previously considered bootstrap which was ", optim_maxit0, " Re-run the bootstrap function by using the same value or by setting the considerPreviousBoots as FALSE!")
  } # end if conditioning on r

  if(!silent) cat('The run of kmcure "bootstrap" function is started at', format(Sys.time(), "%H:%M:%S (%Y-%m-%d)..."), "\n")
  if(!silent) cat("The requested bootstrap repeats is", R, "where", r, "repeats is previously completed!\n")
  if(!silent) cat("Please be patient to complete the bootstrap replications...\n")

  while(r < R){

    n = length(fitObj$data$time)
    rsample = sample(1:n, n, replace=TRUE)
    rtime = fitObj$data$time[rsample]
    revent = fitObj$data$event[rsample]
    if(is.null(fitObj$data$survPreds)){
      rsurvPreds = NULL
    }else{
      rsurvPreds = fitObj$data$survPreds[rsample,]
    }
    if(is.null(fitObj$data$curePreds)){
      rcurePreds = NULL
    }else{
      rcurePreds = fitObj$data$curePreds[rsample,]
    }

    rfit = kmcure (rtime, revent, rsurvPreds, rcurePreds,
                   optim_init=optim_init, silent=TRUE,
                   multiOptim_maxit=multiOptim_maxit, multiOptim_reltol=multiOptim_reltol,
                   optim_reltol=optim_reltol, optim_maxit=optim_maxit)

    if(rfit$exitcode==0){
      r = r + 1
      fitObj$boot$repeats = r
      fitObj$boot$coef = cbind(fitObj$boot$coef, rfit$coef)
      fitObj$boot$pcure = c(fitObj$boot$pcure, rfit$pcure)
      fitObj$boot$pcens = c(fitObj$boot$pcens, rfit$pcens)
      rfit$data = NULL
      fitObj$boot$fitsList[[r]] = rfit
      fitObj$boot$settings = list(useFitEstAsBootInit=useFitEstAsBootInit,
                                  multiOptim_maxit=multiOptim_maxit,
                                  multiOptim_reltol=multiOptim_reltol,
                                  optim_reltol=optim_reltol,
                                  optim_maxit=optim_maxit)
      assign(fitObjName, fitObj, envir = .GlobalEnv) # update the original object with fitObjName name
      if(!silent) cat("The bootstrap replication #", r, "of", R, "is now complate.\n")
    } # end if condition on exitcode
  } # end while loop on r
  if(!silent) cat('The run of kmcure "bootstrap" function is finished at', format(Sys.time(), "%H:%M:%S (%Y-%m-%d)..."), "\n")
}
