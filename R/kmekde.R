#' The computational function that is called by the "kmcure" function to fits AFT Semiparametric Mixture Cure Model using the KME-KDE method
#'
#' Fits AFT (Accelerated Failure Time) Semiparametric Mixture Cure Model using the KME-KDE (Kaplan-Meier Estimation and Kernel Density Estimator) method.
#'
#' @param time a survival time to event variable.
#' @param event a survival status variable: 1 for event and 0 for censoring.
#' @param survPreds a matrix of survival predictor variable(s).
#' @param curePreds an optional matrix of curing predictor variable(s).
#' @param multiOptim_maxit a number showing the maximum of allowed multiple optimization. The program uses multiple optimization if the convergence of "optim" does not meet.
#' @param multiOptim_reltol a number showing the relative tolerance in continuing multiple optimization procedure.
#' @param multiOptim_stopTime an optional number showing time limit per minutes to stop multi-optimization based on calculation time per minutes.
#' @param multiOptim_stopLLp an optional proportion defining a stopping rule based on changes in log-likelihood in successive multi-optimization. 0 disable this stopping rule, and 0.1 stop multi-optimization when the difference in the latest loglik runs becomes less than or equal to 0.1 of difference between loglik values in the first and second "optim" runs.
#' @param optim_reltol a number showing the relative tolerance in continuing of each optimization.
#' @param optim_maxit a number showing the maximum of allowed iterations in each optimization.
#' @param silent a Boolean value which if set to TRUE it prevent from showing output messages.
#' @param conditional a Boolean value which, if set to TRUE it uses an iterative procedure that estimate parameters survival and cure sub-models conditionally on the last estimation of the other one.
#' @param cond_reltol a number showing the relative tolerance in continuing the conditional algorithm optimization.
#' @param cond_maxit a number showing the maximum of allowed iterations in the conditional algorithm optimization.
#' @param cond_reltol_beta a number showing the relative tolerance in continuing optimization of the beta part in the conditional algorithm.
#' @param cond_maxit_beta a number showing the maximum of allowed iterations in the optimization of the beta part in the conditional algorithm.
#' @param cond_reltol_gamma a number showing the relative tolerance in continuing optimization of the gamma part in the conditional algorithm.
#' @param cond_maxit_gamma a number showing the maximum of allowed iterations in the optimization of the gamma part in the conditional algorithm.
#' @param fix_gammacoef an optional numeric vector of fix gamma coefficients that could be used to estimate beta coefficients based of them. This can be used for conditional optimization.
#' @param fix_betacoef an optional numeric vector of fix beta coefficients that could be used to estimate gamma coefficients based of them. This can be used for conditional optimization.
#' @param bandcoef an optional coefficient to multiply the optimal kernel smoothing band-width (Please note the loglik values resulted by applying different "bandcoef" are not comparable. So, changing the default value of this option is Not recommended).
#' @param try_hessian a Boolean value with default value of FALSE. If this set to TRUE, the Hessian matrix will be evaluated after last optimization by applying a final optimization that also try to estimate the hessian matrix.
#' @param optim_method a string value showing the method of optimization: "Nelder-Mead" and "SANN" are supported.
#' @param optim_init an optional numeric vector of initial values. For example, it could be the estimated coefficients of a previous fit to be used in continuing of optimization.
#'
#' @return a list of "exitcode" (0: no warning or error, 1: warning, 2 or 3: error), "coef" (estimated coefficients), "AIC", etc. check names(fit) for more information.
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
#' # fit1 = kmekde (time, event, survPreds, curePreds, multiOptim_maxit = 10)
#'
#' # names(fit1)
#'
#' # fit1$exitcode==TRUE # if TRUE the fit completed without any warning/error
#' # fit1$loglik # the loglik of the fitted model
#' # fit1$timeD # the calculation time
#' # fit1$coef # the estimated coefficients
#'
#' @import survival stats
#' @importFrom survival Surv survfit
#' @importFrom stats binomial dnorm glm lm optim var
#' @export
kmekde <- function(time, event, survPreds, curePreds=NULL,
                   multiOptim_maxit = 1, multiOptim_reltol = 0.001,
                   multiOptim_stopTime = NULL, multiOptim_stopLLp = 0,
                   optim_reltol = 1e-8, optim_maxit = 500,
                   silent = FALSE, conditional = FALSE,
                   cond_reltol = 1e-8, cond_maxit = 250,
                   cond_reltol_beta = 1e-8, cond_maxit_beta = 500,
                   cond_reltol_gamma = 1e-5, cond_maxit_gamma = 50,
                   fix_gammacoef = NULL, fix_betacoef = NULL,
                   bandcoef=1, try_hessian=FALSE,
                   optim_method = "Nelder-Mead", optim_init = NULL)
  {

  call <- match.call()

  if(multiOptim_maxit > 1000) multiOptim_maxit = 1000
  if(multiOptim_stopLLp < 0) multiOptim_stopLLp = 0
  if(multiOptim_stopLLp > 1) multiOptim_stopLLp = 1

  methodName = "KMEKDE"

  tryCatch({

    Y = as.matrix(time) # to be used with exactly this name
    delta = as.matrix(event) # to be used with exactly this name
    Xnames = colnames(survPreds)
    X = as.matrix(survPreds) # to be used with exactly this name
    if(is.null(Xnames)) Xnames = paste0("X", 1:ncol(X))
    if(!is.null(curePreds)){
      Znames = colnames(curePreds)
      Z = as.matrix(curePreds) # to be used with exactly this name
      if(is.null(Znames)) Znames = paste0("Z", 1:ncol(Z))
      oneZ = cbind(1,Z) # to be used with exactly this name
    } else{
      Znames = ""
      oneZ = matrix(1, nrow=length(Y)) # to be used with exactly this name
    }



    ## claculate logY and n
    logY = log(Y) # to be used with exactly this name
    n = length(logY) # to be used with exactly this name

    if(!silent) cat("The program run is started at", format(Sys.time(), "%H:%M:%S (%Y-%m-%d)."), "Please be patient...\n")

    ## calculate cure fraction and censoring percent just to be used in report
    sfit = survival::survfit(survival::Surv(Y, delta)~1)
    pcure = min(sfit$surv)
    pcens = mean(delta==0)

    timeS = Sys.time() # start recording time

    ### Initial values

    ## Initial Values for Logistic model (Uncure Prob.)
    logit.fit = glm( delta ~ oneZ - 1 , family = binomial(link="logit"), control = list(maxit = 100)) # the number of maxit increased from the default value of 25 to 100 to avoid possible non-convergence warnings
    b.ini = coef(logit.fit)

    ## Initial Values for AFT Model of Uncured Individuals
    lm.fit = lm(logY[delta==1]~X[delta==1,])
    beta.ini = lm.fit$coef[-1]

    ## Vector of total initial values for all parameters
    theta = c(b.ini,beta.ini)

    ### define functions

    ## -----------------KERNEL BANDWIDTH ------------------------------------------------------------
    sigma = sqrt(var(logY[delta==1]-as.matrix(X[delta==1,])%*%as.matrix(beta.ini)))
    sigma = as.numeric(sigma)
    optimband = (8*sqrt(2)/3)^(1/5)*sigma*n^(-1/5)
    h = bandcoef*optimband # to be used with exactly this name

    ###------------------------------------------------------------------------------

    ##---------------------- KAPLAN-MEIER --------------------------------------------------------------------------
    s0 = function(eps,delta) # the input delta here is not the delta variable
    {
      if(length(eps)!=length(delta)) return(NA)
      ss = survival::survfit(survival::Surv(eps,delta)~1)
      ss = summary(ss,times=eps)$surv
      return((ss-ss[n])/(1-ss[n]))
    }

    ##------------------- KERNEL -----------------------------------------------------------------------------
    f0 = function(eps,delta,h) # the input delta here is not the delta variable
    {
      # Other input: n
      if(length(eps)!=length(delta)) return(NA)
      ss = s0(eps,delta)
      ss1 = c(1,ss[1:n-1])
      a = ss1-ss
      epss = matrix(eps,nrow=n,ncol=n,byrow=TRUE)
      u = (eps-epss)/h
      # epsmatrix = matrix(eps,nrow=n,ncol=n,byrow=FALSE)
      # u= (epsmatrix-epss)/h
      K = dnorm(u)/h
      dens = K%*%a
      dens = replace(dens,dens==0,0.001)
      return(as.numeric(dens))
    }

    ##----------------------------------------------------------------------------------------------------

    fixInf = function(xvec){
      maxNonInf = 1.797 * 10^308
      xvec[xvec == Inf] = maxNonInf
      xvec[xvec == -Inf] = -maxNonInf
      return(xvec)
    }


    loglik = function(theta)
    {
      # Other inputs: logY , delta , X , oneZ

      tryCatch({

        h = as.numeric(h)[1] # h is a value

        bvec = as.matrix(theta[1:length(b.ini)])

        betavec = as.matrix(theta[(length(b.ini)+1):(length(b.ini)+length(beta.ini))])

        eps = as.numeric(logY - as.matrix(X)%*%betavec)

        De = delta[order(eps)]

        Ze = oneZ[order(eps),]

        Ze_b = as.numeric(as.matrix(Ze)%*%bvec)

        eps = sort(eps) # needs to be sorted
        ##--------------------------------------------------------------------------------------------------------------------------
        px = exp(Ze_b)/( 1 + exp(Ze_b))
        ##---------------------------------------------------------------------------------------------------------------------------

        part1 = log(px*f0(eps,De,h))
        part2 = log(1-px+px*s0(eps,De))
        return(-sum(De*fixInf(part1) + (1-De)*fixInf(part2)))

      }, warning = function(cond) {
        return(NA)
      }, error = function(cond) {
        return(NA)
      }
      ) # end tryCatch
    }

    #---------------------------------------------------------------------------------------------------------------------------


    ### Model Fitting

    #loglik(theta)


    ##### Start the code of the main optimization section (multi-optimization #1)

    ## define conditional likelihood functions
    loglik_fix_b <- function(beta.par, b.par){
      theta = c(b.par, beta.par)
      return( loglik(theta) )
    }

    loglik_fix_beta <- function(b.par, beta.par){
      theta = c(b.par, beta.par)
      return( loglik(theta) )
    }


    ## now we are going to estimate the requested parameters
    conv = NA

    theta = c(b.ini,beta.ini)
    if(length(optim_init)==length(theta)) theta = as.numeric(optim_init)

    if(length(fix_gammacoef)==length(b.ini)){ # then use fix_gammacoef to estimate just the beta part
      methodName = "KMEKDE-fixGamma"
      b.par = fix_gammacoef
      beta.par = theta[(length(b.ini)+1) : length(theta)]

      #Beta-Part
      if(length(beta.par>1)){
        beta.optim = optim(par=beta.par, fn=loglik_fix_b, b.par=b.par, hessian = FALSE, control = list(reltol=cond_reltol_beta, maxit = cond_maxit_beta), method = optim_method)
      }else{
        beta.optim = optim(par=beta.par, fn=loglik_fix_b, b.par=b.par, hessian = FALSE, control = list(reltol=cond_reltol_beta, maxit = cond_maxit_beta), method = "Brent")
      }
      beta.par = beta.optim$par
      conv = beta.optim$convergence
      thetaUpdate = c(b.par, beta.par)

    }else if(length(fix_betacoef)==length(beta.ini)){ # then use fix_betacoef to estimate just the b part
      methodName = "KMEKDE-fixBeta"
      b.par = theta[1 : length(b.ini)]
      beta.par = fix_betacoef

      #Gamma-Part
      if(length(b.par>1)){
        b.optim = optim(par=b.par, fn=loglik_fix_beta, beta.par=beta.par, hessian = FALSE, control = list(reltol=cond_reltol_gamma, maxit = cond_maxit_gamma), method = optim_method)
      }else{
        b.optim = optim(par=b.par, fn=loglik_fix_beta, beta.par=beta.par, hessian = FALSE, control = list(reltol=cond_reltol_gamma, maxit = cond_maxit_gamma), method = "Brent")
      }
      b.par = b.optim$par
      conv = b.optim$convergence
      thetaUpdate = c(b.par, beta.par)

    }else if(conditional==TRUE){ # then use complete conditional optimization to estimate both of the b and beta parts
      methodName = "KMEKDE-conditional"

      iter = 1
      repeat{

        b.par = theta[1 : length(b.ini)]
        beta.par = theta[(length(b.ini)+1) : length(theta)]


        #Gamma-Part
        if(length(b.par>1)){
          b.optim = optim(par=b.par, fn=loglik_fix_beta, beta.par=beta.par, hessian = FALSE, control = list(reltol=cond_reltol_gamma, maxit = cond_maxit_gamma), method = optim_method)
        }else{
          b.optim = optim(par=b.par, fn=loglik_fix_beta, beta.par=beta.par, hessian = FALSE, control = list(reltol=cond_reltol_gamma, maxit = cond_maxit_gamma), method = "Brent")
        }
        b.par = b.optim$par

        #Beta-Part
        if(length(beta.par>1)){
          beta.optim = optim(par=beta.par, fn=loglik_fix_b, b.par=b.par, hessian = FALSE, control = list(reltol=cond_reltol_beta, maxit = cond_maxit_beta), method = optim_method)
        }else{
          beta.optim = optim(par=beta.par, fn=loglik_fix_b, b.par=b.par, hessian = FALSE, control = list(reltol=cond_reltol_beta, maxit = cond_maxit_beta), method = "Brent")
        }
        beta.par = beta.optim$par


        thetaUpdate = c(b.par, beta.par)
        reltolcondition = all(abs((thetaUpdate-theta)) < (cond_reltol * (abs(theta) + cond_reltol)))
        if(iter>cond_maxit) conv = 1
        if(reltolcondition) conv = 0
        if(reltolcondition | (iter>cond_maxit)) (break)

        iter = iter+1
        theta = thetaUpdate
      }

    }else{# then use concurrent optimization to estimate both of the b and beta parts
      methodName = "KMEKDE-concurrent"

      theta.optim = optim(par=theta, fn=loglik, hessian = FALSE, control = list(reltol=optim_reltol, maxit = optim_maxit), method = optim_method)
      thetaUpdate = theta.optim$par
      conv = theta.optim$convergence
    }

    timeE = Sys.time()
    timeD = difftime(timeE, timeS, units = "mins")
    vectimeD = timeD

    ##### End the code of the main optimization section (multi-optimization #1)
    vecLL = -1*loglik(thetaUpdate)
    vecconv = conv
    if(!silent) cat("Optimization # 1 is completed. The optim log-likelihood is", vecLL[1], "\n")
    mat_coef = matrix(thetaUpdate, nrow=length(theta), ncol = 1)

    ### Start applying further optimization if multiOptim_maxit > 1
    if(multiOptim_maxit >= 2){
      for(k in 2:round(multiOptim_maxit)){
        if(vecconv[length(vecconv )]!=0){ # do only if the previous optimization was not converged
        tryCatch({
          if(!is.null(multiOptim_stopTime)){
            timeE = Sys.time()
            timeD = difftime(timeE, timeS, units = "mins")
            if (timeD > multiOptim_stopTime) (break)
          }
          theta = thetaUpdate
          timeSloop = Sys.time()
          theta.optim = optim(par=theta, fn=loglik, hessian = FALSE, control = list(reltol=optim_reltol, maxit = optim_maxit), method = optim_method)
          timeE = Sys.time()
          timeD = difftime(timeE, timeSloop, units = "mins")
          vectimeD = c(vectimeD, timeD)
          thetaUpdate = theta.optim$par
          mat_coef = cbind(mat_coef, thetaUpdate)
          LL = theta.optim$value
          vecLL = c(vecLL,-LL)
          conv = theta.optim$convergence
          vecconv = c(vecconv,conv)
          if(!silent) cat("Optimization #", length(vecLL) ,"is completed. The optim log-likelihood is", vecLL[length(vecLL)], "\n")
          LL_exitCondition = FALSE
          if( abs(vecLL[length(vecLL)-1]-vecLL[length(vecLL)]) <= abs(vecLL[2]-vecLL[1]) * multiOptim_stopLLp ) LL_exitCondition = TRUE
          reltolcondition = all(abs((thetaUpdate-theta)) < (multiOptim_reltol * (abs(theta) + multiOptim_reltol)))
          if(reltolcondition | LL_exitCondition) (break)
        }, warning = function(cond) {
          if(!silent) cat("Skip optimization #", k, " duo to ocurance of a warning.\n")
        }, error = function(cond) {
          if(!silent) cat("Skip optimization #", k, " duo to ocurance of an error.\n")
        }
        ) # end tryCatch
      } # end if of vecconv
    } # end for loop of k
    } # end if of multiOptim_maxit

    ### End applying further optimization if multiOptim_maxit > 1

    ### start calculation hessian if try_hessian is TRUE
    estHessian = NA
    if(try_hessian){
      if(!silent) cat("Optimization #", length(vecLL)+1 ,"is started that try to find an estimate for the Hessian matrix.\n")
      tryCatch({
        theta = thetaUpdate
        theta.optim = optim(par=theta, fn=loglik, hessian = TRUE, control = list(reltol=optim_reltol, maxit = optim_maxit), method = optim_method)
        thetaUpdate = theta.optim$par
        mat_coef = cbind(mat_coef, thetaUpdate)
        LL = theta.optim$value
        vecLL = c(vecLL,-LL)
        estHessian = theta.optim$hessian
        if(!silent) cat("Optimization #", length(vecLL) ,"and claculating the Hessian matrix is completed. The optim log-likelihood is", vecLL[length(vecLL)], "\n")
      }, warning = function(condHessian) {
        if(!silent) cat("As expected, a warning is occurred in trying to find the Hessian matrix!\n")
      }, error = function(condHessian) {
        if(!silent) cat("As expected, an error is occurred in trying to find the Hessian matrix!\n")
      }
      ) # end tryCatch
    } # end if condition on try_hessian
    ### end calculation of hessian if try_hessian is TRUE

    timeE = Sys.time() # end recording time
    timeD = difftime(timeE, timeS, units = "mins")

    coef = thetaUpdate # use the latest estimation of theta that is result of the optimization
    coef = as.matrix(coef)
    if(all(Znames=="")){
      rownames(coef) = c( "Gamma(Intercept)", paste("Beta(",Xnames,")",sep="") )
    }else{
      rownames(coef) = c( "Gamma(Intercept)", paste("Gamma(",Znames,")",sep=""), paste("Beta(",Xnames,")",sep="") )
    }

    colnames(coef) = ("Estimate")
    #print(coef)
    rownames(mat_coef) = rownames(coef)

    vecAIC = -2*vecLL+2*length(coef)
    vecBIC = -2*vecLL+log(n)*length(coef)

    indexMaxLL = which.max(vecLL)
    coef = as.matrix(mat_coef [,indexMaxLL])
    loglik = vecLL[indexMaxLL]
    AIC=vecAIC[indexMaxLL]
    BIC=vecBIC[indexMaxLL]

    if(!silent) cat("The program run is successfully finished at", format(Sys.time(), "%H:%M:%S (%Y-%m-%d)."), "\n")
    return( list(exitcode = 0, call=call, coef=coef, method = methodName, timeD=timeD, AIC=AIC, BIC=BIC, loglik=loglik, pcure=pcure, pcens=pcens, optimband=optimband, bandcoef=bandcoef, bandwidth=h, mat_coef=mat_coef, vecLL=vecLL, vecAIC=vecAIC, hessian = estHessian, vecconv=vecconv, vectimeD=vectimeD))
  }, warning = function(w) {
    if(!silent){
      message("A warning in the KMEKDE fitting is occured!")
      message(w)
    }
    return( list(exitcode = 1, call=call, coef=coef, method = methodName, timeD=timeD, AIC=AIC, BIC=BIC, loglik=loglik, pcure=pcure, pcens=pcens, optimband=optimband, bandcoef=bandcoef, bandwidth=h, mat_coef=mat_coef, vecLL=vecLL, vecAIC=vecAIC, hessian = estHessian, vecconv=vecconv, vectimeD=vectimeD))
  }, error = function(e) {
    if(!silent){
      message("An error in the KMEKDE fitting is occured!")
      message(e)
    }
    return(list(exitcode = 2, call=call, coef=NA, method = methodName, timeD = NA, AIC=NA, BIC=NA, loglik=NA, pcure=pcure, pcens=pcens, optimband=optimband, bandcoef=bandcoef, bandwidth=h, mat_coef=NULL, vecLL=NULL, vecAIC=NULL, hessian = NA, vecconv=NA, vectimeD=NA))
  }, finally = {
    if(all(typeof(coef)!="closure")) if(any(is.na(coef))){
      if(!silent) {
        message("NA coefficients in the KMEKDE fitting is occured!")
      }
      return(list(exitcode = 3, call=call, coef=NA, method = methodName, timeD = NA, AIC=NA, BIC=NA, loglik=NA, pcure=pcure, pcens=pcens, optimband=optimband, bandcoef=bandcoef, bandwidth=h, mat_coef=NULL, vecLL=NULL, vecAIC=NULL, hessian = NA, vecconv=NA, vectimeD=NA))
    }
  }
  ) # end tryCatch

} # end KMEKDE
