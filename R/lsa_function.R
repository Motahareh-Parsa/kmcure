# This modified version is based on the lsa_function.R in the RobMixReg package
# This modified lars.lsa function also exports the used lambda values, and
# the lsa function also exports best.bic, best.aic, lambda.bic, lambda.aic values.

#require(lars)

## Least square approximation. Based on the Oct 19, 2006, version of the R RobMixReg package. This modified version also provide lambda values similar to lars::lars function.

## Reference Wang, H. and Leng, C. (2006) and Efron et al. (2004).

##

## Written by Chenlei Leng

## Comments and suggestions are welcome

##

## Input

## obj: lm/glm/coxph or other object

##

## Output

## beta.ols: the MLE estimate

## beta.bic: the LSA-BIC estimate

## beta.aic: the LSA-AIC estimate



#' Ready to use function to apply the lars.lsa function to a fitted object of lm/glm/coxph.
#'
#' Ready to use function to apply the lars.lsa function to a fitted object of lm/glm/coxph/etc. This modified version is based on the version of Oct 19, 2006 of the lsa function in the RobMixReg package. Now it also reports lambda, BIC, and AIC.
#' @author Reference Wang, H. and Leng, C. (2006) and Efron et al. (2004).
#' @param obj lm/glm/coxph or other object.
#' @return beta.ols: the MLE estimate ; beta.bic: the LSA-BIC estimate ; beta.aic: the LSA-AIC estimate.
#' @export
lsa <- function(obj)

{

  if(class(obj)[1]=='kmcure'){
    if(is.null(obj$boot$coef)){
      stop("There is no bootstrap replication in the inputted kmcure object. See help(bootstrap)")
    }else if(ncol(obj$boot$coef)<50){
      stop("There must be at least 50 bootstrap replications in the inputted kmcure object")
    }else{
      coef = obj$coef
      n = length(obj$data$time)
      covmat = cov(t(obj$boot$coef))
      output = suppressWarnings(LSAkmcure(coef, covmat, n, type="lasso"))
      return(output)
    }
  } # end block of the kmcure object


  intercept <- attr(obj$terms,'intercept')

  if(class(obj)[1]=='coxph') intercept <- 0



  n <- length(obj$residuals)



  Sigma <- vcov(obj)

  SI <- solve(Sigma)

  beta.ols <- coef(obj)

  l.fit <- suppressWarnings(lars.lsa(SI, beta.ols, intercept, n))



  t1 <- sort(l.fit$BIC, ind=T)

  t2 <- sort(l.fit$AIC, ind=T)

  beta <- l.fit$beta

  if(intercept) {

    beta0 <- l.fit$beta0+beta.ols[1]

    beta.bic <- c(beta0[t1$ix[1]],beta[t1$ix[1],])

    beta.aic <- c(beta0[t2$ix[1]],beta[t2$ix[1],])

  }

  else {

    beta0 <- l.fit$beta0

    beta.bic <- beta[t1$ix[1],]

    beta.aic <- beta[t2$ix[1],]

  }

  names(beta.bic) = names(beta.ols) #MM: use beta.ols names to name the beta.bic
  names(beta.aic) = names(beta.ols) #MM: use beta.ols names to name the beta.aic

  best.bic = l.fit$BIC[t1$ix[1]] #MM: best BIC value among the fitted models
  best.aic = l.fit$AIC[t2$ix[1]] #MM: best AIC value among the fitted models

  lambda.bic = l.fit$lambda[t1$ix[1]] #MM: best lambda based on the BIC criteria
  lambda.aic = l.fit$lambda[t2$ix[1]] #MM: best lambda based on the AIC criteria

  obj <- list(beta.ols = beta.ols,
              beta.bic = beta.bic, beta.aic = beta.aic,
              best.bic = best.bic, best.aic = best.aic,
              lambda.bic = lambda.bic, lambda.aic = lambda.aic)

  obj

}





###################################

## lars variant for LSA


#' lars variant for LSA (modifed version that also reports lambda)
#'
#' Solve penalized Least square approximation (LSA) using the Least Angle Regression (LAR) algorithm.
#' @author Reference Wang, H. and Leng, C. (2006) and Efron et al. (2004).
#' @param Sigma0 The parameter.
#' @param b0 The intercept of the regression line.
#' @param intercept The bool variable of whether consider the intercept situation
#' @param n The number of data point.
#' @param type Regression options, choose form "lasso" or "lar".
#' @param eps The converge threshold defined by the machine.
#' @param max.steps The maximum iteration times to stop.
#' @return object.
#' @export
lars.lsa <- function (Sigma0, b0, intercept,  n,

                      type = c("lasso", "lar"),

                      eps = .Machine$double.eps, max.steps)

{

  type <- match.arg(type)

  TYPE <- switch(type, lasso = "LASSO", lar = "LAR")



  n1 <- dim(Sigma0)[1]



  ## handle intercept

  if (intercept) {

    a11 <- Sigma0[1,1]

    a12 <- Sigma0[2:n1,1]

    a22 <- Sigma0[2:n1,2:n1]

    Sigma <- a22-outer(a12,a12)/a11

    b <- b0[2:n1]

    beta0 <- crossprod(a12,b)/a11

  }

  else {

    Sigma <- Sigma0

    b <- b0

  }



  Sigma <- diag(abs(b))%*%Sigma%*%diag(abs(b))

  b <- sign(b)



  nm <- dim(Sigma)

  m <- nm[2]

  im <- inactive <- seq(m)



  Cvec <- drop(t(b)%*%Sigma)

  ssy <- sum(Cvec*b)

  if (missing(max.steps))

    max.steps <- 8 * m

  beta <- matrix(0, max.steps + 1, m)

  lambda=double(max.steps) #MM: initiate the lambda vector similar to the lars::lars function

  Gamrat <- NULL

  arc.length <- NULL

  R2 <- 1

  RSS <- ssy

  first.in <- integer(m)

  active <- NULL

  actions <- as.list(seq(max.steps))

  drops <- FALSE

  Sign <- NULL

  R <- NULL

  k <- 0

  ignores <- NULL



  while ((k < max.steps) & (length(active) < m)) {

    action <- NULL

    k <- k + 1

    C <- Cvec[inactive]

    Cmax <- max(abs(C))

    lambda[k]=Cmax #MM: fill the lambda vector similar to the lars::lars function

    if (!any(drops)) {

      new <- abs(C) >= Cmax - eps

      C <- C[!new]

      new <- inactive[new]

      for (inew in new) {

        R <- lars::updateR(Sigma[inew, inew], R, drop(Sigma[inew, active]),

                     Gram = TRUE,eps=eps)

        if(attr(R, "rank") == length(active)) {

          ##singularity; back out

          nR <- seq(length(active))

          R <- R[nR, nR, drop = FALSE]

          attr(R, "rank") <- length(active)

          ignores <- c(ignores, inew)

          action <- c(action,  - inew)

        }

        else {

          if(first.in[inew] == 0)

            first.in[inew] <- k

          active <- c(active, inew)

          Sign <- c(Sign, sign(Cvec[inew]))

          action <- c(action, inew)

        }

      }

    }

    else action <- -dropid

    Gi1 <- backsolve(R, lars::backsolvet(R, Sign))

    dropouts <- NULL

    A <- 1/sqrt(sum(Gi1 * Sign))

    w <- A * Gi1



    if (length(active) >= m) {

      gamhat <- Cmax/A

    }

    else {

      a <- drop(w %*% Sigma[active, -c(active,ignores), drop = FALSE])

      gam <- c((Cmax - C)/(A - a), (Cmax + C)/(A + a))

      gamhat <- min(gam[gam > eps], Cmax/A)

    }

    if (type == "lasso") {

      dropid <- NULL

      b1 <- beta[k, active]

      z1 <- -b1/w

      zmin <- min(z1[z1 > eps], gamhat)

      # cat('zmin ',zmin, ' gamhat ',gamhat,'\n')

      if (zmin < gamhat) {

        gamhat <- zmin

        drops <- z1 == zmin

      }

      else drops <- FALSE

    }

    beta[k + 1, ] <- beta[k, ]

    beta[k + 1, active] <- beta[k + 1, active] + gamhat * w



    Cvec <- Cvec - gamhat * Sigma[, active, drop = FALSE] %*% w

    Gamrat <- c(Gamrat, gamhat/(Cmax/A))



    arc.length <- c(arc.length, gamhat)

    if (type == "lasso" && any(drops)) {

      dropid <- seq(drops)[drops]

      for (id in rev(dropid)) {

        R <- lars::downdateR(R,id)

      }

      dropid <- active[drops]

      beta[k + 1, dropid] <- 0

      active <- active[!drops]

      Sign <- Sign[!drops]

    }



    actions[[k]] <- action

    inactive <- im[-c(active)]

  }

  beta <- beta[seq(k + 1), ]

  lambda=lambda[seq(k)] #MM: trim the lambda vector similar to the lars::lars function

  dff <- b-t(beta)



  RSS <- diag(t(dff)%*%Sigma%*%dff)



  if(intercept)

    beta <- t(abs(b0[2:n1])*t(beta))

  else

    beta <- t(abs(b0)*t(beta))



  if (intercept) {

    beta0 <- beta0-drop(t(a12)%*%t(beta))/a11

  }

  else {

    beta0 <- rep(0,k+1)

  }

  dof <- apply(abs(beta)>eps,1,sum)

  BIC <- RSS+log(n)*dof

  AIC <- RSS+2*dof

  object <- list(AIC = AIC, BIC = BIC, lambda=lambda, #MM: export the lambda vector similar to the lars::lars function

                 beta = beta, beta0 = beta0)

  object

}



##Note that rq() object implemented the coef()

##but without vcov() implementation. We provide

##a rather simple implementation here.

##This part is written by Hansheng Wang.

# vcov.rq <- function(object,...)
#
# {
#
#   q=object$tau
#
#   x=as.matrix(object$x)
#
#   resid=object$residuals
#
#   f0=density(resid,n=1,from=0,to=0)$y
#
#   COV=q*(1-q)*solve(t(x)%*%x)/f0^2
#
#   COV
#
# }

