#' LSA of the kmcure fit
#'
#' LSA of the kmcure fit that need coef, covmat, n, and type
#' @param coef The estimated coef of a kmcure fit.
#' @param covmat The estimated covariance matrix of a kmcure fit
#' @param n The number of events.
#' @param type variable selection type, choose form "lasso" or "lar".
#' @return object.
#' @export
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

  # ---- shrinkage the b part using our modified lars.lsa which export lambda
  bPart_lsa <- suppressWarnings(lars.lsa(bSI, bcoef, bintercept, n, type))

  bExtract = extract_larslsa(l.fit = bPart_lsa, beta.ols=bcoef, intercept=bintercept, n=n)

  # ---- shrinkage the beta part using our modified lars.lsa which export lambda
  betaPart_lsa <- suppressWarnings(lars.lsa(betaSI, betacoef, betaintercept, n, type))

  betaExtract = extract_larslsa(l.fit = betaPart_lsa, beta.ols=betacoef, intercept=betaintercept, n=n)

  # prepare objects for output
  theta.ols = as.matrix(c(bExtract$beta.ols, betaExtract$beta.ols))
  theta.bic = as.matrix(c(bExtract$beta.bic, betaExtract$beta.bic))
  theta.aic = as.matrix(c(bExtract$beta.aic, betaExtract$beta.aic))
  rownames(theta.ols) = rownames(coef); colnames(theta.ols) = "theta.ols"
  rownames(theta.bic) = rownames(coef); colnames(theta.bic) = "theta.bic"
  rownames(theta.aic) = rownames(coef); colnames(theta.aic) = "theta.aic"

  best.bic = c(bExtract$best.bic, betaExtract$best.bic)
  best.aic = c(bExtract$best.aic, betaExtract$best.aic)
  names(best.bic)=c("best.bic(Gamma)","best.bic(Beta)")
  names(best.aic)=c("best.aic(Gamma)","best.aic(Beta)")

  lambda.bic = c(bExtract$lambda.bic, betaExtract$lambda.bic)
  lambda.aic = c(bExtract$lambda.aic, betaExtract$lambda.aic)
  names(lambda.bic)=c("lambda.bic(Gamma)","lambda.bic(Beta)")
  names(lambda.aic)=c("lambda.aic(Gamma)","lambda.aic(Beta)")

  LSAkmcure_outlist = list(theta.ols=theta.ols,
                           theta.bic=theta.bic, theta.aic=theta.aic,
                           best.bic=best.bic, best.aic=best.aic,
                           lambda.bic=lambda.bic, lambda.aic=lambda.aic)
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

  names(beta.bic) = names(beta.ols)
  names(beta.aic) = names(beta.ols)

  best.bic = l.fit$BIC[t1$ix[1]]
  best.aic = l.fit$AIC[t2$ix[1]]

  lambda.bic = l.fit$lambda[t1$ix[1]] #MM: best lambda based on the BIC criteria
  lambda.aic = l.fit$lambda[t2$ix[1]] #MM: best lambda based on the AIC criteria

  obj <- list(beta.ols = beta.ols,
              beta.bic = beta.bic, beta.aic = beta.aic,
              best.bic = best.bic, best.aic = best.aic,
              lambda.bic = lambda.bic, lambda.aic = lambda.aic)

  return(obj)
}

