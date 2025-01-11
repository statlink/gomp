cv.gomp <- function(y, x, kfolds = 10, folds = NULL, tol = seq(4, 9, by = 1), task = "C", metric = NULL,
                    metricbbc = NULL, modeler = NULL, test = NULL, method = "ar2", B = 1) {

    if ( is.null(tol) )   tol <- seq(4, 9, by = 1)
    tol <- sort(tol)
    ntol <- length(tol);
    sel.vars <- list()
    cv_results_all <- list()
    nama <- paste("tol=", tol, sep = "")

    if ( is.null(folds) ) {
      if (task == "R" ) {
        folds <- gomp::makefolds(y, nfolds = kfolds, stratified = FALSE, seed = FALSE)
      } else if (task == "S") {
        folds <- gomp::makefolds(y[, 1], nfolds = kfolds, stratified = FALSE, seed = FALSE)
      } else if (task == "C")  folds <- gomp::makefolds(y, nfolds = kfolds, stratified = TRUE, seed = FALSE)
    } else  kfolds <- length( folds );

    if ( is.null(task) ) {
      stop("Please provide a valid task argument 'C'-classification, 'R'-regression, 'S'-survival.")
      #to do: select automatically the appropriate task due to the data, target
    } else if(task == 'C') {

      #Classification task (logistic regression)
      if ( is.null(metric) ) {
        metricFunction <- .auc.gomp;
      } else   metricFunction <- metric;

      if ( is.null(modeler) ) {
        modelerFunction <- .glm.gomp;
      } else   modelerFunction <- modeler;

      if ( is.null(test) ) {
        test <- 'logistic';
      } else  test <- test;

    } else if (task == 'R') {

      #Regression task (linear regression)
      if ( is.null(metric) ) {
        metricFunction <- .mse.gomp;
      } else  metricFunction <- metric;

      if ( is.null(modeler) ) {
        modelerFunction <- .lm.gomp;
      } else  modelerFunction <- modeler;

      if ( is.null(test) ) {
        test <- 'normal';
      } else  test <- test;

    } else if(task == 'S') {

      #cox survival analysis (cox regression)
      if ( is.null(metric) ) {
        metricFunction <- .ci.gomp;
      } else  metricFunction <- metric;

      if ( is.null(modeler) ) {
        modelerFunction <- .coxph.gomp;
      } else  modelerFunction <- modeler;

      if ( is.null(test) ) {
        test <- "cox";
      } else  test <- test;

    } else  stop("Please provide a valid task argument 'C'-classification, 'R'-regression, 'S'-survival.")

    tic <- proc.time()

    for (k in 1:kfolds) {
      #print(paste('CV: Fold', k, 'of', kfolds));
      train_samples <- c();
      for ( i in which(c(1:kfolds) != k) )  train_samples = c( train_samples, folds[[ i ]] )
      #leave one fold out each time as a test set and the rest as train set
      xtrain <- x[train_samples, ] #Set the training set
      ytrain <- y[ train_samples ]
      xtest <- x[ folds[[ k ]], ] #Set the validation set
      ytest <- y[ folds[[ k ]] ]
      dm <- dim(ytest)

      results <- gomp::gomp.path(y = ytrain, x = xtrain, tol = tol, test = test, method = method)
      sel.vars[[ k ]] <- results$res[-1, -c(ntol + 1), drop = FALSE]
      cv_results_all[[ k ]] <- list()
      cv_results_all[[ k ]]$preds <- matrix(0, dm[1], ntol)
      colnames(cv_results_all[[ k ]]$preds) <- nama
      cv_results_all[[ k ]]$performances <- numeric(ntol)
      names(cv_results_all[[ k ]]$performances) <- nama
      cv_results_all[[ k ]]$selectedVars <- sel.vars[[ k ]]
      colnames(cv_results_all[[ k ]]$selectedVars) <- nama

      for ( j in 1:ntol ) {

        variables <- sel.vars[[ k ]][, j]
        sign_data <- xtrain[, variables, drop = FALSE]
        sign_test <- xtest[, variables, drop = FALSE]

        if ( sum( variables > 0 ) > 0 ) {
          #generate a model due to the task and find the performance
          #logistic model for a classification task, linear model for the regression task and a cox model for the survival task
          moda <- modelerFunction(ytrain, sign_data, sign_test, wei = NULL)
          preds <- moda$preds
          theta <- moda$theta
        } else  {
          moda <- modelerFunction(ytrain, rep(1, nrow(sign_data)), rep(1, nrow(sign_test)), wei = NULL)
          preds <- moda$preds
          theta <- moda$theta
        }
        if ( is.null(preds) ) {
          cv_results_all[[ k  ]]$preds[, j] <- NULL
          cv_results_all[[ k ]]$performances[j] <- NA
        } else {
          performance <- metricFunction(preds, ytest, theta)
          cv_results_all[[ k ]]$preds[, j] <- preds
          cv_results_all[[ k ]]$performances[j] <- performance
        }
      }  ##  end for ( i in 1:ntol ) {

    }  ## end for (k in 1:kfolds) {

    bbc_best_performance <- NULL

    if (B > 1) {
       n <- dim(x)[1]
       predictions <- cv_results_all[[ 1 ]]$preds
       for ( i in 2:kfolds )  predictions <- rbind(predictions, cv_results_all[[ i ]]$preds )
       bbc_best_performance <- gomp::bbc(predictions, y[unlist(folds)], metric = metricbbc, B = B )$bbc.perf
    }

    runtime <- proc.time() - tic
    perf <- matrix(0, nrow = kfolds, ncol = ntol)
    for (i in 1:kfolds)  perf[i, ] <- cv_results_all[[ i ]]$performances
    perf <- colMeans(perf, na.rm = TRUE)
    best_performance <- max(perf)
    best_configuration <- tol[ which.max(perf) ]
    list(cv_results_all = cv_results_all, best_performance = best_performance, best_configuration = best_configuration,
         bbc_best_performance = bbc_best_performance, runtime = runtime)
}







# auc (binary)
.auc.gomp <- function(predictions, ytest, theta = NULL) {
  ytest <- as.numeric( as.factor(ytest) )
  ri <- rank(predictions)
  up <- max(ytest)
  n <- length(predictions)
  n1 <- sum( ytest == up )
  n0 <- n - n1
  s1 <- sum( ri[ytest == up ] )
  ( s1 - 0.5 * ( n1 * (n1 + 1) ) ) / (n0 * n1)
}

# F score (binary)
.fscore.gomp <- function(predictions, ytest, theta = NULL) {
  ytest <- as.numeric( as.factor(ytest) ) - 1
  predictions <- round(predictions)
  tab <- table(ytest, predictions)
  prec <- tab[2, 2]/(tab[2, 2] + tab[1, 2])
  rec <- tab[2, 2] / (tab[2, 2] + tab[2, 1])
  2 * prec * rec / (prec + rec)
}

# precision (binary)
.prec.gomp <- function(predictions, ytest, theta = NULL) {
  ytest <- as.numeric( as.factor(ytest) ) - 1
  predictions <- round(predictions)
  tab <- table(ytest, predictions)
  tab[2, 2]/(tab[2, 2] + tab[1, 2])
}

# euclid_sens.spec score (binary)
.euclid_sens.spec.gomp <- function(predictions, ytest, theta = NULL) {
  ytest <- as.numeric( as.factor(ytest) ) - 1
  predictions <- round(predictions)
  tab <- table(ytest, predictions)
  spec <- tab[1, 1]/(tab[1, 1] + tab[1, 2])
  sens <- tab[2, 2] / (tab[2, 2] + tab[2, 1])
  sqrt( (1 - sens)^2 + (1 - spec)^2 )
}

# specificity (binary)
.spec.gomp <- function(predictions, ytest, theta = NULL) {
  ytest <- as.numeric( as.factor(ytest) ) - 1
  predictions <- round(predictions)
  tab <- table(ytest, predictions)
  tab[1, 1]/(tab[1, 1] + tab[1, 2])
}

# sensitivity or recall (binary)
.sens.gomp <- function(predictions, ytest, theta = NULL) {
  ytest <- as.numeric( as.factor(ytest) ) - 1
  predictions <- round(predictions)
  tab <- table(ytest, predictions)
  tab[2, 2] / (tab[2, 2] + tab[2, 1])
}

# accuracy (binary)
.acc.gomp <- function(predictions, ytest, theta = NULL) {
  ytest <- as.numeric( as.factor(ytest) ) - 1
  sum( (predictions > 0.5) == ytest ) / length(ytest)
}

#accuracy
.acc_multinom.gomp <- function(predictions, ytest, theta = NULL) {
  sum( predictions == ytest ) / length(ytest)
}

#mse lower values indicate better performance so we multiply with -1 in order to have higher values for better performances
.mse.gomp <- function(predictions, ytest, theta = NULL) {
  - sum( (predictions - ytest)^2 ) / length(ytest)
}

#mse lower values indicate better performance so we multiply with -1 in order to have higher values for better performances
.pve.gomp <- function(predictions, ytest, theta = NULL) {
  co <- length(ytest) - 1
  1 - sum( (predictions - ytest)^2 ) / ( co * Rfast::Var(ytest) )
}

#mean absolute error lower values indicate better performance so we multiply with -1 in order to have higher values for better performances
.ord_mae.gomp <- function(predictions, ytest, theta = NULL) {
  - sum( abs(as.numeric(predictions) - as.numeric(ytest)) ) / length(ytest)
}

#mean absolute error lower values indicate better performance so we multiply with -1 in order to have higher values for better performances
.mae.gomp <- function(predictions, ytest, theta = NULL) {
  - sum( abs(predictions - ytest) ) / length(ytest)
}

#cindex
.ci.gomp <- function(predictions, ytest, theta = NULL) {
  1 - Hmisc::rcorr.cens(ytest, predictions)[1]
  #survival::survConcordance(ytest ~ predictions)$concordance
}

#cindex for weibull, exponential and log-logistic regession
.ciwr.gomp <- function(predictions, ytest, theta = NULL) {
  Hmisc::rcorr.cens(ytest, predictions)[1]
  #1 - survival::survConcordance(ytest ~ predictions)$concordance
}

#Poisson deviance. Lower values indicate better performance so we multiply with -1 in order to have higher values for better performances
.poisdev.gomp <- function(predictions, ytest, theta = NULL) {
  - 2 * sum( ytest * log(ytest / predictions), na.rm = TRUE )
}

#Negative binomial deviance. Lower values indicate better performance so we multiply with -1 in order to have higher values for better performances
.nbdev.gomp <- function(predictions, ytest, theta) {
  - 2 * sum( ytest * log(ytest / predictions), na.rm = TRUE ) +
    2 * sum( ( ytest + theta ) * log( (ytest + theta) / (predictions + theta) ) )
}


#Binomial deviance. Lower values indicate better performance so we multiply with -1 in order to have higher values for better performances
.binomdev.gomp <- function(predictions, ytest, theta = NULL) {
  ya = ytest[, 1]     ;    N = ytest[, 2]
  yb = N - ya
  esta = predictions     ;    estb = N - esta
  - 2 * sum( ya * log(ya / esta), na.rm = TRUE ) - 2 * sum( yb * log(yb / estb), na.rm = TRUE )
}








# Modeling Functions

# input : ytrain, sign_data, sign_test
# output : preds


## binary logistic regression
.glm.gomp <- function(ytrain, sign_data, sign_test) {
  #using this variable x to overcome the structure naming problems when we have just one variable as a sign_data. For more on this contact athineo ;)
  x <- sign_data
  sign_model <- glm( ytrain ~ ., data = data.frame(x), family = binomial() );
  x <- sign_test
  preds <- predict( sign_model, newdata = data.frame(x), type = 'response' )
  #   preds[ preds>=0.5 ] = 1
  #   preds[ preds<0.5 ] = 0
  list(preds = preds, theta = NULL)
  #  }
}

## poisson regression
.pois.gomp <- function(ytrain, sign_data, sign_test) {
  #using this variable x to overcome the structure naming problems when we have just one variable as a sign_data. For more on this contact athineou ;)
  x <- sign_data
  sign_model <- glm( ytrain ~ ., data = data.frame(x), family = poisson() );
  x <- sign_test
  preds <-predict( sign_model, newdata = data.frame(x), type = 'response' )
  list(preds = preds, theta = NULL)
}

## binomial regression
.binom.gomp <-  function(ytrain, sign_data, sign_test){
  #using this variable x to overcome the structure naming problems when we have just one variable as a sign_data. For more on this contact athineou ;)
  y1 <- ytrain[, 1]
  N1 <- ytrain[, 2]
  x <- sign_data
  sign_model <- glm( y1 / N1 ~ ., data = data.frame(x), weights = N1, family = binomial );
  x <- sign_test
  preds <- predict( sign_model, newdata = data.frame(x), type = 'response' ) * N1
  list(preds = preds, theta = NULL)
}

## negative binomial regression
.nb.gomp <- function(ytrain, sign_data, sign_test) {
  #using this variable x to overcome the structure naming problems when we have just one variable as a sign_data. For more on this contact athineou ;)
  x <- sign_data
  sign_model <- MASS::glm.nb( ytrain ~ ., data = data.frame(x) );
  x <- sign_test
  preds <- predict( sign_model, newdata = data.frame(x), type = 'response' )
  list(preds = preds, theta = sign_model$theta)
}

## multinomial regression
.multinom.gomp <- function(ytrain, sign_data, sign_test) {
  #using this variable x to overcome the structure naming problems when we have just one variable as a sign_data. For more on this contact athineou ;)
  x <- sign_data
  sign_model <- nnet::multinom( ytrain ~ ., data = data.frame(x), trace = FALSE );
  x <- sign_test
  preds <- predict( sign_model, newdata = data.frame(x) )
  list(preds = preds, theta = NULL)
}

## oridnal regression
.ordinal.gomp <- function(ytrain, sign_data, sign_test) {
  x <- sign_data
  sign_model <- ordinal::clm( ytrain ~ ., data = data.frame(x), trace = FALSE );
  x <- sign_test
  preds <- predict( sign_model, newdata = data.frame(x), type = "class" )$fit
  list(preds = preds, theta = NULL)
}

## linear regression
.lm.gomp <- function(ytrain, sign_data, sign_test) { ## used for univariate and multivariate target in classical regression
  x <- sign_data
  sign_model <- lm( ytrain ~ ., data = data.frame(x) );
  x <- sign_test
  preds <- predict( sign_model, newdata = data.frame(x) )
  preds <- list(preds = preds, theta = NULL)
}

## quantile (median) regression
.rq.gomp <- function(ytrain, sign_data, sign_test) { ## used for univariate and multivariate target in classical regression
  x <- sign_data
  sign_model <- quantreg::rq( ytrain ~ ., data = data.frame(x));
  x <- sign_test
  preds <- predict( sign_model, newdata = data.frame(x) )
  list(preds = preds, theta = NULL)
}

## robust linear regression
.lmrob.gomp <- function(ytrain, sign_data, sign_test) { ## used for univariate and multivariate target in classical regression
  x <- sign_data
  sign_model <- MASS::rlm( ytrain ~ ., data = data.frame(x), maxit = 2000, method = "MM" );
  x <- sign_test
  preds <- predict( sign_model, newdata = data.frame(x) )
  list(preds = preds, theta = NULL)
}

## beta regression
.beta.gomp <- function(ytrain, sign_data, sign_test) { ## used for univariate and multivariate target in classical regression
  preds <- .beta.mod( ytrain, x = sign_data, xnew = sign_test )$est
  preds <- log( preds / (1 - preds) )  ## logit transformation to make it comparable with the normal regression
  list(preds = preds, theta = NULL)
}

## cox regression
.coxph.gomp <- function(ytrain, sign_data, sign_test) {
  x <- sign_data
  sign_model <- survival::coxph(ytrain~., data = data.frame(x))
  x <- sign_test
  preds <- predict(sign_model, newdata = data.frame(x), type = "risk")
  list(preds = preds, theta = NULL)
}

## weibull regression
.weibreg.gomp <- function(ytrain, sign_data, sign_test) {
  x <- sign_data
  sign_model <- survival::survreg(ytrain~., data = data.frame(x))
  x <- sign_test
  preds <- predict(sign_model, newdata = data.frame(x) )
  list(preds = preds, theta = NULL)
}

## exponential regression
.exporeg.gomp <- function(ytrain, sign_data, sign_test) {
  x <- sign_data
  sign_model <- survreg(ytrain~., data = data.frame(x), dist = "exponential")
  x <- sign_test
  x <- model.frame
  preds <- predict(sign_model, newdata = data.frame(x) )
  list(preds = preds, theta = NULL)
}

## log-logistic regression
.llrreg.gomp <- function(ytrain, sign_data, sign_test) {
  x <- sign_data
  sign_model <- survival::survreg(ytrain~., data = data.frame(x), dist = "loglogistic")
  x <- sign_test
  preds <- predict(sign_model, newdata = data.frame(x) )
  list(preds = preds, theta = NULL)
}



.beta.mod <- function(y, x, xnew = NULL) {

  regbeta <- function(pa, ly, sly1, x, n) {
    phi <- exp(pa[1])    ;    b <- pa[-1]
    m <- exp( tcrossprod(b, x) )
    m <- m / ( 1 + m )
    a1 <- m * phi   ;   a2 <- phi - a1
    - n * lgamma(phi) + sum( lgamma(a1) ) + sum( lgamma(a2) ) - sum(a1 * ly) - phi * sly1
  }

  n <- length(y)
  if ( is.null(x) ) {
    mod <- Rfast::beta.mle(y)
    param <- mod$param
    res <- list(be = param, phi = sum(param), loglik = mod$loglik, est = param[1]/param[2])

  } else {
    x <- model.matrix(y~ ., data.frame(x) )
    iniphi <- log( sum( y * (1 - y) ) / Rfast::Var(y) / n )

    ly1 <- log(1 - y)     ;    ly <- log(y) - ly1
    sly1 <- sum(ly1)           ;    sly2 <- sum( log(y) ) + sly1
    mod1 <- nlm(regbeta, c( iniphi, numeric(dim(x)[2]) ), ly = ly, sly1 = sly1, x = x, n = n, iterlim = 10000 )
    mod2 <- optim(mod1$estimate, regbeta, ly = ly, sly1 = sly1, x = x, n = n, control = list(maxit = 10000), hessian = TRUE  )

    seb <- sqrt( diag( solve(mod2$hessian) ) )
    be <- cbind(mod2$par[-1], seb[-1] )
    be <- cbind(be, (be[, 1] / be[, 2] )^2 )
    pv <- pchisq(be[, 3], 1, lower.tail = FALSE)
    be <- cbind(be, pv)
    rownames(be) <- colnames(x)
    colnames(be) <- c("coefficient", "std error", "chisq-statistic", "p-value")

    if ( !is.null(xnew) ) {
      xnew <- model.matrix(~., data.frame(xnew) )
      est <- exp( as.vector( xnew %*% be[, 1] ) )
      est <- est / (1 + est)
    } else  est <- NULL
    res <- list(be = be, phi = exp(mod2$par[1]), loglik = - mod2$value - sly2, est = est)
  }
  res
}












