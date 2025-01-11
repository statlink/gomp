gomp <- function(y, x, xstand = TRUE, tol = qchisq(0.95, 1), test = "logistic", method = "ar2") {

  if ( xstand )  x <- Rfast::standardise(x)

  if ( test == "normal" | test == "cor" ) {
    tic <- proc.time()
    res <- Rfast::ompr(y, x, xstand = FALSE, method = method, tol = tol)
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = NULL, res = res$info)

  } else if ( test == "logistic" ) {
    tic <- proc.time()
    res <- Rfast::omp(y, x, xstand = FALSE, tol, type = "logistic")
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = NULL, res = res$info)

  } else if ( test == "poisson" ) {
    tic <- proc.time()
    res <- Rfast::omp(y, x, xstand = FALSE, tol = tol, type = "poisson")
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = NULL, res = res$info)

  } else if ( test == "qpoisson" ) {
    tic <- proc.time()
    res <- Rfast::omp(y, x, xstand = FALSE, tol = tol, type = "quasipoisson")
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = res$phi, res = res$info)

  } else if (test == "mvreg") {
    tic <- proc.time()
    res <- Rfast::omp(y, x, xstand = FALSE, tol = tol, type = "mv")
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = NULL, res = res$info)

  } else if ( test == "qlogistic" ) {
    tic <- proc.time()
    res <- Rfast::omp(y, x, xstand = FALSE, tol = tol, type = "quasibinomial")
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = res$phi, res = res$info)

  } else if ( test == "normlog" ) {
    tic <- proc.time()
    res <- Rfast::omp(y, x, xstand = FALSE, tol = tol, type = "normlog")
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = res$phi, res = res$info)

  } else if ( test == "gamma" ) {
    tic <- proc.time()
    res <- Rfast::omp(y, x, xstand = FALSE, tol = tol, type = "gamma")
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = res$phi, res = res$info)

  } else if ( test == "multinomial" ) {
    tic <- proc.time()
    res <- Rfast::omp(y, x, xstand = FALSE, tol = tol, type = "multinomial")
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = res$phi, res = res$info)

  } else {
    d <- dim(x)[2]
    ind <- 1:d
    can <- which( is.na( Rfast::colsums(x) ) )
    ind[can] <- 0

    if ( test == "negbin" ) {
      tic <- proc.time()
     	mod <- MASS::glm.nb(y ~ 1)
    	rho <-  - 2 * as.numeric( logLik(mod) )
    	res <-  y - fitted(mod)
	    ela <- as.vector( cov(res, x) )
	    sel <- which.max( abs(ela) )
      sela <- sel
      names(sela) <- NULL
      mod <- MASS::glm.nb( y ~ x[, sela], control = list(epsilon = 1e-07, maxit = 50, trace = FALSE) )
      res <-  mod$residuals
      rho[2] <-  - 2 * as.numeric( logLik(mod) )
      ind[sel] <- 0
      i <- 2
      while ( (rho[i - 1] - rho[i]) > tol ) {
        r <- rep(NA, d)
        i <- i + 1
        ## r[ind] <- Rfast::colsums( x[, ind, drop = FALSE] * res )
        r[ind] <- Rfast::eachcol.apply(x, res, indices = ind[ind > 0 ], oper = "*", apply = "sum")
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- MASS::glm.nb( y ~ x[, sela], control = list(epsilon = 1e-07, maxit = 50, trace = FALSE) )
        res <- y - fitted(mod)
        rho[i] <-  - 2 * as.numeric( logLik(mod) )
        ind[sela] <- 0
      } ## end while ( (rho[i - 1] - rho[i]) > tol )
      runtime <- proc.time() - tic
      len <- length(sela)
      res <- cbind(c(0, sela[-len]), rho[1:len])
      colnames(res) <- c("Selected Vars", "Deviance")
      result <- list(runtime = runtime, phi = NULL, res = res)

    } else if ( test == "beta" ) {
      tic <- proc.time()
    	mod <- Rfast::beta.mle(y)
      rho <-  - 2 * mod$loglik
	    res <- y - mod$param[1]/sum(mod$param)
      ela <- as.vector( cov(res, x) )
      sel <- which.max( abs(ela) )
      sela <- sel
      names(sela) <- NULL
      mod <- try( .beta.reg(y, x[, sela]), silent = TRUE )
      if ( identical( class(mod), "try-error" ) ) {
        rho[2] <- rho[1]
      } else {
        est <- exp( - mod$be[1] - x[, sela] * mod$be[2] )
        res <- y - 1 / (1 + est)
        rho[2] <-  - 2 * mod$loglik
        ind[sel] <- 0
      }
      i <- 2
      while ( (rho[i - 1] - rho[i]) > tol ) {
        r <- rep(NA, d)
        i <- i + 1
        ##r[ind] <- Rfast::colsums( x[, ind, drop = FALSE] * res )
        r[ind] <- Rfast::eachcol.apply(x, res, indices = ind[ind > 0 ], oper = "*", apply = "sum")
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- try( .beta.reg(y, x[, sela]), silent = TRUE )
        if ( identical( class(mod), "try-error" ) ) {
          rho[i] <- rho[i - 1]
        } else {
          est <- exp( - mod$be[1] - x[, sela] %*% mod$be[-1] )
          res <- y - 1 / (1 + est)
          rho[i] <-  - 2 * mod$loglik
          ind[sela] <- 0
        }
      } ## end while ( (rho[i - 1] - rho[i]) > tol )
      runtime <- proc.time() - tic
      len <- length(sela)
      res <- cbind(c(0, sela[-len]), rho[1:len])
      colnames(res) <- c("Selected Vars", "Deviance")
      result <- list(runtime = runtime, phi = NULL, res = res)

    } else if ( test == "mm" ) {
      tic <- proc.time()
      mod <- MASS::rlm(y ~ 1, method = "MM", maxit = 2000)
      rho <-  - 2 * as.numeric( logLik(mod) )
      res <- mod$residuals
      ela <- as.vector( cov(res, x) )
      sel <- which.max( abs(ela) )
      sela <- sel
      names(sela) <- NULL
      mod <- MASS::rlm(y ~ x[, sela], method = "MM", maxit = 2000 )
      res <- mod$residuals
      rho[2] <-  - 2 * as.numeric( logLik(mod) )
      ind[sel] <- 0
      i <- 2
      while ( (rho[i - 1] - rho[i]) > tol ) {
        r <- rep(NA, d)
        i <- i + 1
        ##r[ind] <- Rfast::colsums( x[, ind, drop = FALSE] * res )
        r[ind] <- Rfast::eachcol.apply(x, res, indices = ind[ind > 0 ], oper = "*", apply = "sum")
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- try( MASS::rlm(y ~ x[, sela], method = "MM", maxit = 2000 ), silent = TRUE )
        if ( identical( class(mod), "try-error" ) ) {
          rho[i] <- rho[i - 1]
        } else {
          res <- mod$residuals
          rho[i] <-  - 2 * as.numeric( logLik(mod) )
          ind[sela] <- 0
        }
      } ## end while ( (rho[i - 1] - rho[i]) > tol )
      runtime <- proc.time() - tic
      len <- length(sela)
      res <- cbind(c(0, sela[-len]), rho[1:len])
      colnames(res) <- c("Selected Vars", "Deviance")
      result <- list(runtime = runtime, phi = NULL, res = res)

    } else if ( test == "quantreg") {
      tic <- proc.time()
      mod <- quantreg::rq(y ~ 1)
      rho <-  - 2 * as.numeric( logLik(mod) )
      res <- mod$residuals
      ela <- as.vector( cov(res, x) )
      sel <- which.max( abs(ela) )
      sela <- sel
      names(sela) <- NULL
      mod <- quantreg::rq(y ~ x[, sela])
      res <- mod$residuals
      rho[2] <-  - 2 * as.numeric( logLik(mod) )
      ind[sel] <- 0
      i <- 2
      while ( (rho[i - 1] - rho[i]) > tol ) {
        r <- rep(NA, d)
        i <- i + 1
        ##r[ind] <- Rfast::colsums( x[, ind, drop = FALSE] * res )
        r[ind] <- Rfast::eachcol.apply(x, res, indices = ind[ind > 0 ], oper = "*", apply = "sum")
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- try( quantreg::rq(y ~ x[, sela]), silent = TRUE )
        if ( identical( class(mod), "try-error" ) ) {
          rho[i] <- rho[i - 1]
        } else {
          res <- mod$residuals
          rho[i] <-  - 2 * as.numeric( logLik(mod) )
          ind[sela] <- 0
        }
      } ## end while ( (rho[i - 1] - rho[i]) > tol )
      runtime <- proc.time() - tic
      len <- length(sela)
      res <- cbind(c(0, sela[-len]), rho[1:len])
      colnames(res) <- c("Selected Vars", "Deviance")
      result <- list(runtime = runtime, phi = NULL, res = res)

    } else if ( test == "ordinal" ) {
      tic <- proc.time()
      mod <- MASS::polr(y ~ 1)
      rho <-  - 2 * as.numeric( logLik(mod) )
      res <- .ord.resid(y, mod$fitted.values)
      ela <- as.vector( cov(res, x) )
      sel <- which.max( abs(ela) )
      sela <- sel
      names(sela) <- NULL
      mod <- MASS::polr(y ~ x[, sela])
      res <- .ord.resid(y, mod$fitted.values)
      rho[2] <-  - 2 * as.numeric( logLik(mod) )
      ind[sel] <- 0
      i <- 2
      while ( (rho[i - 1] - rho[i]) > tol ) {
        r <- rep(NA, d)
        i <- i + 1
        ##r[ind] <- Rfast::colsums( x[, ind, drop = FALSE] * res )
        r[ind] <- Rfast::eachcol.apply(x, res, indices = ind[ind > 0 ], oper = "*", apply = "sum")
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- try( MASS::polr(y ~ x[, sela]), silent = TRUE)
        if ( identical( class(mod), "try-error" ) ) {
          rho[i] <- rho[i - 1]
        } else {
          res <- .ord.resid(y, mod$fitted.values)
          rho[i] <-  - 2 * as.numeric( logLik(mod) )
          ind[sela] <- 0
        }
      } ## end while ( (rho[i - 1] - rho[i]) > tol )
      runtime <- proc.time() - tic
      len <- length(sela)
      res <- cbind(c(0, sela[-len]), rho[1:len])
      colnames(res) <- c("Selected Vars", "Deviance")
      result <- list(runtime = runtime, phi = NULL, res = res)

    } else if ( test == "tobit" ) {
      tic <- proc.time()
      mod <- survival::survreg(y ~ 1, dist = "gaussian")
      rho <-  - 2 * as.numeric( logLik(mod) )
	    res <- resid(mod)
      ela <- as.vector( cov(res, x) )
      sel <- which.max( abs(ela) )
      sela <- sel
      names(sela) <- NULL
      mod <- survival::survreg(y ~ x[, sela], dist = "gaussian" )
	    res <- resid(mod)
      rho[2] <-  - 2 * as.numeric( logLik(mod) )
      ind[sel] <- 0
      i <- 2
      while ( (rho[i - 1] - rho[i]) > tol ) {
        r <- rep(NA, d)
        i <- i + 1
        ##r[ind] <- Rfast::colsums( x[, ind, drop = FALSE] * res )
        r[ind] <- Rfast::eachcol.apply(x, res, indices = ind[ind > 0 ], oper = "*", apply = "sum")
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- try( survival::survreg(y ~ x[, sela], dist = "gaussian" ), silent = TRUE )
        if ( identical( class(mod), "try-error" ) ) {
          rho[i] <- rho[i - 1]
        } else {
          res <- resid(mod)
          rho[i] <-  - 2 * as.numeric( logLik(mod) )
          ind[sela] <- 0
        }
      } ## end while ( (rho[i - 1] - rho[i]) > tol )
      runtime <- proc.time() - tic
      len <- length(sela)
      res <- cbind(c(0, sela[-len]), rho[1:len])
      colnames(res) <- c("Selected Vars", "Deviance")
      result <- list(runtime = runtime, phi = NULL, res = res)

    } else if ( test == "cox" ) {
      tic <- proc.time()
      mod <- survival::coxph(y ~ 1)
      rho <-  - 2 * summary( mod) [[1]]
      res <- mod$residuals   ## martingale residuals
      ela <- as.vector( cov(res, x) )
      sel <- which.max( abs(ela) )
      sela <- sel
      names(sela) <- NULL
      mod <- survival::coxph(y ~ x[, sela] )
      res <-  mod$residuals
      rho[2] <-  - 2 * as.numeric( logLik(mod) )
      ind[sel] <- 0
      i <- 2
      while ( (rho[i - 1] - rho[i]) > tol ) {
        r <- rep(NA, d)
        i <- i + 1
        ##r[ind] <- Rfast::colsums( x[, ind, drop = FALSE] * res )
        r[ind] <- Rfast::eachcol.apply(x, res, indices = ind[ind > 0 ], oper = "*", apply = "sum")
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- try( survival::coxph(y ~ x[, sela] ), silent = TRUE )
        if ( identical( class(mod), "try-error" ) ) {
          rho[i] <- rho[i - 1]
        } else {
          res <- mod$residuals   ## martingale residuals
          rho[i] <-  - 2 * as.numeric( logLik(mod) )
          ind[sela] <- 0
        }
      } ## end while ( (rho[i - 1] - rho[i]) > tol )
      runtime <- proc.time() - tic
      len <- length(sela)
      res <- cbind(c(0, sela[-len]), rho[1:len])
      colnames(res) <- c("Selected Vars", "Deviance")
      result <- list(runtime = runtime, phi = NULL, res = res)

    } else if ( test == "weibull" ) {
      tic <- proc.time()
      mod <- survival::survreg(y ~ 1)
      rho <-  - 2 * as.numeric( logLik(mod) )
      res <- resid(mod)
      ela <- as.vector( cov(res, x) )
      sel <- which.max( abs(ela) )
      sela <- sel
      names(sela) <- NULL
      mod <- survival::survreg(y ~ x[, sela], control = list(iter.max = 5000) )
      res <- resid(mod)
      rho[2] <-  - 2 * as.numeric( logLik(mod) )
      if ( is.na(rho[2]) )  rho[2] <- rho[1]
      ind[sel] <- 0
      i <- 2
      while ( (rho[i - 1] - rho[i]) > tol ) {
        r <- rep(NA, d)
        i <- i + 1
        ##r[ind] <- Rfast::colsums( x[, ind, drop = FALSE] * res )
        r[ind] <- Rfast::eachcol.apply(x, res, indices = ind[ind > 0 ], oper = "*", apply = "sum")
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- try( survival::survreg(y ~ x[, sela], control = list(iter.max = 5000) ), silent = TRUE )
        if ( identical( class(mod), "try-error" ) ) {
          rho[i] <- rho[i - 1]
        } else {
          res <- resid(mod)
          rho[i] <-  - 2 * as.numeric( logLik(mod) )
          ind[sela] <- 0
        }
      } ## end while ( (rho[i - 1] - rho[i]) > tol )
      runtime <- proc.time() - tic
      len <- length(sela)
      res <- cbind(c(0, sela[-len]), rho[1:len])
      colnames(res) <- c("Selected Vars", "Deviance")
      result <- list(runtime = runtime, phi = NULL, res = res)

    } else if ( test == "log-logistic" ) {
      tic <- proc.time()
      mod <- survival::survreg(y ~ 1, dist = "loglogistic")
      rho <-  - 2 * as.numeric( logLik(mod) )
      res <- resid(mod)
      ela <- as.vector( cov(res, x) )
      sel <- which.max( abs(ela) )
      sela <- sel
      names(sela) <- NULL
      mod <- survival::survreg(y ~ x[, sela], control = list(iter.max = 5000), dist = "loglogistic" )
      res <- resid(mod)
      rho[2] <-  - 2 * as.numeric( logLik(mod) )
      if ( is.na(rho[2]) )  rho[2] <- rho[1]
      ind[sel] <- 0
      i <- 2
      while ( (rho[i - 1] - rho[i]) > tol ) {
        r <- rep(NA, d)
        i <- i + 1
        ##r[ind] <- Rfast::colsums( x[, ind, drop = FALSE] * res )
        r[ind] <- Rfast::eachcol.apply(x, res, indices = ind[ind > 0 ], oper = "*", apply = "sum")
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- try( survival::survreg(y ~ x[, sela], control = list(iter.max = 5000), dist = "loglogistic" ), silent = TRUE )
        if ( identical( class(mod), "try-error" ) ) {
          rho[i] <- rho[i - 1]
        } else {
          res <- resid(mod)
          rho[i] <-  - 2 * as.numeric( logLik(mod) )
          ind[sela] <- 0
        }
      } ## end while ( (rho[i - 1] - rho[i]) > tol )
      runtime <- proc.time() - tic
      len <- length(sela)
      res <- cbind(c(0, sela[-len]), rho[1:len])
      colnames(res) <- c("Selected Vars", "Deviance")
      result <- list(runtime = runtime, phi = NULL, res = res)
    } ##  end if (test == "censIndWR")

  }  ##  end if ( test == "testIndReg" | test == "testIndfisher" )

  result
}


.beta.reg <- function(y, x) {

  regbeta <- function(pa, ly, sly1, x, n) {
    phi <- exp(pa[1])    ;    b <- pa[-1]
    m <- exp( tcrossprod(b, x) )
    m <- m / ( 1 + m )
    a1 <- m * phi   ;   a2 <- phi - a1
    - n * lgamma(phi) + sum( lgamma(a1) ) + sum( lgamma(a2) ) - sum(a1 * ly) - phi * sly1
  }

  n <- length(y)
  if ( NCOL(x) == 0 ) {
    mod <- Rfast::beta.mle(y)
    res <- list(phi = sum(mod$param), be = mod$param[1]/mod$param[2], loglik = mod$loglik)

  } else {
    x <- model.matrix(y ~ ., data.frame(x) )
    iniphi <- log( sum( y * (1 - y) ) / Rfast::Var(y) / n )

    ly1 <- log(1 - y)     ;    ly <- log(y) - ly1
    sly1 <- sum(ly1)      ;    sly2 <- sum( log(y) ) + sly1
    mod1 <- nlm(regbeta, c( iniphi, numeric(dim(x)[2]) ), ly = ly, sly1 = sly1, x = x, n = n, iterlim = 10000 )
    mod2 <- nlm(regbeta, mod1$estimate, ly = ly, sly1 = sly1, x = x, n = n, iterlim = 10000 )
   res <-list(phi = exp(mod2$estimate[1]), be = mod2$estimate[-1], loglik = - mod2$minimum - sly2)
  }

  res
}


.ord.resid <- function(y, est) {
  y <- as.numeric(y)
  res <- y
  dm <- dim(est)
  n <- dm[1]
  d <- dm[2]
  ind1 <- which(y == 1)
  indd <- which(y == d)
  ind <- setdiff( 1:n, c(ind1, indd) )
  res[ind1] <- est[ind1, 1] - 1
  res[indd] <- 1 - est[indd, d]
  w <- 1 - Rfast::design_matrix(y, ones = FALSE)
  w <- (w * est)[ind, ]
  a <- numeric( dim(w)[1] )
  for ( i in 1:dim(w)[1] ) {
    a[i] <- sum(w[i, 1:(y[ind[i]] - 1)]) - sum(w[i, (y[ind[i]] + 1):d])
  }
  res[ind] <- a
  res
}






