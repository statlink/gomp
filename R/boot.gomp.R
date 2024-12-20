boot.gomp <- function(y, x, tol = qchisq(0.95, 1), test = "logistic", method = "ar2", B = 500, ncores = 1) {

  runtime <- proc.time()

  sel <- NULL
  n <- dim(x)[1]
  if ( !is.matrix(y) )  dim(y) <- c(n, 1)

  if ( ncores <= 1 ) {

    for (i in 1:B) {
      ina <- Rfast2::Sample.int(n, n, replace = TRUE)
      sel <- c(sel, gomp::gomp(y[ina, ], x[ina, ], tol = tol, test = test, method = method)$res[-1, 1])
    }

  } else {
    requireNamespace("doParallel", quietly = TRUE, warn.conflicts = FALSE)
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)

    mod <- foreach::foreach(i = 1:B, .packages = c("gomp", "Rfast2"), .export = c("gomp", "Sample.int") )  %dopar% {
      ina <- Rfast2::Sample.int(n, n, replace = TRUE)
      sel <- gomp::gomp(y[ina, ], x[ina, ], tol = tol, test = test, method = method)$res[-1, 1]
      return( sel )
    }
    parallel::stopCluster(cl)
    sel <- unlist(mod)
  }

  res <- cbind(unique(sel), Rfast::Table(sel)/B )
  colnames(res) <- c("Variable", "Selection proportion")
  runtime <- proc.time() - runtime
  list(runtime = runtime, res = res)
}

