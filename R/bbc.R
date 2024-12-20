bbc <- function(predictions, y, metric = "auc.gomp", conf = 0.95, B = 1000) {

  dm <- dim(predictions)
  n <- dm[1]    ;    p <- dm[2]
  out.perf <- numeric(B)

  if (metric == "auc.gomp") {
    y <- as.numeric( as.factor(y) ) - 1
    for (i in 1:B) {
      ind <- Rfast2::Sample.int(n, n, replace = TRUE)
      in.perf <- Rfast::colaucs(y[ind], predictions[ind, ])
      poio <- which.max(in.perf)
      out.perf[i] <- Rfast::auc(y[-ind], predictions[-ind, poio])
    }

  } else if (metric == "fscore.gomp") {
    y <- as.numeric( as.factor(y) ) - 1
    in.perf <- numeric(p)
    for (i in 1:B) {
      ind <- Rfast2::Sample.int(n, n, replace = TRUE)
      for (j in 1:p) {
        tab <- table(y[ind], predictions[ind, j])
	  if ( sum( dim(tab) ) < 4 ) {
	    in.perf[j] <- 0.5
	  } else {
          prec <- tab[2, 2]/(tab[2, 2] + tab[1, 2])
          rec <- tab[2, 2] / (tab[2, 2] + tab[2, 1])
          in.perf[j] <- prec * rec / (prec + rec)
        }
      }
      poio <- which.max(in.perf)
      tab <- table(y[-ind], predictions[-ind, poio])
      if ( sum( dim(tab) ) < 4 ) {
        out.perf[i] <- NA
      } else {
        prec <- tab[2, 2]/(tab[2, 2] + tab[1, 2])
        rec <- tab[2, 2] / (tab[2, 2] + tab[2, 1])
        out.perf[i] <- 2 * prec * rec / (prec + rec)
      }
    }

  } else if (metric == "prec.gomp") {
    y <- as.numeric( as.factor(y) ) - 1
    in.perf <- numeric(p)
    for (i in 1:B) {
      ind <- Rfast2::Sample.int(n, n, replace = TRUE)
      for (j in 1:p) {
        tab <- table(y[ind], predictions[ind, j])
	  if ( sum( dim(tab) ) < 4 ) {
	    in.perf[j] <- 0.5
	  } else  in.perf[j] <- tab[2, 2]/(tab[2, 2] + tab[1, 2])
      }
      poio <- which.max(in.perf)
      tab <- table(y[-ind], predictions[-ind, poio])
      if ( sum( dim(tab) ) < 4 ) {
        out.perf[i] <- NA
      } else {
        out.perf[i] <- tab[2, 2]/(tab[2, 2] + tab[1, 2])
      }
    }

  } else if (metric == "euclid_sens.spec.gomp") {
    y <- as.numeric( as.factor(y) ) - 1
    in.perf <- numeric(p)
    for (i in 1:B) {
      ind <- Rfast2::Sample.int(n, n, replace = TRUE)
      for (j in 1:p) {
        tab <- table(y[ind], predictions[ind, j])
	  if ( sum( dim(tab) ) < 4 ) {
	    in.perf[j] <- 0.5
	  } else {
          spec <- tab[1, 1]/(tab[1, 1] + tab[1, 2])
          sens <- tab[2, 2] / (tab[2, 2] + tab[2, 1])
          in.perf[j] <- sqrt( (1 - sens)^2 + (1 - spec)^2 )
	  }
      }
      poio <- which.max(in.perf)
      tab <- table(y[-ind], predictions[-ind, poio])
      if ( sum( dim(tab) ) < 4 ) {
        out.perf[i] <- NA
      } else {
        spec <- tab[1, 1]/(tab[1, 1] + tab[1, 2])
        sens <- tab[2, 2] / (tab[2, 2] + tab[2, 1])
        out.perf[i] <- sqrt( (1 - sens)^2 + (1 - spec)^2 )
      }
    }

  } else if (metric == "spec.gomp") {
    y <- as.numeric( as.factor(y) ) - 1
    in.perf <- numeric(p)
    for (i in 1:B) {
      ind <- Rfast2::Sample.int(n, n, replace = TRUE)
      for (j in 1:p) {
        tab <- table(y[ind], predictions[ind, j])
	  if ( sum( dim(tab) ) < 4 ) {
	    in.perf[j] <- 0.5
	  } else  in.perf[j] <- tab[1, 1]/(tab[1, 1] + tab[1, 2])
      }
      poio <- which.max(in.perf)
      tab <- table(y[-ind], predictions[-ind, poio])
      if ( sum( dim(tab) ) < 4 ) {
        out.perf[i] <- NA
      } else {
        out.perf[i] <- tab[1, 1]/(tab[1, 1] + tab[1, 2])
      }
    }

  } else if (metric == "sens.gomp") {
    y <- as.numeric( as.factor(y) ) - 1
    in.perf <- numeric(p)
    for (i in 1:B) {
      ind <- Rfast2::Sample.int(n, n, replace = TRUE)
      for (j in 1:p) {
        tab <- table(y[ind], predictions[ind, j])
	  if ( sum( dim(tab) ) < 4 ) {
	    in.perf[j] <- 0.5
	  } else  in.perf[j] <- tab[2, 2] / (tab[2, 2] + tab[2, 1])
      }
      poio <- which.max(in.perf)
      tab <- table(y[-ind], predictions[-ind, poio])
      if ( sum( dim(tab) ) < 4 ) {
        out.perf[i] <- NA
      } else {
        out.perf[i] <- tab[2, 2] / (tab[2, 2] + tab[2, 1])
      }
    }

  } else if (metric == "acc.gomp") {
    y <- as.numeric( as.factor(y) ) - 1
    for (i in 1:B) {
      ind <- Rfast2::Sample.int(n, n, replace = TRUE)
      in.perf <- Rfast::colmeans( predictions[ind, ] == y[ind] )
      poio <- which.max(in.perf)
      out.perf[i] <- mean( predictions[-ind, poio] == y[-ind] )
    }

  } else if (metric == "acc_multinom.gomp") {
    for (i in 1:B) {
      ind <- Rfast2::Sample.int(n, n, replace = TRUE)
      in.perf <- Rfast::colmeans( predictions[ind, ] == y[ind] )
      poio <- which.max(in.perf)
      out.perf[i] <- mean( predictions[-ind, poio] == y[-ind] )
    }

  } else if (metric == "mse.gomp") {
    for (i in 1:B) {
      ind <- Rfast2::Sample.int(n, n, replace = TRUE)
      in.perf <-  - Rfast::colmeans( (predictions[ind, ] - y[ind])^2 )
      poio <- which.max(in.perf)
      out.perf[i] <-  - mean( (predictions[-ind, poio] - y[-ind])^2 )
    }

  } else if (metric == "pve.gomp") {
    co <- length(y) - 1
    for (i in 1:B) {
      ind <- Rfast2::Sample.int(n, n, replace = TRUE)
      in.perf <- 1 - Rfast::colsums( (predictions[ind, ] - y[ind])^2 ) / ( co * Rfast::Var(y[ind]) )
      poio <- which.max(in.perf)
      out.perf[i] <-  1 - sum( (predictions[-ind, poio] - y[-ind])^2 ) / ( co * Rfast::Var(y[-ind]) )
    }

  } else if (metric == "ord_mae.gomp") {
    y <- as.numeric(y)
  	for (i in 1:dim(predictions)[2] ) 	predictions[, i] <- as.numeric(predictions[, i])
	  predictions <- as.matrix(predictions)
    for (i in 1:B) {
      ind <- Rfast2::Sample.int(n, n, replace = TRUE)
      in.perf <-  - Rfast::colmeans( abs(predictions[ind, ] - y[ind]) )
      poio <- which.max(in.perf)
      out.perf[i] <-  - mean( abs(predictions[-ind, poio] - y[-ind]) )
    }

  } else if (metric == "mae.gomp") {
    for (i in 1:B) {
      ind <- Rfast2::Sample.int(n, n, replace = TRUE)
      in.perf <-  - Rfast::colmeans( abs(predictions[ind, ] - y[ind]) )
      poio <- which.max(in.perf)
      out.perf[i] <-  - mean( abs(predictions[-ind, poio] - y[-ind]) )
    }

  } else if (metric == "ci.gomp") {
    in.perf <- numeric(p)
    for (i in 1:B) {
      ind <- Rfast2::Sample.int(n, n, replace = TRUE)
      #for (j in 1:p)   in.perf[j] <- survival::survConcordance(y[ind, ] ~ predictions[ind, j])$concordance
      for (j in 1:p)   in.perf[j] <- 1 - Hmisc::rcorr.cens(y[ind, ], predictions[ind, j])[1]
      poio <- which.max(in.perf)
      #out.perf[i] <- survival::survConcordance(y[-ind, ] ~ predictions[-ind, poio])$concordance
      out.perf[i] <- Hmisc::rcorr.cens(y[-ind, ], predictions[-ind, poio])[1]
    }

  } else if (metric == "ciwr.gomp") {
    in.perf <- numeric(p)
    for (i in 1:B) {
      ind <- Rfast2::Sample.int(n, n, replace = TRUE)
      #for (j in 1:p)   in.perf[j] <- 1 - survival::survConcordance(y[ind, ] ~ predictions[ind, j])$concordance
      for (j in 1:p)   in.perf[j] <- Hmisc::rcorr.cens(y[ind, ], predictions[ind, j])[1]
      poio <- which.max(in.perf)
      #out.perf[i] <- 1 - survival::survConcordance(y[-ind, ] ~ predictions[-ind, poio])$concordance
      out.perf[i] <- Hmisc::rcorr.cens(y[-ind, ],predictions[-ind, poio])[1]
    }

  } else if (metric == "poisdev.gomp") {
    in.perf <- numeric(p)
    for (i in 1:B) {
      ind <- Rfast2::Sample.int(n, n, replace = TRUE)
      in.perf <-  - colSums( y[ind] * log(y[ind] / predictions[ind, ]), na.rm = TRUE )
      poio <- which.max(in.perf)
      out.perf[i] <-  - 2 * sum( y[-ind] * log(y[-ind] / predictions[-ind, poio]), na.rm = TRUE )
    }

  } else if (metric == "binomdev.gomp") {
    in.perf <- numeric(p)
    ya <- y[, 1]     ;    N <- y[, 2]
    yb <- N - ya
    esta <- predictions     ;    estb <- N - esta
    for (i in 1:B) {
      ind <- Rfast2::Sample.int(n, n, replace = TRUE)
      in.perf <-  - colSums( ya[ind] * log(ya[ind] / esta[ind, ]), na.rm = TRUE ) - colSums( yb[ind] * log(yb[ind] / estb[ind, ]), na.rm = TRUE )
      poio <- which.max(in.perf)
      out.perf[i] <-  - 2 * sum( ya[-ind] * log(ya[-ind] / esta[-ind, poio]), na.rm = TRUE ) - 2 * sum( yb[-ind] * log(yb[-ind] / estb[-ind, poio]), na.rm = TRUE )
    }

  } ## end  if (metric == " ") {
  a <- 0.5 * (1 - conf )
  ci <- quantile(out.perf, probs = c(a, 1 - a), na.rm = TRUE )
  list(out.perf = out.perf, bbc.perf = mean(out.perf, na.rm = TRUE), ci = ci )

}

