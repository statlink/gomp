\name{Cross-validation for gamma-Orthogonal Matching Pursuit (gamma-OMP) algorithm}
\alias{cv.gomp}

\title{
Cross-validation for the gamma-Orthogonal Matching Pursuit (gamma-OMP) algorithm
}

\description{
The function performs a k-fold cross-validation for identifying the best tolerance values for the \eqn{\gamma}-Orthogonal Matching Pursuit (\eqn{\gamma}-OMP) algorithm.
}

\usage{
cv.gomp(y, x, kfolds = 10, folds = NULL, tol = seq(4, 9, by = 1),
task = "C", metric = NULL, metricbbc = NULL, modeler = NULL, test = NULL,
method = "ar2", B = 1)
}

\arguments{
\item{y}{
The response variable.
}
\item{x}{
A matrix with the predictor variables.
}
\item{kfolds}{
The number of the folds in the k-fold Cross Validation (integer).
}
\item{folds}{
The folds of the data to use. If NULL the folds are created internally with the same function.
}
\item{tol}{
A vector of tolerance values.
}
\item{task}{
A character ("C", "R" or "S"). It can be "C" for classification (logistic, multinomial or ordinal regression), "R" for regression (robust and non robust linear regression, median regression,
(zero inflated) poisson and negative binomial regression, beta regression), "S" for survival regresion (Cox, Weibull or exponential regression).
}
\item{metric}{
A metric function provided by the user. If NULL the following functions will be used: auc.gomp, mse.gomp, ci.gomp for classification, regression and survival analysis tasks, respectively. See details for more. If you know what you have put it here to avoid the function choosing somehting else.
\bold{Note} that you put these words as they are, without "".
}
\item{metricbbc}{
This is the same argument as "metric" with the difference that " " must be placed. If for example, metric = auc.mxm, here metricbbc = "auc.mxm". The same value must be given here. This argument is to be used with the function \code{\link{bbc}} which does bootstrap bias correction of the estimated performance (Tsamardinos, Greasidou and Borboudakis, 2018). This argument is valid if the last argument (B) is more than 1.
}
\item{modeler}{
A modeling function provided by the user. If NULL the following functions will be used: glm.gomp, lm.gomp, coxph.gomp for classification, regression and survival analysis tasks, respectively. See details for more. If you know what you have put it here to avoid the function choosing somehting else. \bold{Note} that you put these words as they are, without "".
}
\item{test}{
A function object that defines the conditional independence test used in the SES function (see also SES help page). If NULL, "cor", "logistic" and "cox" are used for classification, regression and survival analysis tasks, respectively. If you know what you have put it here to avoid the function choosing sometting else. Not all tests can be included here. "mv", "gamma", and "tobit" are anot available.
}
\item{method}{
This is only for the "cor". You can either specify, "ar2" for the adjusted R-square or "sse" for the sum of squares of errors. The tolerance value in both cases must a number between 0 and 1. That will denote a percentage. If the percentage increase or decrease is less than the nubmer the algorithm stops. An alternative is "BIC" for BIC and the tolerance values are like in all other regression models.
}
\item{B}{
How many bootstrap re-samples to draw. This argument is to be used with the function \code{\link{bbc}} which does bootstrap bias correction of the estimated performance (Tsamardinos, Greasidou and Borboudakis, 2018). If you have thousands of observations then this might not be necessary, as there is no optimistic bias to be corrected. What is the lower limit cannot be told beforehand however. SES and MMPC however were designed for the low sample cases, hence, bootstrap bias correction is perhaps a must thing to do.
}
}

\details{
For more details see also \code{\link{gomp}}.
}

\value{
A list including:
\item{cv_results_all}{
A list with predictions, performances and selected variables for each fold and each tolerance value. The elements are called
"preds", "performances" and "selectedVars".
}
\item{best_performance}{
A numeric value that represents the best average performance.
}
\item{best_configuration}{
A numeric value that represents the best tolerance value.
}
\item{bbc_best_performance}{
The bootstrap bias corrected best performance if B was more than 1, othwerwise this is NULL.
}
\item{runtime}{
The runtime of the cross-validation procedure.
}

Bear in mind that the values can be extracted with the $ symbol, i.e. this is an S3 class output.
}

\references{
Tsamardinos I., Greasidou E. and Borboudakis G. (2018).
Bootstrapping the out-of-sample predictions for efficient and accurate cross-validation.
Machine Learning, 107: 1895--1922.
\url{ https://link.springer.com/article/10.1007/s10994-018-5714-4 }

Tsagris M., Papadovasilakis Z., Lakiotaki K. and Tsamardinos, I. (2022).
The \eqn{\gamma}-OMP algorithm for feature selection with application to gene expression data.
IEEE/ACM Transactions on Computational Biology and Bioinformatics, 19(2): 1214--1224.
}

\author{
R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{gomp}, \link{gomp.path}, \link{bbc} }
}

\examples{
# simulate a dataset with continuous data
x <- matrix( rnorm(200 * 50), ncol = 50 )
# the target feature is the last column of the dataset as a vector
y <- x[, 50]
x <- x[, -50]
# run a 10 fold CV for the regression task
best_model <- cv.gomp(y, x, kfolds = 5, task = "R",
tol = seq(0.001, 0.01,by = 0.001), method = "ar2" )
}

