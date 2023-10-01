\name{minb}
\alias{minb}
\title{Multiple-Inflated Negative Binomial Model}
\description{
The minb is the main function to achieve data-driven inflated values selection and identification of important predictors for multiple inflated count data modeling.
}
\usage{
minb(X,y,pars_init=NULL,lambda1_set=NULL,lambda2_set=NULL,ntune=10,
maxiter=200,tol=1e-03,vrbs=FALSE)
}
\arguments{
\item{X}{The design matrix, without an intercept, minb includes an intercept by default.}
\item{y}{The response count vector.}
\item{pars_init}{an optional list containing the initial values for the corresponding component. See details.}
\item{lambda1_set}{A user supplied tuning sequence for inflated values. Typical usage is to have the program compute its own sequence based on \code{ntune}. Supplying a value of \code{lambda1_set} overrides this.}
\item{lambda2_set}{A user supplied tuning sequence for regression predictors. Typical usage is to have the program compute its own sequence based on \code{ntune}. Supplying a value of \code{lambda2_set} overrides this.}
\item{ntune}{The number of the tuning parameter values, default to \code{10}.}
\item{maxiter}{defines the maximum number of iterations, default to \code{200}.}
\item{tol}{specifies the iteration convergency tolerance, default to \code{1e-03}.}
\item{vrbs}{a boolean variable which defines whether to print the iteration details, default to \code{FALSE}.}
}

\details{
Initialization values can be supplied or estimated by glm.fit (the default). To supply initial values of the parameters, start should be a list with elements \code{"kappa"}, \code{"omega"},\code{"beta"}, and \code{"phi"} containing the starting values for the coefficients of the corresponding component of the model.
}
\value{
The \code{minb} returns  a list containing the following components:
\item{omega}{The vector of the estimated mixing proportions of the selected inflated values in the multiple-inflated negative binomial model.}
\item{kappa}{The vector of selected inflated values.}
\item{phi}{The value of the dispersion parameter of Negative Binomial distribution.}
\item{beta}{The vector of estimated non-zero regression coefficients of the Negative Binomial distribution.}

}
\examples{
# This is an example of minb with simulated data
data(minb_SimuData)
X = minb_SimuData$X
y = minb_SimuData$y
result = minb(X=X,y=y,ntune=2)
result$beta
result$omega
}
