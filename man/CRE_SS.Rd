% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/CRE_SS.R
\name{CRE_SS}
\alias{CRE_SS}
\title{A Sample Selection Model with Correlated Random Effects}
\usage{
CRE_SS(sel_form, out_form, id, data = NULL, par = NULL, killed_par = NULL,
  par_files = NULL, delta = 1, sigma = 1, gamma = 1, max_delta = 3,
  max_sigma = 3, max_gamma = 5, rho = 0, tau = 0, lower = c(rho = -1,
  tau = -1), upper = c(rho = 1, tau = 1), method = "L-BFGS-B", H = c(10,
  10), psnH = 20, prbH = 20, plnreH = 20, accu = 10000,
  reltol = sqrt(.Machine$double.eps), verbose = 0, tol_gtHg = Inf)
}
\arguments{
\item{sel_form}{Formula for selection equation, a probit model with random effects}

\item{out_form}{Formula for outcome equation, a Poisson model with random effects}

\item{id}{A vector that represents the identity of individuals, numeric or character}

\item{data}{Input data, a data frame}

\item{par}{Starting values for estimates}

\item{killed_par}{correlation parameters to swtich off}

\item{par_files}{Loading initial values from saved ProbitRE and PoissonRE estimates}

\item{delta}{Variance of random effects on the individual level for ProbitRE}

\item{sigma}{Variance of random effects on the individual level for PLN_RE}

\item{gamma}{Variance of random effects on the <individual, time> level for PLN_RE}

\item{max_delta}{Largest allowed initial delta}

\item{max_sigma}{Largest allowed initial sigma}

\item{max_gamma}{Largest allowed initial gamma}

\item{rho}{Correlation between random effects on the individual level}

\item{tau}{Correlation between error terms on the <individual, time> level}

\item{lower}{Lower bound for estiamtes}

\item{upper}{Upper bound for estimates}

\item{method}{Searching algorithm, don't change default unless you know what you are doing}

\item{H}{A vector of length 2, specifying the number of points for inner and outer Quadratures}

\item{psnH}{Number of Quadrature points for Poisson RE model}

\item{prbH}{Number of Quddrature points for Probit RE model}

\item{plnreH}{Number of Quddrature points for PLN_RE model}

\item{accu}{L-BFGS-B only, 1e12 for low accuracy; 1e7 for moderate accuracy; 10.0 for extremely high accuracy. See optim}

\item{reltol}{Relative convergence tolerance. default typically 1e-8}

\item{verbose}{Level of output during estimation. Lowest is 0.}

\item{tol_gtHg}{tolerance on gtHg, not informative for L-BFGS-B}
}
\value{
A list containing the results of the estimated model
}
\description{
Estimate a sample selection model in panel counting data, in which the selection equation is a Probit model with random effects on individuals, and the outcome equation is a Poisson Lognormal model with random effects on the same individuals. The random effects on the same individual and the error terms on the same <individual, time> dyad are both correlated across two equations.
}
\examples{
\donttest{
data(rt)
# Note: estimation may take up 10~15 minutes
est = CRE_SS(isRetweet~fans+tweets+as.factor(tweet.id),
                       num.words~fans+tweets+as.factor(tweet.id),
                       id=rt$user.id, data=rt)
}
}
\references{
1. Jing Peng and Christophe Van den Bulte. Participation vs. Effectiveness of Paid Endorsers in Social Advertising Campaigns: A Field Experiment. Working Paper.

2. Jing Peng and Christophe Van den Bulte. How to Better Target and Incent Paid Endorsers in Social Advertising Campaigns: A Field Experiment. In Proceedings of the 2015 International Conference on Information Systems.
}
\seealso{
Other PanelCount: \code{\link{CRE}}; \code{\link{PLN_RE}};
  \code{\link{PoissonRE}}; \code{\link{ProbitRE}}
}

