\name{OptimDesControl}
\alias{OptimDesControl}

\title{Set parameters controlling numerical methods for OptimDes and np.OptimDes}
\description{
  Set parameters controlling numerical methods for \code{\link{OptimDes}} and \code{\link{np.OptimDes}}}

\usage{
OptimDesControl(trace=FALSE,tol=0.01,conv=1e-6,rho.int=c(0,1),aboveMin=1.05)
}

\arguments{
  \item{trace}{A logical value indicating if a
    trace of the iteration
    progress should be printed.  Default is \code{FALSE}.  If
    \code{TRUE} the sample size \code{n}, the corresponding optimal
    correlation \code{rho} and minimized \code{EDA}(\code{ETSL}) are
    printed at the conclusion of each iteration.}
  \item{tol}{\code{tol} is the
    desired accuracy for \code{optimize}. Default is 0.01.}
  \item{conv}{The desired convergence accuracy for \code{C1} and
    \code{C2}. Default is 1e-6.}
  \item{rho.int}{A vector containing
    the end-points of the interval to
    be searched for the optimal \code{rho}. Default is (0, 1).}
  \item{aboveMin}{The minimization method searches by increasing order
   from the single-stage
   sample size until the criteria exceeds the current minimum by a 
   multiplicative factor of \code{aboveMin}.  The search will also 
   terminate if \code{sum(m.init)} is reached.}
}

\author{Bo Huang \email{<bo.huang@pfizer.com>} and Neal Thomas \email{<neal.thomas@pfizer.com>}}

\references{

  Huang B., Talukder E. and Thomas N. Optimal two-stage Phase II designs with long-term
  endpoints. \emph{Statistics in Biopharmaceutical Research}, \bold{2(1)}, 51--61. 
  
  Case M. D. and Morgan T. M. (2003) Design of Phase II cancer trials
  evaluating survival probabilities. \emph{BMC Medical Research
  Methodology}, \bold{3}, 7.

  Lin D. Y., Shen L., Ying Z. and Breslow N. E. (1996) Group seqential
  designs for monitoring survival probabilities. \emph{Biometrics},
  \bold{52}, 1033--1042.

  Simon R. (1989) Optimal two-stage designs for phase II clinical
  trials. \emph{Controlled Clinical Trials}, \bold{10}, 1--10.
  
}


\seealso{\code{\link{OptimDes}}, \code{\link{np.OptimDes}}}

\keyword{design}
\keyword{optimize}
