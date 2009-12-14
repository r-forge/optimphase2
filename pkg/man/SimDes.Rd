\name{SimDes}
\alias{SimDes}

\title{Simulation Studies for Two-Stage Designs with Time to Event
  Endpoints from function OptimDes}

\description{
  Simulation experiments to compare the
  alpha level and power of two-stage designs from function
\code{\link{OptimDes}} with the targetted values.
}

\usage{
SimDes(object,B.init,m.init,weib0,weib1,interimRule='e1',
       sim.n=1000,e1conv=1/365,CMadj=F,attainI=1,attainT=1,FixDes="F", 
       Rseed)
}

\arguments{
  \item{object}{Output object of function \code{\link{OptimDes}}.}
  \item{B.init}{A vector of user-specified time points (B1, ..., Bb) that determine a
    set of time intervals with uniform accrual.  This vector 
needs to be specified only if its values differ from the call to
\code{OptimDes}.}
  \item{m.init}{The projected
    number of patients that can be accrued within the time intervals
determined by \code{B.init}.  This vector 
needs to be specified only if its values differ from the call to
\code{OptimDes}.}
  \item{weib0}{A vector with the shape and scale for the Weibull distribution under the
null hypothesis.  These need to be specified only if they differ from
the input to \code{OptimDes}.}
\item{weib1}{A vector with the shape and scale for the Weibull distribution under the
alternative hypothesis.  These need to be specified only if they differ from
the input to \code{OptimDes}.}
  \item{interimRule}{The interim analysis is performed when the planned
  \code{n1} patients are accrued regardless of the time required when
interimRule='n1'.  The interim analysis is performed at the planned
time \code{t1} regardless of the number of patients accrued when
interimRule='t1'.  The interim analysis is performed when the
truncated (by x) total exposure matches the total expected exposure 
when interimRule='e1'.  The default is 'e1'.}
  \item{sim.n}{The number of simulation replications.}
  \item{e1conv}{Convergence criteria for matching the truncated
exposure when \code{interimRule='e1'}.  The default is 1/365, which
is appropriate provided \code{B.init} is specified in years}
  \item{CMadj}{If true, the \code{n}, \code{n1}, and \code{t1} are adjusted by
the ratio of the exact binomial to asymptotic normal sample size for
the single stage design, as in Case and Morgan (2003).  Proportional
adjustment of times and sample sizes are made even if the accrual rates
are not constant. The adjustment to the \code{mda} is made through the
adjustment to \code{n} rather than by multiplication to ensure consistency
with accrual boundaries.  The truncated exposure time is matched to the
adjusted time of the interim analysis. Default is false.
}
 \item{attainI}{Samples sizes and times of the interim and final
analyses often differ from the exact targetted values for operational
reasons. The \code{attainI} permits simulations with a
different interim time or sample size (depending on
\code{interimRule}) by a specified fraction.}
 \item{attainT}{Simulations with a
 total sample size (assuming the trial does not stop based on
the interim analysis) that differs from the planned total  by a specified fraction.}
\item{FixDes}{If FixDes="E" or "N", a fixed design is simulated with the 
sample size determined by the Exact or Normal approximation.  
All other options for modifying the simuations are ignored.  
The alpha level and power based on an exact test and the normal 
approximation are returned.  
All other output variables are 0.  The default is "F"}
 \item{Rseed}{Optional integer for input to function \code{set.seed}. If unspecified, the random seed status at the time of the function call is used.}
}
\details{
  \code{sim.n}(default is \code{1000}) simulation experiments are
  conducted to assess how close the empirical type I and II errors come
  to the target values. 

  Simulation studies can also be used to assess the performance of the
  optimal design under mis-specification of the design parameters. For
  example,  if the Weibull
  shape and scale parameters of the time to event distributions are
changed, or if  the accrual rates are changed. (see Case and Morgan,
2003, for discussion of this topic).
  
The function \code{weibPmatch} can be used to select 
Weibull parameters that yield a target event-free rate at a 
specified time.
}

\value{
  A vector with:

\item{alphaExact}{Estimated alpha level using an exact test for the second stage test.}
\item{alphaNorm}{Estimated alpha level using an approximate normal test at the second stage test.}
\item{powerExact}{Estimated power using an exact test for the second stage test.}
\item{powerNorm}{Estimated power using an approximate normal test at the second stage test.} 
\item{eda}{Estimated mean duration of accrual under the null hypothesis.}
\item{etsl}{Estimated mean total study length under the null hypothesis.}
\item{es}{Estimated mean total sample size under the null hypothesis.}
\item{pstopNull}{The proportion of trials stopped at the interim analysis under the null hypothesis.}
\item{pstopAlt}{The proportion of trials stopped at the interim analysis under the alternative hypothesis.}
\item{aveE}{Average total (truncated at x) exposure at time of interim
analysis.}
\item{pinfoNull}{The proportion of the total information obtained at
the interim analysis under the null hypothesis.}
\item{pinfoAlt}{The proportion of the total information obtained at
the interim analysis under the althernative hypothesis.}
\item{n1}{Average sample size at interim.}
\item{t1}{Average time at interim.}
\item{phatKl}{Lowest interim survival rate continued to final analysis 
             based on the normal approximation.} 
\item{phatKh}{Highest interim survival rate resulting in terimination              based on the normal approximation.}
\item{phatRl}{Lowest final survival rate rejecting null 
              based on the normal approximation.}
\item{phatRh}{Highest final survival rate without rejecting null
              based on the normal approximation.} 
}

\author{Bo Huang \email{<bo.huang@pfizer.com>} and Neal Thomas \email{<neal.thomas@pfizer.com>}}

\references{

  Huang B., Talukder E. and Thomas N. Optimal two-stage Phase II designs with long-term
  endpoints. \emph{In press}. 
  
  Case M. D. and Morgan T. M. (2003) Design of Phase II cancer trials
  evaluating survival probabilities. \emph{BMC Medical Research
  Methodology}, \bold{3}, 7.

  Lin D. Y., Shen L., Ying Z. and Breslow N. E. (1996) Group seqential
  designs for monitoring survival probabilities. \emph{Biometrics},
  \bold{52}, 1033--1042.

  Simon R. (1989) Optimal two-stage designs for phase II clinical
  trials. \emph{Controlled Clinical Trials}, \bold{10}, 1--10.
  
}

\seealso{\code{\link{OptimDes}}, \code{\link{Test2stage}},
         \code{\link{weibPmatch}}}
\examples{
B.init <- c(1, 2, 3, 4, 5)
m.init <- c(15, 20, 25, 20, 15)
alpha <- 0.05
beta <- 0.1
param <- c(1, 1.09, 2, 1.40)
x <- 1

# H0: S0=0.40 H1: S1=0.60

object1 <- OptimDes(B.init,m.init,alpha,beta,param,x,target="EDA")

SimDes(object1)

SimDes(object1,interimRule='t1')

### accrual rates differ from planned

SimDes(object1,m.init=c(5,5,25,25,25))
}

\keyword{iteration}
\keyword{optimize}
