\name{plot.OptimDes}
\alias{plot.OptimDes}

\title{Plot efficiency of optimal two-stage designs as a function of
total sample size}

\description{
Output from function \code{OptimDes} is used to display the ETSL, ES and
EDA for a two-stage design relative to a single-stage design as a
function of the combined stage 1 and 2 sample size.
}
\usage{
\method{plot}{OptimDes}(x, xscale = "t", l.type = 1:4, l.col =
c("blue", "green", "purple", "red"), CMadj=F,...)
}

\arguments{
  \item{x}{Output from function \code{OptimDes}.}
  \item{xscale}{Scale of the x-axis. \code{"t"} for combined stage 1
and 2 study length.  \code{"n"} for combined sample size. Default is \code{t}.}
  \item{l.type}{Line types for the plot. Default is \code{1-4}.}
  \item{l.col}{Line colors for the plot. Default is \code{"blue"} for
    \code{ETSL}, \code{"green"} for \code{EDA}, \code{"purple"} for \code{ES} and \code{"red"} for
    \code{t1}.}
  \item{CMadj}{If true, the sample sizes and times are adjusted by
the ratio of the exact binomial to asymptotic normal sample size for
the single stage design, as in Case and Morgan 2003.  Proportional
adjustment of times and sample sizes are made even if the accrual rates
are not constant.  Default is
false.  
}
  \item{...}{Additional graphical parameters passed to function \code{plot}.}
}

\details{
The plot displays the tradeoff between ETSL, EDA and ES as a function of
the combined sample size.  Robustness of the optimal two-stage design to
deviations from the target sample size can be explored.  The plots
often suggest compromised designs achieving near-optimal results for
both EDA and ETSL may be a better choice.  Test boundary values (C1, C2), and numerical
values of other design parameters, can be obtained for a design
selected from the plots using function \code{\link{np.OptimDes}}.

The plot also includes the time of the interim analysis (\code{t1}) as a ratio to
the time for a corresonding single-stage analysis.
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

\seealso{\code{\link{print.OptimDes}}, \code{\link{OptimDes}}, \code{\link{np.OptimDes}}}

\keyword{hplot}
\keyword{optimize}
