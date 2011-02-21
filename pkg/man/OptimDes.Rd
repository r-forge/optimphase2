\name{OptimDes}
\alias{OptimDes}

\title{Construct an Optimal Two-stage Phase II Design with Time-to-event
Endpoints Evaluated at a Pre-specified Time}

\description{
  Construct an optimal two-stage phase II design with a time-to-event endpoint evaluated 
at a pre-specified time (e.g., 6 months) comparing treatment versus either a historical control 
rate with possible stopping for futility (single-arm), or an active control arm with possible 
stopping for both futility and superiority (two-arm), after the end of Stage I utilizing 
time to event data. The design  minimizes either the expected
duration of accrual (EDA), expected sample size (ES), or the expected total study length
(ETSL). Requires 2-3 minutes on a 2GHz laptop.
}

\usage{
OptimDes(
B.init, m.init, alpha, beta, param, x, target = c("EDA", "ETSL","ES"), 
sf=c("futility","OF","Pocock"), num.arm,r=0.5, recover=TRUE, 
control = OptimDesControl(),...)

}

\arguments{
  \item{B.init}{A vector of user-specified time points (B1, ..., Bb) that determine a
    set of time intervals with uniform accrual.}
  \item{m.init}{The projected
    number of patients that can be accrued within the time intervals
determined by \code{B.init}.  A large number of potential patients
results in long execution times for \code{OptimDes}, so
unrealistically large values should not be entered.}
  \item{alpha}{Type I error.}
  \item{beta}{Type II error.}
  \item{param}{Events should be defined as poor outcomes (e.g. death, progression).  Computations 
   and reporting are based on the proportion without an event at a 
   pre-specified time, \code{x}.  For constructing an optimal design, 
   complete event-free distributions at all times must be specified for 
   the control condition (Null), and for
   the alternative "effective" treatment.  Weibull distributions are currently
   implemented.  \code{param} is a vector of length 4:  (shape null, 
   scale null, shape alternative, scale alternative).  The R parameterization of the Weibull
distribution is used.}
  \item{x}{Pre-specified time for the event-free endpoint (e.g., 1 year).}
  \item{target}{The expected duration of
    accrual (EDA) is minimized with \code{target="EDA"}, the
    expected total study length is minimized with 
    \code{target="ETSL"}, or the expected sample size with \code{target="ES"}.}
  \item{sf}{Spending function for \code{alpha} at the end of Stage 1. There are three
    types of spending functions: no efficacy stopping with \code{sf="futility"}, 
    O'Brien-Fleming boundaries with \code{sf="OF"}, and Pocock boundaries with 
    \code{sf="Pocock"}.}  
  \item{num.arm}{Number of arms: a single-arm design with \code{num.arm=1}, or a 
    randomized two-arm design with \code{num.arm=2}.}
  \item{r}{Proportion of patients randomized to the treatment arm when \code{num.arm=2}.
    By default, r=0.5.}
  \item{recover}{The Simon Two-Stage design recovers alpha from the Type I errors
that would have occured if the study had not been stopped early for
futility.  This recovery allows for a slightly lower boundary for the
final test, thereby increasing power.  The default with partial data
at the interim analysis is to recover the alpha as with the standard
Simon design.  When recover=F, the optimal design is found assuming
that the boundary for the final test will not be adjusted downward
ignoring this small potential gain to ensure protection of the alpha
level even when the accrual rates and survivor functions are mis-specified. 
}
  \item{control}{An optional list of control settings.  See
\code{\link{OptimDesControl}.}
for the parameters that can be set and their default values.}
   \item{\dots}{No additional optional parameters are currently implemented.}   
}

\details{

\code{OptimDes} finds a two-stage phase II design with a time to event endpoint
with potential stopping after the first stage. 

For single arm designs, it implements the Case and Morgan (2003) and Huang, Talukder
and Thomas (2010) generalizaton of the Simon
(1989) two-stage design for comparing a treatment to a known standard rate with 
possible stopping for futility at the interim.

For randomized two-arm comparative designs, it allows an early stopping for both 
futility and superiority. The spending function for superiority can be chosen with
argument \code{sf}.

The design minimizes either the expected duration of accrual (EDA), expected
sample size (ES), or
the expected total study length (ETSL).

The design calculations assume Weibull distributions for the event-free
endpoint in the treatment group, and for the (assumed known, "Null") control
distribution.  The function \code{weibPmatch} can be used to select 
Weibull parameters that yield a target event-free rate at a 
specified time. Estimation is  based on the Kaplan-Meier or
  Nelson-Aalen estimators evaluated at a target time (e.g., 1 year).
The full treatment and control distributions and the accrual
distribution  affect power (and alpha level in some settings), see Case and Morgan (2003)).

Accrual rates are specified by the user.  These rates can differ
across time intervals specified by the user (this generalizes the
results in Case and Morgan).

A package vignette as user manual can be found in
the /doc subdirectory of the \code{OptimPhase2} package.  It can be accessed from the
HTML help page for the package.

}

\value{
  A list with components:
  \item{target}{The optimizaton target ("EDA" or "ETSL" or "ES").}
  \item{sf}{The alpha spending function ("futility" or "OF" or "Pocock").}
  \item{test}{A vector giving the type I error \code{alpha}, type II
    error \code{beta}, \code{recover=T/F}, Weibull parameters \code{param} and survival time
    of interest \code{x}.}
  \item{design}{A vector giving the number of study arms \code{num.arm} and 
    treatment randomization rate \code{r}.}
  \item{accrual}{A list containing the input vectors \code{B.init} and
   \code{m.init}.}
  \item{result}{A three-element vector containing the expected duration of
    accrual (EDA), the expected total study length (ETSL), and the expected sample size of the optimal
    design (ES).}
  \item{n}{A two-element vector containing the sample size for the
    interim analysis and the sample size if both stages are completed.}
  \item{stageTime}{A two-element vector giving the times for the interim and
    final analysis.}
  \item{boundary}{A three-element vector giving the rejection cutpoints \code{C1L} and \code{C1U}
    for the interim test statistic (see \code{Test2stage}), and  \code{C2}
    for the final test statistic.}
  \item{se}{A vector of length 4 with the asymptotic standard errors
    at the iterim and final analysis under the null hypothesis, followed
    by the corresponding SEs under the alternative hypothesis.  These SEs
    must be divided by the square root of sample size.}
  \item{u}{A two-element vector giving means of Z1 and Z2 under H1.  See detailed description.  
  It is used to compute conditional power.}
  \item{exposure}{The expected total exposure of patients at the time
of the planned interim analysis.  Patient exposure is truncated by both
the interim analysis time and the target surival time (i.e., no
exposure after \code{x}).  Exposure is a vector of length 2.  The
first value is the expected exposure without the Case-Morgan finite sample
adjustment to sample sizes.  The second value is the expected exposure with the adjustment applied.}
  \item{all.info}{A data frame containing the results
   for all of the evaluated sample sizes. }
  \item{single.stage}{A six-element vector giving the sample
    size \code{fix.n}, duration of accrual \code{DA},  study length
    \code{SL}, and corresponding values based on the exact binomial test for a one-arm single-stage design 
    and the Fisher exact test for a two-arm single-stage design with the design distributional assumptions. }

}

\author{Bo Huang \email{<bo.huang@pfizer.com>} and Neal Thomas \email{<neal.thomas@pfizer.com>}}


\references{

  Huang B., Talukder E. and Thomas N. (2010). Optimal two-stage Phase II designs with long-term
  endpoints. \emph{Statistics in Biopharmaceutical Research}, \bold{2}, 51--61. 
  
  Case M. D. and Morgan T. M. (2003). Design of Phase II cancer trials
  evaluating survival probabilities. \emph{BMC Medical Research
  Methodology}, \bold{3}, 7.

  Lin D. Y., Shen L., Ying Z. and Breslow N. E. (1996). Group seqential
  designs for monitoring survival probabilities. \emph{Biometrics},
  \bold{52}, 1033--1042.

  Simon R. (1989). Optimal two-stage designs for phase II clinical
  trials. \emph{Controlled Clinical Trials}, \bold{10}, 1--10.
  
}

\note{
  The algorithm will search for the optimal \code{n} between the 
  sample size for a single-stage design and the
  user specified maximum sample size \code{sum(m.init)}.

  When the length of \code{B.init} or \code{m.init} is 1, the accrual
  rate is constant as  in Lin et al. (1996), Case and Morgan
  (2003).
  
}
\seealso{\code{\link{np.OptimDes}}, \code{\link{print.OptimDes}},
  \code{\link{plot.OptimDes}}, \code{\link{weibPmatch}}}
\examples{
B.init <- c(1, 2, 3, 4, 5)
m.init <- c(15, 20, 25, 20, 15)
alpha <- 0.05
beta <- 0.1
param <- c(1, 1.09, 2, 1.40)
x <- 1

# H0: S0=0.40 H1: S1=0.60

object1 <- OptimDes(B.init,m.init,alpha,beta,param,x,target="EDA",sf="futility",num.arm=1)
print(object1)


m.init <- 6*c(15, 20, 25, 20, 15)
object2 <- OptimDes(B.init,m.init,alpha,beta,param,x,target="EDA",sf="futility",num.arm=2)
print(object2)

object3 <- OptimDes(B.init,m.init,alpha,beta,param,x,target="EDA",sf="OF",num.arm=2)
print(object3)


}

\keyword{design}
\keyword{optimize}

