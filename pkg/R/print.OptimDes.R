"print.OptimDes" <-
function(x,dig=3,all=FALSE,condPow=F,CMadj=F,...)
{
  alpha <- x$test[1]
  beta <- x$test[2]
  shape0 <- x$test[3]
  scale0 <- x$test[4]
  shape1 <- x$test[5]
  scale1 <- x$test[6]
  xx <- x$test[7]
  r <- x$design[2]
  num.arm <- x$design[1]
  t1<-  x$stageTime["t1"]
  C1<-  x$boundary["C1L"]
  mtsl<- x$stageTime["MTSL"]
  se<- x$se
  Stime<-x$single.stageTime
  exposure<-x$exposure[1]

  l0<-lambda(shape0,scale0,xx)
  l1<-lambda(shape1,scale1,xx)

  if(CMadj==T & num.arm!=1)stop("The CM adjustment is for single-arm trials only.")

  ### apply adjustment to exact binomial fixed design size
  ### per CM 2003
  if (CMadj){
        CMadjfac=Stime["n0E"]/Stime["n0"]
        ### handle etsl and  mtsl separately
        etsl<-CMadjfac*(t1+(1-pnorm(C1))*(mtsl-xx-t1)) + (1-pnorm(C1))*xx
        mtsl<-CMadjfac*(mtsl-xx) + xx
        exposure<-x$exposure[2]
  }else CMadjfac<-1
  
  cat("             Optimal Design Results \n\n")
  if(num.arm==1)
  cat("    One-Arm Study \n\n")
  if(num.arm==2)
  cat("    Two-Arm Study:", r/(1-r),": 1 randomization \n\n")
  
  cat("    H0: S0=S1=",round(s(shape0,scale0,xx),dig),"H1: S1=",round(s(shape1,scale1,xx),dig),"\n\n")
  cat("    Type I error(1-sided upper):",alpha,"type II error:",beta,"\n")
  cat(paste("    Recover alpha: ",as.logical(x$test["recover"]),"\n",sep=""))
  cat("    Event-free endpoint time:",xx," \n\n")
  
  cat("                  target:",x$target,"\n")
  prt1 <- CMadjfac*x$result
  if(CMadj)prt1[2]<-etsl   ### special adjustment required
  prt1<-round(prt1,dig)
  names(prt1) <- c("EDA","ETSL","           ES")
  print(prt1,...)

  cat("\n")
  cat("          Sample Size at Each Stage \n")
  prt2 <- ceiling(CMadjfac*x$n)
  names(prt2) <- c("n1","          n.last")
  print(prt2,...)

  cat("\n")
  cat("          Study time at Each Stage \n")
  prt3 <- round(c(CMadjfac*x$stageTime[1],mtsl),dig)
  names(prt3) <- c("t1","             MTSL")
  print(prt3,...)

  cat("\n")
  cat("    Projected patient exposure at interim analysis: ",
      round(exposure,2))
  cat("\n\n")
  cat("    Proportion of the total information at the interim analysis:\n")
  if(num.arm==1)
  prt3b<-round(c((se[2]/se[1])^2,(se[4]/se[3])^2),3)
  if(num.arm==2)
  {
    sig10 <- se[1]
    sig20 <- se[2]
    sig11 <- se[3]
    sig21 <- se[4]
    v10 <- (sig10^2)/((1-r)*(l0)^2)
    v11 <- (sig11^2)/(r*(l1)^2)
    v20 <- (sig20^2)/((1-r)*(l0)^2)
    v21 <- (sig21^2)/(r*(l1)^2)
    
    prt3b<-round(c((sig20/sig10)^2,(v20+v21)/(v10+v11)),3)
  }  
  
  names(prt3b)<-c("Under NuLL","    Under Alternative")
  print(prt3b,...)

  cat("\n\n")
  cat("          Hypothesis Test Boundaries \n")
  prt4 <- round(x$boundary,dig)
  names(prt4) <- c("C1L","               C1U","               C2")
  print(prt4,...)

  cat("\n\n")
  cat("    Approximate Rates Corresponding to Test Boundaries* \n\n")
  ### compute boundaries using the null and alternative SE's
  if(num.arm==1) {
  se10<-x$se[1]/sqrt(CMadjfac*x$n[2])  
  se20<-x$se[2]/sqrt(CMadjfac*x$n[2])
  se11<-x$se[3]/sqrt(CMadjfac*x$n[2])
  se21<-x$se[4]/sqrt(CMadjfac*x$n[2])
  C1<- x$boundary[1]
  C2<- x$boundary[3]
  p10<-exp(log(l0)-C1*se10/l0)
  p10<-exp(-p10)
  p11<-exp(log(l0)-C1*se11/l1)
  p11<-exp(-p11)
  p1<-(p10+p11)/2

  p20<-exp(log(l0)-C2*se20/l0)
  p20<-exp(-p20)
  p21<-exp(log(l0)-C2*se21/l1)
  p21<-exp(-p21)
  p2<-(p20+p21)/2
  cat("          Event-free rate for C1: ", 
      round(p1,dig),"\n")

  cat("          Event-free rate for C2: ", 
      round(p2,dig),"\n\n")
  }
  
  if(num.arm==2) {
  vnull1 <- (sig10^2)/((1-r)*r*(l0^2))
  vnull2 <- (sig20^2)/((1-r)*r*(l0^2))
  valt1  <- v10+v11
  valt2  <- v20+v21
  C1<- x$boundary[1] 
  C1U<- x$boundary[2]
  C2<- x$boundary[3]
  
  p10<-exp(log(l0)-C1*sqrt(vnull1/(x$n[2])))
  d10<-exp(-p10) - exp(-l0)
  p11<-exp(log(l0)-C1*sqrt(valt1/(x$n[2])))
  d11<-exp(-p11) - exp(-l0)
  d1<-(d10+d11)/2
  
  pu0<-exp(log(l0)-C1U*sqrt(vnull1/(x$n[2])))
  du0<-exp(-pu0) - exp(-l0)
  pu1<-exp(log(l0)-C1U*sqrt(valt1/(x$n[2])))
  du1<-exp(-pu1) - exp(-l0)
  d1u<-(du0+du1)/2
  
  p20<-exp(log(l0)-C2*sqrt(vnull2/(x$n[2])))
  d20<-exp(-p20) - exp(-l0)
  p21<-exp(log(l0)-C2*sqrt(valt2/(x$n[2])))
  d21<-exp(-p21) - exp(-l0)
  d2<-(d20+d21)/2
  cat("          Difference in Event-free rate for C1L: ", 
      round(d1,dig),"\n")
  cat("          Difference in Event-free rate for C1U: ", 
      round(d1u,dig),"\n")
  cat("          Difference in Event-free rate for C2: ", 
      round(d2,dig),"\n\n")
  }

  if(condPow){
     cat("\n")
     rho<-x$se[4]/x$se[3]
     # P(Z2>C2|Z1=C1,H1) based on bivariate normal distributions
     if(num.arm==1){
          condpow<-1-pnorm(C2,mean=x$u[2] + rho*(C1-x$u[1]),sd=sqrt(1-rho^2))
          cat("          Conditional power at interim test boundary under H1: ",
         round(condpow,dig), "\n\n")
          }
     else {
          condpow1<-1-pnorm(C2,mean=x$u[2] + rho*(C1-x$u[1]),sd=sqrt(1-rho^2))
          condpow1u<-1-pnorm(C2,mean=x$u[2] + rho*(C1U-x$u[1]),sd=sqrt(1-rho^2))
          
         cat("          Conditional power at interim test boundary C1L under H1: ",
         round(condpow1,dig), "\n\n")
         cat("          Conditional power at interim test boundary C1U under H1: ",
         round(condpow1u,dig), "\n\n")
         }
     ### compute boundaries using the null and alternative SE's
  }


  cat("\n")
  if(CMadj)comment="(Exact binomial calculation)"
  else comment="(Asymptotic normal calculation)"
  cat(paste("    Single-stage Design ",comment,sep=""), "\n")
  prt5 <- round(Stime[c(1:3)],dig)
  if(CMadj)prt5 <- round(Stime[c(4:6)],dig)
  names(prt5) <- c("      Single stage N","DA","          SL")
  print(prt5,...)



  if(all==TRUE )
  {
     
    all.info<-x$all.info
    etsl<-CMadjfac*(all.info[,"t1"]+(pnorm(all.info[,"C1U"])-
          pnorm(all.info[,"C1L"]))*(all.info[,"MTSL"]-xx-all.info[,"t1"])) 
    etsl<- etsl + (pnorm(all.info[,"C1U"])-
          pnorm(all.info[,"C1L"]))*xx
    mtsl<-CMadjfac*(all.info[,"MTSL"]-xx) + xx
    all.info[,7]<-mtsl
    all.info[,9]<-etsl
    all.info[,c(1:2,6,8,10)]<-CMadjfac*all.info[,c(1:2,6,8,10)]
    all.info[,1]<-ceiling(all.info[,1])
    all.info<-all.info[,c(1:2,6:10)]
    cat("\n")
    cat("   Design Parameters for Each n \n")
    print(round(all.info,dig))
    
  }

  cat("\n\n *Note:  Rates corresponding to test boundaries are a function \n",
      "of the non-parametric SE computed at the time of the analyses. \n", 
      "The approximate rates are based on the asymptotic SE computed  \n",
      "under the null and alternative hypotheses.\n")

  if(CMadj){
     cat("\n Note:  All sample sizes and times are adjusted by the exact ",
            "binomial correction \n","factor: ",Stime["n0E"],"/",
            Stime["n0"],"\n",sep="")
  }

  return(invisible())
}  
