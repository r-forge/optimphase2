"SimDes" <- 
function(object,B.init,m.init,weib0,weib1,interimRule='e1',
         sim.n=1000,e1conv=1/365,CMadj=F,attainI=1,attainT=1,FixDes="F",
         Rseed)
{  

  interimRule<-match.arg(interimRule,c('e1','n1','t1'))
  if(interimRule!='n1' & interimRule!='t1' & interimRule!='e1')stop('Invalid interimRule')


  if(missing(B.init))B.init<-object$accrual$B.init
  if(missing(m.init))m.init<-object$accrual$m.init
  if(length(B.init)!=length(m.init))stop(
     "The length of the accrual rates and intervals must be equal")

  if(missing(weib0))weib0<-c(object$test["param1"],object$test["param2"])
  if(missing(weib1))weib1<-c(object$test["param3"],object$test["param4"])

  if(!missing(Rseed))set.seed(Rseed)

  if(FixDes=="E" | FixDes=="N")interimRule='t1'

  cat("\n The interimRule is ",interimRule,"\n\n")

  x<-object$test["x"]
  n1<-object$n[1]
  n<- object$n[2]
  t1 <- object$stageTime[1]
  C1 <- object$boundary[1]
  C1U <- object$boundary[2]
  C2 <- object$boundary[3]
  num.arm<-object$design["num.arm"]
  r<-object$design["r"]
  if(num.arm==2)ntrt<-ceiling(r*n)
  else ntrt<- n
  ncntl<- n-ntrt

  if(CMadj & num.arm!=1)stop("CM adjustment available with num.arm=1 only")

  if(length(B.init)!=length(m.init))
  {
    stop("projected patient times and numbers should be of equal length")
  }

  if(x>=max(B.init))
  {
    stop("the survival time of interest is beyond the projected accrual time")
  } 

  ### apply CM correction factor

  Stime<-object$single.stageTime

  ### apply adjustment to exact fixed design size
  ### per CM 2003
  if (CMadj){
        CMadjfac=Stime["n0E"]/Stime["n0"]
        exposure<-object$exposure[2]
  }else{
        CMadjfac<-1
        exposure<-object$exposure[1]
  }

  if(FixDes=="E"){n<-Stime["n0E"]
  }else if(FixDes=="N"){n<-Stime["n0"]
  }else{
      n<-ceiling(CMadjfac*attainT*n)
      ### the attained levels are applied at the beginning of each
      ### iteration because they may change based on interimRule
      n1<-ceiling(CMadjfac*n1)
      t1<-CMadjfac*t1
      if(ceiling(CMadjfac*n1*attainI)>=n)stop("Single-stage design is optimal")
  }


  shape0 <- weib0[1]
  scale0 <- weib0[2]
  shape1 <- weib1[1]
  scale1 <- weib1[2]  

  ### event rate under null hypothesis (not event-free as in OptimDes)
  p0<-pweibull(x,shape=shape0,scale=scale0)
  
  ### cutoff for number of events if exact binomial test is used at second stage for one-arm design
  if(num.arm==1) {
      cexact<-qbinom(object$test["alpha"],n,p0)       # one-arm
      if(pbinom(cexact,n,p0)>object$test["alpha"])cexact<-cexact-1
  }

  if(n>sum(m.init))stop("sample size n greater than total potential accrual")

  ### compute mda and distribution restricted to (0,mda)
  ### the CMadjfac is applied to mda through the adjustment to n
  mdaout<-compMDA(B.init,m.init,n)
  mda<-mdaout[1]
  l<-mdaout[2]

  B <- B.init[1:l]  # "l" not "one"
  B[l]<-mda
  M <- m.init[1:l]
  M[l]<- ifelse(l>1, n-sum(m.init[1:(l-1)]) , n )
  Pi<-M/sum(M)  # proportion in each accrual interval
  B<-c(0,B)     # add B0 for diff calculation


  ### compute truncated exposure conditional on Y<=t1

  EW<-exposure
  tlow<-0
  thigh<-max(B.init)+x

  
  alphaNorm <- 0
  alphaExact<- 0
  powNorm   <- 0
  powExact  <- 0
  eda<-0
  etsl<-0
  es<-0
  aveE<-0
  intvarNull<-0
  finvarNull<-0
  intvarAlt<-0
  finvarAlt<-0
  n1int  <- 0
  t1int <- 0
  difKl<- Inf
  difKh<- -Inf
  difRl<- Inf
  difRh<- -Inf
  lambda0<-lambda(shape0,scale0,x)
  pstopNull<-0
  pstopAlt<-0

  for(i in 1:sim.n)
  {
   
    randI <- sample(1:l,ntrt,replace=TRUE,prob=Pi) #sample accrual intervals
    u<-runif(ntrt)
    sample.Y1 <- sort( B[randI] + u*(B[randI+1]-B[randI]) )
    if(num.arm==1) {
        sample.Y0<-NULL
    }else{
        randI <- sample(1:l,ncntl,replace=TRUE,prob=Pi) #sample accrual intervals
        u<-runif(ncntl)
        sample.Y0 <- sort( B[randI] + u*(B[randI+1]-B[randI]) )
    }
    sample.Y<-c(sample.Y1,sample.Y0)
   
    
    ### compute time for interim to match planned exposure
    if(interimRule=='e1'){
        
       t1.root <- uniroot(ft,c(tlow,thigh),tol=e1conv,y=sample.Y,x=x,EW=EW)
       if(suppressWarnings(warning()!=""))
          stop("The algorithm of uniroot() did not converge for exposure time calculation")
       else{
          t1<-t1.root$root
          }
       t1<-t1*attainI
       if(t1<x+sample.Y1[1]){
          warning("Matched exposure yielded interim time less than x from enrollment")
          t1<-x+sample.Y1[1]
       }
       if(num.arm==2){
           if(t1<x+sample.Y0[1]){
              warning("Matched exposure yielded interim time less than x from enrollment")
              t1<-x+sample.Y0[1]
           }
       }
    }else if(interimRule=='t1'){
       t1<-object$stageTime[1]*CMadjfac*attainI
       if(t1<x+sample.Y1[1]){
          warning("Adjusted t1 occurred before time x from enrollment")
          t1<=x+sample.Y1[1]
       }
       if(num.arm==2){
           if(t1<x+sample.Y0[1]){
              warning("Adjusted t1 occurred before time x from enrollment")
              t1<=x+sample.Y0[1]
           }
       }
    }else n1<-ceiling(object$n[1]*CMadjfac*attainI)
       

    ### alpha level
    sample.T1 <- rweibull(ntrt,shape=shape0,scale=scale0) #null survival times
    if(num.arm==2)sample.T0 <- rweibull(ncntl,shape=shape0,scale=scale0) 

    if(interimRule=='t1' | interimRule=='e1'){
       ntrt1<-sum(sample.Y1<=t1)
       ncntl1<-0
       if(num.arm==2)ncntl1 <- sum(sample.Y0<=t1)
       n1<-ntrt1+ncntl1

       Y11 <- sample.Y1[sample.Y1<=t1]
       T11 <- sample.T1[sample.Y1<=t1]
       Y21 <- sample.Y1[sample.Y1>t1]
       T21 <- sample.T1[sample.Y1>t1]

       if(num.arm==2) {
           Y10 <- sample.Y0[sample.Y0<=t1]
           T10 <- sample.T0[sample.Y0<=t1]
           Y20 <- sample.Y0[sample.Y0>t1]
           T20 <- sample.T0[sample.Y0>t1]
       }else{
           Y10 <- NULL
           T10 <- NULL
           Y20 <- NULL
           T20 <- NULL
       }
    }else{
       ntrt1<-ceiling(ntrt*n1/n)
       ncntl1<-0
       if(num.arm==2)ncntl1<-ceiling(ncntl*n1/n)
       n1<-ntrt1+ncntl1

       if(x>sample.Y1[ntrt]-sample.Y1[1]){
          t1<-x+sample.Y1[1]
          ntrt1<-min(ntrt,sum(sample.Y1<t1)+1)
          warning("ntrt1 patients were accrued before time x from enrollment")
       }else  t1<-sample.Y1[ntrt1]
       if(num.arm==2){
           if(x>sample.Y0[ncntl]-sample.Y0[1]){
              t1<-max(t1,x+sample.Y0[1])
              ncntl1<-min(ncntl,sum(sample.Y0<t1)+1)
              warning("ncntl1 patients were accrued before time x from enrollment")
           }else  t1<-sample.Y0[ncntl1]
       }
       
       Y11 <- sample.Y1[1:ntrt1]
       T11 <- sample.T1[1:ntrt1]
       if(ntrt1<ntrt){
          Y21 <- sample.Y1[(ntrt1+1):ntrt]
          T21 <- sample.T1[(ntrt1+1):ntrt]
       }else{
          Y21<-NULL
          T21<-NULL
       }
       
       if(num.arm==2) {
           Y10 <- sample.Y0[1:ncntl1]
           T10 <- sample.T0[1:ncntl1]
           
           if(ncntl1<ncntl) {
               Y20 <- sample.Y0[(ncntl1+1):ncntl]
               T20 <- sample.T0[(ncntl1+1):ncntl]
           }else{
              Y20 <- NULL 
              T20 <- NULL
           }
       }else{
           Y10 <- NULL
           T10 <- NULL
           Y20 <- NULL
           T20 <- NULL
       }
    }

    
    ### compute interim and z-statistics
    ### the final z-statistic is computed even if the interim
    ### is less than c1.  the stopping rule is imposed below

     intout <- Test2stage(x, C1, C1U, C2, t1, num.arm, Y11, T11, p0, 
               Y21=NULL, T21=NULL, Y10, T10, Y20=NULL, T20=NULL, 
               printTest=FALSE)

     finout <- Test2stage(x, C1, C1U, C2, mda+x, num.arm, Y11, T11, p0, 
               Y21, T21,Y10, T10, Y20, T20,
                        printTest=FALSE)


    r.int<-intout["z"]
    r.fin<-finout["z"]
    
    
    
    yobs1<-sum(c(sample.T1)<=x)
    if(num.arm==2){
       yobs0<-sum(c(sample.T0)<=x)
       pv <- fisher.test(matrix(c(yobs1,yobs0,ntrt-yobs1,ncntl-yobs0),nr=2),
               alternative="less",conf.int=FALSE)$p.value
    }

    if(FixDes=="E" | FixDes=="N"){
       if(num.arm==1)alphaExact<- alphaExact + 1*(yobs1<=cexact)
       else alphaExact<- alphaExact + 1*(pv<object$test["alpha"])
       alphaNorm<-  alphaNorm  + 1*(r.fin>C2)
    }else{
        intY<-(r.int<=C1|r.int>=C1U)
        eda<-eda + (!intY)*max(sample.Y1,sample.Y0) + intY*max(Y11,Y10)
        etsl<-etsl + (!intY)*(max(sample.Y1,sample.Y0)+x) + intY*t1
        es<-es + intY*n1 + (!intY)*n
        aveE<-aveE + sum(pmin(x,pmax(0,t1-c(Y11,Y10))))
        pstopNull<-pstopNull + 1*intY*(n1<n)
  
        if(r.int>C1) {
           if(num.arm==1)alphaExact<- alphaExact + 1*(yobs1<=cexact)
           else alphaExact<- alphaExact + 1*(pv<object$test["alpha"]|r.int>=C1U)
        }


        dif1<-exp(-exp(intout["Log(cumL1)"])) - exp(-exp(intout["Log(cumL0)"]))
        dif2<-exp(-exp(finout["Log(cumL1)"])) - exp(-exp(finout["Log(cumL0)"]))
  
        if(r.int>C1&(r.fin>C2)){
          alphaNorm <- alphaNorm + 1
          difRl<-min(difRl,dif2)
          difKl<-min(difKl,dif1)
        }else if(r.int>C1){
          difRh<-max(difRh,dif2)
          difKl<-min(difKl,dif1)
        }else{
           difKh<-max(difKh,dif1)
        }
        ### proportion of info at interim
        intvarNull<-intvarNull+(intout["se"])^2
        finvarNull<-finvarNull+(finout["se"])^2
  
        n1int<-n1int+n1
        t1int<-t1int+t1
     } 
     
     
     

    ### power
    ### re-use patient start times and associated calculations
    sample.T1 <- rweibull(ntrt,shape=shape1,scale=scale1)

    if(interimRule=='t1' | interimRule=='e1'){
       ntrt1<-sum(sample.Y1<=t1)

       T11 <- sample.T1[sample.Y1<=t1]
       T21 <- sample.T1[sample.Y1>t1]

    }else{
       T11 <- sample.T1[1:ntrt1]
       if(ntrt1<ntrt) T21 <- sample.T1[(ntrt1+1):ntrt]
    }


    intout <- Test2stage(x, C1, C1U, C2, t1, num.arm, Y11, T11, p0, 
              Y21=NULL, T21=NULL, Y10, T10, Y20=NULL, T20=NULL, 
              printTest=FALSE)

    finout <- Test2stage(x, C1, C1U, C2, mda+x, num.arm, Y11, T11, p0, 
              Y21, T21,Y10, T10, Y20, T20,
                       printTest=FALSE)

    r.int<-intout["z"]
    r.fin<-finout["z"]
 
   yobs1<-sum(sample.T1<=x)
   if(num.arm==2){
       yobs0<-sum(sample.T0<=x)
       pv <- fisher.test(matrix(c(yobs1,yobs0,ntrt-yobs1,ncntl-yobs0),nr=2),
               alternative="less",conf.int=FALSE)$p.value
   }


    if(FixDes=="E" | FixDes=="N"){
       if(num.arm==1)powExact<- powExact + 1*(yobs1<=cexact)
       else powExact<- powExact + 1*(pv<object$test["alpha"])
       powNorm<-  powNorm  + 1*(r.fin>C2)
    }else{
        intY<-(r.int<=C1|r.int>=C1U)
        pstopAlt<-pstopAlt + 1*intY*(n1<n)

        if(r.int>C1) {
           if(num.arm==1)powExact<- powExact + 1*(yobs1<=cexact)
           else powExact<- powExact + 1*(pv<object$test["alpha"]|r.int>=C1U)
        }
  
        dif1<-exp(-exp(intout["Log(cumL1)"])) - exp(-exp(intout["Log(cumL0)"]))
        dif2<-exp(-exp(finout["Log(cumL1)"])) - exp(-exp(finout["Log(cumL0)"]))
  
        if(r.int>C1&(r.fin>C2)){
          powNorm <- powNorm + 1
          difRl<-min(difRl,dif2)
          difKl<-min(difKl,dif1)
        }else if(r.int>C1){
          difRh<-max(difRh,dif2)
          difKl<-min(difKl,dif1)
        }else{
           difKh<-max(difKh,dif1)
        }
        ### proportion of info at interim
        intvarAlt<-intvarAlt+(intout["se"])^2
        finvarAlt<-finvarAlt+(finout["se"])^2
     } 

  }


  alphaExact <- alphaExact/sim.n
  alphaNorm  <- alphaNorm/sim.n
  powExact   <- powExact/sim.n
  powNorm    <- powNorm/sim.n
  eda   <- eda/sim.n
  etsl  <- etsl/sim.n
  es    <- es/sim.n
  aveE  <- aveE/sim.n
  n1int <- n1int/sim.n
  t1int <- t1int/sim.n
  intvarNull<-intvarNull/sim.n
  finvarNull<-finvarNull/sim.n
  pinfoNull <- finvarNull/intvarNull
  pstopNull <- pstopNull/sim.n
  intvarAlt<-intvarAlt/sim.n
  finvarAlt<-finvarAlt/sim.n
  pinfoAlt  <- finvarAlt/intvarAlt
  pstopAlt  <- pstopAlt/sim.n
  names(t1int)<-""
  names(n1int)<-""
  names(eda)<-""
  names(etsl)<-""
  names(pinfoNull)<-""
  names(pinfoAlt)<-""
  names(pstopNull)<-""
  names(pstopAlt)<-""
  names(es)<-""
  names(alphaExact)<-""
  names(alphaNorm)<-""
  names(powExact)<-""
  names(powNorm)<-""


  return(c(alphaExact=alphaExact,alphaNorm=alphaNorm,
           powerExact=powExact,powerNorm=powNorm,
           eda=eda,etsl=etsl,
           es=es,pstopNull=pstopNull,pstopAlt=pstopAlt,aveE=aveE,
           pinfoNull=pinfoNull,pinfoAlt=pinfoAlt,
           n1=n1int,t1=t1int,
           difKl=difKl,difKh=difKh,
           difRl=difRl,difRh=difRh))

}


ft<-function(ti,y,x,EW){
  tymin<-pmax(0,ti-y)
  return( sum(pmin(x,tymin))-EW )
}

