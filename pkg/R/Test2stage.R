"Test2stage" <-
function(x,C1,C1U,C2,tan,num.arm,Y11,T11,p0=NULL,Y21=NULL,T21=NULL,Y10=NULL,T10=NULL,Y20=NULL,T20=NULL,printTest=TRUE,
cen11=rep(1,length(T11)),cen21=rep(1,length(T21)),cen10=rep(1,length(T10)),cen20=rep(1,length(T20)))
{
  if(length(Y11)!=length(T11)|length(Y21)!=length(T21))
    stop("study times and failure times should be of equal length")

  if(length(Y10)!=length(T10)|length(Y20)!=length(T20))
    stop("study times and failure times should be of equal length")

  if(x>=tan)
    stop("the survival time of interest should not exceed the interim study time")

  if(num.arm==1 & is.null(p0))
    stop("The null rate must be specified in single arm trials.")

  if(num.arm==2 & is.null(Y10))
    stop("A control group must be specified in two arm trials.")

  if(printTest!=FALSE&printTest!=TRUE)
    stop("printTest should be a logical value.")

  if(length(Y21)==0)
  {
    z1out <- tst(tan,x,num.arm,Y11,T11,p0,Y10,T10,cen11,cen10)         # p0: event rate under the null hypothesis
    z1<-z1out["z"]
    names(z1)<-"z1"
    if(printTest==TRUE)
    {  
      cat(paste("interim test statistic Z1 =",z1,"C1L =",C1),"\n")
      if(z1<C1)
        cat("Z1 < C1L, stop the study and accept the null hypothesis \n") 
         
      if(z1>C1U)  
        cat("Z1 > C1U, stop the study and reject the null hypothesis \n") 
      
      if(z1>=C1&z1<=C1U)  
        cat("Z1 >= C1L and Z1 <= C1U, continue to the second stage \n")
        
    }  
    return(z1out)
    
  }
  else
  {
    Y1=c(Y11,Y21)
    Y0=c(Y10,Y20)
    T1=c(T11,T21)
    T0=c(T10,T20)
    cen1<-c(cen11,cen21)
    cen0<-c(cen10,cen20)
    z2out <- tst(tan,x,num.arm,Y1,T1,p0,Y0,T0,cen1,cen0)
    z2<-z2out["z"]
    if(printTest==TRUE)
    {  
      cat(paste("final test statistic Z2 =",z2,"C2 =",C2),"\n")    
      if(z2>C2)
        cat("since Z2 > C2, reject the null hypothesis \n")
      if(z2<=C2)
        cat("since Z2 <= C2, cannot reject the null hypothesis \n")
    }  
    return(z2out)
     
  }

}  


  # Nelson-Aalen estimate of the cumulative hazard function
  cum <- function(Y,T,t,x,cen)
  {
    X <- pmin(T,x,pmax(0,t-Y))
    delta <- as.numeric(T<=pmin(x,pmax(0,t-Y)) & cen)
    ### set to censored status if explicitly specified
    n <- length(X)
    R <- as.numeric(n)
    for(i in 1:n)
      R[i] <- sum(X>=X[i])
    cum<-sum(delta[X<=x]/R[X<=x])
    varcum<-sum(delta[X<=x]/(R[X<=x]^2))
    return(c(cum=cum,varcum=varcum))
  }

  # test statistic Z(x;t) of equation 3
  tst <- function(t,x,num.arm,Y1,T1,p0,Y0,T0,cen1,cen0)
  {
    cumout1<-cum(Y1,T1,t,x,cen1)
    if(num.arm==1){cumout0<- c(cum=-log(1-p0),varcum=0)
                   names(cumout0)<-c("cum","varcum")
    }
    else cumout0<-cum(Y0,T0,t,x,cen0)
    se<-  sqrt( cumout0["varcum"]/(max(cumout0["cum"],1e-4)^2)+cumout1["varcum"]/(max(cumout1["cum"],1e-4)^2) )
    z<-(log(max(cumout0["cum"],1e-4))-log(max(cumout1["cum"],1e-4)))/se

	names(z)<-"z"
    names(se)<-"se"
    cumL<-log(c(max(cumout1["cum"],1e-4), max(cumout0["cum"],1e-4)))
    names(cumL)<-c("Log(cumL1)","Log(cumL0)")
    return(c(z,se,cumL))
  }
