"Test2stage" <-
function(Y1,T1,Y2=NULL,T2=NULL,p0,x,C1,C2,t1,MTSL=NULL,printTest=TRUE,
cen1=rep(1,length(T1)),cen2=rep(1,length(T2)))
{
  if(length(Y1)!=length(T1)|length(Y2)!=length(T2))
    stop("study times and failure times should be of equal length")

  if(length(T2)>0 & length(MTSL)==0) 
    stop("MTSL must be specified when Y2 is specified")
  
  if(x>=t1)
    stop("the survival time of interest should not exceed the interim study time")

  if(printTest!=FALSE&printTest!=TRUE)
    stop("printTest should be a logical value.")

  if(length(Y2)==0)
  {
    z1out <- tst(Y1,T1,t1,p0,x,cen1)
    z1<-z1out["z"]
    names(z1)<-"z1"
    if(printTest==TRUE)
    {  
      cat(paste("interim test statistic Z1 =",z1,"C1 =",C1),"\n")
      if(z1<C1)
        cat("Z1 < C1, stop the study and accept the null hypothesis \n")
      else
        cat("Z1 >= C1, continue to the second stage \n")
    }  
    return(z1out)
    
  }
  else
  {
    Y=c(Y1,Y2)
    T=c(T1,T2)
    cen<-c(cen1,cen2)
    z2out <- tst(Y,T,MTSL,p0,x,cen)
    z2<-z2out["z"]
    if(printTest==TRUE)
    {  
      cat(paste("final test statistic Z2 =",z2,"C2 =",C2),"\n")    
      if(z2>C2)
        cat("since Z2 > C2, reject the null hypothesis \n")
      else
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

  # test statistic Z(x;t) of equation 2
  tst <- function(Y,T,t,p0,x,cen)
  {
    cumout<-cum(Y,T,t,x,cen)
    z<-(log(-log(1-p0))-log(max(cumout["cum"],1e-4)))*max(cumout["cum"],
        1e-4)/(sqrt(cumout["varcum"]))
	names(z)<-"z"
    se<-sqrt(cumout["varcum"])
    names(se)<-"se"
    cumL<-cumout["cum"]
    names(cumL)<-"cumL"
    return(c(z,se,cumL))
  }
