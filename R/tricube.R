###################################
# returns the tricube of x
###################################

tricube<-function(x) 
{ yy<-x-median(x)
  mi<-min(yy)
  mx<-max(yy)
  if(abs(mx)>abs(mi))
     {a<-1/mx}
  else {a<--1/mi}
  xn<-a*yy	
  (1-(abs(xn))^3)^3
}
##################################