####################################
# scales a vector to a 0-1 bounded interval
##################################


norm01<-function(X) 
{stopifnot(all(!is.na(X)))
 dims<-dim(X)[2]
 mins<-maxs<-array(0,dims);
 for (i in 1:dims) 
  { 
  mx<-max(X[,i])
  mi<-min(X[,i])
  X[,i]<-((X[,i]-mi)/(mx-mi)); 
  mins[i]<-mi;
  maxs[i]<-mx;	
  }
list(Xn=X,mins=mins,maxs=maxs)

}
###########################################3333