#########################################
# scales a vector to a 0.2 0.8 interval
########################################


norm28<-function(X) 
{stopifnot(all(!is.na(X)))
 dims<-dim(X)[2]
 mins<-maxs<-array(0,dims);
 for (i in 1:dims) 
  { 
  mx<-max(X[,i])
  mi<-min(X[,i])
  X[,i]<-((X[,i]-mi)/(mx-mi))*0.6+0.2; 
  mins[i]<-mi;
  maxs[i]<-mx;	

  }
list(Xn=X,mins=mins,maxs=maxs)

}

#######################################