#############################################
# denormalizes a matrix normalized with norm01
###########################################


denorm01<-function(Xn,mins,maxs) 
{Xn<-as.data.frame(Xn);

 for (i in 1:length(Xn)) 
  { 
   Xn[,i]<-(Xn[,i]*(maxs[i]-mins[i])+mins[i]); 
	
  }
Xn

}

###########################################