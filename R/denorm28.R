#############################################
# denormalizes a vector normalized with norm28
########################################

denorm28<-function(Xn,mins,maxs) 
{Xn<-as.data.frame(Xn);

 for (i in 1:length(Xn)) 
  { 
   a<-Xn[,i]*1.0/0.6-1.0/3.0
   Xn[,i]<-(a*(maxs[i]-mins[i])+mins[i]); 
   	
  }
Xn

}

##############################################3333