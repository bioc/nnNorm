###############################################
# Predicts the bias estimate provided a neural network distortion color model
##############################################

predictBias<-function(Avals,Models,slide,pT) 
{ A<-Avals;
  models<-Models;

 if (class(models)!="list") {
   stop("models must be a list produced by maNormNN function!")
 }

 if (slide>dim(Models$par)[3]) {
   stop(paste("The Modles object has only",dim(Models$par)[3]," slides"))
 }
 

 lims<-models$lims[,pT,slide]
 minA<-lims[1]; minM<-lims[2]; maxA<-lims[3]; maxM<-lims[4];

 if((min(A)<minA)||(max(A)>maxA)){
   stop(paste("A must not contain NaNs or have the minimum inferior to ",minA," or have the maximum superior to ",maxA))
 }
 
 I<-1;J<-models$nonlins;
 ws<-models$par[,pT,slide]

 wij<-array(ws[1:((I+1)*J)],c(I+1,J))
 ord<-c(2:(I+1),1)
 wij<-wij[ord,] #ok

 wj<-ws[((I+1)*J+1):length(ws)]
 ord<-c(2:length(wj),1)
 wj<-wj[ord]

 sig<-function(x){1/(1+exp(-x))}
 
 #scale A in interval 0 1
 As<-(A-minA)/(maxA-minA) 


 est<-NULL
 for (i in 1:length(A)){
  x<-c(As[i],1)
  est[i]<-sig(c(sig(x%*%wij),1)%*%wj)
 }
 
 denorm28(cbind(est),minM,maxM)$est

}
##################################