########################################################################
# Compares the distribution of several data sets
#   
#######################################################################



compNorm<-function(x,...,bw="AUTO",xlim=c(-3,3),titles="AUTO",type="d") 
{ 
 a<-list(...)
 lna<-length(a)
 l <- as.list(substitute(list(...)))[-1]
 if(length(titles)==1)
 {nms<-deparse(substitute(x))
  for (i in 1:lna){
  nms<-c(nms,l[[i]])
    }
  }
  else {nms<-titles}
 
  # check validity of bw
  qx<-quantile(x,na.rm=TRUE);bws<-(qx[4]-qx[2])/2
 if (bw=="AUTO"){bww<-bws} 
 else 
  {bww<-bws}

 # check validity of xlim
 xlims<-range(xlim)

  #Create a boxplot with the vlaues of vectors
 if (type=="b"){ 
   boxplot(x,...,range=5,names=nms,cex=0.75) 
   lines(c(0,lna+2),c(0,0),col="red",lwd=2)
 } 
 if(type=="d") {
   # density distributions
   cols<-1:(1+lna);
   pchs<-rep("*",(1+lna))
   lin<-NULL;

   lin<-density(na.omit(x),bw=bww); 
   maxy<-max(lin$y)
   lst<-list(lin)
   for (i in 1:length(a)){
    lst[[1+i]]<-density(na.omit(a[[i]]),bw=bww)
    if (maxy<=max(lst[[1+i]]$y)) {maxy<-max(lst[[1+i]]$y) }
   }

 
   # do the plots
   plot(c(0,0),c(0,maxy),xlim=xlims,ylim=range(c(0,maxy)),xlab="M",ylab="Density",type="l")
   for (i in 1:(lna+1)){
     lines(lst[[i]],col=cols[i],lwd=2)
   }
   legend(xlims[1], maxy, legend = nms, col = cols, fill=cols, cex = 0.7)
 }
 madval<-NULL
 #compute mad values
 madval[1]=mad(x,na.rm=TRUE)
 cat("MAD ",as.character(nms[1]),":",madval[1],"\n")
 for (i in 1:length(a)){
   madval[1+i]=mad(a[[i]],na.rm=TRUE)
 cat("MAD ",as.character(nms[i+1]),":",madval[i+1],"\n")
 }
}
# end with compNorm############################################################3#