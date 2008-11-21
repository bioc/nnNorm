########################################################################
# Intensity and spatial normalization using robust neural networks fitting
#######################################################################
require(nnet)
require(marray)


maNormNN<-function(mbatch,w=NULL,binWidth=3,binHeight=3,model.nonlins=3,iterations=100,nFolds=10,maplots=FALSE,verbose=FALSE) 
{
 #see what kind of object is mbatch and initialize a new normalized one "mbatchn"
 if (class(mbatch)=="marrayRaw") {
     mbatchn<-as(mbatch, "marrayNorm")
 }else if(class(mbatch)=="marrayNorm"){
         mbatchn<-mbatch
 }else stop("Mbatch must be an object of type marrayRaw or marrayNorm")

 # check validity of binWidth
 if(is.null(binWidth)){nbx=1} else if ((binWidth<=maNsc(mbatch))&&(binWidth>=1))
    {nbx<-floor(maNsc(mbatch)/binWidth)
 }else stop(paste("binWidth must be either NULL or an integer superior to 1 and inferior to ",maNsc(mbatch),"\n")) 

 # check validity of binHeight
 if(is.null(binHeight)){nby=1}else if ((binHeight<=maNsr(mbatch))&&(binHeight>=1))
    {nby<-floor(maNsr(mbatch)/binHeight)
 }else stop(paste("binHeight must be either NULL or an integer superior to 1 and inferior to ",maNsr(mbatch),"\n")) 

 #decide whether or not to use the spatial coordinates, function of nbx and nby vals.
 vectd<-NULL
 if (nbx>=2) {vectd<-c(vectd,1)
 }else cat(paste("binWidth being too large or NULL, the space coordinate X was not used in normalization","\n"))

 if (nby>=2){vectd<-c(vectd,2)
 }else cat(paste("binHeight being too large or NULL, the space coordinate Y was not used in normalization","\n"))
 vectd<-c(vectd,3)
 
 #check validity of model.nonlins
 if ((model.nonlins>=1)&&(model.nonlins<=20))
    {nodes<-model.nonlins
 }else stop(paste("model.nonlins must be an integer superior to 0 and inferior to ",20)) 

 #check validity of ncv
 if ((nFolds>=3)&&(nFolds<=50))
    {ncv=nFolds
 }else stop(paste("model.nonlins must be an integer superior to 0 and inferior to ",50)) 


 #check validity of iterations
 if ((iterations>=50)&&(iterations<=500))
    {ite<-iterations
 }else stop(paste("iterations must be an integer superior to 50 and inferior to ",500)) 
 
#check validity of w
 if(is.null(w)){
 w=rep(1,dim(maM(mbatch))[1])
 } else if(length(w)!=maNspots(mbatch)){stop("the weights vector w should have the length maNspots(mbatch)!!!")}


 nv<-length(vectd)
 nreps=5

 set.seed(1234) #require a fixed seed to obtain always the same normalized results 
 
#define some functions
 norm01<-function(x){(x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE))} 
 norm28<-function(x){((x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE)))*0.6+0.2} 
 denorm28<-function(x,mn,mx){
   a<-x*1.0/0.6-1.0/3.0;
   (a*(mx-mn)+mn); 
 }

 
 #define bin limits on X and Y direction on the print-tip group
 bXc<-seq(0,maNsc(mbatch),length=(nbx+1))
 bYr<-seq(0,maNsr(mbatch),length=(nby+1))
 breaksXc<-as.integer(bXc); XcLev<-seq(0,1,length=nbx);
 breaksYr<-as.integer(bYr); YrLev<-seq(0,1,length=nby);

 NSlides<-maNsamples(mbatch)           #get the number of slides in batch
 Npt<-maNgr(mbatch)*maNgc(mbatch)      #get the number of print tips

 for (s in 1:NSlides){                                        #for each slide
   cat(paste("\n","Processing array ",s," of ",NSlides,"\n",sep=""));
   for (pt in 1:Npt){                                         #for each print Tip
        cat("*");

        #get positions in each print tip
        ind<-maPrintTip(mbatch[,s])==pt 
        A<-maA(mbatch[ind,s])	                             #get M (log ratios)
        M<-maM(mbatch[ind,s])                              #get A (log intensity)
        xc<-factor(cut(maSpotCol(mbatch[ind,s]), breaks = breaksXc))
        yr<-factor(cut(maSpotRow(mbatch[ind,s]), breaks = breaksYr))
        levels(xc)<-XcLev;
        levels(yr)<-YrLev;
        #buid the predictors-target matrix
        Xy<-cbind(xc,yr,A,M);colnames(Xy)<-c("xc","yr","A","M")
        Xyn<-cbind(apply(Xy[,c("xc","yr","A")],2,norm01),norm28(Xy[,"M"]))
        colnames(Xyn)<-c("rxc","ryr","rA","rM")

        sbs<-sample(rep(1:ncv,dim(Xyn)[1]/(ncv-1))[1:(dim(Xyn)[1])])  
        prd<-matrix(NA,dim(Xyn)[1],nreps);

        for(rps in 1:nreps){
        for(cv in 1:ncv){        
         ssetT<-(1:(dim(Xyn)[1]))[sbs!=cv]
         ssetG<-(1:(dim(Xyn)[1]))[sbs==cv]
         if(nv==3){        
         nety1 <- nnet(rM~rxc+ryr+rA,data=Xyn, subset=ssetT,weights=w[ind], size = nodes, rang = 0.5,
               maxit = ite,reltol=0.75e-7,trace=verbose)
         }else if(nv==1){
          nety1 <- nnet(rM~rA,data=Xyn, subset=ssetT,weights=w[ind], size = nodes, rang = 0.5,
               maxit = ite,reltol=0.75e-7,trace=verbose)
         }else if(sum(vectd==c(1,3))==2){
          nety1 <- nnet(rM~rxc+rA,data=Xyn, subset=ssetT,weights=w[ind], size = nodes, rang = 0.5,
               maxit = ite,reltol=0.75e-7,trace=verbose)       
         }else {
          nety1 <- nnet(rM~ryr+rA,data=Xyn, subset=ssetT,weights=w[ind], size = nodes, rang = 0.5,
               maxit = ite,reltol=0.75e-7,trace=verbose)       
         }
         yen<-predict(nety1,Xyn[ssetG,])
         prd[sbs==cv,rps]<-denorm28(yen,mx=max(M,na.rm=TRUE),mn=min(M,na.rm=TRUE))

        }
       }
       
       
       maM(mbatchn)[ind,s]<- M-apply(prd,1,median,na.rm=TRUE)


   } #end with print tip 
   if(maplots) {
     #ploting the original M vs A for each slide 
     if(!interactive()){ 
       par(mfrow=c(2,1))
       par(mar=c(2, 4, 2, 2))

     } else {x11()}
     
     plot(maA(mbatch[,s]),maM(mbatch[,s]),xlab="A",ylab="M",main=paste("MA-plot before ANN normalization. Slide #",s,sep=""),pch=20,cex.main=0.7,cex.lab=0.7);
     lines(c(min(na.omit(maA(mbatch[,s]))),max(na.omit(maA(mbatch[,s])))),c(0,0),col="grey",lwd=2); 
     
     if (nv==1){
       cola<-rep(c(2,3,4,5,6,7,8),10);  
       for(ppt in 1:Npt){
        a<-maA(mbatch[maPrintTip(mbatch)==ppt,s])
        be<-maM(mbatch[maPrintTip(mbatch)==ppt,s])-maM(mbatchn[maPrintTip(mbatch)==ppt,s])
        or<-order(a)
        xx<-a[or];yy<-be[or];
         points(xx,yy,col=cola[ppt],cex=0.5);
        }
       legend(min(na.omit(maA(mbatch[,s]))),max(na.omit(maM(mbatch[,s]))),legend=paste(rep("Pt.",Npt),1:Npt,sep=""),text.col=cola[1:Npt],cex=0.5,ncol=Npt/2)
     }
     if(interactive()){x11()}
     plot(maA(mbatchn[,s]),maM(mbatchn[,s]),xlab="A",ylab="M",main=paste("MA-plot after ANN normalization. Slide #",s,sep=""),pch=20,cex.main=0.7,cex.lab=0.7);
     lines(c(min(na.omit(maA(mbatch[,s]))),max(na.omit(maA(mbatch[,s])))),c(0,0),col="red",lwd=2); 

    }

 }#end with slide  

 #return normalized object
 slot(mbatchn, "maNormCall")<-call("maNormNN")
 return(mbatchn)
}



