########################################################################
# Intensity and spatial normalization using robust neural networks fitting
#######################################################################
require(nnet)
require(marray)


maNormNN<-function(mbatch,binWidth=3,binHeight=3,model.nonlins=3,iterations=200,save.models=TRUE,robust=TRUE,maplots=FALSE) 
{
 #see what kind of object is mbatch and initialize a new normalized one "mbatchn"
 if (class(mbatch)=="marrayRaw") {
     mbatchn<-as(mbatch, "marrayNorm")
 }else if(class(mbatch)=="marrayNorm"){
         mbatchn<-mbatch
 }else stop("Mbatch must be an object of type marrayRaw or marrayNorm")

 # check validity of binWidth
 if ((binWidth<=maNsc(mbatch))&&(binWidth>=1))
    {nbx<-floor(maNsc(mbatch)/binWidth)
 }else stop(paste("binWidth must be an integer superior to 1 and inferior to ",maNsc(mbatch))) 

 # check validity of binHeight
 if ((binHeight<=maNsr(mbatch))&&(binHeight>=1))
    {nby<-floor(maNsr(mbatch)/binHeight)
 }else stop(paste("binHeight must be an integer superior to 1 and inferior to ",maNsr(mbatch))) 

 #decide whether or not to use the spatial coordinates, function of nbx and nby vals.
 vectd<-NULL
 if (nbx>=2) {vectd<-c(vectd,1)
 }else paste("binWidth being too large, the space coordinate X was not used in normalization")

 if (nby>=2) {vectd<-c(vectd,2)
 }else paste("binHeight being too large, the space coordinate Y was not used in normalization")
 vectd<-c(vectd,3)
 
 #check validity of model.nonlins
 if ((model.nonlins>=1)&&(model.nonlins<=30))
    {nodes<-model.nonlins
 }else stop(paste("model.nonlins must be an integer superior to 0 and inferior to ",30)) 

 #check validity of iterations
 if ((iterations>=50)&&(iterations<=500))
    {ite<-iterations
 }else stop(paste("iterations must be an integer superior to 50 and inferior to ",500)) 


 #define bins limits on X and Y direction on the print-tip group
 bXc<-seq(0,maNsc(mbatch),length=(nbx+1))
 bYr<-seq(0,maNsr(mbatch),length=(nby+1))
 breaksXc<-as.integer(bXc); XcLev<-seq(0,1,length=nbx);
 breaksYr<-as.integer(bYr); YrLev<-seq(0,1,length=nby);

 NSlides<-maNsamples(mbatch)           #get the number of slides in batch
 Npt<-maNgr(mbatch)*maNgc(mbatch)	   #get the number of print tips
 Nspots<-maNspots(mbatch)		   #get the number of spots per slide
 Models.w<-NULL				   #here will be put the weigths of the models
 Models.lims<-NULL				   #here will be put the ranges of aplicability

 for (s in 1:NSlides){                                        #for each slide
   for (pt in 1:Npt){                                         #for each print Tip
     
        #get positions in each sector/tip group
        ind<-maPrintTip(mbatch[,s])==pt 
        A<-maA(mbatch[ind,s])	                             #get M (log ratios)
        M<-maM(mbatch[ind,s])                              #get A (log intensity)
        xc<-factor(cut(maSpotCol(mbatch[ind,s]), breaks = breaksXc))
        yr<-factor(cut(maSpotRow(mbatch[ind,s]), breaks = breaksYr))
        levels(xc)<-XcLev;
        levels(yr)<-YrLev;
        #buid the predictors-target matrix
        Xy<-cbind(xc,yr,A,M);
        Xyclean<-na.omit(Xy);	#remove NaN
        XLst<-norm01(Xyclean[,1:3]); #normalize predictors to be in interval 0_1 
	  yLst<-norm28(cbind(Xyclean[,4])); #normalize target to be in interval 0.2_0.8
        yn<-yLst$Xn;
	  Xn<-XLst$Xn;
        #compute weigths for samples
        bins<-20;  nv<-length(vectd)
        breaksA <-seq(0, dim(Xn)[1],length=(bins+1))
        ax<-factor(cut(1:dim(Xn)[1],breaks=breaksA))
        levels(ax)<-1:bins
        ord<-order(Xn[,3])
	  Mat<-cbind(Xn[ord,],yn[ord],ax)
        weigthsamp<-NULL;
	  for (i in 1:bins)
	    { Minbini<-Mat[Mat[,5]==i,4]
            weigthi<-tricube(Minbini)
 	      weigthsamp<-c(weigthsamp, weigthi)
	     }
	  #end compute weights

        trials<-1; mm<-1;
	  while((trials<=2)&(mm>0.05)){
         if(robust){	  
          nety <- nnet(cbind(Mat[,vectd]), Mat[,4], weights=weigthsamp,size = nodes, rang = 0.5,
                decay = 0, maxit = ite,reltol=0.75e-7)
         }else {nety <- nnet(cbind(Mat[,vectd]), Mat[,4], size = nodes, rang = 0.5,
                decay = 0, maxit = ite,reltol=0.75e-7)
		       }  
          # get the fitted values 
          yen<-predict(nety,cbind(Xn[,vectd]))
          #and denormalize them back
          Mfit<-denorm28(yen,yLst$mins,yLst$maxs)
	    Mfitted<-cbind(Mfit[[1]])
          mm<-abs(mean(Xyclean[,4]-Mfitted)) #compute the mean of normalized M values, to know if
							     #by chance, training was not done well	
	    trials<-trials+1
        }
         # get the fitted values for all points  
          yen<-predict(nety,cbind(Xn[,vectd]))
          #and denormalize them back
          Mfit<-denorm28(yen,yLst$mins,yLst$maxs)
	    Mfit<-cbind(Mfit[[1]])
        
        # insert NaNs where they where in the original M vector if is the case
         NAloc<-attr(Xyclean,"na.action")
        if (length(NAloc)>0){
          yelong<-matrix(NaN,dim(M)[1],1)
          yelong[-NAloc,1]<-Mfit
        }else{yelong<-Mfit }

      #normalizing by substracting the predicted values
       maM(mbatchn)[ind,s]<-cbind(M)-yelong
      #store the neural network parameters if requested
      if((save.models)&&(nv==1)) {
        Models.w<-cbind(Models.w,nety$wts)
        Models.lims<-cbind(Models.lims,c(XLst$mins[vectd],yLst$mins,XLst$maxs[vectd],yLst$maxs))
      }  
	
   } #end with print tip 
   if(maplots) {
     #ploting the original M vs A for each slide 
     x11()
     plot(maA(mbatch[,s]),maM(mbatch[,s]),xlab="A",ylab="M",main=c("MA-plot before ANN   normalization for slide No.  ",s),pch=20);
     lines(c(min(na.omit(maA(mbatch[,s]))),max(na.omit(maA(mbatch[,s])))),c(0,0),col="red",lwd=2); 
  
     #ploting the normalized M vs A for each slide 
     x11()
     plot(maA(mbatchn[,s]),maM(mbatchn[,s]),xlab="A",ylab="M",main=c("MA-plot after ANN normalization for slide  No.  ",s),pch=20);
     lines(c(min(na.omit(maA(mbatch[,s]))),max(na.omit(maA(mbatch[,s])))),c(0,0),col="red",lwd=2); 
   }

 }#end with slide  

 #return normalized object
 slot(mbatchn, "maNormCall")<-call("maNormNN")
 if((save.models)&&(nv==1)){
   dim(Models.w)<-c(length(nety$wts),Npt,NSlides)
   dim(Models.lims)<-c(4,Npt,NSlides)
   mList<-list(par=Models.w,lims=Models.lims,nonlins=model.nonlins)
   } else{mList<-NULL}
 return(list(batchn=mbatchn,models=mList))
}

#end maNormNN#######################################################################333