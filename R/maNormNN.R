########################################################################
# Intensity and spatial normalization using robust neural networks fitting
#######################################################################
require(nnet)
require(marray)


maNormNN<-function(mbatch,binWidth=3,binHeight=3,model.nonlins=3,iterations=200,robust=TRUE,maplots=FALSE) 
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
 if ((model.nonlins>=1)&&(model.nonlins<=20))
    {nodes<-model.nonlins
 }else stop(paste("model.nonlins must be an integer superior to 0 and inferior to ",20)) 

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
 set.seed(100) #require a fixed seed to obtain always the same normalized results 

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
        len<-dim(Mat)[1]
        aS<-1:len;
        smp1<-sample(aS,round(len/4)) 
        smp2<-sample(aS[-smp1],round(len/4)) 
        smp3<-sample(aS[-c(smp1,smp2)],round(len/4)) 
        smp4<-aS[-c(smp1,smp2,smp3)] 

        s234<-c(smp2,smp3,smp4);
        s134<-c(smp1,smp3,smp4);
        s124<-c(smp1,smp2,smp4); 
        s123<-c(smp1,smp2,smp3); 

        trials<-1; mm<-1;
        if(!robust){weigthsamp<-rep(1,length(weigthsamp))}
	  while((trials<=2)&(mm>0.1)){
          nety1 <- nnet(cbind(Mat[s234,vectd]), Mat[s234,4], weights=weigthsamp[s234],size = nodes, rang = 0.5,
                decay = 0, maxit = ite,reltol=0.75e-7)
          nety2 <- nnet(cbind(Mat[s134,vectd]), Mat[s134,4], weights=weigthsamp[s134],size = nodes, rang = 0.5,
                decay = 0, maxit = ite,reltol=0.75e-7)
          nety3 <- nnet(cbind(Mat[s124,vectd]), Mat[s124,4], weights=weigthsamp[s124],size = nodes, rang = 0.5,
                decay = 0, maxit = ite,reltol=0.75e-7)
          nety4 <- nnet(cbind(Mat[s123,vectd]), Mat[s123,4], weights=weigthsamp[s123],size = nodes, rang = 0.5,
                decay = 0, maxit = ite,reltol=0.75e-7)

          # get the fitted values 
          yen<-array(dim=len) 
          yen[-s234]<-predict(nety1,cbind(Xn[-s234,vectd]))
          yen[-s134]<-predict(nety2,cbind(Xn[-s134,vectd]))
          yen[-s124]<-predict(nety3,cbind(Xn[-s124,vectd]))
          yen[-s123]<-predict(nety4,cbind(Xn[-s123,vectd]))
          #and denormalize them back
          Mfit<-denorm28(yen,yLst$mins,yLst$maxs)
	    Mfitted<-cbind(Mfit[[1]])
          mm<-abs(mean(Xyclean[,4]-Mfitted)) #compute the mean of normalized M values, to know if all is ok
	    trials<-trials+1
        }

	  Mfit<-Mfitted;
        
        # insert NaNs where they where in the original M vector if is the case
         NAloc<-attr(Xyclean,"na.action")
        if (length(NAloc)>0){
          yelong<-matrix(NaN,dim(M)[1],1)
          yelong[-NAloc,1]<-Mfit
        }else{yelong<-Mfit }

      #normalizing by substracting the predicted values
       maM(mbatchn)[ind,s]<-cbind(M)-yelong

   } #end with print tip 
   if(maplots) {
     #ploting the original M vs A for each slide 
     if(interactive()){x11()}
    
     plot(maA(mbatch[,s]),maM(mbatch[,s]),xlab="A",ylab="M",main=paste("MA-plot before ANN normalization. Slide #",s,sep=""),pch=20,cex.main=0.5,cex.lab=0.5);
     lines(c(min(na.omit(maA(mbatch[,s]))),max(na.omit(maA(mbatch[,s])))),c(0,0),col="black",lwd=2); 
     if (nv==1){
      cola<-rep(c(2,3,4,5,6,7,8),10) 

       for(ppt in 1:Npt){
        a<-maA(mbatch[maPrintTip(mbatch)==ppt,s])
        be<-maM(mbatch[maPrintTip(mbatch)==ppt,s])-maM(mbatchn[maPrintTip(mbatch)==ppt,s])
        or<-order(a)
        points(a[or],be[or],col=cola[ppt],pch=20,cex=0.5);}
     }
     legend(min(na.omit(maA(mbatch[,s]))),max(na.omit(maM(mbatch[,s]))),legend=paste(rep("Pt.",Npt),1:Npt,sep=""),text.col=cola[1:Npt],cex=0.5,ncol=Npt/2)
     if(interactive()){x11()}
     plot(maA(mbatchn[,s]),maM(mbatchn[,s]),xlab="A",ylab="M",main=paste("MA-plot after ANN normalization. Slide #",s,sep=""),pch=20,cex.main=0.5,cex.lab=0.75);
     lines(c(min(na.omit(maA(mbatch[,s]))),max(na.omit(maA(mbatch[,s])))),c(0,0),col="red",lwd=2); 

    }

 }#end with slide  

 #return normalized object
 slot(mbatchn, "maNormCall")<-call("maNormNN")
 return(mbatchn)
}