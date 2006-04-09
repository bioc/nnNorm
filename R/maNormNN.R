########################################################################
# Intensity and spatial normalization using robust neural networks fitting
#######################################################################
require(nnet)
require(marray)


maNormNN<-function(mbatch,binWidth=3,binHeight=3,model.nonlins=3,iterations=200,robust=TRUE,maplots=FALSE, verbose=FALSE) 
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
 Npt<-maNgr(mbatch)*maNgc(mbatch)      #get the number of print tips
 Nspots<-maNspots(mbatch)          #get the number of spots per slide
 set.seed(200) #require a fixed seed to obtain always the same normalized results 

 for (s in 1:NSlides){                                        #for each slide
   cat(paste("\n","Processing array ",s," of ",NSlides,"\n",sep=""));
   for (pt in 1:Npt){                                         #for each print Tip
        cat("*");
        #get positions in each sector/tip group
        ind<-maPrintTip(mbatch[,s])==pt 
        A<-maA(mbatch[ind,s])                                #get M (log ratios)
        M<-maM(mbatch[ind,s])                              #get A (log intensity)
        xc<-factor(cut(maSpotCol(mbatch[ind,s]), breaks = breaksXc))
        yr<-factor(cut(maSpotRow(mbatch[ind,s]), breaks = breaksYr))
        levels(xc)<-XcLev;
        levels(yr)<-YrLev;
        #buid the predictors-target matrix
        Xy<-cbind(xc,yr,A,M);
        Xyclean<-na.omit(Xy);   #remove NaN
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
        smp1<-sample(aS,round(len/5)) 
        smp2<-sample(aS[-smp1],round(len/5)) 
        smp3<-sample(aS[-c(smp1,smp2)],round(len/5)) 
        smp4<-sample(aS[-c(smp1,smp2,smp3)],round(len/5))
        smp5<-aS[-c(smp1,smp2,smp3,smp4)] 

        s2345<-c(smp2,smp3,smp4,smp5);
        s1345<-c(smp1,smp3,smp4,smp5);
        s1245<-c(smp1,smp2,smp4,smp5); 
        s1235<-c(smp1,smp2,smp3,smp5); 
        s1234<-c(smp1,smp2,smp3,smp4);
         
      if(!robust){weigthsamp<-rep(1,length(weigthsamp))}
      estm<-NULL;
     for (rr in 1:4){ #get 4 estimates
        trials<-1; mm<-1;
        while((trials<=20)&(mm>0.5)){
          nety1 <- nnet(cbind(Mat[s2345,vectd]), Mat[s2345,4], weights=weigthsamp[s2345],size = nodes, rang = 0.5,
                decay=1e-4, maxit = ite,reltol=0.75e-7,trace=verbose)
          yen<-NULL; 
          yen<-predict(nety1,cbind(Xn[-s2345,vectd]))
          Mfit<-denorm28(yen,yLst$mins,yLst$maxs)
         Mfitted<-cbind(Mfit[[1]])
          mm<-abs(median(Xyclean[-s2345,4]-Mfitted)) #compute the mean of normalized M values, to know if all is ok
          trials<-trials+1;
         } 
        if(trials>=20){warning(paste("At slide #",s," and subarray #",pt," there was a convergence problem.",sep=""))}   

        trials<-1; mm<-1;
        while((trials<=20)&(mm>0.5)){
          nety2 <- nnet(cbind(Mat[s1345,vectd]), Mat[s1345,4], weights=weigthsamp[s1345],size = nodes, rang = 0.5,
                decay=1e-4, maxit = ite,reltol=0.75e-7,trace=verbose)
          yen<-NULL; 
          yen<-predict(nety2,cbind(Xn[-s1345,vectd]))
          #and denormalize them back
          Mfit<-denorm28(yen,yLst$mins,yLst$maxs)
         Mfitted<-cbind(Mfit[[1]])
          mm<-abs(median(Xyclean[-s1345,4]-Mfitted)) #compute the mean of normalized M values, to know if all is ok
          trials<-trials+1;
         } 
         if(trials>=20){warning(paste("At slide #",s," and subarray #",pt," there was a convergence problem.",sep=""))}   

        trials<-1; mm<-1;
        while((trials<=20)&(mm>0.5)){
          nety3 <- nnet(cbind(Mat[s1245,vectd]), Mat[s1245,4], weights=weigthsamp[s1245],size = nodes, rang = 0.5,
                decay=1e-4, maxit = ite,reltol=0.75e-7,trace=verbose)
          yen<-NULL; 
          yen<-predict(nety3,cbind(Xn[-s1245,vectd]))
          #and denormalize them back
          Mfit<-denorm28(yen,yLst$mins,yLst$maxs)
         Mfitted<-cbind(Mfit[[1]])
          mm<-abs(median(Xyclean[-s1245,4]-Mfitted)) #compute the mean of normalized M values, to know if all is ok
          trials<-trials+1;
        }
        if(trials>=20){warning(paste("At slide #",s," and subarray #",pt," there was a convergence problem.",sep=""))}   

        trials<-1; mm<-1;
        while((trials<=20)&(mm>0.5)){
          nety4 <- nnet(cbind(Mat[s1235,vectd]), Mat[s1235,4], weights=weigthsamp[s1235],size = nodes, rang = 0.5,
                decay=1e-4, maxit = ite,reltol=0.75e-7,trace=verbose)
          yen<-NULL;
          yen<-predict(nety4,cbind(Xn[-s1235,vectd]))
          #and denormalize them back
          Mfit<-denorm28(yen,yLst$mins,yLst$maxs)
         Mfitted<-cbind(Mfit[[1]])
          mm<-abs(median(Xyclean[-s1235,4]-Mfitted)) #compute the mean of normalized M values, to know if all is ok
          trials<-trials+1;
        }
        if(trials>=20){warning(paste("At slide #",s," and subarray #",pt," there was a convergence problem.",sep=""))}   
       
         trials<-1; mm<-1;
         while((trials<=20)&(mm>0.5)){
          nety5 <- nnet(cbind(Mat[s1234,vectd]), Mat[s1234,4], weights=weigthsamp[s1234],size = nodes, rang = 0.5,
                decay=1e-4, maxit = ite,reltol=0.75e-7,trace=verbose)
          yen<-NULL;
          yen<-predict(nety5,cbind(Xn[-s1234,vectd]))
          #and denormalize them back
          Mfit<-denorm28(yen,yLst$mins,yLst$maxs)
         Mfitted<-cbind(Mfit[[1]])
          mm<-abs(median(Xyclean[-s1234,4]-Mfitted)) #compute the mean of normalized M values, to know if all is ok
          trials<-trials+1;
       }
        if(trials>=20){warning(paste("At slide #",s," and subarray #",pt," there was a convergence problem.",sep=""))}   

       
       
        # get the fitted values 
        yen<-array(dim=len) 
        yen[-s2345]<-predict(nety1,cbind(Xn[-s2345,vectd]))
        yen[-s1345]<-predict(nety2,cbind(Xn[-s1345,vectd]))
        yen[-s1245]<-predict(nety3,cbind(Xn[-s1245,vectd]))
        yen[-s1235]<-predict(nety4,cbind(Xn[-s1235,vectd]))
        yen[-s1234]<-predict(nety5,cbind(Xn[-s1234,vectd]))
        #and denormalize them back
        Mfit<-denorm28(yen,yLst$mins,yLst$maxs)
        Mfitted<-cbind(Mfit[[1]])
        Mfit<-Mfitted;
        
        # insert NaNs where they where in the original M vector if is the case
         NAloc<-attr(Xyclean,"na.action")
        if (length(NAloc)>0){
          yelong<-matrix(NA,dim(M)[1],1)
          yelong[-NAloc,1]<-Mfit
        }else{yelong<-Mfit }
        estm<-cbind(estm,yelong);
      } #end get 3 estimates
      meann<-function(x){if(sum(is.na(x))>0){NA}else{mean(x,trim=0.25)}}
      yelong<-apply(estm,1,mean);
      #normalizing by substracting the predicted values
       maM(mbatchn)[ind,s]<-cbind(M)-yelong

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



