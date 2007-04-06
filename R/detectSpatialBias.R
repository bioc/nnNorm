########################################################################
# Intensity and spatial normalization using robust neural networks fitting
#######################################################################
require(nnet)
require(marray)


detectSpatialBias<-function(mbatch, corThreshold=0.6) 
{
 #see what kind of object is mbatch and initialize a new normalized one "mbatchn"
 if (class(mbatch)=="marrayRaw") {
     mbatchn<-as(mbatch, "marrayNorm")
 }else if(class(mbatch)=="marrayNorm"){
         mbatchn<-mbatch
 }else stop("Mbatch must be an object of type marrayRaw or marrayNorm")


 NSlides<-maNsamples(mbatch)           #get the number of slides in batch
 Npt<-maNgr(mbatch)*maNgc(mbatch)      #get the number of print tips
 Nspots<-maNspots(mbatch)          #get the number of spots per slide
 lrow<-maNsc(mbatch)
 lcol<-maNsr(mbatch)
  resCorlCol<-array(NA,dim=c(lcol,Npt,NSlides))
 resCorlRow<-array(NA,dim=c(lrow,Npt,NSlides))


 for (s in 1:NSlides){ #for each slide
   cat(paste("\n","Processing array ",s," of ",NSlides,"\n",sep=""));
   for (pt in 1:Npt){      #for each print Tip
        cat("*");
    for(coli in 1:lrow){
        #get positions in each sector
        ccorl<-NULL  
        ind<-(maPrintTip(mbatch)==pt)&maSpotCol(mbatch)==coli 
        #M<-maM(mbatch[ind,s])[1:round(lcol/2)]  
        #resCorlRow[1,coli,pt,s]<-cor(M[!is.na(M)],(1:round(lcol/2))[!is.na(M)])      

        #M<-maM(mbatch[ind,s])[(round(lcol/2)+1):lcol]  
        #resCorlRow[2,coli,pt,s]<-cor(M[!is.na(M)],((round(lcol/2)+1):lcol)[!is.na(M)])   

        M<-maM(mbatch[ind,s])  
        resCorlRow[coli,pt,s]<-cor(M[!is.na(M)],(1:lcol)[!is.na(M)])                           
                      
     }

    for(rowi in 1:lcol){
        #get positions in each sector/tip group
        ind<-(maPrintTip(mbatch)==pt)&maSpotRow(mbatch)==rowi 
        #M<-maM(mbatch[ind,s])[1:round(lrow/2)]  
        #resCorlCol[rowi,1,pt,s]<-cor(M[!is.na(M)],(1:round(lrow/2))[!is.na(M)])      

        #M<-maM(mbatch[ind,s])[(round(lrow/2)+1):lrow]  
        #resCorlCol[rowi,2,pt,s]<-cor(M[!is.na(M)],((round(lrow/2)+1):lrow)[!is.na(M)])    

        M<-maM(mbatch[ind,s])  
        resCorlCol[rowi,pt,s]<-cor(M[!is.na(M)],(1:lrow)[!is.na(M)])                            
                        
                      
     }

   } #end with print tip 

 }#end with slide  

 biasCol<-matrix(NA,Npt,NSlides)
 biasRow<-matrix(NA,Npt,NSlides)
 sCol<-NULL;ptCol<-NULL;
 sRow<-NULL;ptRow<-NULL;
 
 for (s in 1:NSlides){                                        #for each slide
      for (pt in 1:Npt){                                         #for each print Tip
        
   
        biasRow[pt,s]<-round(max(c(sum(resCorlRow[,pt,s]>corThreshold),sum(resCorlRow[,pt,s]<(-corThreshold))))/lrow*100,1)
        biasCol[pt,s]<-round(max(c(sum(resCorlCol[,pt,s]>corThreshold),sum(resCorlCol[,pt,s]<(-corThreshold))))/lcol*100,1)
              
     }#end with print tip 


  } #end with slide  
colnames(biasRow)<-colnames(biasCol)<-colnames(maM(mbatch))
rownames(biasRow)<-rownames(biasCol)<-paste("PrintTip",1:Npt,sep="")
list(biasRow=biasRow,biasCol=biasCol)
}










