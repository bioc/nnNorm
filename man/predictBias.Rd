\name{predictBias}
\alias{predictBias}
\title{
  Computes the bias estimate provided a neural network color distortion model
}
\description{Given a neural network color distorsion model and an intensity value, this function computes
    de estimate of the bias of log ratios, \eqn{M}. The color distorsion model should be obtained with
 \code{\link{maNormNN}}.
}
\usage{
predictBias(Avals,Models,slide,pT) 
}

\arguments{
  \item{Avals}{An array of average log intensities for which the bias will be estimated. All the values
   in array must be inside the validity ranges of the model, i.e. form \code{Models$lims[1,pT,slide]} to    
   \code{Models$lims[3,pT,slide]}
  }
  \item{Models}{The \code{models} component of the list returned by the function \code{\link{maNormNN}}. 
  }	
  \item{slide}{The slide number for which prediction is desired. 
  }	
  \item{pT}{The print tip number for which prediction is desired.
  }
  
}

\details{This function is used to make interpolations usings the color distortion models saved during the 
 normalization effectuated with the function \code{\link{maNormNN}}.
}


\value{An array of the same length as \code{Avals} with the coresponding bias estimates.
}


\examples{ 
data(swirl)
#normalize swirl
sNNr<-maNormNN(swirl,binWidth=maNsc(swirl),binHeight=maNsr(swirl),save.models=TRUE,robust=TRUE)

#retrive original M-A values for slide 4 and print tip 3
s<-4;pt<-3;
MM<-maM(swirl[maPrintTip(swirl)==pt,s])
AA<-maA(swirl[maPrintTip(swirl)==pt,s])

#generate a series of A values in the validity range of the model for slide s and print tip pt
A<-seq(sNNr$models$lims[1,pt,s],sNNr$models$lims[3,pt,s],length=100)

#do the plots
x11()
plot(AA,MM,pch=20,xlab="A",ylab="M",main=paste("Slide=",s," Print tip=",pt)) #raw data
lines(A,predictBias(A,sNNr$models,s,pt),col="red",lwd=2) #robust NN

#for ilustration add the loess normalization curve as computed in marray package
lo<-loess(MM~AA,span=0.4,degree=1,family="symmetric",control=loess.control(trace.hat="approximate",iterations=5,surface="direct"))
lines(A,predict(lo,A),col= "green",lty="longdash",lwd=2) #loess

}
  
\author{Tarca, A.L.}
\references{
 A. L. Tarca, J. E. K. Cooke, and J. Mackay. Robust neural networks approach for spatial and
intensity dependent normalization of cDNA data. Bioinformatics. 2004,submitted.\cr

}	

\seealso{\code{\link{maNormNN}}}



\keyword{smooth}



