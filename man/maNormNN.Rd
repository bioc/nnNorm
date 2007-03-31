\name{maNormNN}
\alias{maNormNN}
\title{Intensity and spatial normalization using robust neural networks fitting}
\description{This function normalizes a batch of cDNA arrays by removing the intensity and spatial dependent bias.
}
\usage{
maNormNN(mbatch,w=NULL,binWidth=3,binHeight=3,model.nonlins=3,iterations=80,nFolds=5,maplots=FALSE,verbose=FALSE) 
}

\arguments{
  \item{mbatch}{A \code{marrayRaw} or \code{marrayNorm} batch of arrays.
  }
  \item{w}{Weights to be assigned to each spot. If provided, it should be a vector with the same length as maNspots(mbatch).
  }
  \item{binWidth}{Width of the bins in the \eqn{X} direction (spot column) in which the 
    print tip will be divided in order to account for spatial variation. Max value 
    is \code{maNsc(mbatch)}, Min value is 1. However if it is set to a number larger than 
    \code{maNsc(mbatch)/2} (so less than two bins in \eqn{X} direction) the variable \eqn{X} will not 
    be used as predictor to estimate the bias.
  }	
  \item{binHeight}{Height of the bins in the \eqn{Y} direction (spot row)in which the 
    print tip will be divided in order to account for spatial variation. Max value 
    is \code{maNsr(mbatch)}, Min value is 1. However if it is set to a number larger than 
    \code{maNsr(mbatch)/2} (so less than two bins in \eqn{Y} direction) the variable \eqn{Y} will not 
    be used as predictor to estimate the bias.
  }
  \item{model.nonlins}{Number of nodes in the hidden layer of the neural network model. 
  }
  \item{iterations}{The number of iterations at which (if not converged) the training of the neural net will be 
                   stopped.
  }
  \item{nFolds}{Number of cross-validation folds. It represents the number of equal parts in which the data from a 
     print tip is divided into: the model is trained on nFolds-1 parts and the bias is estimated for one part at the 
    time. Higher values improve the results but increase the computation time. Ideal values are between 5 and 10. 
  }
  \item{maplots}{If set to \code{"TRUE"} will produce a \eqn{M-A} plot for each slide before and after normalization. 
  }
  \item{verbose}{If set to \code{"TRUE"} will show the output of the nnet function which is training the neural
    network models. 
  }

}

\details{This function uses neural networks to model the bias in cDNA data sets.
}


\value{
  A \code{marrayNorm} object containing the normalized log ratios. See \code{marrayNorm} 
    class for details
  }
}

  


\examples{
# Normalization of swirl data
data(swirl)
# print-tip, intensity and spatial normalization of the first slide in swirl data set
swirlNN<-maNormNN(swirl[,1])   

#do not consider spatial variations, and display M-A plots before and after normalization
swirlNN<-maNormNN(swirl[,1],binWidth=maNsc(swirl),binHeight=maNsr(swirl),maplots=TRUE)  

}

\author{Tarca, A.L.}
\references{
 A. L. Tarca, J. E. K. Cooke, and J. Mackay. Robust neural networks approach for spatial and
intensity dependent normalization of cDNA data. Bioinformatics. 2004,submitted.\cr
}	

\seealso{\code{\link{compNorm}},\code{nnet}}

\keyword{smooth}
\keyword{robust}


