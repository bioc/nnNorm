\name{maNormNN}
\alias{maNormNN}
\title{Intensity and spatial normalization using robust neural networks fitting}
\description{This function normalizes a batch of cDNA arrays by removing the intensity and spatial dependent bias.
}
\usage{
maNormNN(mbatch,binWidth=3,binHeight=3,model.nonlins=3,iterations=200,save.models=TRUE,robust=TRUE,maplots=FALSE) 
}

\arguments{
  \item{mbatch}{A \code{marrayRaw} or \code{marrayNorm} batch of arrays.
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
  \item{save.models}{If set to \code{"TRUE"} will enable storage of the models (parameters and ranges of aplicability) in the
    component \code{models} of the list that this function returns.                    
  }
  \item{robust}{If set to \code{"TRUE"}, each spot will be assigned a weight in the model identification, providing resistance to
          outliers.
  }
  \item{maplots}{If set to \code{"TRUE"} will produce a \eqn{M-A} plot for each slide before and after normalization. 
  }
}

\details{This function uses neural networks to model the bias in cDNA data sets.
}


\value{
  A list with components:
  \item{batchn}{A \code{marrayNorm} object containing the normalized log ratios. See \code{marrayNorm} 
    class for details
  }
  \item{models}{A list containing the parameters and ranges of aplicability for the models. This component
     is to be used only as a argument to the function \code{predictBias}. 
  }
}

  


\examples{
# Normalization of swirl data
data(swirl)
# print-tip, intensity and spatial normalization of the first two slides in swirl data set
swirlNN<-maNormNN(swirl[,1:2])$batchn    

#do not consider spatial variations, and display M-A plots before and after normalization
swirlNN<-maNormNN(swirl[,1:2],binWidth=maNsc(swirl),binHeight=maNsr(swirl),maplots=TRUE)$batchn    

}

\author{Tarca, A.L.}
\references{
 A. L. Tarca, J. E. K. Cooke, and J. Mackay. Robust neural networks approach for spatial and
intensity dependent normalization of cDNA data. Bioinformatics. 2004,submitted.\cr
}	

\seealso{\code{\link{predictBias}},\code{\link{compNorm}},\code{nnet}}

\keyword{smooth}
\keyword{robust}


