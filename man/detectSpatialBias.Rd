\name{detectSpatialBias}
\alias{detectSpatialBias}
\title{Detecting spatial bias within the print-tips of a two channel array}
\description{This function allows to identify in two channel batch of arrays, which are the print-tips where 
spatial bias is present.}
\usage{
detectSpatialBias(mbatch, corThreshold=0.6) 
}

\arguments{
  \item{mbatch}{A \code{marrayRaw} or \code{marrayNorm} batch of two channel arrays.
  }
  \item{corThreshold}{The correlation treshold to be used. 
  }	
 }

\details{This function computes two matrices: \code{biasRow} and \code{biasCol}. The elements of these 
matrices represent the fraction of rows (columns) for which the correlation coefficient between log-ratios, M, and 
column index (row index) is higher than a user specified treshold (default corThreshold=0.6).
The idea here is to see in which print-tip a important fraction of the rows (columns) are highly correlated with the 
column (row) index. Since some rows (columns) will show positive correlation while the other negative correlation, we 
are only interested in a sigle direction of the correlation, i.e. either positive or negative. 
}


\value{
  This function returns a list with two matrices. \code{biasRow} and \code{biasCol}. The rows of these matrices correspond
to the print tips counted metaRow wise, and the columns correspond to arrays. Values in these matrices superior to 33 
point to print-tips that have more tha a third of the rows (columns) with important spatial bias.  
}


 
\examples{
# detecting spatial bias in swirl data
data(swirl)
# print-tip, intensity and spatial normalization of the first slide in swirl data set
myres<-detectSpatialBias(swirl)   
}

\author{Tarca, A.L.}
\references{
A robust neural networks approach for spatial and intensity dependent normalization of cDNA microarray data,
 Adi. L. Tarca , Janice. E. K. Cooke, and John Mackay, Bioinformatics, 21, 2005, 2674 - 2683.\cr
}	

\seealso{\code{\link{maNormNN}}}
\keyword{spatial}


