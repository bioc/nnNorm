\name{compNorm}
\alias{compNorm}
\title{Compares the distribution of several vectors at a time using either boxplots or density curves}
\description{This function was concieved to easily compare several normalization methods in terms 
	of variability of log-ratios, M. Basically it produces two plots: The first is a the density 
	plot of the several matrices passed as arguments, while the second is a box plot. 
	Median of absolute deviations for each method is printed on screen.
}
\usage{
compNorm(x,...,bw="AUTO",xlim=c(-3,3),titles="AUTO",type="d") 
}
\arguments{
  \item{x}{A vector of numerical values, e.q. the \eqn{M} values of a data set: \code{as.vector(maM(swirl))}.
  }
  \item{...}{An undefined number of objects similar with \code{x}.
  }	
  \item{bw}{Band width required to compute the density distribution. \code{"AUTO"} 
     will adjust \code{bw} to a suitable value.
  }
  \item{xlim}{The range for abscissa of the density plots. 
  }
  \item{titles}{Names to be displayed the charts legend. \code{"AUTO"} will use the matrices names passed as arguments.              .
  }
  \item{type}{If set to \code{"d"}, density plot will be shown; if set to \code{"d"} box plot will be shown.
  }

  
}

\details{This function is used to compare the normalized log ratios \eqn{M} obtained with several normalization methods.
}


\value{NULL, this function only displays charts and prints on the screen some statistics.
}

  
\examples{
 # Normalize swirl data with two methods
 data(swirl)
 swirlNN<-maNormNN(swirl[,1])
 swirlLoess<-maNormMain(swirl[,1]) 
 nms<-c("None","Loess","NNets")
 #compare distributions: density plot
 compNorm(as.vector(maM(swirl[,1])),as.vector(maM(swirlLoess)),as.vector(maM(swirlNN)),xlim=c(- 2,2),bw="AUTO",titles=nms,type="d")
 #compare distributions: box plot
 compNorm(as.vector(maM(swirl[,1])),as.vector(maM(swirlLoess)),as.vector(maM(swirlNN)),xlim=c(- 2,2),bw="AUTO",titles=nms,type="b")

}
\author{Tarca, A.L.}
\references{
 A. L. Tarca, J. E. K. Cooke, and J. Mackay. Robust neural networks approach for spatial and
intensity dependent normalization of cDNA data. Bioinformatics. 2004,submitted.\cr

}	

\seealso{\code{\link{maNormNN}}}
\keyword{univar}



