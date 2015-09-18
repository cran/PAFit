\name{summary.PAFitData}
\alias{summary.PAFitData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to summarize a PAFitData object 
}
\description{
This function summarizes information of a PAFitData object. 
}
\usage{
\method{summary}{PAFitData}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class "PAFitData", containing the summary statistics needed in estimation.
}
\item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\value{
Outputs summary information.
}

\author{
Thong Pham \email{thongpham@thongpham.net}
}
\references{
1. Pham, T. and Sheridan, P. and Shimodaira, H. (2015). Nonparametric estimation of the preferential attachment function in complex networks: evidence of deviations from log linearity, in press. Proceedings of ECCS 2014: European Conference on Complex Systems.

2. Pham T, Sheridan P, Shimodaira H (2015) PAFit: A Statistical Method for Measuring Preferential Attachment in Temporal Complex Networks. PLoS ONE 10(9): e0137796. doi:10.1371/journal.pone.0137796 (\url{http://dx.doi.org/10.1371/journal.pone.0137796})
}

\examples{
library("PAFit")
data   <- GenerateNet(N = 1000,m = 1,mode = 1, alpha = 1, shape = 5, rate = 5)
stats  <- GetStatistics(data$graph)
summary(stats)
}
