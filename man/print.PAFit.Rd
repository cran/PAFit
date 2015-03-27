\name{print.PAFit}
\alias{print.PAFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to print a PAFit object
}
\description{
This function print a PAFit object's information.
}
\usage{
\method{print}{PAFit}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class "PAFit", containing the result
}
   \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}

\value{
Outputs some simple information. 
}

\author{
Thong Pham \email{thongpham@thongpham.net}
}
\references{
Pham, T. and Sheridan, P. and Shimodaira, H. (2015) Nonparametric estimation of the preferential attachment function in complex networks: evidence of deviations from log linearity (submitted).
}

\examples{
library("PAFit")
data   <- GenerateNet(N = 1000,m = 1, mode = 1, alpha = 1, shape = 5, rate = 5)
stats  <- GetStatistics(data$graph)
result <- PAFit(stats)
print(result)
}
