\name{Jeong}
\alias{Jeong}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Jeong's method for estimating  the preferential attachment function
}
\description{
  This function estimates the preferential attachment function by Jeong's method. 
}
\usage{
Jeong(raw_net       , 
      net_stat      , 
      T_0_start = 0                        ,
      T_0_end   = round(net_stat$T * 0.75) ,
      T_1_start = T_0_end + 1              ,
      T_1_end   = net_stat$T               ,
      interpolate = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{raw_net}{
    a three-column matrix that contains the network.
  }
  \item{net_stat}{
    An object of class \code{PAFit_data} which contains summerized statistics needed in estimation. This object is created by the function \code{\link{GetStatistics}}.
  }
  \item{T_0_start}{Positive integer. The starting time-step of the \code{T_0_interval}. Default value is \code{0}.}
  \item{T_0_end}{Positive integer. The ending time-step of \code{T_0_interval}. Default value is \code{round(net_stat$T * 0.75)}.}
  
  \item{T_1_start}{Positive integer. The starting time-step of the \code{T_1_interval}. Default value is \code{T_0_end + 1}.}
  \item{T_1_end}{Positive integer. The ending time-step of \code{T_1_interval}. Default value is \code{net_stat$T}.}
   
  \item{interpolate}{
    Logical. If \code{TRUE} then all the gaps in the estimated PA function are interpolated by linear interpolating in logarithm scale. Default value is \code{FALSE}.
  }
}
\value{
  Outputs an \code{PA_result} object which contains the estimated attachment function. It also includes the estimated attachment exponenent \eqn{\alpha} (the field \code{alpha}) and the confidence interval of \eqn{\alpha} (the field \code{ci}) when possible.
}
\author{
  Thong Pham \email{thongpham@thongpham.net}
}
\references{
  1. Jeong, H., \enc{Néda}{Neda}, Z. & \enc{Barabási}{Barabasi}, A. . Measuring preferential attachment in evolving networks. Europhysics Letters. 2003;61(61):567–572. doi: 10.1209/epl/i2003-00166-9 (\url{http://iopscience.iop.org/article/10.1209/epl/i2003-00166-9/fulltext/}) .
}
\seealso{

 See \code{\link{GetStatistics}} for how to create summerized statistics needed in this function.

See \code{\link{Newman}} and \code{\link{OnlyA_Estimate}} for other methods to estimate the attachment function in isolation.

}

\examples{
  library("PAFit")
  net        <- GenerateNet(N = 1000 , m = 1 , mode = 1 , alpha = 1 , shape = 0)
  net_stats  <- GetStatistics(net$graph)
  result     <- Jeong(net$graph       , net_stats, 
                      T_0_start = 0   , T_0_end = 700, 
                      T_1_start = 800 , T_1_end = 900)
  # true function
  true_A     <- result$center_k
  #plot the estimated PA function
  plot(result , net_stats)
  lines(result$center_k, true_A, col = "red") # true line
  legend("topleft" , legend = "True function" , col = "red" , lty = 1 , bty = "n")
}
\concept{preferential attachment}
\concept{attachment function}