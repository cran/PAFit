\name{plot.PAFit}
\alias{plot.PAFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to plot the estimated attachment function and node fitness
}
\description{
This function plots the estimated attachment function \eqn{A_k} and node fitness \eqn{f_i}, together with additional information such as their confidence intervals. 
}
\usage{
\method{plot}{PAFit}(x, data, true_f = NULL, plot = c("A","f","true_f"),
plot_bin = TRUE, line = TRUE, confidence = TRUE, high_deg = NULL, shade = 0.5, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class "PAFit", containing the result
}
  \item{data}{
An object of class "PAFitData", containing the summerized statistics.
}
  \item{true_f}{
Vector. Optional parameter for the true value of node fitnesses (only available in simulated datasets). If this parameter is specified and plot == "true_f", a plot of estimated \eqn{f} versus true \eqn{f} is produced (after a suitable rescaling of the estimated \eqn{f}).  
}
  \item{plot}{
String. Indicates which plot is produced. Default value is "A".
}
\item{plot_bin}{Logical. If TRUE then only the center of each bin is plotted. Default is TRUE.}
\item{line}{
Logical. Indicates whether to plot the line fitted from the log-linear model or not. Default value is TRUE.
}
\item{confidence}{
Logical. Indicates whether to plot the confidence intervals of \eqn{A_k} and \eqn{f_i} or not. If confidence == TRUE, a 2-sigma confidence interval will be plotted at each \eqn{A_k} and \eqn{f_i}.
}
\item{high_deg}{Integer. If this parameter is specified, only nodes whose number of edges acquired is not less than \eqn{high_deg} is plotted.}
\item{shade}{
Numeric. Value between 0 and 1. This is the transparency level of the confidence intervals. Default value is 0.5.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\value{
Outputs the desired plot.
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
result <- PAFit(stats)
#plot A
plot(result,stats,plot = "A")
#plot f
plot(result,stats,plot = "f")
#plot true_f
plot(result,stats,data$fitness, plot = "true_f")
}
