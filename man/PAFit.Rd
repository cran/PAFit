\name{PAFit}
\alias{PAFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to estimate the attachment function and node fitness in temporal complex networks
}
\description{
From a PAFitData object, which contains summary statistics summerized from the dataset, PAFit estimates the attachment function \eqn{A_k} and node fitness \eqn{f_i} by penalized log-likelihood maximization.  It also calculates confidence intervals for \eqn{A_k} and \eqn{f_i}.  Estimation of the attachment function in isolation (while fixing \eqn{f_i = 1}) can be specified by setting only_PA = TRUE.  Estimation of node fitness in isolation (fixing either \eqn{A_k = k} or \eqn{A_k = 1}) can be specified by setting only_f = TRUE.
}
\usage{
PAFit(data, only_PA , only_f, mode_f, true_A, true_f,
shape, rate, auto_lambda, ratio, lambda, weight_PA_mode, auto_stop, stop_cond,
iteration , max_iter, debug , step_size, q , normalized_f, interpolate,...)
}
\arguments{
  \item{data}{
An object of class "PAFitData" containing all the summary statistics summerized from the data by the function GetStatistics
}
\item{only_PA}{
Logical. TRUE means that the attachment function \eqn{A_k} is estimated in isolation(fixing \eqn{f_i = 1}). Default is FALSE.
}
\item{only_f}{
Logical. TRUE means that the fitness function is estimated in isolation. Default is FALSE.
}
\item{mode_f}{
String. Only effective when only_f == TRUE. If mode_f == "Linear_PA" then \eqn{A_k = k}. If mode_f == "Constant_PA" then \eqn{A_k = 1}. Default values is "Linear_PA".
}
\item{true_A}{

}
\item{true_f}{

}
\item{shape}{
Numeric. The shape parameter in the Gamma prior for \eqn{f_i}. Default value is \eqn{0.1}.
}
  \item{rate}{
Numeric. The rate parameter in the Gamma prior for \eqn{f_i}. Default value is \eqn{0.1}. Either Shape == 0 or Rate == 0 means no regularization for \eqn{f_i}.
}
\item{auto_lambda}{
Logical. If auto_lambda == TRUE, lambda will be determined automatically from the data. Default is TRUE.
}
\item{ratio}{
Numeric. The ratio between strength of the regularization term and the total number of data. Default is \eqn{0}.
}

\item{lambda}{
Numeric. The strength of the regularization for PA function. Ignored when auto_lambda == TRUE. Default value is \eqn{1}. lambda == 0 means no regularization for \eqn{A}. 
}
\item{weight_PA_mode}{
Integer. Indicates how the regularization terms for \eqn{A_k} are weighted. If weight_PA_mode == 0, the regularization term for \eqn{A_k} is weighted by the total number of edges connected to degree \eqn{k} nodes. If weight_PA_mode == 1, the regularization terms have uniform weights. Default value is \eqn{0}.
}
\item{auto_stop}{Logical. Indicates whether the algorithm stop automatically or not. Default is TRUE}
\item{stop_cond}{Numeric. If \eqn{auto_stop = TRUE}, the iterative algorithm stops when \eqn{abs(l(i) - l(i+1)/l(i)) < stop_cond}. Default is \eqn{10^-4}.}
\item{iteration}{
Numeric. The number of iterations. Ignored if auto_stop == TRUE. Default value is \eqn{20}. 
}
\item{max_iter}{
Numeric. The maximum number of iterations. Regardless of other settings, the algorithm will stop once the number of  iterations reaches this threshold. Default value is \eqn{1000}.
}
\item{debug}{Logical. if debug == TRUE, the value of the log likelihood function is printed out at each step. Defaule is FALSE.}
\item{step_size}{Numeric. A numeric value between \eqn{0} and \eqn{1}. Determines the convex combination of the default direction given by the MM update and the quasi-Newton direction. Only effective if \eqn{q \ge 2}. \eqn{step_size = 1} means we ignore the default MM direction. Default value is \eqn{1}.}
\item{q}{Integer. The number of previous values used in the quasi-Newton acceleration of Zhou et al. A quasi-newton acceleration for high-dimensional optimization algorithms (2011). The acceleration is activated only when \eqn{q \ge 2}. Default value is \eqn{1}.}
\item{normalized_f}{Logical. Indicates whether we should normalize the estimated value of \eqn{f} after estimation. Default value is FALSE.}
\item{interpolate}{Logical. Indicates whether we should perform interpolation for the missing values of \eqn{A_k}. Default value is FALSE.}
\item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\value{
an object of class "PAFit", which is a list. Some important fields of this object:
\item{A }{The estimated attachment function}
\item{k}{The corresponding degree}
\item{var_A}{Variances of the estimated A}
\item{linear_fit}{Result of fitting the log-linear model \eqn{log A_k = \alpha log k} to the estimated \eqn{A_k}}
\item{alpha}{The estimated attachment exponent of the log-linear model \eqn{A_k =k^\alpha}}
\item{weight_of_A}{The number of A in each bin}
\item{var_logA}{Variances of log(A)}        
\item{upper_A}{The upper value of the 2-sigma confidence interval of A}
\item{lower_A}{The lower value of the 2-sigma confidence interval of A}
\item{center_k}{The logarithmic center of the bins}
\item{theta}{Attachment value of the bins (before converting back to \eqn{A_k})}
\item{upper_bin}{The upper value of the 2-sigma confidence interval of theta}
\item{lower_bin}{The lower value of the 2-sigma confidence interval of theta}
\item{f}{The estimated node fitness}
\item{var_f}{Variances of the estimated f}
\item{upper_f}{The upper value of the 2-sigma confidence interval of \eqn{f}}
\item{lower_f}{The lower value of the 2-sigma confidence interval of \eqn{f}}
\item{log_likelihood}{Values of the (penalized) log-likelihood function recorded at each iteration}
}

\author{
Thong Pham
\email{thongpham@thongpham.net}
}

\references{
1. Pham, T. and Sheridan, P. and Shimodaira, H. (2015). Nonparametric estimation of the preferential attachment function in complex networks: evidence of deviations from log linearity, in press. Proceedings of ECCS 2014: European Conference on Complex Systems.

2. Pham T, Sheridan P, Shimodaira H (2015) PAFit: A Statistical Method for Measuring Preferential Attachment in Temporal Complex Networks. PLoS ONE 10(9): e0137796. doi:10.1371/journal.pone.0137796 (\url{http://dx.doi.org/10.1371/journal.pone.0137796})
}


\examples{
library("PAFit")
data   <- GenerateNet(N = 1000,m = 1, mode = 1, alpha = 1, shape = 5, rate = 5)
stats  <- GetStatistics(data$graph)
result <- PAFit(stats)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ package }
\keyword{fitness model}
\keyword{Preferential Attachment function}
\keyword{Minorize-Maximization algorithms}
\keyword{Barabasi-Albert model}
\keyword{MM algorithims}
\keyword{Bianconi-Barabasi model}
\keyword{fitness}
\keyword{maximum likelihood function}% __ONLY ONE__ keyword per line
