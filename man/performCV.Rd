\name{performCV}
\alias{performCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to perform cross validation.
}
\description{
A function to perform cross validation.
}
\usage{performCV(cv_data,r = 10^c(-2,-1,0,1,2),s = 10^c(-2,-1,0,1,2), 
                      stop_cond = 10^-7, only_PAFit = TRUE, 
                      silent = FALSE, only_linear = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cv_data}{
Cross-validation data.}
\item{r}{
Numeric vector. Specifying the values of the regularization parameter \eqn{r} of the PA function. Default value is the vector \eqn{c(0.01,0.1,1,10,100)}.
}
\item{s}{
Numeric vector. Specifying the values of the regularization parameter \eqn{s} of node fitness. Default value is the vector \eqn{c(0.01,0.1,1,10,100)}.
}
\item{stop_cond}{
Numeric. The iterative algorithm will halt if the relative difference of the objective function falls below this value.
}
\item{only_PAFit}{
Logical. If TRUE then only perform the CV for PAFit full model. Default is TRUE.
}
\item{silent}{
Logical. If TRUE then the progress is not printed out. Default is FALSE.
}
\item{only_linear}{
Logical. If TRUE then only perform the CV assuming the linear functional form \eqn{A_k = k}. Default is FALSE. 
}
}
\value{An object of class "CV_Result" containing the cross validation result.
}

\author{
Thong Pham \email{thongpham@thongpham.net}
}
\references{
1. Pham, T., Sheridan, P. & Shimodaira, H. (2016). Nonparametric Estimation of the Preferential Attachment Function in Complex Networks: Evidence of Deviations from Log Linearity, Proceedings of ECCS 2014, 141-153 (Springer International Publishing) (\url{http://dx.doi.org/10.1007/978-3-319-29228-1_13}).

2. Pham, T., Sheridan, P. & Shimodaira, H. (2015). PAFit: A Statistical Method for Measuring Preferential Attachment in Temporal Complex Networks. PLoS ONE 10(9): e0137796. doi:10.1371/journal.pone.0137796 (\url{http://dx.doi.org/10.1371/journal.pone.0137796}).

3. Pham, T., Sheridan, P. & Shimodaira, H. (2016). Joint Estimation of Preferential Attachment and Node Fitness in Growing Complex Networks. Scientific Reports 6, Article number: 32558. doi:10.1038/srep32558   (\url{www.nature.com/articles/srep32558}).
}

\examples{
library("PAFit")
data      <- GenerateNet(N = 100,m = 5,mode = 1, alpha = 1, shape = 5, rate = 5)
data_cv   <- CreateDataCV(data$graph)
cv_result <- performCV(data_cv,r = c(0.1), s = c(5),only_PAFit = TRUE,stop_cond = 10^-3)
}

\keyword{fitness model}
\keyword{Preferential Attachment function}
\keyword{Minorize-Maximization algorithms}
\keyword{Bianconi-Barabasi model}
\keyword{Barabasi-Albert model}% __ONLY ONE__ keyword per line
