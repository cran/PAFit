\name{CreateDataCV}
\alias{CreateDataCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to create cross-validation data.
}
\description{
A function to create cross-validation data.
}
\usage{CreateDataCV(data, ratio = 0.75, G = 50, 
                        net_type = "directed",deg_thresh = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Matrix. A 3-column matrix whose each row contains information of one edge in the form (from_node id, to_node id, time_stamp). from_node id is the id of the source node. to_node id is the id of the destination node. time_stamp is the arrival time of the edge. from_node id and to_node id are assumed to be integers starting from \eqn{0}. time_stamp can be either numeric or string. We assume that a smaller time_stamp represents an earlier arrival time.}
\item{ratio}{
Numeric. Indicates the ratio of number of new edges in the learning data to that of the testing data. Default is \eqn{0.75}.
}
\item{G}{
Integer. Number of bins. Default value is \eqn{50}.
}
 \item{net_type}{
String. The type of the network: "directed" or "undirected". Default is "directed". 
}
\item{deg_thresh}{
Integer. We only consider nodes with degree at least this threshold. Default value is \eqn{1}.
}

}

\value{An object of class "CV_Data" containing the data needed for cross validation.
}


\author{
Thong Pham \email{thongpham@thongpham.net}
}
\references{
1. Pham, T., Sheridan, P. & Shimodaira, H (2016). Nonparametric Estimation of the Preferential Attachment Function in Complex Networks: Evidence of Deviations from Log Linearity, Proceedings of ECCS 2014, 141-153 (Springer International Publishing) (\url{http://dx.doi.org/10.1007/978-3-319-29228-1_13}).

2. Pham T, Sheridan P, Shimodaira H (2015) PAFit: A Statistical Method for Measuring Preferential Attachment in Temporal Complex Networks. PLoS ONE 10(9): e0137796. doi:10.1371/journal.pone.0137796 (\url{http://dx.doi.org/10.1371/journal.pone.0137796})
}

\examples{
library("PAFit")
data     <- GenerateNet(N = 100,m = 1,mode = 1, alpha = 1, shape = 5, rate = 5)
data_cv  <- CreateDataCV(data$graph)
}

\keyword{fitness model}
\keyword{Preferential Attachment function}
\keyword{Minorize-Maximization algorithms}
\keyword{Bianconi-Barabasi model}
\keyword{Barabasi-Albert model}% __ONLY ONE__ keyword per line
