\name{UCIrvine}
\alias{UCIrvine.data}
\docType{data}
\title{
A message sending-receiving network
}
\description{
Message sending-receiving network of an online community in which users are students at University of California, Irvine. 
}
\usage{data(UCIrvine)}
\format{
  The format is:
 int [1:59835, 1:3] 
}
\details{
 There are total 59,835 messages of 1899 users.  The data format is a 59,835x3 matrix, which each row describes a triple (sending user ID, receiving user ID, time stamp) of a message.
}
\source{
http://toreopsahl.com/datasets/
}
\references{
Opsahl, T. and Panzarasa, P. (2009) Clustering in weighted networks. \emph{Social Networks}, \bold{31}(\bold{2}), 155--163, doi: 10.1016/j.socnet.2009.02.002
}
\examples{
data(UCIrvine)
}
\keyword{datasets}
\keyword{scale_free}
