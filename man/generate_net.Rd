\name{generate_net}
\alias{generate_net}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulating networks from preferential attachment and fitness mechanisms}
\description{
This function generates networks from the General Temporal model, a generative temporal network model that includes many well-known models such as the \enc{Erdős–Rényi}{Erdos-Renyi} model, the \enc{Barabási}{Barabasi}-Albert model or the Bianconi-\enc{Barabási}{Barabasi} model as special cases. This function also includes some flexible mechanisms to vary the number of new nodes and new edges at each time-step in order to generate realistic networks. 
}
\usage{
generate_net (N, 
             num_seed           = 2      , 
             multiple_node      = 1      , 
             specific_start     = NULL   ,
             m                  = 1      ,
             prob_m             = FALSE  ,
             increase           = FALSE  , 
             log                = FALSE  , 
             no_new_node_step   = NULL   ,
             m_no_new_node_step = m      ,
             custom_PA          = NULL   ,
             mode               = 1      , 
             alpha              = 1      , 
             beta               = 2      , 
             sat_at             = 100    ,
             offset             = 1      ,
             mode_f             = "gamma", 
             rate               = 0      , 
             shape              = 0      , 
             meanlog            = 0      , 
             sdlog              = 1      ,
             scale_pareto       = 2      ,
             shape_pareto       = 2      )
}

%- maybe also 'usage' for other objects documented here.
\arguments{
The parameters can be divided into four groups. 

The first group specifies basic properties of the network:
  \item{N}{
Integer. Total number of nodes in the network (including the nodes in the seed graph). Default value is \code{1000}.
}
\item{num_seed}{
Integer. The number of nodes of the seed graph (the initial state of the network). The seed graph is a cycle. Default value is \code{2}.
}
\item{multiple_node}{
Positive integer. The number of new nodes at each time-step. Default value is \code{1}.
}
\item{specific_start}{
Positive Integer. If \code{specific_start} is specified, then all the time-steps from time-step \code{1} to \code{specific_start} are grouped to become the initial time-step in the final output. This option is usefull when we want to create a network with a large initial network that follows a scale-free degree distribution. Default value is \code{NULL}.
}
The second group specifies the number of new edges at each time-step:
  \item{m}{
Positive integer. The number of edges of each new node. Default value is \code{1}.
}
 \item{prob_m}{
Logical. Indicates whether we fix the number of edges of each new node as a constant, or let it follows a Poisson distribution. If \code{prob_m == TRUE}, the number of edges of each new node follows a Poisson distribution. The mean of this distribution depends on the value of \code{increase} and \code{log}. Default value is \code{FALSE}.
}
\item{increase}{
Logical. Indicates whether we increase the mean of the Poisson distribution over time. If \code{increase == FALSE}, the mean is fixed at \code{m}. If \code{increase == TRUE}, the way the mean increases depends on the value of \code{log}. Default value is \code{FALSE}.
}
\item{log}{
Logical. Indicates how to increase the mean of the Poisson distribution. If \code{log == TRUE}, the mean increases logarithmically with the number of current nodes. If \code{log == FALSE}, the mean increases linearly with the number of current nodes. Default value is \code{FALSE}.
}
\item{no_new_node_step}{
Positive integer. The number of time-steps in which no new node is added, while new edges are added between existing nodes. Default value is \code{NULL}, i.e. new nodes are added at each time-step.
}
\item{m_no_new_node_step}{
Positive integer. The number of new edges in the no-new-node steps. Default value is equal to \code{m}. 
}


The third group of parameters specifies the preferential attachment function:
\item{custom_PA}{
Numeric vector. This is the user-input PA function: \eqn{A_0, A_1,..., A_K}. If \code{custom_PA} is specified, then \code{mode} is ignored, and we grow the network using the PA function \code{custom_PA}. Degree greater than \eqn{K} will have attachment value \eqn{A_k}. Default value is \code{NULL}.
}
\item{mode}{
Integer. Indicates the attachment function to be used in generating the network. If \code{mode == 1}, the attachment function is \eqn{A_k = k^\alpha}. If \code{mode == 2}, the attachment function is \eqn{A_k = min(k,sat.at)^\alpha}. If \code{mode == 3}, the attachment function is \eqn{A_k = \alpha log (k)^\beta}. Default value is \code{1}.
} 
  \item{alpha}{
Numeric. If \code{mode == 1}, this is the attachment exponent in the attachment function \eqn{A_k = k^\alpha}. If \code{mode == 2}, this is the attachment exponenet in the attachment function \eqn{A_k = min(k,sat.at)^\alpha}. If \code{mode == 3}, this is the \eqn{\alpha} in the attachment function \eqn{A_k = \alpha log (k)^\beta + 1}.
}
\item{beta}{
Numeric. This is the beta in the attachment function \eqn{A_k = \alpha log (k)^\beta + 1}.
}
\item{sat_at}{
Integer. This is the saturation position \eqn{sat.at} in the attachment function \eqn{A_k = min(k,sat.at)^\alpha}.
}

\item{offset}{
Numeric. The attachment value of degree \code{0}. Default value is \code{1}.
}
The final group of parameters specifies the distribution from which node fitnesses are generated:
  \item{mode_f}{
String. Possible values:\code{"gamma"}, \code{"log_normal"} or \code{"power_law"}. This parameter indicates the true distribution for node fitness. \code{"gamma"} = gamma distribution, \code{"log_normal"} = log-normal distribution. \code{"power_law"} = power-law (pareto) distribution. Default value is "gamma".
}
\item{rate}{
Positive numeric. The rate parameter in the Gamma prior for node fitness. If either rate or shape is \code{0}, all node fitnesses \eqn{\eta} are fixed at \code{1} (i.e. Barabasi-Albert model)
}
\item{shape}{
Positive numeric. The shape parameter in the Gamma prior for node fitness. If either rate or shape is \code{0}, all node fitnesses \eqn{\eta} are fixed at \code{1} (i.e. Barabasi-Albert model)
}
  \item{meanlog}{
Numeric. Mean of the log-normal distribution in log scale. Default value is \code{0}.
}
\item{sdlog}{
Positive numeric. Standard deviation of the log-normal distribution in log scale. Default value is \code{1}.
}
\item{scale_pareto}{
Numeric. The scale parameter of the Pareto distribution. Default value is \code{2}.
}
\item{shape_pareto}{
Numeric. The shape parameter of the Pareto distribution. Default value is \code{2}.
}
}

\value{
The output is a List contains the following two fields:
\item{graph}{a three-column matrix, where each row contains information of one edge, in the form of \code{(from_id, to_id, time_stamp)}. \code{from_id} is the id of the source, \code{to_id} is the id of the destination.}
\item{fitness}{fitness values of nodes in the network. The name of each value is the ID of the node.}
}
\author{
Thong Pham \email{thongpham@thongpham.net}
}
\seealso{
For subsequent estimation procedures, see \code{\link{get_statistics}}.

For simpler functions to generate networks from well-known models, see \code{\link{generate_BA}}, \code{\link{generate_ER}}, \code{\link{generate_BB}} and \code{\link{generate_fit_only}}. }

\examples{
library("PAFit")
#Generate a network from the original BA model with alpha = 1, N = 100, m = 1
net <- generate_net(N = 100,m = 1,mode = 1, alpha = 1, shape = 0)
str(net)
}

\concept{fitness model}
\concept{\enc{Barabási}{Barabasi}-Albert model}
\concept{Bianconi-\enc{Barabási}{Barabasi} model}