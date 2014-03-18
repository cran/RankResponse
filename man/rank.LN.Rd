\name{rank.LN}
\alias{rank.LN}
\title{
Rank responses under the Bayesian framework according to the loss function in Method 1 of Wang and Huang (2004).
}
\description{
Rank responses of a single response question or a multiple response question under the Bayesian framework according to the loss function in Method 1 of Wang and Huang (2004).
}
\usage{
rank.LN(data,response.number,prior.parameter,c)
}
\arguments{
\item{data}{A m x n matrix dij, where dij = 0 or 1. If the 
ith respondent selects the jth response, then dij = 1, otherwise
dij = 0.}
\item{response.number}{The number of the responses.}
\item{prior.parameter}{The parameter vector of the Dirichlet prior distribution,
where the vector dimension is 2^response.number.}
\item{c}{The value of c in the loss function}
}
\value{
The rank.LN returns the estimated probabilities of the responses being selected in the first line and the ranks of the responses in the second line.
}
\author{
Hsiuying Wang\email{wang@stat.nctu.edu.tw},Yu-Jun Lin\email{restart79610@hotmail.com}
}
\references{
Wang, H. and Huang, W. H. (2014). Bayesian Ranking Responses in Multiple Response
Questions. Journal of the Royal Statistical Society: Series A (Statistics in Society), 177, 191-208.
}
\seealso{
\code{\link{rank.btmm}},\code{\link{rank.btnr}},\code{\link{rank.btqn}},\code{\link{rank.L2R}},\code{\link{rank.LR}},\code{\link{rank.gs}},\code{\link{rank.wald}}
}
\examples{
## This is an example to rank three responses in a multiple response
## question when the number of respondents is 1000 and the value e2R 
## is 0.15. In this example, we do not use a real data, but generate 
## data in the first three lines.
A <-sample.int(2,1000,replace=TRUE,prob=c(0.37,0.63))-1
B <-sample.int(2,1000,replace=TRUE,prob=c(0.71,0.29))-1
C <-sample.int(2,1000,replace=TRUE,prob=c(0.22,0.78))-1
D <-cbind(A,B,C)
data <-matrix(D,nrow=1000,ncol=3)
## or upload the true data
response.number <-3
prior.parameter <- c(5,98,63,7,42,7,7,7)
c <- 0.05
rank.LN(data,response.number,prior.parameter,c)
}