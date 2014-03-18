\name{rank.btnr}
\alias{rank.btnr}
\title{
Rank responses based on the Bradley-Terry model with Newton-Raphson method
}
\description{
Adopt the Bradley-Terry model to rank responses in a single response question or in a multiple response question with the Newton-Raphson method. This method associates
each response with a value 'gamma', and use this value  to rank responses.
}
\usage{
rank.btnr(data)
}
\arguments{
\item{data}{A m x n matrix dij, where dij = 0 or 1. If the 
ith respondent selects the jth response, then dij = 1, otherwise
dij = 0.}
}
\value{
The rank.btnr returns the associated  values in the first line and the ranks of the responses in the second line.
}
\author{
Hsiuying Wang\email{wang@stat.nctu.edu.tw},Yu-Jun Lin\email{restart79610@hotmail.com}
}
\references{
Hunter DR (2004). MM algorithms for generalized Bradley-Terry models. The Annals
of Statistics, 32, 384-406.
}
\seealso{
\code{\link{rank.btmm}},\code{\link{rank.btqn}},\code{\link{rank.L2R}},\code{\link{rank.LN}},\code{\link{rank.LR}},\code{\link{rank.gs}},\code{\link{rank.wald}}
}
\examples{
## This is an example to rank three responses in a multiple response question when
## the number of respondents is 1000. In this example, we do not use a real data,
## but generate data in the first three lines.
A <-sample.int(2,1000,replace=TRUE,prob=c(0.37,0.63))-1
B <-sample.int(2,1000,replace=TRUE,prob=c(0.65,0.35))-1
C <-sample.int(2,1000,replace=TRUE,prob=c(0.5,0.5))-1
D <-cbind(A,B,C)
data <-matrix(D,nrow=1000,ncol=3)
## or upload the true data
rank.btnr(data)
}