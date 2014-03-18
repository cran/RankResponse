\name{rank.wald}
\alias{rank.wald}
\title{
Rank responses based on the Wald test
}
\description{
Rank responses of a single response question or a multiple response question by the Wald test procedure.
}
\usage{
rank.wald(data, alpha, type=2)
}
\arguments{
\item{data}{A m x n matrix dij, where dij = 0 or 1. If the 
ith respondent selects the jth response, then dij = 1, otherwise
dij = 0.}
\item{alpha}{The significance level used in the Wald test.}
\item{type}{type=1 for a single response question ;type=2 for a multiple response question .}
}
\value{
The rank.wald returns the estimated probabilities of the responses being selected and the ranks of the responses by the Wald trest procedure.
}
\author{
Hsiuying Wang\email{wang@stat.nctu.edu.tw},Yu-Jun Lin\email{restart79610@hotmail.com}
}
\references{
Wang, H. (2008). Ranking Responses in Multiple-Choice Questions. Journal of Applied Statistics, 35, 465-474.
}
\seealso{
\code{\link{rank.btmm}},\code{\link{rank.btnr}},\code{\link{rank.btqn}},\code{\link{rank.L2R}},\code{\link{rank.LN}},\code{\link{rank.LR}},\code{\link{rank.gs}}
}
\examples{
## This is an example to rank three responses in a multiple response question when
## the number of respondents is 1000 and the significance level is 0.05. In this 
## example,we do not use a real data, but generate data in the first three lines.
A <-sample.int(2,1000,replace=TRUE,prob=c(0.21,0.79))-1
B <-sample.int(2,1000,replace=TRUE,prob=c(0.81,0.19))-1
C <-sample.int(2,1000,replace=TRUE,prob=c(0.62,0.28))-1
D <-cbind(A,B,C)
data <-matrix(D,nrow=1000,ncol=3)
## or upload the true data
alpha<-0.05
rank.wald(data,alpha,2)
}
