% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rank.wald.R
\name{rank.wald}
\alias{rank.wald}
\title{Rank Responses based on the Wald Test}
\usage{
rank.wald(data, alpha = 0.05, ranktype = 1)
}
\arguments{
\item{data}{A m by n matrix \eqn{d_{ij}}, where \eqn{d_{ij}} = 0 or 1.
If the ith respondent selects the jth
response, then \eqn{d_{ij}} = 1, otherwise \eqn{d_{ij}} = 0.}

\item{alpha}{The significance level is used to control the type I error rate.
The default is 0.05.}

\item{ranktype}{A numerical value specifies which type of ranking
method is used. The default is 1 (see 'Details').}
}
\value{
rank.wald returns a table contains the estimated probabilities of the responses being selected in the first line and
        the ranks of the responses in the second line.
}
\description{
Rank responses of a single response question
or a multiple response question by the wald test procedure.
}
\details{
Suppose that the question has k responses.
Let \eqn{\pi_{j}} denote the probability that the jth response is selected.
Using the survey data, \eqn{\pi_{j}}  can be estimated.

If \code{ranktype} is 1, the ranking rule is the following steps.
Let \eqn{\pi_{(j)}} denote the order statistic.
If the hypothesis \eqn{\pi_{(k)}} = \eqn{\pi_{(k-1)}} is rejected,
 we rank the response corresponding to \eqn{\pi_{(k)}} first.
If it is not rejected, we compare \eqn{\pi_{(k)}} with \eqn{\pi_{(j)}}
, \eqn{j \le k-2} sequentially.

If \code{ranktype} is 2, the rank of the ith response can be defined as
\deqn{
R_{i} = k - \sum_{j=1, j\ne i}^{k} I(\pi_{i} > \pi_{j})
}
}
\examples{
set.seed(12345)
# This is an example to rank k responses in a multiple response question
# when the number of respondents is 1000.
# In this example, we do not use a real data, but generate data in the first six lines.
k <- 5
data <- matrix(NA, nrow = 1000, ncol = k)
for(i in 1:k){
  p <- runif(1)
  data[, i] <- sample(c(0, 1), 1000, p = c(p, 1-p), replace = TRUE)
}
## or upload the true data
rank.wald(data)

}
\references{
Wang, H. (2008). Ranking Responses in Multiple-Choice Questions.
            Journal of Applied Statistics, 35, 465-474.

            Wang, H. and Huang, W. H. (2014). Bayesian Ranking Responses in Multiple Response Questions.
            Journal of the Royal Statistical Society: Series A (Statistics in Society), 177, 191-208.
}
\seealso{
\code{\link{rankL2R}}, \code{\link{rankLN}}, \code{\link{rank.gs}}
}
\author{
Hsiuying Wang \email{wang@stat.nycu.edu.tw}
, Wan-Ting Huang \email{wthuang.sc09@nycu.edu.tw}
, Yu-Chun Lin \email{restart79610@hotmail.com}
}
