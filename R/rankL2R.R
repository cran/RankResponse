#' @title Rank responses under the Bayesian framework according to
#' the loss function in Method 3 of Wang and Huang (2004).
#'
#' @description Rank responses of a single response question or
#' a multiple response question under the Bayesian framework according to the
#' loss function in Method 3 of Wang and Huang (2004).
#'
#' @usage rankL2R(data, response.number, prior.parameter, e)
#'
#' @param data A m by n matrix \eqn{d_{ij}}, where \eqn{d_{ij}} = 0 or 1.
#'             If the ith respondent selects the jth
#'             response, then \eqn{d_{ij}} = 1, otherwise \eqn{d_{ij}} = 0.
#'
#' @param response.number The number of the responses.
#' @param prior.parameter The parameter vector of the Dirichlet prior distribution,
#' where the vector dimension is 2^{response.number}.
#'
#' @param e A cut point used in the loss function which depends on the economic costs.
#'
#' @export rankL2R
#'
#' @return The rankL2R returns the estimated probabilities of the responses
#' being selected in the first line and the ranks of the responses in the second line.
#'
#' @author Hsiuying Wang \email{wang@stat.nycu.edu.tw}
#' , Yu-Chun Lin \email{restart79610@hotmail.com}
#'
#' @references Wang, H. and Huang, W. H. (2014). Bayesian Ranking Responses in Multiple Response Questions.
#'             Journal of the Royal Statistical Society: Series A (Statistics in Society), 177, 191-208.
#'
#' @seealso \code{\link{rankLN}}, \code{\link{rank.wald}}, \code{\link{rank.gs}}
#'
#' @import stats
#' @importFrom stats kmeans pnorm
#'
#' @examples
#' set.seed(12345)
#' # This is an example to rank k responses in a multiple response question
#' # when the number of respondents is 1000 and the value e is 0.15.
#' # In this example, we do not use a real data, but generate data in the first six lines.
#' k <- 3
#' data <- matrix(NA, nrow = 1000, ncol = k)
#' for(i in 1:k){
#'   p <- runif(1)
#'   data[, i] <- sample(c(0, 1), 1000, p = c(p, 1-p), replace = TRUE)
#' }
#' ## or upload the true data
#' response.number <- 3
#' prior.parameter <- c(5, 98, 63, 7, 42, 7, 7, 7)
#' e <- 0.15
#' rankL2R(data, response.number, prior.parameter, e)
#'

rankL2R <- function(data, response.number, prior.parameter, e)
{
      data <- as.matrix(data)
      v <- response.number
      alpha <- prior.parameter
      A <- matrix(0, 2^v, v)
      for(j in 1:v)
      {
		    A[,j] <- rep(c(rep(0, 2^(v-j)), rep(1, 2^(v-j))), 2^(j-1))
	    }

      Temp <- rbind(data, A)
      group <- kmeans(Temp,(2^v))
      obs <- numeric(2^v)
      for(j in 1:(2^v))
      {
           for(i in 1:(2^v))
           {
              if(all(group$center[j,]==A[i,])==T){obs[i]=(group$size[j])-1}
           }
      }

      n <- (v-1)*v
      a1 <- numeric(n)
      a2 <- numeric(n)
      for(i in 1:v)
      {
		    a1[((i-1)*(v-1)+1):((i-1)*(v-1)+(v-1))] <- rep(i,(v-1))
	    }

      a2[1:(v-1)] <- c(2:v)
	    a2[((v-1)*(v-1)+1):((v-1)*(v-1)+(v-1))] <- c(1:(v-1))
	    for(i in 2:(v-1)){
		    a2[((v-1)*(i-1)+1):((v-1)*(i-1)+(v-1))] <- c(1:(i-1), (i+1):v)
	    }

	AA <- cbind(A, alpha, obs)

	V <- matrix(0, 1, n)

	for(i in 1:n){

		alpha0 <- sum(AA[, (v+1)] + AA[, (v+2)])
		alphai <- AA[((AA[, a1[i]]!= AA[, a2[i]])& AA[, a1[i]]==1),][, (v+1)]+AA[((AA[, a1[i]]!=AA[, a2[i]])& AA[, a1[i]]==1),][, (v+2)]
		alphaj <- AA[((AA[, a1[i]]!= AA[, a2[i]])& AA[, a2[i]]==1),][, (v+1)]+AA[((AA[, a1[i]]!=AA[, a2[i]])& AA[, a2[i]]==1),][, (v+2)]

		E <- sum(alphai - alphaj) / alpha0

		Var1 <- sum(alphai * alpha0 - alphai ^ 2)
		Var2 <- sum(alphaj * alpha0 - alphaj ^ 2)
		Var3 <- 2 * sum(alphai %*% t(alphaj))

		Var <- (Var1+Var2+Var3) / (alpha0^3+alpha0)
		V[1, i] <- pnorm(-E/sqrt(Var), 0, 1)
	}

	tv <- round(V[1,], 6)
	pt <- cbind(a1, a2, tv)

      pi <- numeric(v)
      N <- sum(AA[, (v+2)])
      for(i in 1:v){

            pi[i] <- sum(AA[AA[, i]==1, (v+2)]) / N
      }

      names(pi) <- c(1:v)
      pi_temp <- sort(pi, decreasing = TRUE)
      rank_temp <- as.numeric(names(pi_temp))

      u <- numeric(v-1)
      for(i in 1:(v-1)){

            u[i] <- pt[pt[, 1]==rank_temp[i] & pt[, 2]==rank_temp[i+1], 3]
      }


      epsilon <- 0.001
	FD_bar <- numeric(1) ; FDR_bar <- numeric(100)

      for(j in 1:100){
  		  t <- 0.01*j
  		  D <- numeric(1) ; FD_bar <- numeric(1)
  		  for(i in 1:(v-1)){
  			  if(u[i] > t){
  				  FD_bar <- FD_bar + (1-u[i])
  				  D <- D+1
  			  }
  		  }

		    FDR_bar[j] <- FD_bar / (D+epsilon)
	    }

	for(i in 1:100){
		if(FDR_bar[i] <= e){break}
		t <- i+1
	}

      tD_star <- t*0.01

      sumI_temp =! (tv >= tD_star)
      sumI <- matrix(sumI_temp, ncol=v-1, nrow=v, byrow=TRUE)
      rank <- v - apply(sumI, 1, sum)

      probability <- apply(data, 2, mean)
      result <- rbind(probability, rank)

      return(result)
}
