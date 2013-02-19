D <- 151

m <- c(.6,.4)
X1j <- c(0,1,1)
X2j <- c(1,1,0)

Sj <- (m[1]*X1j + m[2]*X2j)*D

n <- 30
h <- .8

phij <- 1/((1-h)/(1+h)*((1+h)^n-1)-1)

alphaj <- Sj*(1+h)^n*phij
betaj <- 1/phij

result <- matrix(NA, nrow=10000, ncol=3)
for (i in 1:10000)
{
    lambdaj <- sapply(alphaj, function(alphajk) rgamma(1, alphajk, scale=betaj) )
  
    pcrj <- sapply(lambdaj, function(lambdajk) rpois(1, lambdajk))
    result[i,] <- pcrj
}

signif(apply(result, 2, mean),3)
signif(apply(result, 2, var), 3)

signif(Sj*(1+h)^n, 3)
signif(Sj*(1+h)^n*(1-h)/(1+h)*((1+h)^n-1), 3)

subresult <- 1/1E4*result

signif(apply(subresult, 2, mean),3)
signif(apply(subresult, 2, var), 3)

