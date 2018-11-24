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
result2 <- sapply(alphaj, function(alphajk) rnbinom(10000, phij*alphajk*betaj, phij/(1+phij)))

signif(apply(result, 2, mean),3)
signif(Sj*(1+h)^n, 3)
signif(apply(result2, 2, mean), 3)

signif(apply(result, 2, var), 3)
signif(Sj*(1+h)^n*(1-h)/(1+h)*((1+h)^n-1), 3)
signif(apply(result2, 2, var), 3)

sumResult <- apply(result, 1, sum)
normResult <- apply(result, 2, "/", sumResult)

apply(normResult, 2, function(x) {windows(); hist(x)})

countResult <- apply(normResult, 1, function(x) rmultinom(1, 1500, prob=x))

apply(countResult, 1, function(x) {windows(); hist(x)})


apply(countResult, 1, mean)
apply(countResult, 1, var)


# not overdispersed to have a multinomial with set probabilities
Z <- rmultinom(10000, 1500, prob=c(0.5, 0.3, 0.2))

apply(Z, 1, function(x) {windows(); hist(x)})

apply(Z, 1, mean)
apply(Z, 1, var)

# overdispersed to sample a dirichlet first
result <- sapply(c(2,3,4), function(x) rgamma(10000, x, 1))
sumResult <- apply(result, 1, sum)
normResult <- apply(result, 2, "/", sumResult)

apply(normResult, 2, function(x) {windows(); hist(x)})

countResult <- apply(normResult, 1, function(x) rmultinom(1, 1500, prob=x))

apply(countResult, 1, function(x) {windows(); hist(x)})

apply(countResult, 1, mean)
apply(countResult, 1, var)




##############################################################3

D <- 151

m <- c(.6,.4)
X1j <- c(0,1,1)
X2j <- c(1,1,0)

Sj <- (m[1]*X1j + m[2]*X2j)*D

n <- 30
h <- .8
g <- 0.00001

phij <- 1/((1-h)/(1+h)*((1+h)^n-1)-1)

alphaj <- Sj*(1+h)^n*phij*g
betaj <- 1/phij

result <- matrix(NA, nrow=10000, ncol=3)
for (i in 1:10000)
{
  lambdaj <- sapply(alphaj, function(alphajk) rgamma(1, alphajk, scale=betaj) )
  pcrj <- sapply(lambdaj, function(lambdajk) rpois(1, lambdajk))
  result[i,] <- pcrj
}

signif(apply(result, 2, mean),3)
signif(Sj*(1+h)^n, 3)
signif(Sj*(1+h)^n*g, 3)

signif(apply(result, 2, var), 3)
signif(Sj*(1+h)^n*(1-h)/(1+h)*((1+h)^n-1), 3)
signif(Sj*(1+h)^n*(1-h)/(1+h)*((1+h)^n-1), 3)
