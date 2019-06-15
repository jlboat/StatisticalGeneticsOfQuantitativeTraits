y    <- c(1339, 154, 151, 1195)
m    <- length(y)
n    <- sum(y)
nsim <- 10000

#-----Log-likelihood function-----
lambda < function(p1, p2, p3, x1, x2, x3){
    (x1*log(p1) + x2*log(p2) + x3*log(p3) +
    (n - x1 - x2 - x3) * log(1 - p1 - p2 - p3)) 
}

p1hat <- (y[1] + y[4]) / (2*n)
p2hat <- (y[2] + y[3]) / (2*n)

Lnull <- lambda(p1hat, p2hat, p2hat, y[1], y[2], y[3])
Lmax  <- lambda(y[1]/n, y[2]/n, y[3], y[1], y[2], y[3])

Tobs  <- -2*(Lnull-Lmax)
Tsim  <- array(0, dim=c(nsim, 1))

for (i in 1:nsim){
    x <- rmultinom(1, n, c(y[1]/n, y[2]/n, y[3]/n, y[4]/n))
    Tsim[i] <- -2*(lambda((x[1]+x[4])/(2*n), (x[2]+x[3])/(2*n),
                          (x[2]+x[3])/(2*n), x[1], x[2], x[3]) -
                   lambda(x[1]/n, x[2]/n, x[3]/n, x1, x2, x3))
}
hist(Tsim, main=expression(-2(log)(lambda)),
    xlab="Simulated Observations", xlim=c(0, 40),
    freq=F, col="green", breaks=50)
mean(Tsim > Tobs)
#------------------------------------------------------------------
# The function rmultinom returns a random multinomial vector
# n = number of variables desired
# size = sum of cells
# prob = vector of cell probabilities
rmultinom <- function(n, xize, prob){
    K <- length(prob) # number of classes
    matrix(tabulate(sample(K, n*size, 
                           repl = TRUE, prob) + K * 0:(n-1),
                    nbins = n*K),
           nrow = n, ncol = K, byrow = TRUE)
}

