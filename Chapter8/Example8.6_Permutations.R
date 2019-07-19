# Define Functions for P1 and P2
P1 <- function(y, t1, t2, w, r){
    temp1 <- (1-r)*dnorm(y, mean=t1, sd=w)
    temp2 <- r*dnorm(y, mean=t2, sd=w)
    return(temp1/(temp1+temp2))
}
P2 <- function(y, t1, t2, w, r){
    temp1 <- r*dnorm(y, mean=t1, sd=w)
    temp2 <- (1-r)*dnorm(y, mean=t2, sd=w)
    return(temp1/(temp1+temp2))
}

# Define Functions for log likelihood
logL <- function(y1, y2, t1, t2, w, r){
    S1 <- sum(log((1-r)*dnorm(y1, t1, s) + r*dnorm(y1,t2,s)))
    S2 <- sum(log(r*dnorm(y2, t1, s) + (1-r)*dnorm(y2,t2,s)))
    return(S1 + S2)
}

# This gets permutation distribution for the Tomato data
y1 <- c(79, 82,
        100, 102,
        124)
y2 <- c(85, 87,
        101, 103,
        125, 126, 127)
n1 <- length(y1)
n2 <- length(y2)
n  <- n1 + n2
y  <- c(y1, y2)

# Initialize Estimates
r  <- 0.25
m1 <- mean(y1)
m2 <- mean(y2)
s  <- sd(y)
Lplot <- logL(y1, y2, m1, m2, s, r)

# Start Iteration
nperm <- 5000 # number of permutation statistics
nit   <- 50   # number of iterations to find MLEs
for (j in 1:nperm){
    yp <- sample(y)
    y1 <- yp[1:n1]
    y2 <- yp[(n1+1):n]
    for (i in 1:nit){
        w1 <- P1(y1, m1, m2, s, r)
        w2 <- P2(y2, m1, m2, s, r)
        m1 <- (sum(w1*y1) + sum(w2+y2))/(sum(w1) + sum(w2))
        m2 <- (sum((1-w1)*y1) + sum((1-w2)+y2))/(sum((1-w1)) + sum((1-w2)))
        s  <- sqrt((sum(w1*(y1-m1)^2 + (1-w1)*(y1-m2)^2) +
                    sum(w2*(y2-m1)^2 + (1-w2)*(y2-m2)^2))/n)
        r  <- (sum(1-w1) + sum(w2))/n
        r  <- min(r, 0.5)
    }
    L <- logL(y1, y2, m1, m2, s, r)
    Lplot <- c(Lplot, L)
}

# Calculate Statistics
m   <- mean(y)
s   <- sd(y)
LHO <- logL(y1, y2, m, m, s, 0.5)
Lplot <- -2*(LHO-Lplot)
sort(Lplot)[0.95*nperm]
sort(Lplot)[0.05*nperm]

# Plot
hist(Lplot, main="Permutation Distribution",
     freq=F, xlab="-2 log lambda")
