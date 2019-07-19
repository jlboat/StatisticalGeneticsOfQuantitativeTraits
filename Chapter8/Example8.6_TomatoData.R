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
# This gets MLEs for the Tomato data
y1 <- c(79, 82,
       100, 102,
       124)
y2 <- c(85, 87, 
       101, 103, 
       125, 126, 127)
n1 <- length(y1)
n2 <- length(y2)
n <- n1 + n2

# Initialize Estimates
r <- 0.25
m1 <- mean(y1)
m2 <- mean(y2)
s <- sqrt(var(c(y1,y2)))
m1plot <- m1
m2plot <- m2
splot <- s
rplot <- r

# Start Iteration
nit <- 100
for (i in 1:nit){
    w1 <- P1(y1, m1, m2, s, r)
    w2 <- P2(y2, m1, m2, s, r)
    m1 <- (sum(w1*y1) + sum(w2*y2))/(sum(w1)+sum(w2))
    m2 <- (sum((1-w1)*y1) + sum((1-w2)*y2))/(sum((1-w1))+sum((1-w2)))
    s <- sqrt((sum(w1*(y1-m1)^2+(1-w1)*(y1-m2)^2) +
               sum(w2*(y2-m1)^2+(1-w2)*(y2-m2)^2))/n)
    r <- (sum(1-w1) + sum(w2))/n
    r <- min(r, 0.5)
    m1plot <- c(m1plot, m1)
    m2plot <- c(m2plot, m2)
    splot <- c(splot, s)
    rplot <- c(rplot, r)
}

# Plot
par(mfrow=c(2,2))
par(mar=c(4,4,1,1))
plot(m1plot, type="l", lwd=2, main="",
     ylab="mean 1", xlab="Iteration")
plot(m2plot, type="l", lwd=2, main="",
     ylab="mean 2", xlab="Iteration")
plot(splot, type="l", lwd=2, main="",
     ylab="std", xlab="Iteration")
plot(rplot, type="l", lwd=2, main="",
     ylab="r", xlab="Iteration")


