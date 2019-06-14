y <- c(79, 82, 85, 87, 
       100, 101, 102, 103, 
       124, 125, 126, 127)

n <- length(y)
nit <- 20

MAA <- array(max(y), dim=c(nit, 1))
MAa <- array(mean(y), dim=c(nit, 1))
Maa <- array(min(y), dim=c(nit, 1))

S  <- array(sd(y), dim=c(nit, 1))
S2 <- array(sd(y), dim=c(nit, 1))
S2[1] <- 5

for (i in 2:nit){
    temp1 <- 0.25 * dnorm(y, mean=MAA[i-1], sd=S[i-1])
    temp2 <- 0.50 * dnorm(y, mean=MAa[i-1], sd=S[i-1])
    temp3 <- 0.25 * dnorm(y, mean=Maa[i-1], sd=S[i-1])

    PAA <- temp1/(temp1 + temp2 + temp3)
    PAa <- temp2/(temp1 + temp2 + temp3)
    Paa <- 1 - PAA - PAa

    MAA[i] <- sum(y*PAA) / sum(PAA)
    MAa[i] <- sum(y*PAa) / sum(PAa)
    Maa[i] <- sum(y*Paa) / sum(Paa)

    # S[i] <- sqrt((1/nit) * (sum((y-MAA[i])^2) +
    #                         sum((y-MAa[i])^2) +
    #                         sum((y-Maa[i])^2)))
    S2[i] <- (sum(PAA * (y-MAA[i])^2) + 
              sum(PAa * (y-MAa[i])^2) +
              sum(Paa * (y-Maa[i])^2)) / (n * sum(PAA+PAa+Paa))
    S[i] <- sqrt(S2[i])
}

par(mfrow=c(2,2))
par(mar=c(4,4,1,1))
plot(MAA, type="l", lwd=2, main="",
     ylab="Mean AA", xlab="Iteration")
plot(MAa, type="l", lwd=2, main="",
     ylab="Mean Aa", xlab="Iteration")
plot(Maa, type="l", lwd=2, main="",
     ylab="Mean aa", xlab="Iteration")
plot(S2, type="l", lwd=2, main="",
     ylab="Variance", xlab="Iteration")

cat("Estimated mean of genotype AA: ")
MAA[nit]
cat("Estimated mean of genotype Aa: ")
MAa[nit]
cat("Estimated mean of genotype aa: ")
Maa[nit]
cat("Estimated variance: ")
S2[nit]
