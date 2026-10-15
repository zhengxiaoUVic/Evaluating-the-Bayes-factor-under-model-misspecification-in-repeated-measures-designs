library(fGarch) # v 4033.92
png(filename="~/Desktop/Figure B.1.png", width=3800, height=1500, res=300)
par(mfrow=c(1,2))

curve(dsstd(x, mean=0, sd=1, nu=4, xi=20), -3, 3,
      lty=2, lwd=3, col="red", xlab=expression(italic(x)), ylab="Density", main="Probability Density Functions")
curve(dsstd(x, mean=0, sd=1, nu=4, xi=1/20),
      lty=3, lwd=3, col="blue", add=T)
curve(dnorm, lwd=3, add=T)

curve(psstd(x, mean=0, sd=1, nu=4, xi=20), -3, 3,
      lty=2, lwd=3, col="red", xlab=expression(italic(x)), ylab="Probability", main="Cumulative Distribution Functions")
curve(psstd(x, mean=0, sd=1, nu=4, xi=1/20),
      lty=3, lwd=3, col="blue", add=T)
curve(pnorm, lwd=3, add=T)
legend("bottomright",c("Right-Skewed","Left-Skewed","Normal"), 
       col=c("red","blue","black"), lty=c(2,3,1), lwd=3, bty="n")
dev.off()

