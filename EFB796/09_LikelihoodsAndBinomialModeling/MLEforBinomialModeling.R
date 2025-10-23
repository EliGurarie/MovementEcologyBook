## ---------------------------------------------------------------------------------------------
Solea <- read.csv("data/Solea.csv")


## ---------------------------------------------------------------------------------------------
head(Solea)


## ----fig.height=4, fig.width = 6, echo=-1-----------------------------------------------------
par(bty="l", cex.lab=1.25, mar=c(5,5,2,2))
Y <- Solea$Solea_solea
X <- Solea$salinity
plot(X, jitter(Y, factor=0.1), col=rgb(0,0,0,.5), pch=19, ylim=c(-.2,1.2), xlab="Salinity", ylab="Presence", yaxt="n")
boxplot(X~Y, horizontal=TRUE, notch=TRUE, add=TRUE, at=c(-.1,1.1), width=c(.1,.1), col="grey", boxwex = 0.1, las=1)


## ---------------------------------------------------------------------------------------------
plotBinary <- function (X, Y, ...) 
 {
     plot(X, jitter(Y, factor = 0.1), col = rgb(0, 0, 0, 0.5), 
         pch = 19, ylim = c(-0.2, 1.2), ...)
     boxplot(X ~ Y, horizontal = TRUE, notch = TRUE, add = TRUE, 
         at = c(-0.1, 1.1), width = c(0.1, 0.1), col = "grey", 
         boxwex = 0.1, yaxt = "n")
 }


## ----fig.height=10, echo = -1-----------------------------------------------------------------
par(mfrow=c(3,1), bty="l", mar=c(4,4,1,1), cex.lab=1.5)
with(Solea, plotBinary(salinity, Solea_solea, xlab="salinity", ylab = "Presence", yaxt = "n"))
with(Solea, plotBinary(temp, Solea_solea, xlab="temperature",   ylab = "Presence", yaxt = "n"))
with(Solea, plotBinary(depth, Solea_solea, xlab="depth", ylab = "Presence", yaxt = "n"))


## ----echo=-1, fig.height=4--------------------------------------------------------------------
par(cex.lab=1.25, bty="l")
curve(exp(x)/(1+exp(x)), xlim=c(-6,6), lwd=2)
curve(exp(1 - x/2)/(1+exp(1-x/2)), add=TRUE, col=2, lwd=2)
curve(exp(-5 + 2*x)/(1+exp(-5+2*x)), add=TRUE, col=3, lwd=2)
legend("left", col=1:3, legend=c("y=x", "y=1-x/2", "y=-5+2x"), title="Predictor", lty=1)


## ---------------------------------------------------------------------------------------------
X <- Solea$salinity
Y <- Solea$Solea_solea


## ---------------------------------------------------------------------------------------------
getP.hat <- function(beta0, beta1, X){
	lP <- (beta0 + beta1 * X)
	exp(lP)/(1 + exp(lP))
}


## ----echo = -1, fig.height = 4----------------------------------------------------------------
par(bty="l")
plotBinary(X,Y)
lines(1:40, getP.hat(beta0 = 20, beta1 = -1, 1:40))
lines(1:40, getP.hat(beta0 = 10, beta1 = -1/2, 1:40), col=2)
lines(1:40, getP.hat(beta0 = 5, beta1 = -1/4, 1:40), col=3)


## ---------------------------------------------------------------------------------------------
logLike.Binomial <- function(coefs, X, Y){
	beta0 <- coefs['beta0']
	beta1 <- coefs['beta1']
	p.hat <- getP.hat(beta0, beta1, X)
	-sum(dbinom(Y, size = 1, prob = p.hat, log = TRUE))
}


## ---------------------------------------------------------------------------------------------
coef0 <- c(beta0 = 0, beta1 = -1)


## ---------------------------------------------------------------------------------------------
(glm.fit <- optim(coef0, fn = logLike.Binomial, X = X, Y = Y, hessian= TRUE))


## ---------------------------------------------------------------------------------------------
(coef.hat <- glm.fit$par)


## ----echo = -1, fig.height = 4----------------------------------------------------------------
par(bty="l")
plotBinary(X,Y)
lines(1:40, getP.hat(beta0 = glm.fit$par['beta0'], beta1 = glm.fit$par['beta1'], 1:40))


## ---------------------------------------------------------------------------------------------
sqrt(diag(solve(glm.fit$hessian)))


## ----cache=TRUE-------------------------------------------------------------------------------
getCoefs <- function(X,Y){
	optim(coef.hat, fn = logLike.Binomial, X = X, Y = Y, hessian= FALSE)$par
}


## ----PerformBootstrap, cache=TRUE-------------------------------------------------------------
nreps <- 1000
Coefs.bs <- matrix(nrow = nreps, ncol = 2)
for(i in 1:nreps){
	sample <- sample(1:length(X), replace = TRUE)
	Coefs.bs[i,] <- getCoefs(X[sample], Y[sample])
}
colnames(Coefs.bs) <- c("beta0", "beta1")


## ----BootStrapFigure, echo = -1, fig.height = 4, fig.width = 8--------------------------------
par(mfrow = c(1,2))
hist(Coefs.bs[,'beta0'], col="grey", bor = "darkgrey", breaks = 50, 
	main = expression("Bootstrap of "~beta[0]))
abline(v = quantile(Coefs.bs[,1], c(0.025, 0.5, 0.975)), col = 2, lwd = 2, lty = c(3,1,3))

hist(Coefs.bs[,'beta1'], col="grey", bor = "darkgrey", breaks = 50, 
	main = expression("Bootstrap of "~beta[1]))
abline(v = quantile(Coefs.bs[,2], c(0.025, 0.5, 0.975)), col = 2, lwd = 2, lty = c(3,1,3))


## ---------------------------------------------------------------------------------------------
apply(Coefs.bs, 2, quantile, p = c(0.025, 0.975))


## ---------------------------------------------------------------------------------------------
summary(glm(Y~X, family = "binomial"))


## ---------------------------------------------------------------------------------------------
logLik(glm(Y~X, family = "binomial"))


## ---------------------------------------------------------------------------------------------
glm.fit$value

