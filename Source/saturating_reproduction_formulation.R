z <- exp(seq(-5,10,length.out=5))
n <- exp(seq(-5,10,length.out=1000))
d <- expand.grid(n=n,z=z)
a <- 1
b <- 2
c <-100
y <- matrix(nc=5,nr=1000)
y[] <- a*d$z/((b+d$n)*(c+d$z))

matplot(log(n),log(y),type="l")
matplot(log(n),log(y),type="l")

y2 <- matrix(nc=5,nr=1000)
y2[] <- exp(log(a)+log(d$z)-log(b+d$n)-log(c+d$z))
matplot(log(n),log(y2),type="l",add=T,col="purple")

persp(log(n),log(z),log(y2),phi=15,theta=45)

      