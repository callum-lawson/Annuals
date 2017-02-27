### SETUP

n <- 10^6

### LOGNORMAL

#mseq <- rlnorm(n,meanlog=0.6,sdlog=0.5)
#Sseq <- exp(-mseq)

### GAMMA

### tried examples:
# 0.17,2

Sbar <- 0.385	
m <- -log(Sbar)
sd <- 1.3
v <- sd^2

shape <- 0.6 # m^2/v
scale <- 1.8 # v/m

mseq <- rgamma(n,shape=shape,scale=scale)
Sseq <- exp(-mseq)

### GRAPHS

par(mfrow=c(3,1))
hist(mseq,prob=T,breaks=10^4)
hist(Sseq,prob=T,breaks=10^4)
hist(qlogis(Sseq),prob=T,breaks=10^4,xaxt="n")

y_labels = c(0.01,0.02,0.05,0.1,0.2,0.5,0.8,0.9,0.98,0.99)
y_at = qlogis(y_labels)
axis(1, labels=y_labels, at=y_at)

### SUMMARY STATS

mean(Sseq)
quantile(Sseq,probs=c(0.01,0.05,0.25,0.318,0.5,0.682,0.75,0.95,0.99))

mu <- mean(qlogis(Sseq))
sig <- sd(qlogis(Sseq))
plogis(mu)
plogis(mu+sig)
plogis(mu-sig)

table(Sseq<0.02 | Sseq>0.8)[2]/10^6



