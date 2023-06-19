
source("thmm1.r")

library(mhsmm)

tmat1 <- matrix( 0.09,5,5 ) + 0.55*diag(5);
tmat2 <- cbind(0,rbind(0.4*diag(4),0))
tmat2 <- tmat2 + t(tmat2)
tmat2[1,5] <- 0.4;
tmat2[5,1] <- 0.4;
tmat2 <- tmat2 + matrix(0.04,5,5)
#tmat2 <- tmat2 + 0.4*diag(5)

crossEnt <- function(
  est, tru
){
  tab = table(est,tru);
  len = nrow(tab);
  ent = 0;
  csum= colSums(tab);
  ftab= t(t(tab)/csum);
  ltab= log2(ftab);
  mx  = apply( ltab, 2, max )
  return(
    -sum( mx )
  )
}


#
# BMWD Simulation
#

###########
# Low Sep #
###########

set.seed(137)
dat.bmwd1med <- sim.thmm(
  n=200,
  prInit=c(1,0,0,0,0),
  prTran=tmat1,
  drift=c(-4,-2,0,2,4),
  type="BMWD",
  len=100
)

out.bmwd1med <- thmm.BaumWelch(
  dat.bmwd1med$data, stateCnt=5,
  verbose=1,tol=1e-8,maxIter=1000
)

table( out.bmwd1med$bestSeq, dat.bmwd1med$ss )
adjustedRandIndex( out.bmwd1med$bestSeq, dat.bmwd1med$ss )
crossEnt( out.bmwd1med$bestSeq, dat.bmwd1med$ss )

## PCA HMM

set.seed(137)
pcs = 2
stcnt = 5;
tmat.init <- matrix(
  -log( runif(stcnt^2) ),stcnt,stcnt
)
tmat.init <- (tmat.init/rowSums(tmat.init))

nr = nrow( dat.bmwd1med$data )
prc.bmwd1med = prcomp( x = dat.bmwd1med$data, retx=T )
prc.data <- prc.bmwd1med$x[,1:pcs]
mu.list = list()
sig.list= list()
for( i in 1:stcnt ){
  mu.list[[i]]  = prc.data[sample(nr,1),]
  sig.list[[i]] = cov(prc.data[sample(nr,nr/2),])
}

hmm.model <- hmmspec( 
  init = rep(1/stcnt,stcnt),
  trans= tmat.init,
  parms.emission=list(
    mu=mu.list, sigma=sig.list 
  ),
  dens.emission=dmvnorm.hsmm
)

train <- simulate(
  hmm.model, nsim=200, seed=137, 
  rand.emis=rmvnorm.hsmm
)

train$s <- dat.bmwd1med$ss;
train$x <- prc.bmwd1med$x[,1:pcs];

fit.hmm1 = hmmfit(
  train,hmm.model,mstep=mstep.mvnorm
)

est.ss <- predict( fit.hmm1, train )

table( est.ss$s, dat.bmwd1med$ss )
adjustedRandIndex( est.ss$s, dat.bmwd1med$ss )
crossEnt( est.ss$s, dat.bmwd1med$ss )

##
set.seed(137)
dat.bmwd2med <- sim.thmm(
  n=200,
  prInit=c(1,0,0,0,0),
  prTran=tmat2,
  drift=c(-4,-2,0,2,4),
  type="BMWD",
  len=100
)

out.bmwd2med <- thmm.BaumWelch(
  dat.bmwd2med$data, stateCnt=5,
  verbose=1,tol=1e-12,maxIter=1000
)

table( out.bmwd2med$bestSeq, dat.bmwd2med$ss )
adjustedRandIndex( out.bmwd2med$bestSeq, dat.bmwd2med$ss )
crossEnt( out.bmwd2med$bestSeq, dat.bmwd2med$ss )

## PCA HMM

set.seed(137)
pcs = 2
stcnt = 5;
tmat.init <- matrix(
  -log( runif(stcnt^2) ),stcnt,stcnt
)
tmat.init <- (tmat.init/rowSums(tmat.init))

nr = nrow( dat.bmwd2med$data )
prc.bmwd2med = prcomp( x = dat.bmwd2med$data, retx=T )
prc.data <- prc.bmwd2med$x[,1:pcs]
mu.list = list()
sig.list= list()
for( i in 1:stcnt ){
  mu.list[[i]]  = prc.data[sample(nr,1),]
  sig.list[[i]] = cov(prc.data[sample(nr,nr/2),])
}

hmm.model <- hmmspec( 
  init = rep(1/stcnt,stcnt),
  trans= tmat.init,
  parms.emission=list(
    mu=mu.list, sigma=sig.list 
  ),
  dens.emission=dmvnorm.hsmm
)

train <- simulate(
  hmm.model, nsim=200, seed=137, 
  rand.emis=rmvnorm.hsmm
)

train$s <- dat.bmwd2med$ss;
train$x <- prc.bmwd2med$x[,1:pcs];

fit.hmm2 = hmmfit(
  train,hmm.model,mstep=mstep.mvnorm
)

est.ss <- predict( fit.hmm2, train )

table( est.ss$s, dat.bmwd2med$ss )
adjustedRandIndex( est.ss$s, dat.bmwd2med$ss )
crossEnt( est.ss$s, dat.bmwd2med$ss )

############
# high Sep #
############

set.seed(137)
dat.bmwd1high <- sim.thmm(
  n=200,
  prInit=c(1,0,0,0,0),
  prTran=tmat1,
  drift=c(-8,-4,0,4,8),
  type="BMWD",
  len=100
)

out.bmwd1high <- thmm.BaumWelch(
  dat.bmwd1high$data, stateCnt=5,
  verbose=1,tol=1e-8,maxIter=1000
)

table( out.bmwd1high$bestSeq, dat.bmwd1high$ss )
adjustedRandIndex( out.bmwd1high$bestSeq, dat.bmwd1high$ss )
crossEnt( out.bmwd1high$bestSeq, dat.bmwd1high$ss )


## PCA HMM

set.seed(137)
pcs = 2
stcnt = 5;
tmat.init <- matrix(
  -log( runif(stcnt^2) ),stcnt,stcnt
)
tmat.init <- (tmat.init/rowSums(tmat.init))

nr = nrow( dat.bmwd1high$data )
prc.bmwd1high = prcomp( x = dat.bmwd1high$data, retx=T )
prc.data <- prc.bmwd1high$x[,1:pcs]
mu.list = list()
sig.list= list()
for( i in 1:stcnt ){
  mu.list[[i]]  = prc.data[sample(nr,1),]
  sig.list[[i]] = cov(prc.data[sample(nr,nr/2),])
}

hmm.model <- hmmspec( 
  init = rep(1/stcnt,stcnt),
  trans= tmat.init,
  parms.emission=list(
    mu=mu.list, sigma=sig.list 
  ),
  dens.emission=dmvnorm.hsmm
)

train <- simulate(
  hmm.model, nsim=200, seed=137, 
  rand.emis=rmvnorm.hsmm
)

train$s <- dat.bmwd1high$ss;
train$x <- prc.bmwd1high$x[,1:pcs];

fit.hmm3 = hmmfit(
  train,hmm.model,mstep=mstep.mvnorm
)

est.ss <- predict( fit.hmm3, train )

table( est.ss$s, dat.bmwd1high$ss )
adjustedRandIndex( est.ss$s, dat.bmwd1high$ss )
crossEnt( est.ss$s, dat.bmwd1high$ss )


###

set.seed(137)
dat.bmwd2high <- sim.thmm(
  n=200,
  prInit=c(1,0,0,0,0),
  prTran=tmat2,
  drift=c(-8,-4,0,4,8),
  type="BMWD",
  len=100
)

out.bmwd2high <- thmm.BaumWelch(
  dat.bmwd2high$data, stateCnt=5,
  verbose=1,tol=1e-8,maxIter=1000
)

table( out.bmwd2high$bestSeq, dat.bmwd2high$ss )
adjustedRandIndex( out.bmwd2high$bestSeq, dat.bmwd2high$ss )
crossEnt( out.bmwd2high$bestSeq, dat.bmwd2high$ss )

## PCA HMM

set.seed(137)
pcs = 2
stcnt = 5;
tmat.init <- matrix(
  -log( runif(stcnt^2) ),stcnt,stcnt
)
tmat.init <- (tmat.init/rowSums(tmat.init))

nr = nrow( dat.bmwd2high$data )
prc.bmwd2high = prcomp( x = dat.bmwd2high$data, retx=T )
prc.data <- prc.bmwd2high$x[,1:pcs]
mu.list = list()
sig.list= list()
for( i in 1:stcnt ){
  mu.list[[i]]  = prc.data[sample(nr,1),]
  sig.list[[i]] = cov(prc.data[sample(nr,nr/2),])
}

hmm.model <- hmmspec( 
  init = rep(1/stcnt,stcnt),
  trans= tmat.init,
  parms.emission=list(
    mu=mu.list, sigma=sig.list 
  ),
  dens.emission=dmvnorm.hsmm
)

train <- simulate(
  hmm.model, nsim=200, seed=137, 
  rand.emis=rmvnorm.hsmm
)

train$s <- dat.bmwd2high$ss;
train$x <- prc.bmwd2high$x[,1:pcs];

fit.hmm4 = hmmfit(
  train,hmm.model,mstep=mstep.mvnorm
)

est.ss <- predict( fit.hmm4, train )

table( est.ss$s, dat.bmwd2high$ss )
adjustedRandIndex( est.ss$s, dat.bmwd2high$ss )
crossEnt( est.ss$s, dat.bmwd2high$ss )

##############
# OU Process #
##############

set.seed(137)
dat.ouproc <- sim.thmm(
  n=200,
  prInit=c(1,0,0,0,0),
  prTran=tmat1,
  drift=c(-2,0,4,2,1),
  theta=c(4,4,8,2,20),
  type="OU",
  len=100
)

out.ouproc <- thmm.BaumWelch(
  dat.ouproc$data, stateCnt=5,
  verbose=1,tol=1e-8,maxIter=1000,
  type='OU'
)

table( out.ouproc$bestSeq, dat.ouproc$ss )
adjustedRandIndex( out.ouproc$bestSeq, dat.ouproc$ss )
crossEnt( out.ouproc$bestSeq, dat.ouproc$ss )

## PCA HMM

set.seed(137)
pcs = 2
stcnt = 5;
tmat.init <- matrix(
  -log( runif(stcnt^2) ),stcnt,stcnt
)
tmat.init <- (tmat.init/rowSums(tmat.init))

nr = nrow( dat.ouproc$data )
prc.ouproc = prcomp( x = dat.ouproc$data, retx=T )
prc.data <- prc.ouproc$x[,1:pcs]
mu.list = list()
sig.list= list()
for( i in 1:stcnt ){
  mu.list[[i]]  = prc.data[sample(nr,1),]
  sig.list[[i]] = cov(prc.data[sample(nr,nr/2),])
}

hmm.model <- hmmspec( 
  init = rep(1/stcnt,stcnt),
  trans= tmat.init,
  parms.emission=list(
    mu=mu.list, sigma=sig.list 
  ),
  dens.emission=dmvnorm.hsmm
)

train <- simulate(
  hmm.model, nsim=200, seed=137, 
  rand.emis=rmvnorm.hsmm
)

train$s <- dat.ouproc$ss;
train$x <- prc.ouproc$x[,1:pcs];

fit.hmm4 = hmmfit(
  train,hmm.model,mstep=mstep.mvnorm
)

est.ss <- predict( fit.hmm4, train )

table( est.ss$s, dat.ouproc$ss )
adjustedRandIndex( est.ss$s, dat.ouproc$ss )
crossEnt( est.ss$s, dat.ouproc$ss )

#############

set.seed(137)
dat.ouproc <- sim.thmm(
  n=200,
  prInit=c(1,0,0,0,0),
  prTran=tmat2,
  drift=c(-2,0,4,2,1),
  theta=c(4,4,8,2,20),
  type="OU",
  len=100
)

out.ouproc <- thmm.BaumWelch(
  dat.ouproc$data, stateCnt=5,
  verbose=1,tol=1e-8,maxIter=1000,
  type='OU'
)

table( out.ouproc$bestSeq, dat.ouproc$ss )
adjustedRandIndex( out.ouproc$bestSeq, dat.ouproc$ss )
crossEnt( out.ouproc$bestSeq, dat.ouproc$ss )


## Four State ##

set.seed(137)
dat.ouproc <- sim.thmm(
  n=200,
  prInit=c(1,0,0,0,0),
  prTran=tmat1,
  drift=c(-2,0,4,2,1),
  theta=c(4,4,8,2,20),
  type="OU",
  len=100
)

out.ouproc <- thmm.BaumWelch(
  dat.ouproc$data, stateCnt=4,
  verbose=1,tol=1e-8,maxIter=1000,
  type='OU'
)

labs <- dat.ouproc$ss
labs[which(labs==5)] <- 4;

table( out.ouproc$bestSeq, labs  )
adjustedRandIndex( out.ouproc$bestSeq, labs )
crossEnt( out.ouproc$bestSeq, labs )

# PCA HMM

set.seed(137)
pcs = 2
stcnt = 4;
tmat.init <- matrix(
  -log( runif(stcnt^2) ),stcnt,stcnt
)
tmat.init <- (tmat.init/rowSums(tmat.init))

nr = nrow( dat.ouproc$data )
prc.ouproc = prcomp( x = dat.ouproc$data, retx=T )
prc.data <- prc.ouproc$x[,1:pcs]
mu.list = list()
sig.list= list()
for( i in 1:stcnt ){
  mu.list[[i]]  = prc.data[sample(nr,1),]
  sig.list[[i]] = cov(prc.data[sample(nr,nr/2),])
}

hmm.model <- hmmspec( 
  init = rep(1/stcnt,stcnt),
  trans= tmat.init,
  parms.emission=list(
    mu=mu.list, sigma=sig.list 
  ),
  dens.emission=dmvnorm.hsmm
)

train <- simulate(
  hmm.model, nsim=200, seed=137, 
  rand.emis=rmvnorm.hsmm
)

train$s <- dat.ouproc$ss;
train$x <- prc.ouproc$x[,1:pcs];

fit.hmm4 = hmmfit(
  train,hmm.model,mstep=mstep.mvnorm
)

est.ss <- predict( fit.hmm4, train )

table( est.ss$s, labs )
adjustedRandIndex( est.ss$s, labs )
crossEnt( est.ss$s, labs )

##############
# Frac Brown #
##############

## Hurst 0.8

set.seed(137)
dat.fracbm <- sim.thmm(
  n=200,
  prInit=c(1,0,0,0,0),
  prTran=tmat1,
  drift=c(-10,-6,-2,0,2),
  type="BMWD",
  hurst=0.8,
  len=100
)

matplot( 
  t(dat.fracbm$data), type='l', 
  lty=dat.fracbm$ss, col=dat.fracbm$ss+1, 
  las=1, xaxt='n', xlab='', ylab="",
  main="Hurst Parameter = 0.8"
); 
axis(side=1,at=(0:5)*20,labels=(0:5)/5)

out.fracbm <- thmm.BaumWelch(
  dat.fracbm$data, stateCnt=5,
  verbose=1,tol=1e-8,maxIter=1000,
  type='BMWD',hurst=0.8
)

table( out.fracbm$bestSeq, dat.fracbm$ss )
adjustedRandIndex( out.fracbm$bestSeq, dat.fracbm$ss )
crossEnt( out.fracbm$bestSeq, dat.fracbm$ss )


## PCA HMM

set.seed(137)
pcs = 2
stcnt = 5;
tmat.init <- matrix(
  -log( runif(stcnt^2) ),stcnt,stcnt
)
tmat.init <- (tmat.init/rowSums(tmat.init))

nr = nrow( dat.fracbm$data )
prc.fracbm = prcomp( x = dat.fracbm$data, retx=T )
prc.data <- prc.fracbm$x[,1:pcs]
mu.list = list()
sig.list= list()
for( i in 1:stcnt ){
  mu.list[[i]]  = prc.data[sample(nr,1),]
  sig.list[[i]] = cov(prc.data[sample(nr,nr/2),])
}

hmm.model <- hmmspec( 
  init = rep(1/stcnt,stcnt),
  trans= tmat.init,
  parms.emission=list(
    mu=mu.list, sigma=sig.list 
  ),
  dens.emission=dmvnorm.hsmm
)

train <- simulate(
  hmm.model, nsim=200, seed=137, 
  rand.emis=rmvnorm.hsmm
)

train$s <- dat.fracbm$ss;
train$x <- prc.fracbm$x[,1:pcs];

fit.hmm4 = hmmfit(
  train,hmm.model,mstep=mstep.mvnorm
)

est.ss <- predict( fit.hmm4, train )

table( est.ss$s, dat.fracbm$ss )
adjustedRandIndex( est.ss$s, dat.fracbm$ss )
crossEnt( est.ss$s, dat.fracbm$ss )



## Hurst 0.25

set.seed(137)
dat.fracbm <- sim.thmm(
  n=200,
  prInit=c(1,0,0,0,0),
  prTran=tmat1,
  drift=c(-16,-8,-4,0,8),
  type="BMWD",
  hurst=0.25,
  len=100
)

matplot( 
  t(dat.fracbm$data), type='l', 
  lty=dat.fracbm$ss, col=dat.fracbm$ss+1, 
  las=1, xaxt='n', xlab='', ylab="",
  main="Hurst Parameter = 0.25"
); 
axis(side=1,at=(0:5)*20,labels=(0:5)/5)

out.fracbm <- thmm.BaumWelch(
  dat.fracbm$data, stateCnt=5,
  verbose=1,tol=1e-8,maxIter=1000,
  type='BMWD',hurst=0.25
)
  

table( out.fracbm$bestSeq, dat.fracbm$ss )
adjustedRandIndex( out.fracbm$bestSeq, dat.fracbm$ss )
crossEnt( out.fracbm$bestSeq, dat.fracbm$ss )


## PCA HMM

set.seed(137)
pcs = 2
stcnt = 5;
tmat.init <- matrix(
  -log( runif(stcnt^2) ),stcnt,stcnt
)
tmat.init <- (tmat.init/rowSums(tmat.init))

nr = nrow( dat.fracbm$data )
prc.fracbm = prcomp( x = dat.fracbm$data, retx=T )
prc.data <- prc.fracbm$x[,1:pcs]
mu.list = list()
sig.list= list()
for( i in 1:stcnt ){
  mu.list[[i]]  = prc.data[sample(nr,1),]
  sig.list[[i]] = cov(prc.data[sample(nr,nr/2),])
}

hmm.model <- hmmspec( 
  init = rep(1/stcnt,stcnt),
  trans= tmat.init,
  parms.emission=list(
    mu=mu.list, sigma=sig.list 
  ),
  dens.emission=dmvnorm.hsmm
)

train <- simulate(
  hmm.model, nsim=200, seed=137, 
  rand.emis=rmvnorm.hsmm
)

train$s <- dat.fracbm$ss;
train$x <- prc.fracbm$x[,1:pcs];

fit.hmm4 = hmmfit(
  train,hmm.model,mstep=mstep.mvnorm
)

est.ss <- predict( fit.hmm4, train )

table( est.ss$s, dat.fracbm$ss )
adjustedRandIndex( est.ss$s, dat.fracbm$ss )
crossEnt( est.ss$s, dat.fracbm$ss )

############
## NonPar ##
############

library(fda.usc)

set.seed(137)
tt = seq(0,1,length.out=100);
dat.nonpar <- sim.thmm(
  n=200,
  prInit=c(1,0,0,0,0),
  prTran=tmat1,
  drift=rbind( 
    sin(2*pi*tt),sin(2*pi*(tt+0.2)),
    sin(2*pi*(tt+0.4)),sin(2*pi*(tt+0.6)),
    sin(2*pi*(tt+0.8))
  ),
  type="nonpar",
  len=100
)

matplot( 
  t(dat.nonpar$data), type='l', 
  lty=dat.nonpar$ss, col=dat.nonpar$ss+1, 
  las=1, xaxt='n', xlab='', ylab=""
); 
axis(side=1,at=(0:5)*20,labels=(0:5)/5)

# L2 #

out.nonpar.l2 <- thmm.BaumWelch(
  dat.nonpar$data, stateCnt=5,
  verbose=1,tol=1e-12,maxIter=1000,
  type='nonpar',norm="L2"
)
  
matplot( 
  tt,t(dat.nonpar$data), type='l', 
  lty=out.nonpar.l2$bestSeq, col=out.nonpar.l2$bestSeq+1, 
  las=1, xlab='', ylab=""
); 
matplot( 
  tt,t(out.nonpar.l2$stPar), type='l', 
  lty=1:5, col="black",add=T,lwd=3
); 

table( out.nonpar.l2$bestSeq, dat.nonpar$ss )
adjustedRandIndex( out.nonpar.l2$bestSeq, dat.nonpar$ss )
crossEnt( out.nonpar.l2$bestSeq, dat.nonpar$ss )

# kmeans #

out.km = kmeans.fd( dat.nonpar$data, ncl=5 )
table( out.km$cluster, dat.nonpar$ss )
adjustedRandIndex( out.km$cluster, dat.nonpar$ss )
crossEnt( out.km$cluster, dat.nonpar$ss )

# kmeans + L2 #

out.nonpar.l2a <- thmm.BaumWelch(
  dat.nonpar$data, stateCnt=5,
  verbose=1,tol=1e-12,maxIter=1000,
  type='nonpar',norm="L2",
  state_mean=out.km$centers$data
)
  
table( out.nonpar.l2a$bestSeq, dat.nonpar$ss )
adjustedRandIndex( out.nonpar.l2a$bestSeq, dat.nonpar$ss )
crossEnt( out.nonpar.l2a$bestSeq, dat.nonpar$ss )

# W21 #

out.nonpar.W21 <- thmm.BaumWelch(
  dat.nonpar$data, stateCnt=5,
  verbose=1,tol=1e-12,maxIter=10000,
  type='nonpar',norm="W21"
)
  
table( out.nonpar.W21$bestSeq, dat.nonpar$ss )
adjustedRandIndex( out.nonpar.W21$bestSeq, dat.nonpar$ss )
crossEnt( out.nonpar.W21$bestSeq, dat.nonpar$ss )

# kmeans + W21 #

out.nonpar.W21a <- thmm.BaumWelch(
  dat.nonpar$data, stateCnt=5,
  verbose=1,tol=1e-12,maxIter=10000,
  type='nonpar',norm="W21",
  state_mean=out.km$centers$data
)
  
table( out.nonpar.W21a$bestSeq, dat.nonpar$ss )
adjustedRandIndex( out.nonpar.W21a$bestSeq, dat.nonpar$ss )
crossEnt( out.nonpar.W21a$bestSeq, dat.nonpar$ss )

## PCA HMM

set.seed(137)
pcs = 4
stcnt = 5;
tmat.init <- matrix(
  -log( runif(stcnt^2) ),stcnt,stcnt
)
tmat.init <- (tmat.init/rowSums(tmat.init))

nr = nrow( dat.nonpar$data )
prc.nonpar = prcomp( x = dat.nonpar$data, retx=T )
prc.data <- prc.nonpar$x[,1:pcs]
mu.list = list()
sig.list= list()
for( i in 1:stcnt ){
  mu.list[[i]]  = prc.data[sample(nr,1),]
  sig.list[[i]] = cov(prc.data[sample(nr,nr/2),])
}

hmm.model <- hmmspec( 
  init = rep(1/stcnt,stcnt),
  trans= tmat.init,
  parms.emission=list(
    mu=mu.list, sigma=sig.list 
  ),
  dens.emission=dmvnorm.hsmm
)

train <- simulate(
  hmm.model, nsim=200, seed=137, 
  rand.emis=rmvnorm.hsmm
)

train$s <- dat.nonpar$ss;
train$x <- prc.nonpar$x[,1:pcs];

fit.hmm4 = hmmfit(
  train,hmm.model,mstep=mstep.mvnorm
)

est.ss <- predict( fit.hmm4, train )

table( est.ss$s, dat.nonpar$ss )
adjustedRandIndex( est.ss$s, dat.nonpar$ss )
crossEnt( est.ss$s, dat.nonpar$ss )

## PCA HMM kmeans

set.seed(137)
pcs = 4
stcnt = 5;
tmat.init <- matrix(
  -log( runif(stcnt^2) ),stcnt,stcnt
)
tmat.init <- (tmat.init/rowSums(tmat.init))

nr = nrow( dat.nonpar$data )
prc.nonpar = prcomp( x = dat.nonpar$data, retx=T )
prc.data <- prc.nonpar$x[,1:pcs]

km <- kmeans( prc.data, 5 )

mu.list = list()
sig.list= list()
for( i in 1:stcnt ){
  mu.list[[i]]  = km$centers[i,]
  sig.list[[i]] = cov(prc.data[sample(nr,nr/2),])
}

hmm.model <- hmmspec( 
  init = rep(1/stcnt,stcnt),
  trans= tmat.init,
  parms.emission=list(
    mu=mu.list, sigma=sig.list 
  ),
  dens.emission=dmvnorm.hsmm
)

train <- simulate(
  hmm.model, nsim=200, seed=137, 
  rand.emis=rmvnorm.hsmm
)

train$s <- dat.nonpar$ss;
train$x <- prc.nonpar$x[,1:pcs];

fit.hmm4 = hmmfit(
  train,hmm.model,mstep=mstep.mvnorm
)

est.ss <- predict( fit.hmm4, train )

table( est.ss$s, dat.nonpar$ss )
adjustedRandIndex( est.ss$s, dat.nonpar$ss )
crossEnt( est.ss$s, dat.nonpar$ss )

###  Double Take ###

set.seed(137)
tt = seq(0,1,length.out=100);
dat.nonpar <- sim.thmm(
  n=200,
  prInit=c(1,0,0,0,0),
  prTran=tmat1,
  drift=rbind( 
    sin(2*pi*tt),sin(2*pi*(tt+0.2)),
    sin(2*pi*(tt+0.4)),sin(2*pi*(tt+0.6)),
    sin(2*pi*(tt+0.8))
  ),
  type="nonpar",
  len=100
)

out.nonpar2 <- thmm.BaumWelch(
  dat.nonpar$data, stateCnt=5,
  verbose=1,tol=1e-8,maxIter=1000,
  type='nonpar',norm="L2"
)
  
matplot( 
  t(dat.nonpar$data), type='l', col=dat.nonpar$ss+1, lty=1, 
  las=1,ylab=""
)
matplot(
  t(out.nonpar2$stPar), type='l', col='black', lty=1:5, add=T,lwd=3 
)

table( out.nonpar2$bestSeq, dat.nonpar$ss )
adjustedRandIndex( out.nonpar2$bestSeq, dat.nonpar$ss )
crossEnt( out.nonpar2$bestSeq, dat.nonpar$ss )

out.nonpar2a <- thmm.BaumWelch(
  dat.nonpar$data, stateCnt=5,
  verbose=1,tol=1e-8,maxIter=1000,
  type='nonpar',norm="W21",
  prob.init=exp(out.nonpar2$prInit),
  prob.trans=exp(out.nonpar2$prTran),
  state_mean=out.nonpar2$stPar
)

table( out.nonpar2a$bestSeq, dat.nonpar$ss )
adjustedRandIndex( out.nonpar2a$bestSeq, dat.nonpar$ss )
crossEnt( out.nonpar2a$bestSeq, dat.nonpar$ss )

############
# Snowfall #
############

edm = c()
for( i in 1940:1990 ){
  fp = paste0("../../../DATA/Weather/EdmontonDaily/en_climate_daily_AB_3012208_",i,"_P1D.csv")
  tmp = read.csv(fp);
  tmp = tmp$Total.Snow..cm.
  #tmp = tmp$Mean.Temp..Â.C.
  if(length(tmp)==366)
    tmp = tmp[-60];
  edm = rbind(edm,(tmp))
}
edm <- edm[-c(53,54,66),]

edm.sm = edm
for( i in 1:nrow(edm) ){
  edm.sm[i,] = ksmooth(1:365,edm[i,],bandwidth = 30)$y
}

edm.h1 = edm[,1:182]
edm.h2 = edm[,183:365]
edm.h  = cbind( edm.h2[-51,], edm.h1[-1,] )
edm.h <- t(apply(edm.h,1,cumsum))

edm.sm = edm.h
for( i in 1:nrow(edm.h) ){
  edm.sm[i,] = ksmooth(1:365,edm.h[i,],bandwidth = 30)$y
}

plot( 
  1940:1989,edm.h[,365],type='l',xlab='Year',las=1,
  ylab="Total Snowfall (cm)",main="Univariate HMM" 
); 

matplot(
  t(edm.sm), type='l' 
)

#########
stcnt = 2  ## Set it here for all!
#########


set.seed(137)

## Rep it!
reps=20
fit.hmm.list = list();
fit.hmm.like = c();
for( rr in 1:reps ){

out.edm.l2 <- thmm.BaumWelch(
  data = edm.sm,
  stateCnt = stcnt,
  maxIter = 10000,
  tol = 1e-13,
  type = "nonpar",
  norm = "L2",
  verbose = 0
)

fit.hmm.list[[rr]] =  out.edm.l2
fit.hmm.like[rr] = out.edm.l2$bestProb

}

indx = which.max( fit.hmm.like )
out.edm.l2 <- fit.hmm.list[[indx]]

#matplot(
#  t(edm.sm), type='l',
#  col = out.edm.l2$bestSeq
#)

set.seed(137)

## Rep it!
reps=20
fit.hmm.list = list();
fit.hmm.like = c();
for( rr in 1:reps ){

out.edm.w21 <- thmm.BaumWelch(
  data = edm.sm,
  stateCnt = stcnt,
  maxIter = 10000,
  tol = 1e-13,
  type = "nonpar",
  norm = "W21",
  verbose = 0
)

fit.hmm.list[[rr]] =  out.edm.w21
fit.hmm.like[rr] = out.edm.w21$bestProb

}

indx = which.max( fit.hmm.like )
out.edm.w21 <- fit.hmm.list[[indx]]

set.seed(137)

## Rep it!
reps=20
fit.hmm.list = list();
fit.hmm.like = c();
for( rr in 1:reps ){

out.edm.w22 <- thmm.BaumWelch(
  data = edm.sm,
  stateCnt = stcnt,
  maxIter = 10000,
  tol = 1e-13,
  type = "nonpar",
  norm = "W22",
  verbose = 0
)

fit.hmm.list[[rr]] =  out.edm.w22
fit.hmm.like[rr] = out.edm.w22$bestProb

}

indx = which.max( fit.hmm.like )
out.edm.w22 <- fit.hmm.list[[indx]]


adjustedRandIndex( out.edm.l2$bestSeq, out.edm.w21$bestSeq )
adjustedRandIndex( out.edm.l2$bestSeq, out.edm.w22$bestSeq )
adjustedRandIndex( out.edm.w22$bestSeq, out.edm.w21$bestSeq )

###  PCA  ###

set.seed(137)

## Rep it!
reps=20
fit.hmm.list = list();
fit.hmm.like = c();
for( rr in 1:reps ){

pcs = 4
#stcnt = 2;
tmat.init <- matrix(
  -log( runif(stcnt^2) ),stcnt,stcnt
)
tmat.init <- (tmat.init/rowSums(tmat.init))

nr = nrow( edm.sm )
prc.nonpar = prcomp( x = edm.sm, retx=T )
prc.data <- prc.nonpar$x[,1:pcs]

km <- kmeans( prc.data, 5 )

mu.list = list()
sig.list= list()
for( i in 1:stcnt ){
  mu.list[[i]]  = km$centers[i,]
  sig.list[[i]] = cov(prc.data[sample(nr,nr/2),])
}

hmm.model <- hmmspec( 
  init = rep(1/stcnt,stcnt),
  trans= tmat.init,
  parms.emission=list(
    mu=mu.list, sigma=sig.list 
  ),
  dens.emission=dmvnorm.hsmm
)

train <- simulate(
  hmm.model, nsim=nr, seed=137, 
  rand.emis=rmvnorm.hsmm
)

train$s <- c(rep(1,nr/2),rep(2,nr/2));
train$x <- prc.nonpar$x[,1:pcs];

fit.hmm4 = hmmfit(
  train,hmm.model,mstep=mstep.mvnorm
)

fit.hmm.list[[rr]] = fit.hmm4
fit.hmm.like[rr] = max(fit.hmm4$loglik)

}

indx = which.max( fit.hmm.like )
est.ss <- predict( fit.hmm.list[[indx]], train )

#est.ss <- predict( fit.hmm4, train )

### Max Snowfall ###

set.seed(137)

## Rep it!
reps=20
fit.hmm.list = list();
fit.hmm.like = c();
for( rr in 1:reps ){

#stcnt = 4
tmat.init <- matrix(
  -log( runif(stcnt^2) ),stcnt,stcnt
)
tmat.init <- (tmat.init/rowSums(tmat.init))

edm.mx <- edm.h[,365]
nr = length(edm.mx)

mu.list = list()
sig.list= list()
for( i in 1:stcnt ){
  mu.list[[i]]  = sample(edm.mx,1)
  sig.list[[i]] = sd(edm.mx[sample(nr,nr/2)])
}

hmm.model <- hmmspec( 
  init = rep(1/stcnt,stcnt),
  trans= tmat.init,
  parms.emission=list(
    mu=sample(edm.mx,nr), 
    sigma=rep(sd(edm.mx[sample(nr,nr/2)]),nr)
  ),
  dens.emission=dnorm.hsmm
)

train <- simulate(
  hmm.model, nsim=nr, seed=137, 
  rand.emis=rnorm.hsmm
)

train$s <- c(rep(1,nr/2),rep(2,nr/2));
train$x <- edm.mx;

fit.hmm4 = hmmfit(
  train,hmm.model,mstep=mstep.norm
)

fit.hmm.list[[rr]] = fit.hmm4
fit.hmm.like[rr] = max(fit.hmm4$loglik)

}

indx = which.max( fit.hmm.like )
est.edm.mx <- predict( fit.hmm.list[[indx]], train )

### ARI Matrix ###

arimat <- matrix( 0,5,5 );

arimat[1,1] = adjustedRandIndex( out.edm.l2$bestSeq, out.edm.l2$bestSeq )
arimat[1,2] = adjustedRandIndex( out.edm.l2$bestSeq, out.edm.w21$bestSeq )
arimat[1,3] = adjustedRandIndex( out.edm.l2$bestSeq, out.edm.w22$bestSeq)
arimat[1,4] = adjustedRandIndex( out.edm.l2$bestSeq, est.ss$s)
arimat[1,5] = adjustedRandIndex( out.edm.l2$bestSeq, est.edm.mx$s)
arimat[2,1] = adjustedRandIndex( out.edm.w21$bestSeq, out.edm.l2$bestSeq )
arimat[2,2] = adjustedRandIndex( out.edm.w21$bestSeq, out.edm.w21$bestSeq )
arimat[2,3] = adjustedRandIndex( out.edm.w21$bestSeq, out.edm.w22$bestSeq)
arimat[2,4] = adjustedRandIndex( out.edm.w21$bestSeq, est.ss$s)
arimat[2,5] = adjustedRandIndex( out.edm.w21$bestSeq, est.edm.mx$s)
arimat[3,1] = adjustedRandIndex( out.edm.w22$bestSeq, out.edm.l2$bestSeq )
arimat[3,2] = adjustedRandIndex( out.edm.w22$bestSeq, out.edm.w21$bestSeq )
arimat[3,3] = adjustedRandIndex( out.edm.w22$bestSeq, out.edm.w22$bestSeq)
arimat[3,4] = adjustedRandIndex( out.edm.w22$bestSeq, est.ss$s)
arimat[3,5] = adjustedRandIndex( out.edm.w22$bestSeq, est.edm.mx$s)
arimat[4,1] = adjustedRandIndex( est.ss$s, out.edm.l2$bestSeq )
arimat[4,2] = adjustedRandIndex( est.ss$s, out.edm.w21$bestSeq )
arimat[4,3] = adjustedRandIndex( est.ss$s, out.edm.w22$bestSeq)
arimat[4,4] = adjustedRandIndex( est.ss$s, est.ss$s)
arimat[4,5] = adjustedRandIndex( est.ss$s, est.edm.mx$s)
arimat[5,1] = adjustedRandIndex( est.edm.mx$s, out.edm.l2$bestSeq )
arimat[5,2] = adjustedRandIndex( est.edm.mx$s, out.edm.w21$bestSeq )
arimat[5,3] = adjustedRandIndex( est.edm.mx$s, out.edm.w22$bestSeq)
arimat[5,4] = adjustedRandIndex( est.edm.mx$s, est.ss$s)
arimat[5,5] = adjustedRandIndex( est.edm.mx$s, est.edm.mx$s)


### Plot some snowfall ###

bestseq = out.edm.l2$bestSeq
meanfal = out.edm.l2$stPar

bestseq = out.edm.w21$bestSeq
meanfal = out.edm.w21$stPar

bestseq = est.ss$s

bestseq = est.edm.mx$s

matplot( 
  t(edm.sm),lty=bestseq,
  col=bestseq+1,type='l',lwd=2,
  las=1,ylab="cumulative snowfall (cm)",
  xlab="Month",xaxt='n',main="Two State THMM"
)
axis( 
  side = 1,
  at = cumsum(c(0,30,31,30,31,30,31,31,28,31,30,31,30)),
  labels = c(
    "Jul","Aug","Sep","Oct","Nov","Dec","Jan",
    "Feb","Mar","Apr","May","Jun","Jul"
  ), las=2 
)
matplot( 
  t(meanfal),
  lwd=4,lty=1:2,type='l',
  add=T,col=1 
)

##
plot(
  1940:1989,edm.mx,type='l',lwd=2,las=1,
  ylab="Total Snowfall (cm)", main="Univariate HMM"
)
points(
  1940:1989,edm.mx,col=est.edm.mx$s+1,
  pch=15+est.edm.mx$s,cex=2
)
##
plot(
  prc.data[,1:2],col=est.ss$s+1,pch=est.ss$s+15,
  xlab="PrComp 1", ylab="PrComp 2"
)


