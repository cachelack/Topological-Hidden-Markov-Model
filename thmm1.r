
#
#  Random Dirichlet Draw
#

rdirichlet <- function( n, alf ){
  len = length(alf)
  if(len<2){
    print("alf must be > 1")
    return(NULL)
  }
  toRet <- matrix(0,nrow = n, ncol = len);
  toRet <- matrix( 
    rgamma(n*len,shape = alf,rate = 1),
    nrow = n, ncol = len, byrow = T
  );
  return(toRet/rowSums(toRet));
}

#
#  emmision function for BMwD
#

emm.bj <- function(thmm, state, time, log.scale=T){
  obs = thmm$data[time,];   # t-th observation
  muj = thmm$stPar[state,]; # s-th drift param
  len = length(obs);        # length of obs
  bj  = abs( diff(obs)-muj )^2
  bj <- -0.5*sum(bj)/len 
  if(log.scale){
    return(bj)
  } else {
    return(exp(bj))
  }
}

thmm.bj.all <- function(thmm, log.scale=T,type="BMWD",norm="W21"){
  obs = thmm$data;  # observations
  len = ncol(obs);  # length of obs
  if(type!="BM2d"){
    if(type=="nonparMulti"){
      obs0 <- (obs[,-1,]+obs[,-len,])/2;
      obs1 <- obs[,-1,]-obs[,-len,] # first diff
    } else {
      obs0 <- (obs[,-1]+obs[,-len])/2;
      obs1 <- obs[,-1]-obs[,-len] # first diff
    }
  } else{
    obs0 <- (obs[,,-1]+obs[,,-len])/2;
    obs1 <- obs[,,-1]-obs[,,-len] # first diff
  }
  if(type=="BMWD"){
    muj = thmm$stPar; # drift params
    bj  = matrix(0,nrow(obs),thmm$stCnt)
    cnst= log(2*pi)
    for( ii in 1:thmm$stCnt ){
      bj[,ii]  = -0.5*apply( 
        abs( obs1-muj[ii]/len ), 1, sum
      )^2 - cnst;
    }
  } else if(type=="BM2d"){
    muj1 = thmm$stPar[1,]; # drift params
    muj2 = thmm$stPar[2,];
    bj  = matrix(0,nrow(obs[1,,]),thmm$stCnt)
    for( ii in 1:thmm$stCnt ){
      bj[,ii]  = -0.5*(
        apply( 
          abs( obs1[1,,]-muj1[ii]/len ), 1, sum
        )^2 + apply( 
          abs( obs1[2,,]-muj2[ii]/len ), 1, sum
        )^2
      );
    }
  } else if(type=="OU"){
    muj = thmm$stPar[1,]; # drift params
    thj = thmm$stPar[2,]; # theta params
    bj  = matrix(0,nrow(obs),thmm$stCnt)
    cnst= log(2*pi)
    for( ii in 1:thmm$stCnt ){
      bj[,ii]  = -0.5*apply( 
        ( obs1-(muj[ii]-thj[ii]*obs0)/len ), 1, sum
      )^2 + thj[ii]/2 - cnst;
    }
  } else if(type=="nonpar"){
    muj = thmm$stPar; # mean curves
    bj  = matrix(0,nrow(obs),thmm$stCnt)
    for( ii in 1:thmm$stCnt ){
      if(norm=="L2"){
        norm.ii = apply( 
          t(t(obs)-muj[ii,])^2,1,sum
        )
      } else if(norm=="W21"){
        norm.ii = apply(
          apply(
            t(t(obs)-muj[ii,]),1,diff
          )^2,2,sum
        )
      } else if(norm=="W22"){
        norm.ii = apply(
          apply(
            apply(
              t(t(obs)-muj[ii,]),1,diff
            ),2,diff
          )^2,
          2,sum
        )
      } else if(norm=="sup"){
        norm.ii = apply( 
          abs(t(t(obs)-muj[ii,])),1,max
        )
      }
      bj[,ii]  = -0.5*( norm.ii )
    }
  } else if(type=="nonparMulti"){
    ncurv <- dim(obs)[3]
    muj = thmm$stPar; # mean curves
    bj  = matrix(0,nrow(obs),thmm$stCnt)
    for( ii in 1:thmm$stCnt ){
      #norm.l2 = apply( 
      #  t(t(obs)-muj[ii,])^2,1,sum
      #)
      norm.W2 = 0
      for( jj in 1:ncurv )
        norm.W2 = norm.W2 +
          apply(
            apply(
              t(t(obs[,,jj])-muj[ii,,jj]),1,diff
            )^2,2,sum
          )
      #norm.W22 = apply(
      #  apply(
      #    apply(
      #      t(t(obs)-muj[ii,]),1,diff
      #    ),2,diff
      #  )^2,
      #  2,sum
      #)
      bj[,ii]  = -0.5*( norm.W2 )
    }
  }
  
  if(log.scale){
    thmm$bj <- bj
  } else {
    thmm$bj <- exp(bj)
  }
  return(thmm)
}


#
#  Initialize the THMM
#

thmm.init <- function(
  data, state_mean, prob.trans=NULL, 
  prob.init=NULL,type="BMWD"
){
  thmm = list();
  thmm$data = data;
  if(type=="BMWD"){
    thmm$stCnt= length(state_mean);
  } else if(type=="BM2d") {
    thmm$stCnt= length(state_mean)/2;
  } else if(type=="OU") {
    thmm$stCnt= length(state_mean)/2;
  } else if(type=="nonpar" || type=="nonparMulti"){
    thmm$stCnt= nrow(state_mean);
  }
  
  thmm$stPar= state_mean;
  if(is.null(prob.trans)){
    thmm$prTran = rdirichlet(
      thmm$stCnt,rep(6,thmm$stCnt)
    )
  } else {
    thmm$prTran = prob.trans
  }
  if(is.null(prob.init)){
    thmm$prInit = rdirichlet(
      1,rep(6,thmm$stCnt)
    )
  } else {
    thmm$prInit = prob.init
  }
  return(thmm)
}

#
# Sum logs in a stable way
#

sumLogs <- function(mat,byRow=T){
  if(!byRow){
    mat <- t(mat);
  }
  nr = nrow(mat);
  nc = ncol(mat);
  ret= mat[1,]
  for( ii in 2:nr ){
    rowMin = pmin(mat[ii,],ret)
    rowMax = pmax(mat[ii,],ret)
    ret = rowMax +
      log( 1 + exp(rowMin-rowMax)  )
  }
  return(ret);
}

sumLogs.vec <- function(vec){
  len = length(vec);
  ret= vec[1]
  for( ii in 2:len ){
    vMax = max(vec[ii],ret)
    vMin = min(vec[ii],ret)
    ret = vMax +
      log( 1 + exp(vMin-vMax)  )
  }
  return(ret);
}

sumLogs.arr <- function(arr){
  dm  = dim(arr)
  len = dm[1];
  ret= arr[1,,];
  for( ii in 2:len ){
    arrMin = pmin(arr[ii,,],ret)
    arrMax = pmax(arr[ii,,],ret)
    ret = arrMax +
      log( 1 + exp(arrMin-arrMax)  )
  }
  return(ret);
}

sumLogs.arr2 <- function(arr){
  dm  = dim(arr)
  ret= arr[,1,1];
  for( ii in 1:dm[2] ){
    for( jj  in 1:dm[3]){
      if((ii==1)&(jj==1))
        next;
      vMin = pmin(arr[,ii,jj],ret)
      vMax = pmax(arr[,ii,jj],ret)
      ret = vMax +
        log( 1 + exp(vMin-vMax)  )
    }
  }
  return(ret);
}

#
# Alpha Pass
#

thmm.alfa <- function(thmm,log.scale=T,type){
  if(type!="BM2d"){
    TT   = nrow(thmm$data); # number of obs
  } else{
    TT   = dim(thmm$data)[2]; # number of obs
  }
  PP   = thmm$stCnt;      # number of states
  alfa = matrix(0,TT,PP); 
  alfa[1,] <- (thmm$prInit) + thmm$bj[1,];
  for( tt in 2:TT ){
    for( jj in 1:PP ){
      tmp = alfa[tt-1,] + (thmm$prTran)[,jj] +
        thmm$bj[tt,jj];
      alfa[tt,jj] = sumLogs.vec(tmp);
      if(is.infinite(alfa[tt,jj])){
        print(cbind(alfa[tt-1,],(thmm$prTran)[,jj],
                      thmm$bj[tt,jj],tmp))
      }
        
    }
    #tmp = t(t(alfa[tt-1,] + log(thmm$prTran)) + thmm$bj[tt,])
    #alfa[tt,] <- sumLogs(tmp)
  }
  thmm$alfa <- alfa;
  return(thmm);
}

#
# beta Pass
#

thmm.beta <- function(thmm,log.scale=T,type){
  if(type!="BM2d"){
    TT   = nrow(thmm$data); # number of obs
  } else{
    TT   = dim(thmm$data)[2]; # number of obs
  }  
  PP   = thmm$stCnt;      # number of states
  beta = matrix(0,TT,PP); 
  beta[TT,] <- 0;  # 0 = log( 1 )
  for( tt in (TT-1):1 ){
    for( jj in 1:PP){
      tmp = beta[tt+1,] + (thmm$prTran)[jj,] +
        thmm$bj[tt+1,]
      beta[tt,jj] = sumLogs.vec(tmp);
      #if(is.infinite(beta[tt,jj]))
      #  beta[tt,jj] <- beta[tt+1,jj]
    }
    #tmp = t(beta[tt+1,] + t(log(thmm$prTran)) + thmm$bj[tt+1,])
    #beta[tt,] <- sumLogs(tmp,byRow = F)
  }
  thmm$beta <- beta;
  return(thmm);
}

#
# Reestimation
#

thmm.rest <- function(thmm,log.scale=T,type="BMWD",hurst=hurst){
  if(type!="BM2d"){
    TT   = nrow(thmm$data); # number of obs
  } else{
    TT   = dim(thmm$data)[2]; # number of obs
  }
  PP   = thmm$stCnt;      # number of states
  xi  = array(0,dim=c(TT-1,PP,PP));
  rng = 1:(TT-1)
  # Compute Gammas
  tmp     = thmm$alfa + thmm$beta;
  tmp.sum = sumLogs(tmp,byRow = F);
  gam = tmp - tmp.sum;

  # Compute Xis
  for( ii in 1:PP ){  # Double Loop, can we do better?!?
    for( jj in 1:PP ){
      tmp = 
        thmm$alfa[rng,ii] + (thmm$prTran[ii,jj]) +
        thmm$bj[rng+1,jj] + thmm$beta[rng+1,jj];
      xi[,ii,jj] = tmp #- sumLogs.vec(tmp);
    }
  }
  xi <- xi - sumLogs.arr2(xi)
  
  # Update init probs
  thmm$prInit <- (gam[1,]);
  # Update tran probs
  thmm$prTran <- ( sumLogs.arr(xi) - sumLogs(gam[rng,]));
  # Update drifts for BMwD
  if(type=="BMWD"){
    nc   = ncol(thmm$data);
    if(hurst==1/2){
      odif = thmm$data[,nc] - thmm$data[,1];
      coef = 1;
    } else {
      diffDat = t(apply(thmm$data,1,diff));
      nu   = 1/2-hurst;
      coef = gamma( 2*nu+2 )/gamma( nu+1 )
      tt   = seq(1:(nc-1))/nc;
      odif = apply(
        t(tt^(nu)*t(diffDat)),1,sum
      )
    }
    alfbet = sumLogs( thmm$alfa + thmm$beta );
    indx.p = which(odif>0);
    indx.n = which(odif<0);
    num.p  = sumLogs((thmm$alfa + thmm$beta + log(abs(odif)))[indx.p,])
    num.n  = sumLogs((thmm$alfa + thmm$beta + log(abs(odif)))[indx.n,])
    thmm$stPar <- coef*( exp(num.p-alfbet) - exp(num.n-alfbet) );
    #print(cbind(thmm$alfa,thmm$beta,thmm$bj))
  } else if(type=="BM2d"){
    nc   = dim(thmm$data)[3];
    odif = thmm$data[,,nc] - thmm$data[,,1];
    alfbet = sumLogs( thmm$alfa + thmm$beta );
    # do dim 1
    indx.p = which(odif[1,]>0);
    indx.n = which(odif[1,]<0);
    num.p  = sumLogs((thmm$alfa + thmm$beta + log(abs(odif[1,])))[indx.p,])
    num.n  = sumLogs((thmm$alfa + thmm$beta + log(abs(odif[1,])))[indx.n,])
    thmm$stPar[1,] <- exp(num.p-alfbet) - exp(num.n-alfbet);
    # do dim 2
    indx.p = which(odif[2,]>0);
    indx.n = which(odif[2,]<0);
    num.p  = sumLogs((thmm$alfa + thmm$beta + log(abs(odif[2,])))[indx.p,])
    num.n  = sumLogs((thmm$alfa + thmm$beta + log(abs(odif[2,])))[indx.n,])
    thmm$stPar[2,] <- exp(num.p-alfbet) - exp(num.n-alfbet);
  } else if(type=='OU'){
    nc   = ncol(thmm$data);
    odif = thmm$data[,-1] - thmm$data[,-nc];
    oave = (thmm$data[,-1] + thmm$data[,-nc])/2;
    alfbet = sumLogs( thmm$alfa + thmm$beta );
    #print(cbind(thmm$alfa,thmm$beta,thmm$bj))
    if(T){
    f2Opt.org <- function(x){
      slop <- x[2];
      intc <- x[1];
      bj <- 0.5*apply( 
        abs( odif-(intc-slop*oave)/nc ), 1, sum
      )^2 - slop/2;
      indx.p = which(bj>0);
      indx.n = which(bj<0);
      ret = 0;
      if(length(indx.n)==0){
        num.p  = sumLogs.vec(
          (thmm$alfa[,ii] + thmm$beta[,ii] + log(abs(bj)))[indx.p]
        )
        ret = exp(num.p-alfbet[ii]) 
      } else if(length(indx.p)==0) {
        num.n = sumLogs.vec(
          (thmm$alfa[,ii] + thmm$beta[,ii] + log(abs(bj)))[indx.n]
        )
        ret = -exp(num.n-alfbet[ii])
      } else {
        num.p  = sumLogs.vec(
          (thmm$alfa[,ii] + thmm$beta[,ii] + log(abs(bj)))[indx.p]
        )
        num.n = sumLogs.vec(
          (thmm$alfa[,ii] + thmm$beta[,ii] + log(abs(bj)))[indx.n]
        )
        ret = exp(num.p-alfbet[ii]) - exp(num.n-alfbet[ii])
      }
      #print(ret)
      return(ret)
    } # OU OM functional 
    f2Opt.dif <- function(x){
      slop <- x[2];
      intc <- x[1];
      bj.h <- apply( 
        ( odif-(intc-slop*oave)/nc )*oave, 1, sum
      ) - 1/2;
      bj.k <- -apply( 
        ( odif-(intc-slop*oave)/nc ), 1, sum
      );
      indx.p.h = which(bj.h>0);
      indx.n.h = which(bj.h<0);
      indx.p.k = which(bj.k>0);
      indx.n.k = which(bj.k<0);
      ret1 = 0;
      ret2 = 0;
      if(length(indx.n.h)==0){
        num.p.h  = sumLogs.vec(
          (thmm$alfa[,ii] + thmm$beta[,ii] + log(abs(bj.h)))[indx.p.h]
        )
        ret1 = exp(num.p.h-alfbet[ii]); 
      } else if(length(indx.p.h)==0){
        num.n.h = sumLogs.vec(
          (thmm$alfa[,ii] + thmm$beta[,ii] + log(abs(bj.h)))[indx.n.h]
        )
        ret1 = - exp(num.n.h-alfbet[ii])
      } else {
        num.p.h  = sumLogs.vec(
          (thmm$alfa[,ii] + thmm$beta[,ii] + log(abs(bj.h)))[indx.p.h]
        )
        num.n.h = sumLogs.vec(
          (thmm$alfa[,ii] + thmm$beta[,ii] + log(abs(bj.h)))[indx.n.h]
        )
        ret1 = exp(num.p.h-alfbet[ii]) - exp(num.n.h-alfbet[ii])
      }
      
      if(length(indx.n.k)==0){
        num.p.k  = sumLogs.vec(
          (thmm$alfa[,ii] + thmm$beta[,ii] + log(abs(bj.k)))[indx.p.k]
        )
        ret2 = exp(num.p.k-alfbet[ii])
      } else if(length(indx.p.k)==0) {
        num.n.k = sumLogs.vec(
          (thmm$alfa[,ii] + thmm$beta[,ii] + log(abs(bj.k)))[indx.n.k]
        )
        ret2 = - exp(num.n.k-alfbet[ii])
      } else {
        num.p.k  = sumLogs.vec(
          (thmm$alfa[,ii] + thmm$beta[,ii] + log(abs(bj.k)))[indx.p.k]
        )
        num.n.k = sumLogs.vec(
          (thmm$alfa[,ii] + thmm$beta[,ii] + log(abs(bj.k)))[indx.n.k]
        )
        ret2 = exp(num.p.k-alfbet[ii]) - exp(num.n.k-alfbet[ii])
      }
      #print(c(ret1,ret2))
      return( c(ret1,ret2) )
    } # OU OM derivative
    #print(thmm$stPar)
    for( ii in 1:thmm$stCnt ){
      #print(f2Opt.org(thmm$stPar[,ii]))
      #print(f2Opt.dif(thmm$stPar[,ii]))
      best = optim( 
        thmm$stPar[,ii], f2Opt.org, f2Opt.dif, method = "L-BFGS-B",
        lower = c(-Inf,0)
      );
      thmm$stPar[,ii] <- best$par;
    }
    }
    if(F){
    tmp.i <- thmm$stPar[2,];
    tmp.s <- thmm$stPar[1,];
    for(jj in 1:10){
      old.i <- tmp.i;
      tmp.i <- solveOU.intercept(thmm$alfa,thmm$beta,thmm$data,tmp.s)
      tmp.s <- solveOU.slope(thmm$alfa,thmm$beta,thmm$data,old.i)
      #print(tmp.i)
      #print(tmp.s)
      #print("")
    }
    thmm$stPar[1,] = tmp.s;
    thmm$stPar[2,] = tmp.i;
    }
  } else if(type=="nonpar"){
    dat.shift = thmm$data - min(thmm$data)+1
    alfbet = sumLogs( thmm$alfa + thmm$beta );
    num.p = matrix(0,nrow=thmm$stCnt,ncol(thmm$data));
    for( ii in 1:thmm$stCnt ){
      num.p[ii,] = sumLogs((thmm$alfa[,ii] + thmm$beta[,ii] + log(dat.shift)))
    }
    
    thmm$stPar <- exp(num.p-alfbet) + min(thmm$data)-1;
  } else if(type=="nonparMulti"){
    ncurv <- dim(thmm$data)[3]
    dat.shift = thmm$data - min(thmm$data)+1
    alfbet = sumLogs( thmm$alfa + thmm$beta );
    num.p = array(0,dim=c(thmm$stCnt,ncol(thmm$data),ncurv));
    for( ii in 1:thmm$stCnt ){
      for( jj in 1:ncurv ){
        num.p[ii,,jj] = sumLogs(
          (thmm$alfa[,ii] + thmm$beta[,ii] + log(dat.shift[,,jj]))
        )
      }
    }
    
    thmm$stPar <- exp(num.p-alfbet) + min(thmm$data)-1;
  }
  
  return(thmm)
}

solveOU.intercept <- function(
  alf, bet, dat, slope
){
  stc  = length(slope)
  len  = ncol(dat);
  obs0 = dat[,-1];
  obs1 = dat[,-1]-dat[,-len];
  toRet <- rep(0,stc);
  for(ii in 1:stc){
    num  = dat[,len] - dat[,1] + slope[ii]*rowSums(dat)/len;
    alfbet = sumLogs.vec( alf[,ii] + bet[,ii] );
    indx.p = which(num>0);
    indx.n = which(num<0);
    if(length(indx.p)==0){
      num.n  = sumLogs.vec((alf[,ii] + bet[,ii] + log(abs(num)))[indx.n])
      num.p  = -10000
    } else if(length(indx.n)==0){
      num.p  = sumLogs.vec((alf[,ii] + bet[,ii] + log(abs(num)))[indx.p])
      num.n  = -10000
    } else{
      num.p  = sumLogs.vec((alf[,ii] + bet[,ii] + log(abs(num)))[indx.p])
      num.n  = sumLogs.vec((alf[,ii] + bet[,ii] + log(abs(num)))[indx.n])
    }
    toRet[ii] = exp(num.p-alfbet) - exp(num.n-alfbet);
  }
  return(toRet);
}

solveOU.slope <- function(
  alf, bet, dat, intercept
){
  stc  = length(intercept)
  len  = ncol(dat);
  obs0 = dat[,-1];
  obs1 = dat[,-1]-dat[,-len];
  toRet <- rep(0,stc);
  for(ii in 1:stc){
    num  = intercept[ii]*rowSums(dat)/len - rowSums( obs0*obs1 )/len + 1/2;
    alfbet = sumLogs.vec( alf[,ii] + bet[,ii] + log(rowSums(dat)^2/len) );
    indx.p = which(num>0);
    indx.n = which(num<0);
    if(length(indx.p)==0){
      num.n  = sumLogs.vec((alf[,ii] + bet[,ii] + log(abs(num)))[indx.n])
      toRet[ii] = -exp(num.n-alfbet);
    } else if(length(indx.n)==0){
      num.p  = sumLogs.vec((alf[,ii] + bet[,ii] + log(abs(num)))[indx.p])
      toRet[ii] = exp(num.p-alfbet);
    } else{
      num.p  = sumLogs.vec((alf[,ii] + bet[,ii] + log(abs(num)))[indx.p])
      num.n  = sumLogs.vec((alf[,ii] + bet[,ii] + log(abs(num)))[indx.n])
      toRet[ii] = exp(num.p-alfbet) - exp(num.n-alfbet);
    }
  }
  return(toRet);
}

#
#  Auto Generate Random Parameters to start
#

autoParam <- function( data, stateCnt, type){
  if(type!="BM2d"){
    nr = nrow(data);
    nc = ncol(data);
  } else {
    dd = dim(data)
    nr = dd[2]
    nc = dd[3]
  }
  
  toRet = c();
  if(type=="BMWD"){
    odif = data[,nc] - data[,1];
    toRet = sample(odif,stateCnt)
  } else if(type=="BM2d"){
    odif1 = data[1,,nc] - data[1,,1];
    odif2 = data[2,,nc] - data[2,,1];
    mag = sqrt(odif1^2 + odif2^2)
    rot = matrix( 
      c(
        cos(2*pi/stateCnt),sin(2*pi/stateCnt),
        -sin(2*pi/stateCnt),cos(2*pi/stateCnt)),
      2,2
    )
    toRet = matrix(0,2,stateCnt)
    toRet[,1] = c(median(mag),0);
    for(ii in 2:stateCnt)
      toRet[,ii] = rot%*%toRet[,ii-1];
  } else if(type=="OU"){
    obs0 = data[,-1];
    obs1 = data[,-1]-data[,-nc];
    intercept = data[,nc] - data[,1];
    slope     = (intercept*rowSums(data)/nc - rowSums( obs0*obs1 )/nc + 1/2)
    #indx = sample(nr,stateCnt);
    slp = 2*median(slope)+rnorm(stateCnt,0,0.5)
    rng = seq(min(intercept),max(intercept),length.out=stateCnt+2)
    toRet = rbind( rng[2:(stateCnt+1)]*slp, slp  )
  } else if(type=="nonpar"){
    indx  = sample(nr,stateCnt)
    toRet = data[ indx, ];
  } else if(type=="nonparMulti"){
    indx  = sample(nr,stateCnt)
    toRet = data[ indx,,];
  }
  return(toRet)
}

#
# Main Algorithm
#

thmm.BaumWelch <- function( 
  data, stateCnt, state_mean=NULL, prob.trans=NULL, prob.init=NULL,
  maxIter = 100, minIter = 10, verbose=F, tol = 1e-3, viterbi=T,type="BMWD",
  PLOT=F,dim=1,hurst=1/2,truth=NULL, norm="W21"
){
  # Initialize
  if(is.null(state_mean))
    state_mean <- autoParam( data, stateCnt, type)
  thmm = thmm.init(
    data, state_mean, prob.trans, prob.init,type=type
  )
  thmm$prInit <- log(thmm$prInit)
  thmm$prTran <- log(thmm$prTran)
  if(verbose==2){
    #TT   = nrow(thmm$alfa);
    #like = sumLogs.vec( thmm$alfa[TT,] )
    #print(paste("Iter",ii,like))
    print("Initalization:")
    print("init prob")
    print(thmm$prInit)
    print("tran prob")
    print(thmm$prTran)
    print("drift")
    print(thmm$stPar)
    print("")
  }
  # Main Loop
  like.old = 0;
  if(verbose==0)
    pb = txtProgressBar(min = 0,max=maxIter,style = 3);
  for( ii in 1:maxIter  ){
    # Update Emission Probs
    thmm <- thmm.bj.all(thmm, log.scale=T,type=type,norm=norm);
    # Alpha & Beta Pass
    thmm <- thmm.alfa(thmm,log.scale=T,type);
    thmm <- thmm.beta(thmm,log.scale=T,type);
    # Reestimation
    thmm <- thmm.rest(thmm,log.scale=T,type=type,hurst=hurst);
    # Plot it
    if(PLOT){
      thmm <- thmm.viterbi(thmm,type)
      if(type!="BM2d"){
        if(type=="nonparMulti"){
          ncurv = dim(data)[3];
          par(mfcol=c(1,ncurv));
          for( jj in 1:ncurv ){
            matplot(t(data[,,jj]),type='l',col=thmm$bestSeq+1,main=paste("Iteration:",ii))
            matplot( t(thmm$stPar[,,jj]),col='black',type='l',lwd=3,add=T  )
          }
        } else {
          matplot(t(data),type='l',col=thmm$bestSeq+1,main=paste("Iteration:",ii))
          if(type=="nonpar")
            matplot( t(thmm$stPar),col='black',type='l',lwd=3,add=T  )
        }
      } else {
        matplot(
          t(data[1,,]),t(data[2,,]),type='l',
          col=thmm$bestSeq,main=paste("Iteration:",ii)
        )
      }
    }
    # Compute "Likelihood"
    if(ii>1){
      like.old = like 
    }
    TT   = nrow(thmm$alfa);
    like = sumLogs.vec( thmm$alfa[TT,] )
    if(verbose==0)
      setTxtProgressBar(pb,ii)
    if(verbose==1){
      if(is.null(truth)){
        if(ii==1){
          print(paste(ii,like))
        } else {
          print(paste(ii,like,(like.old-like)/like))
        }
      } else {
        thmm <- thmm.viterbi(thmm,type)
        ari = adjustedRandIndex(thmm$bestSeq,truth)
        if(ii==1){
          print(paste(ii,like,ari))
        } else {
          print(paste(ii,like,(like.old-like)/like,ari))
        }
      }
    }
    if(verbose==2){
      TT   = nrow(thmm$alfa);
      like = sumLogs.vec( thmm$alfa[TT,] )
      if(ii==1){
        print(paste("Iter",ii,like))
      } else {
        print(paste("Iter",ii,like,(like.old-like)/like))
      }
      print("init prob")
      print(thmm$prInit)
      print("tran prob")
      print(thmm$prTran)
      print("drift")
      print(thmm$stPar)
      print("")
    }
    if(is.infinite(like))
      break;
    if( (ii>minIter) &((like.old-like)/like < tol) )
      break;
  }
  if(verbose==0)
    close(pb)
  if(verbose<2){
    TT   = nrow(thmm$alfa);
    like = sumLogs.vec( thmm$alfa[TT,] )
    print(paste("Iter",ii,like))
    print("init prob")
    print(thmm$prInit)
    print("tran prob")
    print(thmm$prTran)
    if(type!="nonpar" & type!="nonparMulti"){
      print("drift")
      print(thmm$stPar) 
    }
    print("")
  }
  # Run Viterbi to get best state seq
  if(viterbi)
    thmm <- thmm.viterbi(thmm,type)
  return(thmm);
}

#
#  Viterbi Algorithm
#

library(mclust)

thmm.viterbi <- function(thmm,type="BMWD"){
  if(type!="BM2d"){
    TT   = nrow(thmm$data); # number of obs
  } else{
    TT   = dim(thmm$data)[2]; # number of obs
  }  
  PP  = thmm$stCnt;      # number of states
  obs = thmm$data;
  del = array(0,dim = c(TT,PP));
  phi = array(0,dim = c(TT,PP));
  del[1,] <- (thmm$prInit) + thmm$bj[1,];
  for( tt in 2:TT ){ # Forward Pass
    del[tt,] <- apply(
      del[tt-1,]+(thmm$prTran),2,max
    ) + thmm$bj[tt,]
    phi[tt,] <- apply(
      del[tt-1,]+(thmm$prTran),2,which.max
    )
  }
  thmm$bestProb <- max(del[TT,]);
  bestSeq = rep(0,TT);
  bestSeq[TT] <- which.max(del[TT,]);
  for( tt in (TT-1):1  ){ # Backtrack
    bestSeq[tt] <- phi[tt+1,bestSeq[tt+1]];
  }
  thmm$bestSeq <- as.integer(bestSeq);
  return(thmm)
}


#
#  Simulate some data
#

rWiener <- function( n, len, drift=0 ){
  toRet = matrix(
    rnorm(n*len),nrow = n,ncol = len
  ) + drift/sqrt(len);
  if( n == 1){
    toRet <- cumsum(toRet)/sqrt(len)
  } else if( n > 1){
    toRet <- apply(toRet,1,cumsum)/sqrt(len)
  }
  return(toRet);
}

rGaussianProc <- function( n, len, drift=0, sqrtcov ){
  toRet = matrix(
    rnorm(n*len),nrow = n,ncol = len
  );
  toRet = t(sqrtcov%*%t(toRet));
  dterm = (1:len)*drift/len;
  if( n == 1){
    toRet <- toRet + dterm
  } else if( n > 1){
    toRet <- t(t(toRet)+dterm)
  }
  return(toRet);
}

rOrnUhl <- function( n, len, drift=0, theta=0 ){
  if( n == 1){
    wn = rnorm(len,mean = 0,sd = 1/sqrt(len))
    toRet = rep(0,len+1)
    for(i in 1:len){
      toRet[i+1] = toRet[i] + theta*(drift-toRet[i])/len  + wn[i]
    }
    toRet <- toRet[-1]
  } else if( n > 1){
    wn = matrix(
      rnorm(n*len,mean = 0,sd = 1/sqrt(len)),nrow = n,ncol = len
    )
    toRet = matrix(0,nrow=n,ncol=len+1)
    for(i in 1:len){
      toRet[,i+1] = toRet[,i] + theta*(drift-toRet[,i])/len  + wn[,i]
    }
    toRet <- toRet[,-1]
  }
  return(toRet);
}


rMarkov <- function( n, prInit, prTran ){
  toRet = which(rmultinom(1,1,prInit)==1)
  for( ii in 2:n ){
    toRet[ii] = which(
      rmultinom(1,1,prTran[toRet[ii-1],])==1
    );
  }
  return(toRet)
}

sim.thmm <- function( 
  n=100, 
  prInit=c(0.5,0.3,0.2), 
  prTran=matrix( c(0.8,0.1,0.1,0.4,0.2,0.4,0.1,0.4,0.5),3,3,byrow = T ), 
  drift=c(-2,0,2),  
  theta=c(2,10,5),
  type="BMWD",
  smooth=0,
  hurst=1/2,
  len = 1000
){
  dat   = c();
  if(hurst!=1/2){
    # compute square root of covariance matrix
    # for fractional brownian motion
    tt = (1:len)/len;
    tmpmat1= outer( tt^(2*hurst),tt^(2*hurst),FUN="+" );
    tmpmat2= abs(outer( tt,tt,FUN='-' ))^(2*hurst)
    covmat = ( tmpmat1 - tmpmat2 )/2;
    eigcov = eigen(covmat,symmetric = T);
    sqrtcov= (eigcov$vectors)%*%diag(sqrt(eigcov$values))%*%t(eigcov$vectors)
  }
  if(type=="BM2d")
    dat = array(0,dim = c(2,n,len))
  stseq = rMarkov(n,prInit,prTran);
  for( ii in 1:n ){
    if(type=="BMWD"){
      if(hurst==1/2){
        if(smooth==0){
          dat = rbind(
            dat,rWiener(1,len=len,drift = drift[stseq[ii]] )
          )
        } else if(smooth>0){
          dat = rbind(
            dat,lowess(
              1:len,rWiener(1,len=len,drift = drift[stseq[ii]] ),
              f = smooth
            )$y
          )
        }
      } else {
        dat = rbind(
          dat,rGaussianProc(
            1,len=len,drift = drift[stseq[ii]],sqrtcov=sqrtcov
          )
        )
      }
    } else if(type=="BM2d"){
      dat[1,ii,] = rWiener(1,len=len,drift = drift[1,stseq[ii]] )
      dat[2,ii,] = rWiener(1,len=len,drift = drift[2,stseq[ii]] )
    } else if(type=="OU"){
      if(smooth==0){
        dat = rbind(
          dat,rOrnUhl(
            1,len=len,drift = drift[stseq[ii]], 
            theta = theta[stseq[ii]]  
          )
        )
      } else if(smooth>0){
        dat = rbind(
          dat,lowess(
            1:len,rOrnUhl(
              1,len=len,drift = drift[stseq[ii]], 
              theta = theta[stseq[ii]]  
            ),
            f = smooth
          )$y
        )
      }
    } else if(type=="nonpar"){
      dat = rbind(
        dat,rSmoothError(
          drift = drift[stseq[ii],]
        )
      )
    }
  }
  return(list(
    data=dat, ss = stseq
  ))
}

rSmoothError <- function( drift=rep(0,100) ){
  len= length(drift)
  kk = 16;
  tt = seq(0,1,length.out=len);
  sig= 0.4;
  cfs= rnorm( kk, 0, sig );
  kmat= matrix( 1:kk, kk, len );
  tmat= matrix( tt, kk, len, byrow=T )
  zmat= matrix( cfs, kk, len )
  err= colSums( 
    sqrt(2)*zmat*sin( kmat*pi*tmat )/( kmat*pi )
  )
  return( drift + err );
}

normalizeAlfBet <- function( ab ){
  rs  <- sumLogs(ab,byRow=F)
  ret <- ab - rs;
  return(exp(ret)) 
}


