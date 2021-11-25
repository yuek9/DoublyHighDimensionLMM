################################
## functions to estimate fixed effects and construct confidence intervals; and function to compute variance components
## adapted from Sai Li , T. Tony Cai & Hongzhe Li (2021): Inference for high-dimensional linear mixed-effects models: A quasi-likelihood approach, Journal of the American Statistical Association, DOI: 10.1080/01621459.2021.1888740
## Kun Yue
## 11/04/2021
################################


#######################
## main function for fixed effect estimation and inference
#######################
#Inputs-- 
## X: design for fixed effects; 
## y: response; 
## Z: design for random effects;
## grp: a vector of length N indicating group membership; 
## a : tuning parameter in $\Sig_a$
## lm: linear model fitting (T or F, if T will ignore correlations among observations); 
## inf.coord: the coordinates of beta for inference.

#Outputs-- 
## beta.hat: the Lasso estimator based on psedo-likelihood; 
## beta.db: debiased Lasso estimators for the fixed effects in inf.coord;
## beta.db.sd: standard deviation of the debiased Lasso estimators for the fixed effects in inf.coord.

Fixed_effects_est_inf <- function(X, y, Z, grp, a=1, lm=F, inf.coord=NULL){
  
  N=length(y)
  n=length(unique(grp))
  q=ncol(Z)
  p=ncol(X)
  
  ### preprocessing
  y <- y - mean(y)
  X <- scale (X)
  X.sd <- attributes(X)$`scaled:scale` ## we need to divide these sd of X from the coefficients
  Z <- scale (Z)
  Z.sd <- attributes(Z)$`scaled:scale`
  
  lam.a.seq<-c(exp(seq(log(100),log(0.1), length.out = 20)))*sqrt(log(p)/N/a^2/q) ## this corresponds to lambda_a scale under glmnet default objective function
  lam.j.seq <- c(exp(seq(log(100),log(0.1), length.out = 20)))*sqrt(log(p)/N/a^2/q^2)

  if(lm){
    ## linear model fitting, no correlation considered
    ## use cvglmnet to get tuning paramter. it is usually of order sqrt(log(p)/N)
    sig.init<-scalreg(X,y)$hsigma ## scaled lasso
    beta.hat <- as.numeric(glmnet(X, y, lambda = sig.init*sqrt(2*log(p)/N))$beta)
    
    return(list(beta.hat=beta.hat/X.sd))
  
  }else{
    ## mixed effect model fitting
    X.a<-X
    y.a<-y
    tr.siga.inv<-0
    for (i in 1:n){
      cur.mem=which(grp==i)
      mi=length(cur.mem)
      Zi = as.matrix(Z[cur.mem,])
      
      sigmai = a*Zi%*%t(Zi)+diag(1,mi)
      sigmai.svd=svd(sigmai)
      Sig.ai.inv.half <- sigmai.svd$u %*% diag(1/sqrt(sigmai.svd$d)) %*% t(sigmai.svd$u) #compute sqrt of Sigmai
      
      X.a[cur.mem,] =  Sig.ai.inv.half %*% as.matrix(X[cur.mem,]) 
      y.a[cur.mem] =  Sig.ai.inv.half %*% as.matrix(y[cur.mem])
      tr.siga.inv<-tr.siga.inv + sum(1/sigmai.svd$d) #compute the effective sample size, which is the trace of Sigma.a.inv
    }
    
    
    ## use scaled-lasso to compute a tuning parameter
    #sig.init<-scalreg(X.a,y.a)$hsigma #scaled-lasso to get tuning parameter    
    #beta.hat= as.numeric(glmnet(X.a, y.a, lambda =sig.init*sqrt(1.2*log(p)/N))$beta)

    ## use cv.glmnet to select lambda value
    cv.init<-cv.glmnet(X.a, y.a, lam=lam.a.seq)
    ## double check if the tuning parameter range is reasonable
    if(which.lambda(cv.init$lambda.min, lam.a.seq) == length(lam.a.seq)){
      lam.a.seq.tmp = c(exp(seq(log(10),log(1e-4), length.out = 20)))*sqrt(log(p)/N/a^2/q)
      cv.init<-cv.glmnet(X.a, y.a, lam=lam.a.seq.tmp)
      if(which.lambda(cv.init$lambda.min, lam.a.seq.tmp) ==  length(lam.a.seq.tmp)){
        warning(paste0('lambda.a selected on boundary: ', which.lambda(cv.init$lambda.min, lam.a.seq.tmp), ' selected' ))
      }
    }
    if(which.lambda(cv.init$lambda.min, lam.a.seq)== 1){
      lam.a.seq.tmp = c(exp(seq(log(1e6),log(10), length.out = 40)))*sqrt(log(p)/N/a^2/q)
      cv.init<-cv.glmnet(X.a, y.a, lam=lam.a.seq.tmp)
      if(which.lambda(cv.init$lambda.min, lam.a.seq.tmp) ==  1 ){
        warning(paste0('lambda.a selected on boundary: ', which.lambda(cv.init$lambda.min, lam.a.seq.tmp), ' selected' ))
      }
    }
    
    beta.hat<-coef(cv.init, s=cv.init$lambda.min)[-1]
    
    ########
    ## inference for beta (some modification: directly use full-sample estimated beta.hat in this part)
    ########
    ### compute debiased Lasso for LMM, for the specified coordinates
    beta.db.sd<-rep(NA,length=length(inf.coord))
    beta.db<-rep(NA,length=length(inf.coord))
    if(is.null(inf.coord)){
      return( list(beta.hat = beta.hat/X.sd, beta.db=beta.db/X.sd, beta.db.sd=beta.db.sd/X.sd, tr.siga.inv=tr.siga.inv))
    }
    
    
    ## compute for each coordinate of interest
    for(j in 1:length(inf.coord)){
      col.j<-inf.coord[j]
      X.b<-X
      y.b<-y
      tr.sigb.inv<-0
      ## compute the correction score (wj, kappa.j)
      for (i in 1:n){
        cur.mem=which(grp==i)
        mi=length(cur.mem)
        Zi = as.matrix(Z[cur.mem,])
        if(col.j <= q){
          Zi.tmp = Zi[,-col.j]
        }else{
          Zi.tmp = Zi
        }
        sig.bi = a*Zi.tmp%*%t(Zi.tmp)+diag(1,mi)
        sig.bi.svd=svd(sig.bi)
        Sig.bi.inv.half <- sig.bi.svd$u %*% diag(1/sqrt(sig.bi.svd$d)) %*% t(sig.bi.svd$u) #compute sqrt of Sig.bi
        
        X.b[cur.mem,] =  Sig.bi.inv.half %*% as.matrix(X[cur.mem,]) 
        y.b[cur.mem] =  Sig.bi.inv.half %*% as.matrix(y[cur.mem])
        tr.sigb.inv<-tr.sigb.inv + sum(1/sig.bi.svd$d) 
      }
      
      cv.init.j<-cv.glmnet(X.b[, -col.j], X.b[, col.j], lam=lam.j.seq)
      ## double check if the tuning parameter range is reasonable
      if(which.lambda(cv.init.j$lambda.min, lam.j.seq) == length(lam.j.seq)){
        lam.j.seq.tmp = c(exp(seq(log(10),log(1e-4), length.out = 20)))*sqrt(log(p)/N/a^2/q^2)
        cv.init.j<-cv.glmnet(X.b[, -col.j], X.b[, col.j], lam=lam.j.seq.tmp)
        if(which.lambda(cv.init.j$lambda.min, lam.j.seq.tmp) ==  length(lam.j.seq.tmp)){
          warning(paste0('lambda.j for j=', col.j, ' selected on boundary: ', which.lambda(cv.init.j$lambda.min, lam.j.seq.tmp), ' selected' ))
        }
      }
      if(which.lambda(cv.init.j$lambda.min, lam.j.seq)== 1){
        lam.j.seq.tmp = c(exp(seq(log(1e6),log(10), length.out = 40)))*sqrt(log(p)/N/a^2/q^2)
        cv.init.j<-cv.glmnet(X.b[, -col.j], X.b[, col.j], lam=lam.j.seq.tmp)
        if(which.lambda(cv.init.j$lambda.min, lam.j.seq.tmp) ==  1 ){
          warning(paste0('lambda.j for j=', col.j, ' selected on boundary: ', which.lambda(cv.init.j$lambda.min, lam.j.seq.tmp), ' selected' ))
        }
      }
      

      
      kappa.hat.j<-coef(cv.init.j, s=cv.init.j$lambda.min)[-1]
      wj.hat<- X.b[,col.j]-  X.b[,-col.j]%*%kappa.hat.j
      beta.db[j] = beta.hat[col.j] + sum( wj.hat * (y.b - X.b %*% beta.hat ))/sum(wj.hat * X.b[,col.j])
      
      ## and compute sandwich estimates for variance
      num=0
      for(i in 1:n){
        cur.mem=which(grp==i)
        num <- num + sum(wj.hat[cur.mem]*(y.b[cur.mem] - X.b[cur.mem,] %*% beta.hat))^2
      }
      beta.db.sd[j] = sqrt(num)/sum(wj.hat*X.b[,col.j])
    }
    
    return( list(beta.hat = beta.hat/X.sd, beta.db=beta.db/X.sd[inf.coord], beta.db.sd=beta.db.sd/X.sd[inf.coord], tr.siga.inv=tr.siga.inv))
  }
}


#######################
## main function for variance component estimation
#######################

VC_estimation <- function(X, y, Z, grp, a=1, beta.hat){
  
  N=length(y)
  n=length(unique(grp))
  q=ncol(Z)
  p=ncol(X)
  m = max(table(grp))
  
  ### preprocessing
  y <- y - mean(y)
  
  X <- scale (X)
  X.sd <- attributes(X)$`scaled:scale` ## we need to divide these sd of X from the coefficients
  Z <- scale (Z)
  Z.sd <- attributes(Z)$`scaled:scale`

  lam.sig.seq<-c(50,10,exp(seq(log(5),log(0.1), length.out = 20)))*q*sqrt(log(q)/n)/ m 

  
  ## compute VC
  res = lapply(1:n, function(i){

      cur.mem=which(grp==i)
      mi=length(cur.mem)
      Zi = as.matrix(Z[cur.mem,])
      ri = y[cur.mem] - as.matrix(X[cur.mem,]) %*% beta.hat
      index.mat.tmp = lower.tri(diag(1, mi), diag=F) ## extract only off-diagonal entries
      
      A_l = lapply(1:q, function(l) Zi[,l, drop=F] %*% t(Zi[,l, drop=F])) ## these are coefficient matrices for the VC parameters
      
      y_vc = (ri%*%t(ri))[index.mat.tmp]
      X_vc = do.call(cbind, lapply(A_l, function(A)A[index.mat.tmp]))
      return(list(y_vc, X_vc))
    })
  
  y_vc = do.call(c, sapply(res, `[`, 1))  
  X_vc = do.call(rbind, sapply(res, `[`, 2))  
  
  
  
  vc.cv.init = cv.glmnet(X_vc, y_vc, lam = lam.sig.seq)
  ## double check if the tuning parameter range is reasonable
  if(which.lambda(vc.cv.init$lambda.min, lam.sig.seq) == length(lam.sig.seq)){
    lam.sig.seq.tmp = c(exp(seq(log(10),log(1e-4), length.out = 20)))*q*sqrt(log(q)/n)/ m 
    vc.cv.init<-cv.glmnet(X_vc, y_vc, lam = lam.sig.seq.tmp)
    if(which.lambda(vc.cv.init$lambda.min, lam.sig.seq.tmp) ==  length(lam.sig.seq.tmp)){
      warning(paste0('lambda.sig selected on boundary: ', which.lambda(vc.cv.init$lambda.min, lam.sig.seq.tmp), ' selected' ))
    }
  }
  if(which.lambda(vc.cv.init$lambda.min, lam.sig.seq) == 1){
    lam.sig.seq.tmp = c(exp(seq(log(1e4),log(10), length.out = 20)))*q*sqrt(log(q)/n)/ m 
    vc.cv.init<-cv.glmnet(X_vc, y_vc, lam = lam.sig.seq.tmp)
    if(which.lambda(vc.cv.init$lambda.min, lam.sig.seq.tmp) ==  1){
      warning(paste0('lambda.sig selected on boundary: ', which.lambda(vc.cv.init$lambda.min, lam.sig.seq.tmp), ' selected' ))
    }
  }
  
  
  vc.hat<-coef(vc.cv.init, s=vc.cv.init$lambda.min)[-1]
  
  sig2.eps.hat = 1/N*
    sum(sapply(1:n, function(i){
      
      cur.mem=which(grp==i)
      mi=length(cur.mem)
      Zi = as.matrix(Z[cur.mem,])
      ri = y[cur.mem] - as.matrix(X[cur.mem,]) %*% beta.hat
      
      return(sum(ri^2) - sum(diag(Zi %*% diag(vc.hat) %*% t(Zi))) )
    }))
  
  sig2.eps.hat = ifelse(sig2.eps.hat<0, 0, sig2.eps.hat)
  
  ## make estimates non-negative: truncate at 0
  vc.hat[vc.hat<0]<-0
    

    return( list(vc.hat = vc.hat/Z.sd^2 , sig2.eps.hat = sig2.eps.hat))
  
}






#########################
## other utility functions
#########################

## function to quickly generate a covariance matrix
generate_cov_matrix = function(type, param){
  ## type = 'diag': param is the diagonal vector
  ## type = 'sparse-random': param = c(matrix dimension, off-diag unif range, off-diag nonzero prob)
  ## type = 'AR1': param = c(matrix dimension, rho)
  ## type = 'const': param = c(matrix dimension, constant for off-diagonal)
  if(type == 'diag'){ 
    Sig = diag(param)
  }else if(type == 'sparse-random'){
    Sig = diag(0, param[1])
    counter=0
    while(any(eigen(Sig)$values<=0) & counter < 10){
      counter = counter+1
      Sig = matrix(runif(param[1]^2),param[1],param[1])
      Sig = 1*((Sig+t(Sig))/2 < param[3])
      Sig = Sig * matrix(runif(param[1]^2, min = -param[2], max=param[2]), param[1], param[1])
      Sig = (Sig + t(Sig))/2
      diag(Sig)<-1
    }
    
    if(any(eigen(Sig)$values<=0)) stop('Cov matrix not pd')
  }else if(type == 'AR1'){
    Sig = toeplitz(param[2]^(0:(param[1]-1)))
  }else if(type == 'const'){
    Sig = matrix(param[2], param[1], param[1])
    diag(Sig)<-1
  }else{
    stop('specify Cov type among diag, sparse-random, AR1, const')
  }
  
  return(Sig)
}


## function to select best value a as a tuning parameter
select_a<- function(X, y, Z, grp, a.seq=seq(0, 10, 0.5)){
  n<-length(unique(grp))
  n.tr=round(0.8*n) # four fold to train and one fold to test (not a cross validation); subset by whole clusters
  N.tr<- sum(as.numeric(grp) <= n.tr)
  best.err<- sum(y^2)
  best.a <- a.seq[1]
  for(i in 1:length(a.seq)){
    est.re<-Fixed_effects_est_inf(X[1:N.tr,], y[1:N.tr], Z[1:N.tr,], grp[1:N.tr], a=a.seq[i])
    pred.err<-sum((y[-(1:N.tr)]-X[-(1:N.tr),]%*%matrix(est.re$beta.hat, ncol=1))^2)
    if(pred.err<best.err){
      best.err<-pred.err
      best.a=a.seq[i]
    }
  }
  return(best.a=best.a)
}

## function to determine which tuning parameter lambda is selected (due to rounding error in cv.glmnet, the selected lambda is different from supplied lambda)
which.lambda = function(lambda, lambda_list){
  tmp = abs(lambda_list - lambda)
  which(tmp == min(tmp))
}
