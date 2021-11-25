####################################
## simulation of true LMM setting
## explore fixed effect estimation and inference

library(MASS)
library(glmnet)
library(scalreg)
library(tictoc)
library(MLmetrics)

#setwd('\\\\biostat-fs2\\users\\yuek\\Desktop\\LMM_testing')
setwd('~/Desktop/LMM_testing')

source('lib.R')
source('VC_est_func.R')
source('fix_est_func.R')


###########Experiments#####
Niter=200

m = 15
n = 30
p = 20
q = 20

m_list = m
n_list = n
p_list = p
q_list = q





inputs = commandArgs(T)
p_list = as.integer(inputs[1])
q_list = as.integer(inputs[2])
m_list = as.integer(inputs[3])
n_list = as.integer(inputs[4])

filename = paste0('res_true_LMM/p', p_list, 'q', q_list,'m', m_list, 'n', n_list, 'Niter', Niter, '.RData')
print(filename)


if(file.exists(filename)){
  stop('results obtained')
}

if(p_list < q_list){
  stop('need p >=q')
}

set.seed(234)
## the position and values of non-zero betas
beta_supp= c(1,2,6,7,9)
beta_nonzero = c(1,0.5,0.2,0.1,0.05)

# ## eta is for the random effect cov diagonal; not longer than p
eta_supp = c(1,4,7,9,10,12,16,20)
eta_nonzero = c(2, 2, 0.1, 0.1, 4, 0.1, 2, 0.1)
Psi.cov.type = 'diag' ## this give a sparse diagonal cov

## sigma.e is for error sd 
sigma.e = 1
Ri.cov.type = 'diag'

## Sig.X type
Sig.X.type = 'sparse-random'


inf.coord<-c(beta_supp, # for 1: strong signal+mid var; 2: strong sig+0 var; 6: mid sig+0 var; 7: weak sig+small var; 9: very weak sig+small var
             10, # 0 sig+large var 
             11, # 0 sig+0 var
             12  # 0 sig+small var
             )###coordinates to compute p-values
a.seq=seq(1, 10, length.out=10) ## a cannot be zero here in our algorithm



ci.re<-power.re<-beta.est.re<-vc.re <- replicate(length(p_list), list(), simplify = F)
vc.selection <- replicate(length(p_list), list(), simplify = F)



## try to simulate with X=Z, and the precision matrix of X being sparse
counter=0
for(p in p_list){
  
  ## define the values for beta
  beta<-rep(0,p)
  beta[beta_supp]=beta_nonzero

  ## define the covariance matrix for X (and Z)
  set.seed(11)
  Sig.x = generate_cov_matrix(type = Sig.X.type, param = c(p, 0.5, 0.2))
  
  for(q in q_list){ ## size for random effect, we want it as large as p
    counter = counter+1
    ## define the diagonal values of Psi
    set.seed(11)
    eta = rep(0, q)
    eta[eta_supp]<-eta_nonzero
    
    for(n in n_list){
      for(m in m_list){ ## cluster size;
        
        ## create G_list to use Li's method for VC
        G.list = lapply(1:p, function(s){
          tmp = matrix(0, p, p)
          tmp[s,s]<-1
          tmp})
        
        ## create Ri cov for noise, and Psi cov for random effect
        Ri = generate_cov_matrix(type = Ri.cov.type, param = rep(sigma.e^2, m))                
        Psi = generate_cov_matrix(type= Psi.cov.type, param=eta)
        
         
          N = n*m        
          cat('\np', p, 'q', q, 'm', m, 'n', n, 'N', N, '\n')
          

          rej<-matrix(0,ncol=length(inf.coord), nrow=length(a.seq), dimnames = list(as.character(a.seq), as.character(inf.coord)))
          ci.cov <-matrix(0,ncol=length(inf.coord), nrow=length(a.seq), dimnames = list(as.character(a.seq), as.character(inf.coord)))
          tr.a.avg = matrix(0, ncol=length(a.seq), nrow=1, dimnames = list(NULL, as.character(a.seq)))
          fit.time <- matrix(0, ncol=length(a.seq), nrow=1, dimnames = list(NULL, as.character(a.seq)))
          a.choice = 0
          theta.time= matrix(0, ncol=length(a.seq), nrow=1, dimnames = list(NULL, as.character(a.seq)))
          # cur.est.err<- matrix(0, ncol=length(a.seq), nrow=1)
          # ci.sd <-matrix(0,ncol=length(inf.coord), nrow=length(a.seq))
          
          ## also store the empirical bias and sd for each beta
          tmp1 = replicate(length(a.seq), rep(0, p), simplify = F) # for beta est
          tmp2 = replicate(length(a.seq), rep(0, length(inf.coord)), simplify = F) # for beta.db and beta.db.sd
          tmp3 = replicate(length(a.seq), matrix(0, nrow=1, ncol=q+1, dimnames = list(NULL, c(1:q, 'e2'))), simplify = F) # for VC
          names(tmp1) <- names(tmp2) <- names(tmp3) <- as.character(a.seq)
 
          store_est = list(beta.hat= tmp1, 
                           beta.db = tmp2,
                           beta.db.sd = tmp2,
                           ## also store the estimated variance components
                           theta.hat = tmp3,
                           a.opt = 0)
          
          
          my.eval = list(rej = rej,
                         rej.opt = rej[1,, drop=F],
                         ci.cov = ci.cov,
                         ci.cov.opt = ci.cov[1,,drop=F],
                         tr.a.avg = tr.a.avg,
                         fit.time = fit.time,
                         a.choice = a.choice,
                         theta.time = theta.time,
                         store_est = store_est)
          
          Li.eval = my.eval
          
          for (iter in 1:Niter){
            
            cat('\niter=',iter, '\t')
            
            ## first simulate data
            set.seed(Niter*counter+iter)
            {
              X = mvrnorm(N, rep(0, p), Sig.x)
              Z = X[,1:q]
              
              gamma = mvrnorm(n,rep(0,q), Psi)
              
              epsilon = mvrnorm(n, rep(0, m), Ri)
              
              y = rep(0,N) 
              for (i in 1:n){
                cur.mem=((i-1)*m+1):(i*m)
                y[cur.mem]=X[cur.mem,]%*%beta+Z[cur.mem,]%*%gamma[i,]+epsilon[i,]
              }
            }
            
            #########################################################
            ## our method
            cat('\nour method..\n')
            ## treat 'a' as a tuning parameter and select the best 'a' value; compute results based on best a
            a.opt<-select_a(X=X, y=y, Z=Z, grp=factor(rep(1:n,each=m)), a.seq = a.seq)
            res.opta = Fixed_effects_est_inf(X=X, y=y, Z=Z, grp=factor(rep(1:n,each=m)), a=a.opt, inf.coord=inf.coord)
            res.vc.opta = VC_estimation(X=X, y=y, Z=Z, a=a.opt,grp=factor(rep(1:n,each=m)), beta.hat = res.opta$beta.hat)
            
            ## save best a for each iteration
            my.eval$a.choice<- my.eval$a.choice+a.opt
            my.eval$store_est$a.opt = c(my.eval$store_est$a.opt, a.opt)
            ## compute confidence interval coverage and power
            beta.db<- res.opta$beta.db
            beta.db.sd<-res.opta$beta.db.sd
            my.eval$rej.opt<-my.eval$rej.opt+(abs(beta.db/beta.db.sd)>=1.96) # power
            my.eval$ci.cov.opt<- my.eval$ci.cov.opt + 
              sapply(1:length(inf.coord), function(j) (beta[inf.coord[j]]<= beta.db[j]+1.96 * beta.db.sd[j] & beta[inf.coord[j]]>= beta.db[j]- 1.96 * beta.db.sd[j])) # CI coverage
            
            
            ##also study the effect of a on the outcomes
            for(a in a.seq){
              cat('a=', a, '.. ')
              
              tic()
              
              fixed_res_a_i = Fixed_effects_est_inf(X=X, y=y, Z=Z, grp=factor(rep(1:n,each=m)), 
                                                    a=a, inf.coord=inf.coord)
              t1<-toc(quiet=T)
              
              ## store estimates
              my.eval$store_est$beta.hat[[as.character(a)]] <- rbind(my.eval$store_est$beta.hat[[as.character(a)]], fixed_res_a_i$beta.hat)
              my.eval$store_est$beta.db[[as.character(a)]] <- rbind(my.eval$store_est$beta.db[[as.character(a)]], fixed_res_a_i$beta.db)
              my.eval$store_est$beta.db.sd[[as.character(a)]] <- rbind(my.eval$store_est$beta.db.sd[[as.character(a)]], fixed_res_a_i$beta.db.sd)
              
              ## compute confidence interval coverage and power
              beta.db<- fixed_res_a_i$beta.db
              beta.db.sd<-fixed_res_a_i$beta.db.sd
              
              my.eval$rej[as.character(a),]<-my.eval$rej[as.character(a),]+(abs(beta.db/beta.db.sd)>=1.96) # power
              my.eval$ci.cov[as.character(a),]<- my.eval$ci.cov[as.character(a),] + 
                sapply(1:length(inf.coord), function(j) (beta[inf.coord[j]]<= beta.db[j]+1.96 * beta.db.sd[j] & beta[inf.coord[j]]>= beta.db[j]- 1.96 * beta.db.sd[j])) # CI coverage
              
              ## store other relevant quantities
              my.eval$tr.a.avg[1,as.character(a)] <- my.eval$tr.a.avg[1,as.character(a)] + fixed_res_a_i$tr.siga.inv
              my.eval$fit.time[1,as.character(a)]<- my.eval$fit.time[1,as.character(a)] + (t1$toc-t1$tic)
              
              
              ##################
              ## then we estimate the variance components
              tic()
              vc_res_a_i = VC_estimation(X=X, y=y, Z=Z, a=a, grp=factor(rep(1:n,each=m)), beta.hat = res.opta$beta.hat)
              t2 = toc(quiet=T)
              my.eval$theta.time[1,as.character(a)] <- my.eval$theta.time[1,as.character(a)] + t2$toc-t2$tic
              my.eval$store_est$theta.hat[[as.character(a)]] <- rbind(my.eval$store_est$theta.hat[[as.character(a)]], c(vc_res_a_i$vc.hat, vc_res_a_i$sig2.eps.hat))
              
            }
            
            ##########################################################
            ## compare to Li's method  
            if(T){
              cat("\nLi's method...\n")
              ## treat 'a' as a tuning parameter and select the best 'a' value; compute results based on best a
              a.opt = CV_fixed_effects(X=X, y=y, z=Z, grp=factor(rep(1:n,each=m)), a.seq=a.seq)
              res.opta = Fixed_effects_estimation(X=X, y=y, z=Z, grp=factor(rep(1:n,each=m)), a=a.opt, lm=F, inf.coord=inf.coord)
              res.vc.opta = NULL
              tryCatch({
                res.vc.opta = Varcomp.est.Li(r.hat= y-X%*%fixed_res_a_i$beta.hat, z=Z, G.list=G.list, a=a.opt, grp=factor(rep(1:n,each=m)))
              },
              error = function(e)NULL)
              
              ## save best a for each iteration
              Li.eval$a.choice<- Li.eval$a.choice+a.opt
              Li.eval$store_est$a.opt = c(Li.eval$store_est$a.opt, a.opt)
              ## compute confidence interval coverage and power
              beta.db<- res.opta$beta.db
              beta.db.sd<-res.opta$beta.db.sd
              Li.eval$rej.opt<-Li.eval$rej.opt+(abs(beta.db/beta.db.sd)>=1.96) # power
              Li.eval$ci.cov.opt<- Li.eval$ci.cov.opt + 
                sapply(1:length(inf.coord), function(j) (beta[inf.coord[j]]<= beta.db[j]+1.96 * beta.db.sd[j] & beta[inf.coord[j]]>= beta.db[j]- 1.96 * beta.db.sd[j])) # CI coverage
              
              
              ##also study the effect of a on the outcomes
              for(a in a.seq){
                cat('a=', a, '.. ')
                
                tic()
                
                fixed_res_a_i = Fixed_effects_estimation(X=X, y=y, z=Z, grp=factor(rep(1:n,each=m)), a=a, lm=F, inf.coord=inf.coord)
                t1<-toc(quiet=T)
                
                ## store estimates
                Li.eval$store_est$beta.hat[[as.character(a)]] <- rbind(Li.eval$store_est$beta.hat[[as.character(a)]], fixed_res_a_i$beta.hat)
                Li.eval$store_est$beta.db[[as.character(a)]] <- rbind(Li.eval$store_est$beta.db[[as.character(a)]], fixed_res_a_i$beta.db)
                Li.eval$store_est$beta.db.sd[[as.character(a)]] <- rbind(Li.eval$store_est$beta.db.sd[[as.character(a)]], fixed_res_a_i$beta.db.sd)
                
                ## compute confidence interval coverage and power
                beta.db<- fixed_res_a_i$beta.db
                beta.db.sd<-fixed_res_a_i$beta.db.sd
                
                Li.eval$rej[as.character(a),]<-Li.eval$rej[as.character(a),]+(abs(beta.db/beta.db.sd)>=1.96) # power
                Li.eval$ci.cov[as.character(a),]<- Li.eval$ci.cov[as.character(a),] + 
                  sapply(1:length(inf.coord), function(j) (beta[inf.coord[j]]<= beta.db[j]+1.96 * beta.db.sd[j] & beta[inf.coord[j]]>= beta.db[j]- 1.96 * beta.db.sd[j])) # CI coverage
                
                ## store other relevant quantities
                Li.eval$tr.a.avg[1,as.character(a)] <- Li.eval$tr.a.avg[1,as.character(a)] + fixed_res_a_i$tr.a
                Li.eval$fit.time[1,as.character(a)]<- Li.eval$fit.time[1,as.character(a)] + (t1$toc-t1$tic)
                
                
                ##################
                ## then we estimate the variance components
                tic()
                vc_res_a_i = NULL
                tryCatch({
                  vc_res_a_i = Varcomp.est.Li(r.hat= y-X%*%fixed_res_a_i$beta.hat, z=Z, G.list=G.list, a=a, grp=factor(rep(1:n,each=m)))
                },
                error = function(e)NULL)
                t2 = toc(quiet=T)
                Li.eval$theta.time[1,as.character(a)] <- Li.eval$theta.time[1,as.character(a)] + t2$toc-t2$tic
                Li.eval$store_est$theta.hat[[as.character(a)]] <- 
                  rbind(Li.eval$store_est$theta.hat[[as.character(a)]], c(vc_res_a_i$eta.hat, vc_res_a_i$sig.e.sq.hat^2))
                
              }
              
              
              
            }
            
            
          }
          
  
          ########################################################  
          ## summarize the results for this setting (m, n, p, q) combination
          ## separate for different (p,q)
          
          if(T){
            power.re[[counter]]<-rbind(power.re[[counter]], 
                                       ## results for each a
                                       data.frame(N=N, n=n, m=m, p=p, q=q, sigma.e=sigma.e, a=a.seq, 
                                                  rej = my.eval$rej/Niter, fittime = t(my.eval$fit.time/Niter), method = 'my'),
                                       ## result for optimum a
                                       data.frame(N=N, n=n, m=m, p=p, q=q, sigma.e=sigma.e, a=my.eval$a.choice/Niter, 
                                                  rej = my.eval$rej.opt/Niter, fittime = mean(my.eval$fit.time/Niter), method = 'my.opt'),
                                       
                                       ## results for Li's method
                                       data.frame(N=N, n=n, m=m, p=p, q=q, sigma.e=sigma.e, a=a.seq, 
                                                  rej = Li.eval$rej/Niter, fittime = t(Li.eval$fit.time/Niter), method = 'Li'),
                                       data.frame(N=N, n=n, m=m, p=p, q=q, sigma.e=sigma.e, a=Li.eval$a.choice/Niter, 
                                                  rej = Li.eval$rej.opt/Niter, fittime = mean(Li.eval$fit.time/Niter), method = 'Li.opt')
            )
            
            ci.re[[counter]] <- rbind(ci.re[[counter]], 
                                      data.frame(N=N, n=n,m=m, p=p, q=q, sigma.e=sigma.e, a = a.seq,
                                                 CIcov = my.eval$ci.cov/Niter, method = 'my'),
                                      data.frame(N=N, n=n,m=m, p=p, q=q, sigma.e=sigma.e, a = my.eval$a.choice/Niter,
                                                 CIcov = my.eval$ci.cov.opt/Niter, method= 'my.opt'),
                                      ## Li's results
                                      data.frame(N=N, n=n,m=m, p=p, q=q, sigma.e=sigma.e, a = a.seq,
                                                 CIcov = Li.eval$ci.cov/Niter, method = 'Li'),
                                      data.frame(N=N, n=n,m=m, p=p, q=q, sigma.e=sigma.e, a = Li.eval$a.choice/Niter,
                                                 CIcov = Li.eval$ci.cov.opt/Niter, method= 'Li.opt')
            )
            
            
            ## compute MSE, bias and sd for each beta for all the iterations
            
            get_beta_est_summary = function(res.eval){
              beta_hat_MSE = t(sapply(res.eval$store_est$beta.hat, function(x) rowMeans((t(x[-1,]) - beta)^2)))
              colnames(beta_hat_MSE) = paste0('MSE', 1:p)
              
              beta_hat_bias = t(sapply(res.eval$store_est$beta.hat, function(x) colMeans(x[-1,]) - beta))
              colnames(beta_hat_bias) = paste0('bias', 1:p)
              
              beta_hat_sd = t(sapply(res.eval$store_est$beta.hat, function(x) apply(x[-1,],2, sd) ))
              colnames(beta_hat_sd) = paste0('sd', 1:p)
              
              ## for beta.db, only computed for inf.coord
              beta_db_MSE = t(sapply(res.eval$store_est$beta.db, function(x) rowMeans((t(x[-1,]) - beta[inf.coord])^2)))
              colnames(beta_db_MSE) = paste0('dbMSE', 1:length(inf.coord))
              
              beta_db_bias = t(sapply(res.eval$store_est$beta.db, function(x) colMeans(x[-1,]) - beta[inf.coord]))
              colnames(beta_db_bias) = paste0('dbbias', 1:length(inf.coord))
              
              beta_db_sd = t(sapply(res.eval$store_est$beta.db, function(x) apply(x[-1,],2, sd) ))
              colnames(beta_db_sd) = paste0('dbsd', 1:length(inf.coord))
              
              ## also get the summary for optimum a in each iteration
              beta_hat_opt = sapply(1:Niter, function(i)res.eval$store_est$beta.hat[[as.character(res.eval$store_est$a.opt[-1][i])]][i+1,])
              beta_hat_opt_MSE = rowMeans(beta_hat_opt - beta)^2
              beta_hat_opt_bias = rowMeans(beta_hat_opt) - beta
              beta_hat_opt_sd = apply(beta_hat_opt, 1, sd)
              
              beta_db_opt = sapply(1:Niter, function(i)res.eval$store_est$beta.db[[as.character(res.eval$store_est$a.opt[-1][i])]][i+1,])
              beta_db_opt_MSE = rowMeans(beta_db_opt - beta[inf.coord])^2
              beta_db_opt_bias = rowMeans(beta_db_opt) - beta[inf.coord]
              beta_db_opt_sd = apply(beta_db_opt, 1, sd)
              
              return(list(beta_hat_MSE = beta_hat_MSE, beta_hat_bias=beta_hat_bias,beta_hat_sd= beta_hat_sd,
                          beta_hat_opt_MSE = beta_hat_opt_MSE, beta_hat_opt_bias=beta_hat_opt_bias,beta_hat_opt_sd= beta_hat_opt_sd,
                          
                          beta_db_MSE=beta_db_MSE, beta_db_bias=beta_db_bias,beta_db_sd=beta_db_sd,
                          beta_db_opt_MSE=beta_db_opt_MSE, beta_db_opt_bias=beta_db_opt_bias,beta_db_opt_sd=beta_db_opt_sd
                          
              ))
            }
            
            my.tmp = get_beta_est_summary(my.eval)
            Li.tmp = get_beta_est_summary(Li.eval)
            
            beta.est.re[[counter]] <- rbind(beta.est.re[[counter]], 
                                            cbind(data.frame(N=N, n=n,m=m, p=p, q=q, sigma.e=sigma.e, a = a.seq, 
                                                             total_hat_MSE = rowSums(my.tmp$beta_hat_MSE), 
                                                             total_db_MSE = rowSums(my.tmp$beta_db_MSE),
                                                             method = 'my'),
                                                  my.tmp$beta_hat_MSE,
                                                  my.tmp$beta_hat_bias,
                                                  my.tmp$beta_hat_sd),
                                            
                                            cbind(data.frame(N=N, n=n,m=m, p=p, q=q, sigma.e=sigma.e, a = my.eval$a.choice/Niter, 
                                                             total_hat_MSE = sum(my.tmp$beta_hat_opt_MSE), 
                                                             total_db_MSE = sum(my.tmp$beta_db_opt_MSE),
                                                             method = 'my.opt'),
                                                  matrix(my.tmp$beta_hat_opt_MSE, nrow=1, dimnames = list(NULL, c(paste0('MSE', 1:p)))),
                                                  matrix(my.tmp$beta_hat_opt_bias, nrow=1, dimnames = list(NULL, c(paste0('bias', 1:p)))),
                                                  matrix(my.tmp$beta_hat_opt_sd, nrow=1, dimnames = list(NULL, c(paste0('sd', 1:p))))
                                            ),
                                            
                                            cbind(data.frame(N=N, n=n,m=m, p=p, q=q, sigma.e=sigma.e, a = a.seq, 
                                                                 total_hat_MSE = rowSums(Li.tmp$beta_hat_MSE), 
                                                                 total_db_MSE = rowSums(Li.tmp$beta_db_MSE),
                                                                 method = 'Li'),
                                                  Li.tmp$beta_hat_MSE,
                                                  Li.tmp$beta_hat_bias,
                                                  Li.tmp$beta_hat_sd),
                                            cbind(data.frame(N=N, n=n,m=m, p=p, q=q, sigma.e=sigma.e, a = Li.eval$a.choice/Niter, 
                                                             total_hat_MSE = sum(Li.tmp$beta_hat_opt_MSE), 
                                                             total_db_MSE = sum(Li.tmp$beta_db_opt_MSE),
                                                             method = 'Li.opt'),
                                                  matrix(Li.tmp$beta_hat_opt_MSE, nrow=1, dimnames = list(NULL, c(paste0('MSE', 1:p)))),
                                                  matrix(Li.tmp$beta_hat_opt_bias, nrow=1, dimnames = list(NULL, c(paste0('bias', 1:p)))),
                                                  matrix(Li.tmp$beta_hat_opt_sd, nrow=1, dimnames = list(NULL, c(paste0('sd', 1:p))))
                                            )
            )
            
            get_vc_est_summary = function(res.eval){
              
              vc_MSE = t(sapply(res.eval$store_est$theta.hat, function(x) rowMeans((t(x[-1,]) - c(eta, sigma.e^2))^2)))
              colnames(vc_MSE) = paste0('MSE', c(1:q, 'sig2'))
              
              vc_bias = t(sapply(res.eval$store_est$theta.hat, function(x) colMeans(x[-1,]) - c(eta, sigma.e^2)))
              colnames(vc_bias) = paste0('bias', c(1:q, 'sig2'))
              
              vc_sd = t(sapply(res.eval$store_est$theta.hat, function(x) apply(x[-1,],2, sd) ))
              colnames(vc_sd) = paste0('sd', c(1:q, 'sig2'))
              
              
              
              ## also get the summary for optimum a in each iteration
              vc_opt = sapply(1:Niter, function(i)res.eval$store_est$theta.hat[[as.character(res.eval$store_est$a.opt[-1][i])]][i+1,])
              vc_opt_MSE = matrix(rowMeans(vc_opt - c(eta, sigma.e^2))^2, nrow=1, dimnames = list(NULL, c(paste0('MSE',c(1:q, 'sig2')))))
              vc_opt_bias = matrix(rowMeans(vc_opt) - c(eta, sigma.e^2), nrow=1, dimnames = list(NULL, c(paste0('bias',c(1:q, 'sig2')))))
              vc_opt_sd = matrix(apply(vc_opt, 1, sd), nrow=1, dimnames = list(NULL, c(paste0('sd',c(1:q, 'sig2')))))
              
              ## also get the selection consistency measure: 
              ## MCC as a summary of FP, TP, FN, TN, but some problem here to use MCC: for a single vc compute MCC over all replications, then truth is either all 0 or all 1; so TP or FN would be 0 
              ## compute the percentage of replications that select the correct non-zero VCs can also be problematic: seems rarely able to recover exactly the non-zero pattern
              ## F1-Score does not work either. 
              ## So just look at percentage of non-zero for each VC, so for true non-zeros will be TPR; for true zeros will be 1-TNR

              vc_nonzero_perc = t(
                sapply(res.eval$store_est$theta.hat, function(x){
                    colMeans(x[-1,1:q]>0)
              })
              )
              colnames(vc_nonzero_perc) <-paste0('node', 1:q)
                     
              vc_opt_nonzero_perc = rowMeans(vc_opt[1:q,]>0)
              names(vc_opt_nonzero_perc) <- paste0('node', 1:q)

              
              return(list(vc_MSE = vc_MSE, vc_bias=vc_bias,vc_sd= vc_sd,
                          vc_opt_MSE = vc_opt_MSE, vc_opt_bias=vc_opt_bias,vc_opt_sd= vc_opt_sd,
                          vc_nonzero_perc  = vc_nonzero_perc , vc_opt_nonzero_perc =vc_opt_nonzero_perc
              ))
            }
            
            my.tmp = get_vc_est_summary(my.eval)
            vc.re[[counter]] <- rbind(vc.re[[counter]], 
                                      cbind(data.frame(N=N, n=n,m=m, p=p, q=q, sigma.e=sigma.e, a = a.seq, 
                                                       method = 'my', 
                                                       time = t(my.eval$theta.time/Niter), 
                                                       total_vc_MSE = rowSums(my.tmp$vc_MSE[,-(q+1)])), ## this will only be eta, we separate sig.e2
                                            my.tmp$vc_MSE,
                                            my.tmp$vc_bias,
                                            my.tmp$vc_sd),
                                      cbind(data.frame(N=N, n=n,m=m, p=p, q=q, sigma.e=sigma.e, a = my.eval$a.choice/Niter, 
                                                       method = 'my.opt', 
                                                       time = mean(my.eval$theta.time/Niter), 
                                                       total_vc_MSE = sum(my.tmp$vc_opt_MSE[-(q+1)])),
                                            my.tmp$vc_opt_MSE,
                                            my.tmp$vc_opt_bias,
                                            my.tmp$vc_opt_sd)
            )
            
            vc.selection[[counter]] <- rbind(vc.selection[[counter]], 
                                             data.frame(N=N, n=n,m=m, p=p, q=q, sigma.e=sigma.e, a = a.seq, 
                                                        method = 'my',
                                                        my.tmp$vc_nonzero_perc ),
                                             data.frame(N=N, n=n,m=m, p=p, q=q, sigma.e=sigma.e, a = my.eval$a.choice/Niter,
                                                        method = 'my.opt',
                                                        t(as.matrix(my.tmp$vc_opt_nonzero_perc)))
                                             )
            
            if(m>q){
              Li.tmp = get_vc_est_summary(Li.eval)
              vc.re[[counter]] <- rbind(vc.re[[counter]],                          
                                        cbind(data.frame(N=N, n=n,m=m, p=p, q=q, sigma.e=sigma.e, a = a.seq, method = 'Li', 
                                                         time = t(Li.eval$theta.time/Niter), 
                                                         total_vc_MSE = rowSums(Li.tmp$vc_MSE[,-(q+1)])),
                                              Li.tmp$vc_MSE,
                                              Li.tmp$vc_bias,
                                              Li.tmp$vc_sd),
                                        cbind(data.frame(N=N, n=n,m=m, p=p, q=q, sigma.e=sigma.e, a = Li.eval$a.choice/Niter, method = 'Li.opt', 
                                                         time = mean(Li.eval$theta.time/Niter), 
                                                         total_vc_MSE = sum(Li.tmp$vc_opt_MSE[-(q+1)])),
                                              Li.tmp$vc_opt_MSE,
                                              Li.tmp$vc_opt_bias,
                                              Li.tmp$vc_opt_sd)
              )
              
              vc.selection[[counter]] <- rbind(vc.selection[[counter]], 
                                               data.frame(N=N, n=n,m=m, p=p, q=q, sigma.e=sigma.e, a = a.seq, 
                                                          method = 'Li',
                                                          Li.tmp$vc_nonzero_perc ),
                                               data.frame(N=N, n=n,m=m, p=p, q=q, sigma.e=sigma.e, a = Li.eval$a.choice/Niter,
                                                          method = 'Li.opt',
                                                          t(as.matrix(Li.tmp$vc_opt_nonzero_perc)))
              )
              
            }
          }
          
                             
        
      } 
    }
  }
}


save.image(filename)

