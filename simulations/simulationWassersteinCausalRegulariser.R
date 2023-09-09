####################################################################
# Working directory should be the root of the code folder
####################################################################
source("./src/martingaleWDRO.R") #Martingale Estimator
source("./src/wasserstein_causal_regulariser.R") #Iterative scheme WDRORI Estimator

pkgs=c("MASS", "pcalg","ggplot2")
for(package in pkgs){
  if(!require(package,character.only = TRUE)) {install.packages(package)}
  library(package, character.only = TRUE)
}

#Function to generate environment data
env_generator=function(n,d,mu,Sigma,myDAG){
  mu_delta=runif(d,0,1)
  mu_delta=c(mu_delta,0)
  #Scale the mean to get heterogeneous
  mu_delta=sqrt(4) * mu_delta
  G=matrix(runif((d+1)^2)*2-1, ncol=d+1) 
  Sigma_delta=t(G) %*% G
  #Covariance with response is 0
  Sigma_delta[d+1,-(d+1)]=0
  Sigma_delta[-(d+1),d+1]=0
  #Scale the covariance to get heterogeneous
  Sigma_delta=4*Sigma_delta/norm(Sigma_delta,"2")
  eMat=mvrnorm(n,mu=mu_delta,Sigma_delta)
  d.CnormMat=rmvDAG(n,myDAG,errMat = eMat)
  y=d.CnormMat[,11]
  X=d.CnormMat[,-11]
  return(list(y=y,X=X))
}

#Function to generate test data that is mean-shifted 
shift_test_generator=function(n,d,mu,Sigma,myDAG,shift){
  #Shift the mean by amount shift
  mu_delta=runif(d,shift,shift+1)
  mu_delta=c(mu_delta,0)
  #Scale the mean to get heterogeneous
  mu_delta=sqrt(4) * mu_delta
  G=matrix(runif((d+1)^2)*2-1, ncol=d+1) 
  Sigma_delta=t(G) %*% G
  #Covariance with response is 0
  Sigma_delta[d+1,-(d+1)]=0
  Sigma_delta[-(d+1),d+1]=0
  #Scale the covariance to get heterogeneous
  Sigma_delta=4*Sigma_delta/norm(Sigma_delta,"2")
  eMat=mvrnorm(n,mu=mu_delta,Sigma_delta)
  d.CnormMat=rmvDAG(n,myDAG,errMat = eMat)
  y=d.CnormMat[,11]
  X=d.CnormMat[,-11]
  return(list(y=y,X=X))
}

n=100
d=10
test_env=20
n_env=10
k=d+1
I=diag(k)
MSE_causal=rep(0,test_env)
MSE_ols=rep(0,test_env)
MSE_pooled=rep(0,test_env)
MSE_wcr=rep(0,test_env)
MSE_cr=rep(0,test_env)
MSE_m=rep(0,test_env)
MSE_wcr_oracle=rep(0,test_env)
MSE_cr_oracle=rep(0,test_env)
nsim=50
maxmatrix=matrix(0,nsim,8)
meanmatrix=matrix(0,nsim,8)

maxshift=18
mean97=matrix(0,maxshift,8)
mean2=matrix(0,maxshift,8)
meanmean=matrix(0,maxshift,8)

max97=matrix(0,maxshift,8)
max2=matrix(0,maxshift,8)
maxmean=matrix(0,maxshift,8)


gammas_cr=c(seq(0.01,0.91,by=0.05),seq(0.92,0.99,by=0.01),seq(0.991,0.999,by=0.001),1)
set.seed(142)
myDAG=randomDAG(n = 11, prob= 0.2, lB = 0.1, uB = 1)
plot(myDAG)
for(shift in 1:maxshift){
  for(sim in 1:nsim){
    
    #Observational data
    G=matrix(runif((d+1)^2)*2-1, ncol=d+1) 
    Sigma=t(G) %*% G
    Sigma[d+1,-(d+1)]=0
    Sigma[-(d+1),d+1]=0
    Sigma=4*Sigma/norm(Sigma,"2")
    mu=rep(0,d+1)
    eMat=mvrnorm(n,mu=mu,Sigma)
    d.CnormMat=rmvDAG(n,myDAG,errMat = eMat)
    y=d.CnormMat[,11]
    X=d.CnormMat[,-11]
    
    #Generate data for environments
    X_e =list()
    y_e =list() 
    for(i in 1:n_env){
      env=env_generator(n,10,mu,Sigma,myDAG)
      X_e[[i]] = env$X
      y_e[[i]] = env$y
    }
    
    #Construct matrix for pooled OLS, pooled Causal, and pooled Martingale
    pooled_X=do.call(rbind,X_e)
    pooled_X=rbind(X,pooled_X)
    pooled_y=do.call(c,y_e)
    pooled_y=c(y,pooled_y)
    
    Xtest_e =list()
    ytest_e =list() 
    for(i in 1:test_env){
      env=shift_test_generator(n,10,mu,Sigma,myDAG,1.5*shift)
      Xtest_e[[i]] = env$X-colMeans(pooled_X)
      ytest_e[[i]] = env$y-mean(pooled_y)
    }
    X=X-colMeans(pooled_X)
    y=y-mean(pooled_y)
    
    pooled_X=sweep(pooled_X,2,colMeans(pooled_X),"-")
    pooled_y=pooled_y-mean(pooled_y)
    lm(y~X)$coefficients[-1]
    lm(pooled_y~pooled_X)$coefficients[-1]
    coef_ols=lm(y~X)$coefficients[-1]
    coef_pooled=lm(pooled_y~pooled_X)$coefficients[-1]
    
    #Adjacency matrix to find matrix B for Martingale estimator
    Adjacency=matrix(0,k,k)
    B_true=matrix(0,k,k)
    B=matrix(0,k,k)
    o=1
    #Find start and end node to fill adjacency matrix
    for (edge in names(myDAG@edgeData@data)) {
      splitted=strsplit(edge, split = "\\|")[[1]]
      start_node=as.numeric(splitted[1])
      end_node=as.numeric(splitted[2])
      Adjacency[end_node,start_node]=1
      B_true[end_node,start_node]=unlist(myDAG@edgeData@data[[o]])
      o=o+1
    }
    #Regress each variable on its parents to fill B
    for(o in 1:d){
      parent=which(Adjacency[o,]!=0)
      if(sum(Adjacency[o,])>0){
        B[o,parent]=lm(pooled_X[,o]~pooled_X[,parent])$coefficients[-1]
      }
    }
    #Regress on parent of response to fill last row of B
    parent=which(Adjacency[k,]!=0)
    causal=lm(pooled_y~1)
    if(sum(Adjacency[k,])>0){
      causal=lm(pooled_y~pooled_X[,parent])
      B[k,parent]=causal$coefficients[-1]
    }
    #Construct matrix A needed for the Martingale estimator
    A=t(I-B)%*%(I-B)
    
    #Estimate of Martingale estimator
    beta_m=tikhonov_constrained(pooled_X,pooled_y,A,2000)
    #CV on 7 environments, find best value on gamma on left out environments
    samp=sample(n_env,7)
    
    X_etrain=X_e[samp]
    X_eval=X_e[-samp]
    y_etrain=y_e[samp]
    y_eval=y_e[-samp]
    
    
    pooledval_X=do.call(rbind,X_eval)
    pooledval_y=do.call(c,y_eval)
    
    #The value is monotonically increasing beyond a certain value, if value found do a finer search
    low=Inf
    #First search by 5, if the loss increases then the lowest value is around the one with the lowest value
    for(gamma in seq(0,100,by=5)){
      wcr=wasserstein_causal_regulariser(X,y,X_etrain,y_etrain,gamma=gamma,n_iter=50)
      nextvalue=mean((pooledval_X%*%wcr-pooledval_y)^2)
      #print(nextvalue)
      if(nextvalue<=low){
        gamma_tracker=gamma
        low=nextvalue
      }
    }  
    
    low=Inf
    #Finer search
    for(gamma2 in seq(max(gamma_tracker-5,0),gamma_tracker+5,by=1)){
      wcr=wasserstein_causal_regulariser(X,y,X_etrain,y_etrain,gamma=gamma2,n_iter=50)
      nextvalue=mean((pooledval_X%*%wcr-pooledval_y)^2)
      #print(nextvalue)
      if(nextvalue<=low){
        gamma_tracker=gamma2
        low=nextvalue
      }
    }
    
    low=Inf
    for(gamma3 in seq(max(gamma_tracker-1,0),gamma_tracker+1,by=0.1)){
      wcr=wasserstein_causal_regulariser(X,y,X_etrain,y_etrain,gamma=gamma3,n_iter=50)
      nextvalue=mean((pooledval_X%*%wcr-pooledval_y)^2)
      #print(nextvalue)
      if(nextvalue<=low){
        gamma_tracker=gamma3
        low=nextvalue
      }
    }
    
    low=Inf
    #Finest search
    for(gamma4 in seq(max(gamma_tracker-0.1,0.001),gamma_tracker+0.1,by=0.01)){
      wcr=wasserstein_causal_regulariser(X,y,X_etrain,y_etrain,gamma=gamma4,n_iter=50)
      nextvalue=mean((pooledval_X%*%wcr-pooledval_y)^2)
      if(nextvalue<=low){
        gamma_tracker=gamma4
        low=nextvalue
      }
    }
    #Use gamma_tracker to estimate wcr
    beta_wcr=wasserstein_causal_regulariser(X,y,X_e,y_e,gamma=gamma_tracker,n_iter=50)
    
    low=Inf
    #Similar approach for causal regularisation but the value is never greater than 1
    for(gamma in seq(0,1,by=0.1)){
      cr=causal_regulariser(X,y,X_etrain,y_etrain,gamma=gamma)
      nextvalue=mean((pooledval_X%*%cr-pooledval_y)^2)
      if(nextvalue<=low){
        gamma_tracker=gamma
        low=nextvalue
      }
    }
    
    low=Inf
    for(gamma2 in seq(max(gamma_tracker-0.1,0),gamma_tracker+0.1,by=0.01)){
      cr=causal_regulariser(X,y,X_etrain,y_etrain,gamma=gamma2)
      nextvalue=mean((pooledval_X%*%cr-pooledval_y)^2)
      if(nextvalue<=low){
        gamma_tracker=gamma2
        low=nextvalue
      }
    }
    
    low=Inf
    for(gamma3 in seq(max(gamma_tracker-0.01,0),gamma_tracker+0.01,by=0.001)){
      cr=causal_regulariser(X,y,X_etrain,y_etrain,gamma=gamma3)
      nextvalue=mean((pooledval_X%*%cr-pooledval_y)^2)
      if(nextvalue<=low){
        gamma_tracker=gamma3
        low=nextvalue
      }
    }
    
    beta_cr=causal_regulariser(X,y,X_e,y_e,gamma=gamma_tracker)
    
    min_value_cr=Inf
    for(gamma in gammas_cr){
      beta_cr_oracle=causal_regulariser(X,y,X_e,y_e,gamma=gamma)
      for(i in 1:test_env){
        MSE_cr_oracle[i]=mean((Xtest_e[[i]]%*%beta_cr_oracle-ytest_e[[i]])^2)
        
      }
      if(mean(MSE_cr_oracle)<min_value_cr){
        min_gamma=gamma
        min_value_cr=mean(MSE_cr_oracle)
      }
    }
    beta_cr_oracle=causal_regulariser(X,y,X_e,y_e,gamma=min_gamma)
    
    #Do the same search but now for the oracle, look at the test MSE
    low=Inf
    for(gamma in seq(0,100,by=5)){
      wcr=wasserstein_causal_regulariser(X,y,X_e,y_e,gamma=gamma,n_iter=50)
      nextvalue=0
      for(i in 1:test_env){
        nextvalue=nextvalue+mean((Xtest_e[[i]]%*%wcr-ytest_e[[i]])^2)
      }
      if(nextvalue<=low){
        gamma_tracker=gamma
        low=nextvalue
      }
    }  
    low=Inf
    for(gamma2 in seq(max(gamma_tracker-5,0),gamma_tracker+5,by=1)){
      wcr=wasserstein_causal_regulariser(X,y,X_e,y_e,gamma=gamma2,n_iter=50)
      nextvalue=0
      for(i in 1:test_env){
        nextvalue=nextvalue+mean((Xtest_e[[i]]%*%wcr-ytest_e[[i]])^2)
      }
      if(nextvalue<=low){
        gamma_tracker=gamma2
        low=nextvalue
      }
    }  
    
    low=Inf
    for(gamma3 in seq(max(gamma_tracker-1,0),gamma_tracker+1,by=0.1)){
      wcr=wasserstein_causal_regulariser(X,y,X_e,y_e,gamma=gamma3,n_iter=50)  
      nextvalue=0
      for(i in 1:test_env){
        nextvalue=nextvalue+mean((Xtest_e[[i]]%*%wcr-ytest_e[[i]])^2)
      }
      #print(nextvalue)
      if(nextvalue<=low){
        gamma_tracker=gamma3
        low=nextvalue
      }
    }
    
    low=Inf
    for(gamma4 in seq(max(gamma_tracker-0.1,0.001),gamma_tracker+0.1,by=0.01)){
      wcr=wasserstein_causal_regulariser(X,y,X_e,y_e,gamma=gamma4,n_iter=50)
      nextvalue=0
      for(i in 1:test_env){
        nextvalue=nextvalue+mean((Xtest_e[[i]]%*%wcr-ytest_e[[i]])^2)
      }
      if(nextvalue<=low){
        gamma_tracker=gamma4
        low=nextvalue
      }
    }
    beta_wcr_oracle=wasserstein_causal_regulariser(X,y,X_e,y_e,gamma=gamma_tracker,n_iter=50)
    
    #For each test environment find the MSE
    for(i in 1:test_env){
      
      MSE_cr_oracle[i]=mean((Xtest_e[[i]]%*%beta_cr_oracle-ytest_e[[i]])^2)
      MSE_wcr_oracle[i]=mean((Xtest_e[[i]]%*%beta_wcr_oracle-ytest_e[[i]])^2)
      
      MSE_cr[i]=mean((Xtest_e[[i]]%*%beta_cr-ytest_e[[i]])^2)
      MSE_wcr[i]=mean((Xtest_e[[i]]%*%beta_wcr-ytest_e[[i]])^2)
      MSE_m[i]=mean((Xtest_e[[i]]%*%beta_m-ytest_e[[i]])^2)
      MSE_pooled[i]=mean((Xtest_e[[i]]%*%coef_pooled+lm(pooled_y~pooled_X)$coefficients[1]-ytest_e[[i]])^2)
      MSE_ols[i]=mean((Xtest_e[[i]]%*%coef_ols+lm(y~X)$coefficients[1]-ytest_e[[i]])^2)
      if(length(parent)>=1){
        MSE_causal[i]=mean((as.matrix(Xtest_e[[i]][,parent])%*%causal$coefficients[-1]+causal$coefficients[1]-ytest_e[[i]])^2)
      }
      else{MSE_causal[i]=mean((causal$coefficients[1]-ytest_e[[i]])^2)}
    }
    #Find the mean loss for each test environment
    meanmatrix[sim,1]=mean(MSE_ols)
    meanmatrix[sim,2]=mean(MSE_pooled)
    meanmatrix[sim,3]=mean(MSE_causal)
    meanmatrix[sim,4]=mean(MSE_m)
    meanmatrix[sim,5]=mean(MSE_wcr)
    meanmatrix[sim,6]=mean(MSE_cr)
    meanmatrix[sim,7]=mean(MSE_wcr_oracle)
    meanmatrix[sim,8]=mean(MSE_cr_oracle)
    
    #The worst-case test MSE
    maxmatrix[sim,1]=max(MSE_ols)
    maxmatrix[sim,2]=max(MSE_pooled)
    maxmatrix[sim,3]=max(MSE_causal)
    maxmatrix[sim,4]=max(MSE_m)
    maxmatrix[sim,5]=max(MSE_wcr)
    maxmatrix[sim,6]=max(MSE_cr)
    maxmatrix[sim,7]=max(MSE_wcr_oracle)
    maxmatrix[sim,8]=max(MSE_cr_oracle)
    print(sim)
  }
  print(shift)
  #Capture 97.5% quantile of the mean test MSE over the repetitions 
  mean97[shift,1]=quantile(meanmatrix[,1],0.975)
  mean97[shift,2]=quantile(meanmatrix[,2],0.975)
  mean97[shift,3]=quantile(meanmatrix[,3],0.975)
  mean97[shift,4]=quantile(meanmatrix[,4],0.975)
  mean97[shift,5]=quantile(meanmatrix[,5],0.975)
  mean97[shift,6]=quantile(meanmatrix[,6],0.975)
  mean97[shift,7]=quantile(meanmatrix[,7],0.975)
  mean97[shift,8]=quantile(meanmatrix[,8],0.975)
  
  #Capture 2.5% quantile of the mean test MSE over the repetitions 
  mean2[shift,1]=quantile(meanmatrix[,1],0.025)
  mean2[shift,2]=quantile(meanmatrix[,2],0.025)
  mean2[shift,3]=quantile(meanmatrix[,3],0.025)
  mean2[shift,4]=quantile(meanmatrix[,4],0.025)
  mean2[shift,5]=quantile(meanmatrix[,5],0.025)
  mean2[shift,6]=quantile(meanmatrix[,6],0.025)
  mean2[shift,7]=quantile(meanmatrix[,7],0.025)
  mean2[shift,8]=quantile(meanmatrix[,8],0.025)
  
  meanmean[shift,1]=mean(meanmatrix[,1])
  meanmean[shift,2]=mean(meanmatrix[,2])
  meanmean[shift,3]=mean(meanmatrix[,3])
  meanmean[shift,4]=mean(meanmatrix[,4])
  meanmean[shift,5]=mean(meanmatrix[,5])
  meanmean[shift,6]=mean(meanmatrix[,6])
  meanmean[shift,7]=mean(meanmatrix[,7])
  meanmean[shift,8]=mean(meanmatrix[,8])
  
  #Capture 97.5% quantile of the max test MSE over the repetitions 
  max97[shift,1]=quantile(maxmatrix[,1],0.975)
  max97[shift,2]=quantile(maxmatrix[,2],0.975)
  max97[shift,3]=quantile(maxmatrix[,3],0.975)
  max97[shift,4]=quantile(maxmatrix[,4],0.975)
  max97[shift,5]=quantile(maxmatrix[,5],0.975)
  max97[shift,6]=quantile(maxmatrix[,6],0.975)
  max97[shift,7]=quantile(maxmatrix[,7],0.975)
  max97[shift,8]=quantile(maxmatrix[,8],0.975)
  
  #Capture 2.5% quantile of the max test MSE over the repetitions 
  max2[shift,1]=quantile(maxmatrix[,1],0.025)
  max2[shift,2]=quantile(maxmatrix[,2],0.025)
  max2[shift,3]=quantile(maxmatrix[,3],0.025)
  max2[shift,4]=quantile(maxmatrix[,4],0.025)
  max2[shift,5]=quantile(maxmatrix[,5],0.025)
  max2[shift,6]=quantile(maxmatrix[,6],0.025)
  max2[shift,7]=quantile(maxmatrix[,7],0.025)
  max2[shift,8]=quantile(maxmatrix[,8],0.025)
  
  maxmean[shift,1]=mean(maxmatrix[,1])
  maxmean[shift,2]=mean(maxmatrix[,2])
  maxmean[shift,3]=mean(maxmatrix[,3])
  maxmean[shift,4]=mean(maxmatrix[,4])
  maxmean[shift,5]=mean(maxmatrix[,5])
  maxmean[shift,6]=mean(maxmatrix[,6])
  maxmean[shift,7]=mean(maxmatrix[,7])
  maxmean[shift,8]=mean(maxmatrix[,8])
}

meanmean
maxmean

for(b in 1:7){
  loess_fit=loess(smoothedmatrix[,b] ~ shifts,span=1)
  smoothedmatrix[,b]=predict(loess_fit)
}
result=as.data.frame(smoothedmatrix)  # Convert the matrix to a data frame
result$shifter=shifts
colname(result$meanadjusted)="V5"
ggplot(result, aes(x = shifter)) +
  geom_line(aes(y = V1, color = "Pooled OLS")) +
  geom_line(aes(y = V2, color = "Pooled Causal")) +
  geom_line(aes(y = V3, color = "Martingale")) +
  geom_line(aes(y = V4, color = "WCR")) +
  geom_line(aes(y = meanadjusted, color = "CR")) +
  geom_line(aes(y = V6, color = "Oracle WCR")) +
  geom_line(aes(y = V7, color = "Oracle CR")) +
  theme_minimal() +
  xlab("Perturbation strength") +
  ylab("Test MSE") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(
    name = "Method",
    values = c(      "Pooled OLS" = "green", 
                     "Pooled Causal" = "magenta", 
                     "Martingale" = "black", 
                     "WCR" = "#8B9400",
                     "CR" = "#00BFFF",
                     "Oracle WCR" = "#8B0000",
                     "Oracle CR" = "#4169E1" 
    ),
    breaks = c("Pooled OLS", "Pooled Causal", "Martingale", "WCR","CR","Oracle WCR","Oracle CR"),
    labels = c("Pooled OLS", "Pooled Causal", "Martingale", "WCR","CR","Oracle WCR","Oracle CR")
  ) +
  #scale_colour_discrete(labels = c("Pooled OLS", "Pooled Causal", "Martingale", "WCR","CR","Oracle WCR","Oracle CR"))+
  labs(linetype = "Method")+
  theme(legend.position=c(0.1,0.8))