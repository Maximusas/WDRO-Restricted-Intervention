####################################################################
# Working directory should be the root of the code folder
####################################################################

#plot_generator is a function that generates the plot for the Martingale estimator
#of the parent oracle and the Martingale estimator that includes a redundant covariate
#Input:
#seed: seed sets the seed to generate different graphs which are reproducible
library(ggplot2)
source("./src/martingaleWDRO.R") #Martingale estimator

plot_generator=function(seed=302){
  set.seed(seed)
  n=100
  ntest=1000
  d=10
  k=d+1
  B = matrix(0, k, k)
  #The coefficients of the covariates into the response
  B[k,1:10]=0
  B[k,9]=1
  B[k,10]=2
  I=diag(k)
  A=t(I-B)%*%(I-B)
  pa=which(B[k,-k]!=0)
  #Half of the data is from the interventional distribution
  frac=0.5
  #Use setting 4, where a shift is present in the test distribution
  X=matrix(rnorm(n * d,0,1), n, d)
  Xtest=matrix((rnorm(ntest * d,0,1)+rnorm(ntest*d,sqrt(0.5),1)), ntest, d) 
  X[ceiling(frac*n):n,9]=X[ceiling(frac*n):n,9]+rnorm(n-ceiling(frac*n)+1,0,5)
  
  #Normally distributed confounder
  H=rnorm(n,0,1)
  Htest=rnorm(ntest,0,1)
  #Add confounder
  X=X+H
  Xtest=Xtest+Htest
  
  #If propagation then can generate X via this loop:
  for(j in 1:d){
    X[,j]=B[j,1:d]%*%t(X)+X[,j]     
  }
  for(j in 1:d){
    Xtest[,j]=B[j,1:d]%*%t(Xtest)+Xtest[,j]    
  }
  
  #Generate response by 1*X9+2*X10 and a confounder
  y=as.numeric(B[k,-k]%*%t(X))+rnorm(n,0,1)+H
  ytest=as.numeric(B[k,-k]%*%t(Xtest))+rnorm(ntest,0,1)+Htest
  Xtest=sweep(Xtest,2,colMeans(X),FUN="-")
  
  #Standardise
  X=sweep(X,2,colMeans(X),FUN="-")
  ytest=ytest-mean(y)
  y=y-mean(y)

  B_2= matrix(0, k, k)
  B_22= matrix(0, k, k)
  causal=lm(y~X[,pa])
  #Construct B matrix with estimated causal coefficients
  B_2[k,pa]=causal$coefficients[-1]
  #Construct A that is used to solve the estimator
  A_2=t(I-B_2)%*%(I-B_2)
  #Add a redundant covariate to the parents, here covariate 1
  pa2=c(1,pa)
  B_22[k,pa2]=lm(y~X[,pa2])$coefficients[-1]
  A_22=t(I-B_22)%*%(I-B_22)
  
  #Estimate OLS
  ols=lm(y~X)
  beta_estimated_ols=ols$coefficients[-1]
  ypred_ols=Xtest%*%beta_estimated_ols
  MSE_ols=(mean((ypred_ols-ytest)^2))
  
  #Estimate causal estimator
  beta_estimated_causal=causal$coefficients[-1]
  ypred_causal=Xtest[,pa]%*%beta_estimated_causal
  MSE_causal=(mean((ypred_causal-ytest)^2))
  
  lambdas=seq(0.001,5000,by=0.4)
  
  #Save the MSE over different values of the radius (here named lambda)
  MSE_saver=rep(0,length(lambdas))
  MSE_saver2=rep(0,length(lambdas))
  coeffA=matrix(0,nrow=length(lambdas),ncol=d)
  i=1
  for(lambda in lambdas){
    beta_estimated_M2=tikhonov_constrained(X,y,A_2,lambda)
    coeffA[i,]=beta_estimated_M2
    ypred_M2=Xtest%*%beta_estimated_M2
    beta_estimated_M22=tikhonov_constrained(X,y,A_22,lambda)
    ypred_M22=Xtest%*%beta_estimated_M22
    MSE_saver[i]=mean((ypred_M2-ytest)^2)
    MSE_saver2[i]=mean((ypred_M22-ytest)^2)
    i=i+1
  }
  
  #Put results into dataframe to be able to use ggplot
  result=data.frame(Epsilon=lambdas,Test_MSE=MSE_saver,Test_MSE2=MSE_saver2)
  #Make the plot
  p=ggplot(result, aes(x = Epsilon)) +
    geom_line(aes(y=Test_MSE, color = "Martingale MSE "), linewidth = 1) +
    geom_line(aes(y=Test_MSE2, color = "Martingale MSE2 "), linewidth = 1) +
    geom_hline(aes(yintercept = MSE_causal, color = "Causal MSE "), linetype = "dashed", linewidth = 1) +
    geom_hline(aes(yintercept = MSE_ols, color = "OLS MSE "), linetype = "dashed", linewidth = 1) +
    xlab(expression(epsilon)) +
    ylab("Test MSE") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_manual(
      name = "Method",
      values = c("Causal MSE " = "red", "Martingale MSE " = "black", "Martingale MSE2 " = "blue", "OLS MSE "="green"),
      labels = c("Causal estimator","Martingale parent oracle", "Martingale redundant covariate","OLS")
    ) +
    labs(linetype = "Method") +
    guides(
      color=guide_legend(override.aes=list(linetype=c("dashed","solid","solid","dashed"), linewidth=c(0.5,1,1,0.5)))
    ) +
    theme(legend.position=c(0.836,0.196))
  
  return(p)
}

#Both outperform
plot_generator(32)
#Redundant marginally outperforms Oracle:
plot_generator(35001)
plot_generator(35013)
#Redundant covariate has no benefit:
plot_generator(330)
#Redundant covariate outperforms OLS in small region
plot_generator(3409)
#No benefits
plot_generator(37003)