####################################################################
# Working directory should be the root of the code folder
####################################################################
source("./src/ph_anchor.R") #Anchor estimators

simulation_delta=function(nsim,verbose=FALSE){
  set.seed(421)
  ntest=100
  n=100
  d=10
  Id=diag(1,n)
  gamma=5
  deltas=seq(0.05,2.25,0.1)
  lambdas=seq(0.001,0.2,0.005)
  MSE_test_anch_vec=rep(0,nsim)
  MSE_test_anch_ph_vec=matrix(0,nsim,length(deltas))
  epsilon=0.1

  MSE_matrix=rep(0,length(deltas)+1)
  Id=diag(1,n)

  for(sim in 1:nsim){
    #Simulation generation, see paper for details about settings
    A1=rnorm(n,0,1)
    A2=rnorm(n,0,1)
    Atest1=rnorm(ntest,5,1)
    Atest2=rnorm(ntest,5,1)
    M1=rnorm(d)
    M1=(M1-mean(M1))/sd(M1)
    M2=rnorm(d)
    M2=(M2-mean(M2))/sd(M2)
    X=matrix(rnorm(n*d),nrow=n)
    Xtest=matrix(rnorm(ntest*d),nrow=ntest)
    for(i in 1:d){
      X[,i]=X[,i]+M1[i]*A1+M2[i]*A2
      Xtest[,i]=Xtest[,i]+M1[i]*Atest1+M2[i]*Atest2
      Xtest[,i]=Xtest[,i]-mean(X[,i])
      X[,i]=X[,i]-mean(X[,i])
    }
    y=X[,1]+2*X[,2]+rnorm(n)
    ytest=Xtest[,1]+2*Xtest[,2]+rnorm(ntest)
    
    
    #Standardising
    Xtest=sweep(Xtest,2,colMeans(X),FUN="-")
    X=sweep(X,2,colMeans(X),FUN="-")
    ytest=ytest-mean(y)
    y=y-mean(y)
    
    #Transforming the data via the Anchor Transformation
    A=cbind(A1,A2)
    projA=A%*%solve(t(A)%*%A)%*%t(A)
    yanch=(Id+(sqrt(gamma)-1)*projA)%*%y
    Xanch=(Id+(sqrt(gamma)-1)*projA)%*%X
    yanch=as.numeric(unlist(yanch))
    
    #Computing the estimators and the MSE on the test set
    #Anchor
    anchor = Anchor(yanch,Xanch,lambdas = lambdas)
    a0_anch = anchor$a0
    coef_anch = anchor$beta
    ypred_anch = (Xtest %*% coef_anch + a0_anch)
    MSE_test_anch = mean((ypred_anch - ytest)^2)
    #print(MSE_test_anch)
    
    for(i in 1:length(deltas)){
      #PH Anchor without RWP
      ph_anchor=Pseudo_Huber_Anchor(yanch,Xanch,delta=deltas[i],epsilon)
      coef_anch_pha=ph_anchor$coef
      a0_anch_pha=ph_anchor$a0
      ypred_anch_pha=((Xtest%*%coef_anch_pha)+a0_anch_pha)
      MSE_test_anch_pha=mean((ypred_anch_pha-ytest)^2)
      
      
      #Recording the MSE
      MSE_test_anch_ph_vec[sim,i]=MSE_test_anch_pha
      #print(MSE_test_anch_ph_vec[sim,i])
    }
    MSE_test_anch_vec[sim]=MSE_test_anch
    if(verbose==TRUE){
    cat("Iteration:", sim, "\n") #Iteration number progress
    }
    
  }
  MSE_matrix[1:length(deltas)]=apply(MSE_test_anch_ph_vec,2,mean)
  MSE_matrix[24]=mean(MSE_test_anch_vec)
  
  
  
  return(MSE_matrix)
}

delta_sim_matrix=simulation_delta(50)

#Plotting:
deltas=seq(0.05,2.25,0.1)
result=data.frame(Epsilon=deltas,Test_MSE=delta_sim_matrix[1:23])
p=ggplot(result, aes(x = Epsilon, y = Test_MSE)) +
  geom_line(aes(color = "PH-Anchor"), linewidth = 1) +
  geom_hline(aes(yintercept = delta_sim_matrix[24], color = "Anchor"), linetype = "dashed", size = 0.5) +
  xlab(expression(delta)) +
  ylab("Test MSE") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(
    name = "Method",
    values = c("Anchor" = "red", "PH-Anchor" = "black"),
    labels = c("Anchor", "PH-Anchor")
  ) +
  labs(linetype = "Method") +
  guides(
    color=guide_legend(override.aes=list(linetype=c("dashed","solid"),linewidth=c(0.02,1)))
  ) +
  theme(legend.position=c(0.816,0.816))
