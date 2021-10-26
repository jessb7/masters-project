data(loomis)
data=na.omit(loomis)

n=500
mu_beta=c(0, 0)
sigma_beta=1
beta0=3.5
beta1=.5

x2=runif(n)
#x2=data$travel
x2=transpose(x2)
rho=exp(beta0+beta1*x2)

w=rexp(n, rho)
x=rgeom(n, exp(-w))
w=transpose(w)
x=transpose(x)
y=x+1
#y=data$anvisits

nsave=9000
nburn=1000
ntot=nburn+nsave

nMH=2
nrep=1

wold=matrix(w, ncol=ntot, nrow=n)
told=matrix(exp(-w), ncol=ntot, nrow=n)
beta0old=matrix(beta0, ncol=1, nrow=nMH)
beta1old=matrix(beta1, ncol=1, nrow=nMH)
beta0new=matrix(beta0, ncol=nrep, nrow=ntot+1)
beta1new=matrix(beta1, ncol=nrep, nrow=ntot+1)
predictive=rep(0, n)

mean0=0
mean1=0
med0=0
med1=0
CIU0=0
CIL0=0
CIU1=0
CIL1=0
mse0=0
medse0=0
mse1=0
medse1=0

for (rep in 1:nrep){
  for (iter in 1:ntot){
    for (j in 1:n){
      told[j,iter]=rbeta(1, exp(beta0new[iter, 1]+beta1new[iter, 1]*x2[j,])+1, y)
      wold[j,iter]=-log(told[j,iter])
    }
    for (i in 2:nMH){
      mustar_beta=c(beta0old[i-1,1], beta1old[i-1,1])
      tau=0.1
      beta0star=rnorm(1,mustar_beta[1],tau)
      beta1star=rnorm(1,mustar_beta[2],tau)
      betastar=c(beta0star, beta1star)
      betastar=transpose(betastar)
      
      postar=matrix(0, ncol=1, nrow=n)
      for (kk in 1:n){
        postar[kk,1]=exp(beta0star+beta1star*x2[kk,]*wold[kk,iter])
      }
      sustar=-sum(postar)
      lognum2=sustar-log(2*pi)+n*beta0star+beta1star*sum(x2)-.5*(beta0star^2+beta1star^2)/sigma_beta
      
      poold=matrix(0, ncol=1, nrow=n)
      for (kk1 in 1:n){
        poold[kk1,1]=exp(beta0old[i-1,1]+beta1old[i-1,1]*x2[kk1,]*wold[kk1,iter])
      }
      suold=-sum(poold)
      logden2=suold-log(2*pi)+n*beta0old[i-1,1]+beta1old[i-1,1]*sum(x2)-.5*(beta0old[i-1,1]^2+beta1old[i-1,1]^2)/sigma_beta
      
      minbe=lognum2-logden2
      u=runif(1)
      
      if (log(u)<minbe){
        beta0old[i,1]=beta0star
        beta1old[i,1]=beta1star
      }else{
        beta0old[i,1]=beta0old[i-1,1]
        beta1old[i,1]=beta1old[i-1,1]
      }
    }
    beta0new[iter+1,rep]=beta0old[i,1]
    beta1new[iter+1,rep]=beta1old[i,1]

  }
  mean0= mean0+mean(beta0new[nburn:ntot,rep])
  mean1= mean1+mean(beta1new[nburn:ntot,rep])
  
  med0=med0+median(beta0new[nburn:ntot,rep])
  med1=med1+median(beta1new[nburn:ntot,rep])
  
  CIU0=CIU0+mean(beta0new[nburn:ntot,rep])+1.96*sqrt(var(beta0new[nburn:ntot,rep])/n)
  CIL0=CIL0+mean(beta0new[nburn:ntot,rep])-1.96*sqrt(var(beta0new[nburn:ntot,rep])/n)
  CIU1=CIU1+mean(beta1new[nburn:ntot,rep])+1.96*sqrt(var(beta1new[nburn:ntot,rep])/n)
  CIL1=CIL1+mean(beta1new[nburn:ntot,rep])-1.96*sqrt(var(beta1new[nburn:ntot,rep])/n)
  
  mse0=mse0+sqrt(mean((beta0new[nburn:ntot,rep]-beta0)^2))
  medse0=medse0+sqrt(median((beta0new[nburn:ntot,rep]-beta0)^2))
  
  mse1=mse1+sqrt(mean((beta1new[nburn:ntot,rep]-beta1)^2))
  medse1=medse1+sqrt(median((beta1new[nburn:ntot,rep]-beta1)^2))
}


for (loop in 1:500){
rho_new=exp(beta0new[loop+9500,1]+beta1new[loop+9500,1]*x2[loop])
w_new=rexp(1,rho_new)
predictive[loop]=1+rgeom(1,exp(-w_new))
}
hist(predictive, xlim=c(1,350))

for (loop in 1:500){
  rho_new=exp(beta0new[loop+9500,1]+beta1new[loop+9500,1]*x2[loop])
  w_new=rexp(1,rho_new)
  pre[loop]=1+rgeom(1,exp(-w_new))
}
hist(pre, xlim=c(1,350))

mean0/nrep
mean1/nrep
med0/nrep
med1/nrep
CIU0/nrep
CIL0/nrep
CIU1/nrep
CIL1/nrep
mse0/nrep
medse0/nrep
mse1/nrep
medse1/nrep


plot(1000:10000, beta0new[nburn:ntot,nrep], type='l')
plot(1000:10000, beta1new[nburn:ntot,nrep], type='l')
hist(beta0new[nburn:ntot,nrep], main=NULL)
hist(beta1new[nburn:ntot,nrep], main=NULL)
