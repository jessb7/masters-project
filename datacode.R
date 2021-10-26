#This code simulates data from the Yule-Simon distribution. 
#The true parameter is rho, the number of iterations for 
#the Metropolis Hastings algorithm is m, and sample size is n.
library("optimbase") #include library for transpose function

m=5000 #number of MCMC samples
rep=20 #number of replications
n=200 #sample size of the data
b=c(-.5, 5) #true betas
x2=runif(n) #regressor values
rho=exp(b[1]+b[2]*x2)

#generate data from the Yule-Simon distribution
W=rexp(n,rho)
k=rgeom(n,exp(-W))
k=k+1

#initialize variables
mean1=0
mean2=0
med1=0
med2=0
CI1=0
CI11=0
CI2=0
CI21=0

#repeat 5000 iterations 20 times
for (iter in 1:rep) {
  
  #initialize variables
  beta=matrix(ncol=2, nrow=m, byrow=T)
  y=matrix(nrow=m,ncol=2)
  beta[1,]=c(b[1],b[2])
  y[1,]=c(1,1)
  A=vector(length=m)

  #interate through MCMC algorithm
  for (i in 2:m){
    t=rbeta(n, exp(beta[i-1, 1]+beta[i-1, 2]*x2)+1, k)
    w=-log(t)

    tau=0.01
    epsilon1=rnorm(1, 0, 1)
    epsilon2=rnorm(1, 0, 1)
    y[i,1]=beta[i-1, 1]+tau*epsilon1
    y[i,2]=beta[i-1, 2]+tau*epsilon2
    
    #calculate sums to use in the acceptance ratio
    sumy1=0
    sumb1=0
    sumy2=0
    sumb2=0
    for (j in 1:n) {
      sumy1=sumy1+exp(y[i, 1]+y[i, 2]*x2[j])*w[j]
      sumb1=sumb1+exp(beta[i-1, 1]+beta[i-1, 2]*x2[j])*w[j]
      sumy2=sumy2+y[i, 1]+y[i, 2]*x2[j]
      sumb2=sumb2+beta[i-1, 1]+beta[i-1, 2]*x2[j]
    }
  
    #calculate acceptance ratio
    A_n=-sumy1-0.5*(y[i,])%*%transpose(y[i,])+sumy2
    A_d=-sumb1-0.5*(beta[i-1,])%*%transpose(beta[i-1,])+sumb2
    A[i]=A_n-A_d
    
    u=runif(1, 0, 1)

    if (log(u)<= A[i]){
      beta[i,1]=y[i,1]
      beta[i,2]=y[i,2]
    } else{
      beta[i, 1]=beta[i-1, 1]
      beta[i, 2]=beta[i-1, 2]
    }

  }

  #sum the means, medians, and CIs with 95% confidence
  mean1 = mean1+mean(beta[1000:5000,1])
  mean2 = mean2+mean(beta[1000:5000,2])
  med1 = med1+median(beta[1000:5000,1])
  med2 = med2+median(beta[1000:5000,2])
  CI1 = CI1+mean(beta[1000:5000,1])+qnorm(0.975)*sqrt(var(beta[1000:5000,1])/n)
  CI11 = CI11+mean(beta[1000:5000,1])-qnorm(0.975)*sqrt(var(beta[1000:5000,1])/n)
  CI2 = CI2+mean(beta[1000:5000,2])+qnorm(0.975)*sqrt(var(beta[1000:5000,2])/n)
  CI21 = CI21+mean(beta[1000:5000,2])-qnorm(0.975)*sqrt(var(beta[1000:5000,2])/n)

}

#calculate means, medians, and confidence intervals
mean1/rep
mean2/rep
med1/rep
med2/rep
CI1/rep
CI11/rep
CI2/rep
CI21/rep

#plot posterior sample and histogram
plot(1000:5000, beta[1000:5000,1], xlab="Iteration", ylab=expression(beta_0), type='l')
hist(beta[1000:5000,1], xlim= c(-1,.5), main=NULL, xlab=expression(beta_0))
plot(1000:5000, beta[1000:5000,2], xlab="Iteration", ylab=expression(beta_1), type='l')
hist(beta[1000:5000,2], xlim= c(3,5.5), main=NULL, xlab=expression(beta_1))
 