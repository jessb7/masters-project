m=10000
n=100
rho=4
r=3
a=.25
b=.05
theta=.5

W=rexp(n,rho)
lambda=rgamma(n, r, exp(-W)/(1-exp(-W)))
k=rpois(n, lambda)
#k=rgeom(n,exp(-W))
k=k+1

#for (l in 1:20){
  t=vector(length=n)
  p=vector(length=m)
  rnew=as.vector(r)
  A=vector(length=m)
  p[1]=rho
  
  for (j in (2:m)){
    t=rbeta(n, p[j-1]+r[j-1], k)
    w=-log(t)
    
    p[j]=rgamma(1, a+n, b+sum(w))
    
    rnew[j]=r[j-1]
    ran=runif(1)
    if (rnew[j]==1 || ran>0.5) {
      rnew[j]=rnew[j]+1
    }else {
      rnew[j]=rnew[j]-1
    }
    

    
    sumw=0
    prodr=1
    prodrnew=1
    for (i in 1:n){
      sumw=sumw+w[i]
      prodr=prodr*choose(k[i]+r[j-1]-2, k[i]-1)
      prodrnew=prodrnew*choose(k[i]+rnew[j]-2, k[i]-1)
    }
    
    A_num=((1-theta)^(rnew[j]-1)*exp(-rnew[j]*sumw)*prodrnew)
    A_denom=((1-theta)^(r[j-1]-1)*exp(-r[j-1]*sumw)*prodr)
    A[j]=A_num/A_denom
    #A[j]=((1-theta)^rnew[j]*exp(-rnew[j]*sumw)*prodrnew)/((1-theta)^r[j-1]*exp(-r[j-1]*sumw)*prodr)
    
    u=runif(1, 0, 1)
    
    if (u <= A[j]){
      r[j]=rnew[j]
    } else{
      r[j]=r[j-1]
    }
      
  }
  
#}
  