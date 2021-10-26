m=50000
n=50
rho=0.7
a=.25
b=.05
av=0
md=0
msea=0
msem=0
CIU=0
CIL=0

W=rexp(n,rho)
k=rgeom(n,exp(-W))
k=k+1

for (l in 1:20){
t=vector(length=n)
p[1]=1

  for (j in (2:m)){
    for (i in (1:n)) {
     t[i]=rbeta(1, p[j-1]+1, k[i])
    }
    w=-log(t)
    p[j]=rgamma(1, a+n, b+sum(w))
  }

av= av + mean(p[10000:50000])
md = md + median(p[10000:50000])
msea = msea + mean((p[10000:50000]-rho)^2)
msem = msem + median((p[10000:50000]-rho)^2)
CIU=CIU+mean(p[10000:50000])+1.96*sqrt(var(p[10000:50000])/n)
CIL=CIL+mean(p[10000:50000])-1.96*sqrt(var(p[10000:50000])/n)


l=l+1
}

av=av/20
md=md/20
msea=msea/20
msem=msem/20
CIU=CIU/20
CIL=CIL/20

plot(10000:50000, p[10000:50000], xlab="Iteration", ylab=expression(rho), type='l')
hist(p[10000:50000], xlim= c(0,2), main=NULL,xlab=expression(rho))

