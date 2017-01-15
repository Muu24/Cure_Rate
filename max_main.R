##data simulation
library(survival)
library(compoisson)
library(flexsurv)
library(sqldf)
#data(kidney)
#head(kidney)
#str(kidney)
#kidney$disease
##Survival function for suceptible individuals
H_0<-function(gam1,gam2,t){
  return((t/gam2)^(gam1))
}
##group 1, say sex: female
delta_ber<-function(delta,gam1,gam2,p,p0,y_group,beta0,beta1,x_b){
  t=rexp(10000,rate=1)
  dum1=y_group*exp(beta0+beta1*x_b)*H_0(gam1,gam2,t/delta)
  dum2=(sum(exp(-dum1)))/10000
  return(dum2-((p-p0)/(1-p0)))
}
#For group 1 the delta for censor rate of the expential distribution is delta1
p1=0.5
p2=0.4
p01=0.4
p02=0.3
n=20
ni=10
gam1=4
gam2=3
beta0=0
beta1=-log(2)
simulation<-function(gam1,gam2,n,ni,p1,p01,beta0,beta1){
  y_group<-c(rgengamma.orig(10000, shape=1/3, scale=0.001010101,k=9))
  x_b0<-0
  x_b1<-1
  ##group1 p1,p01,x_b=0
  delta1<-uniroot(delta_ber,c(0,100),gam1=gam1,gam2=gam2,p=p1,p0=p01,y_group=y_group,beta0=beta0,beta1=beta1,x_b=0)$root   
  delta2<-uniroot(delta_ber,c(0,100),gam1=gam1,gam2=gam2,p=p1,p0=p01,y_group=y_group,beta0=beta0,beta1=beta1,x_b=1)$root   
  dis1<-rep(NA,n*ni)
  hos1<-rep(NA,n*ni)
  d1<-rep(NA,n*ni)
  t1<-rep(NA,n*ni)
  y_group_sim1<-c(rgengamma.orig(n, shape=1/3, scale=0.001010101,k=9))
  c1<-rexp(n*ni,rate=delta1)
  u1<-runif(n*ni,0,1)
  for(i in 1:n){
    for(j in 1:ni){
      if (u1[10*(i-1)+j]<=p01){
        t1[10*(i-1)+j]=c1[10*(i-1)+j]
        d1[10*(i-1)+j]=0
      }else{
        w=rweibull(1,shape=gam1,scale=(gam2/((y_group_sim1[i]*exp(beta0+beta1*x_b0))^(1/gam1))))
        t1[10*(i-1)+j]=min(w,c1[10*(i-1)+j])
        if(min(w,c1[10*(i-1)+j])==c1[10*(i-1)+j]){
          d1[10*(i-1)+j]=0 
        }else{
          d1[10*(i-1)+j]=1 
        }
      }
      hos1[10*(i-1)+j]=i
    }
  }
  data1<-cbind(t1,d1,x_b0,hos1)
  #group2 p1,p01,x_b=1  
  dis2<-rep(NA,n*ni)
  hos2<-rep(NA,n*ni)
  d2<-rep(NA,n*ni)
  t2<-rep(NA,n*ni)
  y_group_sim2<-c(rgengamma.orig(n,shape=1/3, scale=0.001010101,k=9))
  c2<-rexp(n*ni,rate=delta2)
  u2<-runif(n*ni,0,1)
  for(i in 1:n){
    for(j in 1:ni){
      if (u2[10*(i-1)+j]<=p01){
        t2[10*(i-1)+j]=c2[10*(i-1)+j]
        d2[10*(i-1)+j]=0
      }else{
        w=rweibull(1,shape=gam1,scale=(gam2/((y_group_sim2[i]*exp(beta0+beta1*x_b1))^(1/gam1))))
        t2[10*(i-1)+j]=min(w,c2[10*(i-1)+j])
        if(min(w,c2[10*(i-1)+j])==c2[10*(i-1)+j]){
          d2[10*(i-1)+j]=0 
        }else{
          d2[10*(i-1)+j]=1 
        }
      }
      hos2[10*(i-1)+j]=i
    }
  }
  data2<-cbind(t2,d2,x_b1,hos2)
  data<-rbind(data1,data2)
  return(data)
}
#test cure rate
#cure_rate<-c(rep(NA,1000))
#censor_rate<-c(rep(NA,100))
#for(i in 1:1000){
#  dat<-simulation(gam1=4,gam2=3,n=20,ni=10,p1=0.5,p01=0.4,beta0=0,beta1=-log(2))  
#  dat<-as.data.frame(dat)
#  cure_rate[i]<-length(subset(dat,d1 %in% 0)[,1])/400
  #censor_rate[i]<-length(subset(a,d1 %in% 0.5)[,1])/400
#}
#mean(cure_rate)
#mean(censor_rate)

###Above is the group 1 p0=0.5, p01=0.4




##Max
P_I_1<-function(delta_ij,phi=0,beta0,beta1,beta0_w,beta2_w,beta3_w,x_w,x_b,gam1,gam2,q,sigma,t){
  p0temp<-p0_fct(eta(beta0_w,beta2_w,x_b,beta3_w,x_w))
  Mgftemp<-Mgf_n(N=1000,beta0,beta1,t(matrix(x_b,length(x_b),1000)),gam1,gam2,q,sigma,t(matrix(t,length(x_b),1000)),t(matrix(delta_ij,length(x_b),1000)))$Mgfn
  return(delta_ij+(1-delta_ij)*(1-p0temp)*Mgftemp/(p0temp+(1-p0temp)*Mgftemp))
}
P_I_0<-function(delta_ij,phi=0,beta0,beta1,beta0_w,beta2_w,beta3_w,x_w,x_b,gam1,gam2,q,sigma,t){
  p0temp<-p0_fct(eta(beta0_w,beta2_w,x_b,beta3_w,x_w))
  Mgftemp<-Mgf_n(N=1000,beta0,beta1,t(matrix(x_b,length(x_b),1000)),gam1,gam2,q,sigma,t(matrix(t,length(x_b),1000)),t(matrix(delta_ij,length(x_b),1000)))$Mgfn
  return((1-delta_ij)*(p0temp)/(p0temp+(1-p0temp)*Mgftemp))
}
p0_fct<-function(eta){
  return(1/(1+eta))
}
eta<-function(beta0_w,beta2_w,x_b,beta3_w,x_w){
  eta<-exp(beta0_w+beta2_w*x_b+beta3_w*x_w)
  return(eta)
}
##Mgf numerical for a vector of t
Mgf_n<-function(N,beta0,beta1,x_b,gam1,gam2,q,sigma,t,delta_ij){
  shape=q/sigma
  scale=(gamma(q^(-2)))/(gamma(q^(-2)+sigma/q))
  k=1/(q^2)
  ytemp<-rgengamma.orig(N,shape=shape, scale=scale,k=k)
  #S_s<-exp(-H_0(gam1,gam2,t)%*%ytemp*exp(beta0+beta1*x_b))
  S_s<-exp(-ytemp*exp(beta0+beta1*x_b)*H_0(gam1,gam2,t))
  Mgfn<-colMeans(S_s)
  ##
  h_i<-ytemp*exp(beta0+beta1*x_b)*h_0(gam1,gam2,t)
  f_i<-colMeans(S_s*(h_i^delta_ij))
  return(list(Mgfn=Mgfn, f_i=f_i))
}
#Mgf_n(N=1000,beta0 = 0,beta1=-log(2) ,x_b = 0,gam1 = 4,gam2=1,q=1/3,sigma=10,t=10,delta_ij = 0)
#Mgf_n(N=1000,beta0 = 0,beta1=-log(2) ,x_b = t(matrix(c(1,1,1,0,0),5,1000)),gam1 = 4,gam2=1,q=1/3,sigma=10,t=t(matrix(c(1,2,3,4,10),5,1000)),delta_ij = t(matrix(c(1,1,0,1,0),5,1000)))
##h_0
h_0<-function(gam1,gam2,t){
  return(gam1*(t^(gam1-1))/(gam2^gam1))
}
##
Likelihood<-function(gam1,gam2,beta0,beta1,beta0_w,beta2_w,beta3_w,phi=0,sigma,q,delta_ij,x_w,x_b,t){
  main<-Mgf_n(N=1000,beta0 =0,beta1=beta1 ,x_b = t(matrix(x_b,length(x_b),1000)),gam1 = gam1,gam2=gam2,q=q,sigma=sigma,t=t(matrix(t,length(x_b),1000)),delta_ij = t(matrix(delta_ij,length(x_b),1000)))$f_i
  Mgf<-Mgf_n(N=1000,beta0 = 0,beta1=beta1 ,x_b = t(matrix(x_b,length(x_b),1000)),gam1 = gam1,gam2=gam2,q=q,sigma=sigma,t=t(matrix(t,length(x_b),1000)),delta_ij = t(matrix(delta_ij,length(x_b),1000)))$Mgfn
  prI1<-P_I_1(delta_ij=delta_ij,phi=0,beta0=0,beta1=beta1,beta0_w=beta0_w,beta2_w=0,beta3_w=0,x_w=0,x_b=x_b,gam1=gam1,gam2=gam2,q=q,sigma=sigma,t=t)
  prI0<-1-prI1
  p0temp<-p0_fct(eta(beta0_w,beta2_w,x_b,beta3_w,x_w))
  likelihood1<-(1-p0temp)^delta_ij*main*prI1
  likelihood2<-p0temp*(1-p0temp)^delta_ij*prI0
  likelihood<-likelihood1+likelihood2
  like_final<-prod(likelihood)
  return(like_final)
}
#Likelihood(gam1=4,gam2=3,beta0=0,beta1=-log(2),beta0_w=log(1.5),beta2_w=0,beta3_w=0,phi=0,sigma=1,q=3,delta_ij=c(1,1,0,1,0),x_w=0,x_b=c(1,1,1,0,0),t=c(1,2,3,4,5))
#Likelihood(gam1=4,gam2=3,beta0=0,beta1=-log(2),beta0_w=log(1.5),beta2_w=0,beta3_w=0,phi=0,sigma=1,q=3,delta_ij=delta_ij,x_w=0,x_b=x_b,t=t)
#
dat<-as.data.frame(simulation(gam1=4,gam2=3,n=20,ni=10,p1=0.5,p01=0.4,beta0=0,beta1=-log(2)))  
Max<-function(gam1,gam2,beta1,beta0_w,sigma,q){
  dat<-as.data.frame(simulation(gam1=4,gam2=3,n=20,ni=10,p1=0.5,p01=0.4,beta0=0,beta1=-log(2))) 
  #dat<-dat
  x_b<-c(dat$x_b0)
  t<-c(dat$t1)
  delta_ij<-c(dat$d1)
  beta2_w=0
  beta3_w=0
  x_w=0
  main<-Mgf_n(N=1000,beta0 =0,beta1=beta1 ,x_b = t(matrix(x_b,400,1000)),gam1 = gam1,gam2=gam2,q=q,sigma=sigma,t=t(matrix(t,400,1000)),delta_ij = t(matrix(delta_ij,400,1000)))$f_i
  Mgf<-Mgf_n(N=1000,beta0 = 0,beta1=beta1 ,x_b = t(matrix(x_b,400,1000)),gam1 = gam1,gam2=gam2,q=q,sigma=sigma,t=t(matrix(t,400,1000)),delta_ij = t(matrix(delta_ij,400,1000)))$Mgfn
  prI1<-P_I_1(delta_ij=delta_ij,phi=0,beta0=0,beta1=beta1,beta0_w=beta0_w,beta2_w=0,beta3_w=0,x_w=0,x_b=x_b,gam1=gam1,gam2=gam2,q=q,sigma=sigma,t=t)
  prI0<-1-prI1
  p0temp<-p0_fct(eta(beta0_w,beta2_w,x_b,beta3_w,x_w))
  likelihood1<-((1-p0temp)^delta_ij)*main*prI1
  likelihood2<-p0temp*((1-p0temp)^delta_ij)*prI0
  likelihood<-likelihood1+likelihood2
  like_final<-prod(likelihood)
  return((like_final))
}
Max(4,3,-log(2),log(1.5),1,1/3)
para<-matrix(NA,5,6)
likvalue<-c(rep(NA,5))
for(i in 1:5){
  resultstemp<-optim(c(3,3,-log(2),log(1.5),1,1/3), method = 'L-BFGS-B',lower=c(0,0,-2,0,0,0),function(x) Max(x[1],x[2],x[3],x[4],x[5],x[6]))
  para[i,]<-c(resultstemp$par)
  likvalue[i]<-log(resultstemp$value)
}
colMeans(para)

 