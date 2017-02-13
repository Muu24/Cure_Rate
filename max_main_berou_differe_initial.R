##data simulation
library(survival)
library(compoisson)
library(flexsurv)
#library(e1071)
#library(matrixcalc) # is.positive.definite
#library(expm)
#install.packages('grid')
#install.packages('heR.Misc')
#library(gridSearch)
#install.packages('gridSearch')
#library(sq_1ldf)
#data(kidney)
#head(kidney)
#str(kidney)
#kidney$disease
seed <- 0909
set.seed(seed)
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
ni=2
gam1=4
gam2=3
beta0=0
beta1=-log(2)
beta0_w<-log(1.5)
beta2_w<-log(7/3)-log(1.5)
simulation<-function(gam1,gam2,n,ni,p1,p2,p01,p02,beta0,beta1){
  y_group<-c(rgengamma.orig(10000, shape=0.333333333, scale=0.001010101,k=9.000000000))
  x_b0<-0
  x_b1<-1
  ##group1 p1,p01,x_b=0
  delta1<-uniroot(delta_ber,c(0,100),gam1=gam1,gam2=gam2,p=p1,p0=p01,y_group=y_group,beta0=beta0,beta1=beta1,x_b=0)$root   
  delta2<-uniroot(delta_ber,c(0,100),gam1=gam1,gam2=gam2,p=p2,p0=p02,y_group=y_group,beta0=beta0,beta1=beta1,x_b=1)$root   
  dis1<-rep(NA,n*ni)
  hos1<-rep(NA,n*ni)
  d1<-rep(NA,n*ni)
  t1<-rep(NA,n*ni)
  y_group_sim1<-c(rgengamma.orig(10000, shape=0.333333333, scale=0.001010101,k=9.000000000))
  c1<-rexp(n*ni,rate=delta1)
  u1<-runif(n*ni,0,1)
  for(i in 1:n){
    for(j in 1:ni){
      if (u1[ni*(i-1)+j]<=p01){
        t1[ni*(i-1)+j]=c1[ni*(i-1)+j]
        d1[ni*(i-1)+j]=0
      }else{
        w=rweibull(1,shape=gam1,scale=(gam2/((y_group_sim1[i]*exp(beta0+beta1*x_b0))^(1/gam1))))
        t1[ni*(i-1)+j]=min(w,c1[ni*(i-1)+j])
        if(min(w,c1[ni*(i-1)+j])==c1[ni*(i-1)+j]){
          d1[ni*(i-1)+j]=0 
        }else{
          d1[ni*(i-1)+j]=1 
        }
      }
      hos1[ni*(i-1)+j]=i
    }
  }
  data1<-cbind(t1,d1,x_b0,hos1)
  #group2 p1,p01,x_b=1  
  dis2<-rep(NA,n*ni)
  hos2<-rep(NA,n*ni)
  d2<-rep(NA,n*ni)
  t2<-rep(NA,n*ni)
  y_group_sim2<-c(rgengamma.orig(10000, shape=0.333333333, scale=0.001010101,k=9.000000000))
  c2<-rexp(n*ni,rate=delta2)
  u2<-runif(n*ni,0,1)
  for(i in 1:n){
    for(j in 1:ni){
      if (u2[ni*(i-1)+j]<=p02){
        t2[ni*(i-1)+j]=c2[ni*(i-1)+j]
        d2[ni*(i-1)+j]=0
      }else{
        w=rweibull(1,shape=gam1,scale=(gam2/((y_group_sim2[i]*exp(beta0+beta1*x_b1))^(1/gam1))))
        t2[ni*(i-1)+j]=min(w,c2[ni*(i-1)+j])
        if(min(w,c2[ni*(i-1)+j])==c2[ni*(i-1)+j]){
          d2[ni*(i-1)+j]=0 
        }else{
          d2[ni*(i-1)+j]=1 
        }
      }
      hos2[ni*(i-1)+j]=i
    }
  }
  data2<-cbind(t2,d2,x_b1,hos2)
  data<-rbind(data1,data2)
  return(data)
}
#test cure rate
#cure_rate<-c(rep(NA,1000))
#censor_rate<-c(rep(NA,1000))
#for(i in 1:1000){
#  dat<-simulation(gam1=4,gam2 = 3,n=20,ni=2,p1=0.5,p2=0.4,p01=0.4,p02=0.3,beta0=0,beta1=-log(2))
#  dat<-as.data.frame(dat)
#  cure_rate[i]<-length(subset(dat,d1 %in% 0 &x_b0==0 )[,1])/(n*ni)
#censor_rate[i]<-length(subset(dat,d1 %in% 0 &x_b0==1)[,1])/(n*ni)
##}
#mean(cure_rate)
#mean(censor_rate)

###Above is the group 1 p0=0.5, p01=0.4




##Max
P_I_1<-function(delta_ij,phi=0,beta0,beta1,beta0_w,beta2_w,beta3_w,x_w,x_b,gam1,gam2,q,sigma,t){
  p0temp<-p0_fct(eta(beta0_w,beta2_w,x_b,beta3_w,x_w))
  Mgftemp<-Mgf_n(N=1000,beta0,beta1,t(matrix(x_b,length(x_b),1000)),gam1,gam2,q,sigma,t(matrix(t,length(x_b),1000)),t(matrix(delta_ij,length(x_b),1000)))$Mgfn
  return(delta_ij+(1-delta_ij)*(1-p0temp)*Mgftemp/(p0temp+(1-p0temp)*Mgftemp))
}
#P_I_0<-function(delta_ij,phi=0,beta0,beta1,beta0_w,beta2_w,beta3_w,x_w,x_b,gam1,gam2,q,sigma,t){
#  p0temp<-p0_fct(eta(beta0_w,beta2_w,x_b,beta3_w,x_w))
#  Mgftemp<-Mgf_n(N=1000,beta0,beta1,t(matrix(x_b,length(x_b),1000)),gam1,gam2,q,sigma,t(matrix(t,length(x_b),1000)),t(matrix(delta_ij,length(x_b),1000)))$Mgfn
#  return((1-delta_ij)*(p0temp)/(p0temp+(1-p0temp)*Mgftemp))
#}
p0_fct<-function(eta){
  return(1/(1+eta))
}
eta<-function(beta0_w,beta2_w,x_b,beta3_w,x_w){
  eta<-exp(beta0_w+beta2_w*x_b+beta3_w*x_w)
  return(eta)
}
##Mgf numerical for a vector of t
Mgf_n<-function(N=1000,beta0,beta1,x_b,gam1,gam2,q,sigma,t,delta_ij){
  shape=q/sigma
  scale=(gamma(q^(-2)))/(gamma(q^(-2)+sigma/q))
  k=1/(q^2)
  ytemp<-rgengamma.orig(n=1000,shape=shape, scale=scale,k=k)
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
#Likelihood<-function(gam1,gam2,beta0,beta1,beta0_w,beta2_w,beta3_w,phi=0,sigma,q,delta_ij,x_w,x_b,t){
#  main<-Mgf_n(N=1000,beta0 =0,beta1=beta1 ,x_b = t(matrix(x_b,length(x_b),1000)),gam1 = gam1,gam2=gam2,q=q,sigma=sigma,t=t(matrix(t,length(x_b),1000)),delta_ij = t(matrix(delta_ij,length(x_b),1000)))$f_i
#  Mgf<-Mgf_n(N=1000,beta0 = 0,beta1=beta1 ,x_b = t(matrix(x_b,length(x_b),1000)),gam1 = gam1,gam2=gam2,q=q,sigma=sigma,t=t(matrix(t,length(x_b),1000)),delta_ij = t(matrix(delta_ij,length(x_b),1000)))$Mgfn
#  prI1<-P_I_1(delta_ij=delta_ij,phi=0,beta0=0,beta1=beta1,beta0_w=beta0_w,beta2_w=0,beta3_w=0,x_w=0,x_b=x_b,gam1=gam1,gam2=gam2,q=q,sigma=sigma,t=t)
#  prI0<-1-prI1
#  p0temp<-p0_fct(eta(beta0_w,beta2_w,x_b,beta3_w,x_w))
#  likelihood1<-(1-p0temp)^delta_ij*main*prI1
#  likelihood2<-p0temp*(1-p0temp)^delta_ij*prI0
#  likelihood<-likelihood1+likelihood2
#  like_final<-prod(likelihood)
#  return(like_final)
#}
#Likelihood(gam1=4,gam2=3,beta0=0,beta1=-log(2),beta0_w=log(1.5),beta2_w=0,beta3_w=0,phi=0,sigma=1,q=3,delta_ij=c(1,1,0,1,0),x_w=0,x_b=c(1,1,1,0,0),t=c(1,2,3,4,5))
#Likelihood(gam1=4,gam2=3,beta0=0,beta1=-log(2),beta0_w=log(1.5),beta2_w=0,beta3_w=0,phi=0,sigma=1,q=3,delta_ij=delta_ij,x_w=0,x_b=x_b,t=t)
#
#dat<-as.data.frame(simulation(gam1=4,gam2=3,n=20,ni=10,p1=0.5,p01=0.4,beta0=0,beta1=-log(2)))  

#Max<-function(par=c(gam1,gam2,beta1,beta0_w,sigma,q)){
  #dat<-as.data.frame(simulation(gam1=4,gam2=3,n=20,ni=10,p1=0.5,p01=0.4,beta0=0,beta1=-log(2))) 
  #dat<-dat
  #gam1<-par[1]
  #gam2<-par[2]
  #beta1<-par[3]
  #beta0_w<-par[4]
  #sigma<-par[5]
  #q<-par[6]
  #x_b<-c(dat$x_b0)
  #t<-c(dat$t1)
  #delta_ij<-c(dat$d1)
  #beta2_w=0
  #beta3_w=0
  #x_w=0
  #main<-Mgf_n(N=1000,beta0 =0,beta1=beta1 ,x_b = t(matrix(x_b,400,1000)),gam1 = gam1,gam2=gam2,q=q,sigma=sigma,t=t(matrix(t,400,1000)),delta_ij = t(matrix(delta_ij,400,1000)))$f_i
  #Mgf<-Mgf_n(N=1000,beta0 = 0,beta1=beta1 ,x_b = t(matrix(x_b,400,1000)),gam1 = gam1,gam2=gam2,q=q,sigma=sigma,t=t(matrix(t,400,1000)),delta_ij = t(matrix(delta_ij,400,1000)))$Mgfn
  #prI1<-P_I_1(delta_ij=delta_ij,phi=0,beta0=0,beta1=beta1,beta0_w=beta0_w,beta2_w=0,beta3_w=0,x_w=0,x_b=x_b,gam1=gam1,gam2=gam2,q=q,sigma=sigma,t=t)
  #prI0<-1-prI1
  #p0temp<-p0_fct(eta(beta0_w,beta2_w,x_b,beta3_w,x_w))
  #likelihood1<-((1-p0temp)^delta_ij)*main*prI1
  #likelihood2<-p0temp*prI0
  #likelihood<-likelihood1+likelihood2
  #like_final<-sum(log(likelihood))
  #return(like_final)
#}
#Max(c(gam1=4,gam2=3,beta1=-log(2),beta0_w=log(1.5),sigma=1,q=1/3))
#####search for adequate initial values
#toTest <- expand.grid(
#  alpha = seq(0.1, 1, by = 0.1), 
#  eta = seq(0.1, 1, by = 0.1), 
#  gamma = seq(0.1, 1, by = 0.1))
#ml <- apply(toTest, 1, function(x){
#  exp(sum(x))
#})
#toTest[which.max(ml), ]
grid_initial<-function(gam1,gam2,beta1,beta0_w,beta2_w,sigma,q,sep=0.05){
  v_gam1<-seq(gam1*0.8,gam1*1.2,by=gam1*sep)
  v_gam2<-seq(gam2*0.8,gam2*1.2,by=gam2*sep)
  v_beta1<-seq(beta1*0.8,beta1*1.2,by=beta1*sep)
  v_beta0_w<-seq(beta0_w*0.8,beta0_w*1.2,by=beta0_w*sep)
  v_beta2_w<-seq(beta2_w*0.8,beta2_w*1.2,by=beta2_w*sep)
  v_sigma<-seq(sigma*0.8,sigma*1.2,by=sigma*sep)
  v_q<-seq(q*0.8,q*1.2,by=q*sep)
  #grid_matrix<-rbind(v_gam1,v_gam2,v_beta1,v_beta0_w,v_sigma,v_q)
  #maxtemp<-max(v_gam1[1],v_gam2[1],v_beta1[1],v_beta0_w[1],v_sigma[1],v_q[1])
  #v_initial<-gridSearch(fun=Max,levels=list(v_gam1,v_gam2,v_beta1,v_beta0_w,v_sigma,v_q))
  #v_initial<-optgrid(f=Max,par=c(gam1,gam2,beta1,beta0_w,sigma,q),incr=c(gam1*0.1,gam2*0.1,v_beta1*0.1,v_beta0_w*0.1,sigma*0.1,q*0.1),lower=c(gam1*0.8,gam2*0.8,v_beta1*0.8,v_beta0_w*0.8,sigma*0.8,q*0.8))
  #v_initial<-tune(methods=Max,)
  v_initial<-expand.grid(v_gam1,v_gam2,v_beta1,v_beta0_w,v_beta2_w,v_sigma,v_q)
  v_sample<-v_initial[sample(nrow(v_initial),200,replace=FALSE),]
  v_sample<-as.matrix(v_sample)
  ml <- apply(v_sample,1,Intapp)
  #ml <- apply(v_initial,1,function(x) Max(x[1],x[2],x[3],x[4],x[5],x[6]))
  true_initial<-v_sample[which.max(ml),]
  return(true_initial)
  #return(list(true_initial,ml))
}

Intapp<-function(par=c(gam1,gam2,beta1,beta0_w,beta2_w,sigma=1,q=1/3)){
  ##do the marginal likelihood
  gam1<-par[1]
  gam2<-par[2]
  beta1<-par[3]
  beta0_w<-par[4]
  beta2_w=par[5]
  sigma<-par[6]
  q<-par[7]
  x_b<-c(dat$x_b0)
  t<-c(dat$t1)
  delta_ij<-c(dat$d1)
  beta0=0
  beta3_w=0
  x_w=0
  a<-matrix(rep(NA,3),20,3)
  b<-list(a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a)
  gp=1
  while(gp<=20){ 
    # sqldf('select b from c where a==i') &x_b0 %in% (gp-1)
    b[[gp]]<-subset(dat,hos1 %in% gp)[,c(1,2,3)]
    gp<-gp+1
  }
  shape=q/sigma
  scale=(gamma(q^(-2)))/(gamma(q^(-2)+sigma/q))
  k=1/(q^2)
  #ytemp<-rgengamma.orig(N,shape=shape, scale=scale,k=k)
  main<-rep(NA,20)
  p0temp<-rep(NA,20)
  for(i in 1:20){
    ytemp<-rgengamma.orig(n=1000,shape=shape, scale=scale,k=k)
    ff<-exp(beta0+beta1*b[[i]]$x_b0)*H_0(gam1,gam2,b[[i]]$t1)
    gg<-x_b_cal<-t(matrix(b[[i]]$x_b0,length(b[[i]]$x_b0),1000))
    t_cal<-t(matrix(b[[i]]$t1,length(b[[i]]$t1),1000))
    delta_ij_cal<-t(matrix(b[[i]]$d1,length(b[[i]]$d1),1000))
    #S_s<-exp(-H_0(gam1,gam2,t)%*%ytemp*exp(beta0+beta1*x_b))
    prI1_tol<-P_I_1(delta_ij=b[[i]]$d1,phi=0,beta0=0,beta1=beta1,beta0_w=beta0_w,beta2_w=0,beta3_w=0,x_w=0,x_b=b[[i]]$x_b0,gam1=gam1,gam2=gam2,q=q,sigma=sigma,t=b[[i]]$t1)
    prI0<-1-prI1_tol
    p0temp<-p0_fct(eta(beta0_w,beta2_w,x_b_cal,beta3_w,x_w))
    A<-(exp(-ytemp*(exp(beta0+beta1*x_b_cal)*H_0(gam1,gam2,t_cal)))*(ytemp*(exp(beta0+beta1*x_b_cal)*h_0(gam1,gam2,t_cal)))^(delta_ij_cal)*t(matrix(prI1_tol,length(prI1_tol),1000))+prI0*p0temp)
    main[i]<-mean(apply(A, 1, prod))
    #<-colMeans(S_s)
  }
  Likelihood<-sum(log(main))
  #prI1_tol<-P_I_1(delta_ij=delta_ij,phi=0,beta0=0,beta1=beta1,beta0_w=beta0_w,beta2_w=0,beta3_w=0,x_w=0,x_b=x_b,gam1=gam1,gam2=gam2,q=q,sigma=sigma,t=t)
  #prI0<-1-prI1_tol
  #likelihood1<-c()
  #likelihood2<-p0temp*prI0
  #likelihood<-likelihood1+likelihood2
  #like_final<-sum(log(likelihood))
  return(Likelihood)
}

Averagemax<-function(par=c(gam1,gam2,beta1,beta0_w,beta2_w,sigma,q)){
  gam1<-par[1]
  gam2<-par[2]
  beta1<-par[3]
  beta0_w<-par[4]
  beta2_w<-par[5]
  sigma<-par[6]
  q<-par[7]
  Aver_max<-c(rep(NA,10))
  for(i in 1:10){
    Aver_max[i]<-Intapp(c(gam1,gam2,beta1,beta0_w,beta2_w,sigma,q)) 
  }
  return(mean(Aver_max))
}
#Averagemax(c(4.399999,3.6,-0.7624618,0.4865585,0.9999995,0.2666662 ))
#Averagemax(c(4,3,-0.6,0.4,1,1/3))
#Max(c(4.399999,3.6,-0.7624618,0.4865585,0.9999995,0.2666662 ))
#Max(try_int)
#Averagemax(c(4,3,-log(2),beta0_w=log(1.5),sigma=1,q=1/3))

#aa<-optim(c(try_int),method = 'L-BFGS-B',lower=c(0,0,-5,0,0,0),control=list(fnscale=-1),fn=Max,hessian = TRUE)
#aa<-optim(c(try_int[[1]]),method = 'L-BFGS-B',lower=c(0,0,-5,0,0,0),function(x) Max(x[1],x[2],x[3],x[4],x[5],x[6]),hessian = TRUE)
#para<-matrix(NA,5,6)
#likvalue<-c(rep(NA,5))
dat<-simulation(gam1=4,gam2 = 3,n=20,ni=2,p1=0.5,p2=0.4,p01=0.4,p02=0.3,beta0=0,beta1=-log(2))
dat<-as.data.frame(dat)
Intapp(c(4,3,-log(2),beta0_w=log(1.5),beta2_w=log(7/3)-log(1.5),sigma=1,q=1/3))
#aa<-optim(c(try_int),method = 'L-BFGS-B',lower=c(0,0,-5,0,0,0,0),control=list(fnscale=-1),fn=Intapp)

resultstemp<-vector('list',10)
for(i in 1:10){
  dat<-simulation(gam1=4,gam2 = 3,n=20,ni=2,p1=0.5,p2=0.4,p01=0.4,p02=0.3,beta0=0,beta1=-log(2))
  dat<-as.data.frame(dat)
  try_int<-grid_initial(4,3,-log(2),beta0_w=log(1.5),beta2_w=log(7/3)-log(1.5),sigma=1,q=1/3)
  aa<-try(optim(c(try_int),method = 'L-BFGS-B',lower=c(0,0,-5,0,0,0,0),control=list(fnscale=-1),fn=Intapp),TRUE)
  resultstemp[[i]]<-aa
    #try(optim(c(4,3,-log(2),log(1.5),1,1/3), method = 'L-BFGS-B',lower=c(0,0,-5,0,0,0),function(x) Max(x[1],x[2],x[3],x[4],x[5],x[6]),hessian = TRUE),TRUE)
  #para[i,]<-c(resultstemp$par)
  #likvalue[i]<-log(resultstemp$value)
}
#aa<-try(optim(c(4,3,-log(2),log(1.5),1,1/3), method = 'L-BFGS-B',lower=c(0,0,-5,0,0,0),Max,hessian = TRUE),TRUE)
working <- resultstemp[sapply(resultstemp, function(x) !inherits(x, "try-error"))]

para<-matrix(NA,length(working),8)
for(i in 1:nrow(para)){
  para[i,] <- c(working[[i]]$par,working[[i]]$value)
}
para
#aa<-c(4,3,-log(2),log(1.5),3,4)
#colMeans(para)
#dd<-optim(c(4,3,-log(2),log(1.5),3,4), method = 'L-BFGS-B',lower=c(0,0,-5,0,0,0),function(x) Max(x[1],x[2],x[3],x[4],x[5],x[6]),hessian = TRUE)
#bb<-maxLik(Max,start=c(try_int[[1]]))
#tem<-c(bb$estimate)
#cc<-maxLik(Averagemax,start=c(try_int[[1]]),method = 'NR')
saveRDS(para,paste("results",seed,".RDS",sep=""))

####KM
#po_fit<-survfit(Surv(time=dat$t1[1:200],dat$d1[1:200])~1,conf.type="none")
#po_fit1<-survfit(Surv(time=dat$t1[200:400],dat$d1[200:400])~1,conf.type="none")
#plot(po_fit)
#plot.survfit(po_fit)
#library(survival)
#lines(po_fit1,col = 'red')
#install.packages('KMsurv')
#library(KMsurv)
#data("aids")
#head(aids)
#data(bfeed)
#head(bfeed)
##data(channing)
#head(channing)
