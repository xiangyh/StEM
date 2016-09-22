####generate failure data,N=100####
N <- 50
T <- c(0,25,50)
m <- 1:(length(T)-1)

#### define the censored proportion ####
p <- c(0.2,1)


sem.lambda <- NULL
sem.mu <- NULL
mle.lambda <- NULL
mle.mu <- NULL
em.lambda <- NULL
em.mu <- NULL
semem.lambda <- NULL
semem.mu <- NULL


#### BOOTSTRAP ####
for(B in 1:4000)
{ 
  
  #============================
  #   type-I censoring data
  #============================
  x <- rweibull(N,scale=20, shape=1.5)
  t<-rep(NA,N)
  d <- NULL
  r <- NULL
  for(i in m)
  {
    y1 <- which((x>T[i] & x<=T[i+1]) & is.na(t))
    t[y1] <- 0
    d[i]  <- length(y1)
    r[i]  <- floor(p[i]*(length(which(is.na(t)))))
    z1 <- sample(which(is.na(t)),r[i])
    t[z1] <- m[i]
  }
  t[which(is.na(t))] <- length(T)-1
  r[length(T)-1] <- length(which(t==(length(T)-1)))
  
  ## initial value from sample by interval estimation#
  x1<-rep(NA,N)
  for(q in m)
  {
    y2 <- which((x>T[q] & x<=T[q+1]))
    x1[y2] <- (T[q]+T[q+1])/2
  }
  y2 <- which(is.na(x1))
  x1[y2] <- max(T)
  weibull<-function(theta) {
    -sum(suppressWarnings(dweibull(x1, shape = theta[1], scale = theta[2], log = TRUE)) )
  }
  mt<-suppressWarnings(nlm(weibull, theta <- c(1.5,20), hessian=F)$estimate)
  rm(theta)
  inimu<-mt[1]
  inila<-mt[2]
  
  
  #======================
  #     Stochastic EM
  #======================
  
  #### initial value ####
  mu_sem <- inimu
  lambda_sem <- inila
  y1<-NULL
  
  #### SEM ####
  t<-NULL
  for(i in m)
  {
    t <- c(t,rep(i,d[i]+r[i]))
  }
  v<-  tryCatch({
  for(n in 1:1100)
  { 
    xyz<-rep(0,N)
    pwb <- pweibull(T, shape=mu_sem[n], scale=lambda_sem[n])
    for(i in m)
    {
      #### impute interval censored data ####
      y1 <- which(t==i)
      for(j in 1:d[i])
      { 
        u <- runif(1)
        xyz[y1[j]] <- qweibull(u*(pwb[i+1]-pwb[i])+pwb[i], shape=mu_sem[n], scale=lambda_sem[n])
      }
      if(r[i]>=1)
      {
        #### impute right censored data ####
        for(j in (d[i]+1):(d[i]+r[i]))
        { 
          u <- runif(1)
          xyz[y1[j]] <- qweibull(u*(1-pwb[i+1])+pwb[i+1], shape=mu_sem[n], scale=lambda_sem[n]) 
        }
      }
    }
    
    #### M-step: maximun the Q-function usint the imputed data ####
    fn<-function(theta) {
      -sum(suppressWarnings(dweibull(xyz, shape = theta[1], scale = theta[2], log = TRUE)) )
    }
    m0<-suppressWarnings(nlm(fn, theta <- c(mu_sem[n],lambda_sem[n]), hessian=F)$estimate)
    mu_sem <- c(mu_sem,m0[1])
    lambda_sem <- c(lambda_sem,m0[2])
  }
  mean.lambda <- mean(lambda_sem[-c(1:100)])
  mean.mu <- mean(mu_sem[-c(1:100)])
     v<-1
  },error=function(e){
    print("error")
  })
  rm(theta)
  if(v=="error"){
    sem.lambda<-c(sem.lambda,0)
    sem.mu<-c(sem.mu,0)}else{
      sem.lambda <- c(sem.lambda,mean.lambda)
      sem.mu <- c(sem.mu,mean.mu)}


  #======================
  #        MLE
  #======================
  likef<-0
  mle_weibull <- function(theta){
    for(i in m){
      lf <- ((exp(-(T[i]/theta[2])^theta[1])-exp(-(T[i+1]/theta[2])^theta[1]))^d[i])*(exp(-(T[i+1]/theta[2])^theta[1])^r[i])
      likef<-likef+log(lf)
    }
    return(-likef)
  }
  
  m2<-suppressWarnings(nlm(mle_weibull, theta <- c(1.5,20), hessian=F)$estimate)
  mle.lambda <- c(mle.lambda,m2[2])
  mle.mu <- c(mle.mu,m2[1])
  rm(theta)
  
  
  #######################
  ######### EM ##########
  #######################
  mu_em=inimu
  lambda_em=inila
  
  a<-  tryCatch({
    for(L in 1:100)
    {
      fn1 <- function(y){
        log(y)*dweibull(y,shape=mu_em[L],scale=lambda_em[L])
      }  
      
      en1 <- function(T1,T2){
        integrate(fn1,lower=T1,upper=T2)$value/(pweibull(T2,shape=mu_em[L],scale=lambda_em[L])-pweibull(T1,shape=mu_em[L],scale=lambda_em[L]))
      }
      
      lfem1<-0
      for(i in m){
        lf1<-d[i]*en1(T[i],T[i+1])+r[i]*en1(T[i+1],Inf)
        lfem1 <- lfem1+lf1
      }
      
      em_weibull<-function(theta) {
        fn2<-function(y){
          ((y/theta[2])^theta[1])*dweibull(y,shape=mu_em[L],scale=lambda_em[L])
        }
        en2 <- function(T1,T2){
          integrate(fn2,lower=T1,upper=T2)$value/(pweibull(T2,shape=mu_em[L],scale=lambda_em[L])-pweibull(T1,shape=mu_em[L],scale=lambda_em[L]))
        }
        lfem2<-0
        for(i in m){
          lf<-d[i]*en2(T[i],T[i+1])+r[i]*en2(T[i+1],Inf)
          lfem2 <- lfem2+lf
        }
        
        return(-(N*(log(theta[1])-theta[1]*log(theta[2]))+(theta[1]-1)*lfem1-lfem2))
      }
      m1<-suppressWarnings(nlm(em_weibull, theta <- c(mu_em[L],lambda_em[L]), hessian=F)$estimate)
      
      mu_em<-c(mu_em,m1[1])
      lambda_em<-c(lambda_em,m1[2])
      
      if((abs((mu_em[L+1]-mu_em[L]))<10^(-5))&(abs((lambda_em[L+1]-lambda_em[L]))<10^(-5)))
        break
    }
    a<-1
  },error=function(e){
    print("error")
  })
  if(a=="error"){
    em.lambda<-c(em.lambda,0)
    em.mu<-c(em.mu,0)}else{
      em.lambda<-c(em.lambda,lambda_em[length(lambda_em)])
      em.mu<-c(em.mu,mu_em[length(mu_em)])}
  rm(theta)

  #======================
  #     StEM+EM
  #======================
  mu_semem=inimu
  lambda_semem=inila
  
  b<-  tryCatch({
    for(n in 1:300)
    {
      xyz <- rep(0,N)
      pwb <- pweibull(T, shape=mu_semem[n], scale=lambda_semem[n])
      for(i in m)
      {
        y1 <- which(t==i)
        #### impute interval censored data ####
        for(j in 1:d[i])
        { 
          u <- runif(1)
          xyz[y1[j]] <- qweibull(u*(pwb[i+1]-pwb[i])+pwb[i], shape=mu_semem[n], scale=lambda_semem[n])
        }
        
        if(r[i]>=1)
        { #### impute right censored data ####
          
          for(j in (d[i]+1):(d[i]+r[i]))
          { 
            u <- runif(1)
            xyz[y1[j]] <- qweibull(u*(1-pwb[i+1])+pwb[i+1], shape=mu_semem[n], scale=lambda_semem[n]) 
          }
        }
      }
      #### M-step: maximun the Q-function usint the imputed data ####
      semem<-function(theta)
      {
        -sum(suppressWarnings(dweibull(xyz, shape = theta[1], scale = theta[2], log = TRUE)) )
      }
      
      m3<-suppressWarnings(nlm(semem, theta <- c(mu_semem[n],lambda_semem[n]), hessian=F)$estimate)
      mu_semem <- c(mu_semem,m3[1])
      lambda_semem <- c(lambda_semem,m3[2])
    }
    
    la<-mean(lambda_semem[-c(1:100)])
    mu1<-mean(mu_semem[-c(1:100)])
    rm(theta)

  for(L in 1:100)
  {
    fn1 <- function(y)
    {
      log(y)*dweibull(y,shape=mu1[L],scale=la[L])
    }
    
    en1 <- function(T1,T2)
    {
      integrate(fn1,lower=T1,upper=T2)$value/(pweibull(T2,shape=mu1[L],scale=la[L])-pweibull(T1,shape=mu1[L],scale=la[L]))
    }
    
    lfem1<-0
    for(i in m)
    {
      lf1<-d[i]*en1(T[i],T[i+1])+r[i]*en1(T[i+1],Inf)
      lfem1 <- lfem1+lf1
    }
    
    semem_weibull<-function(theta)
    {
      fn2<-function(y){((y/theta[2])^theta[1])*dweibull(y,shape=mu1[L],scale=la[L])}
      
      en2 <- function(T1,T2){integrate(fn2,lower=T1,upper=T2)$value/(pweibull(T2,shape=mu1[L],scale=la[L])-pweibull(T1,shape=mu1[L],scale=la[L]))}
      
      lfem2<-0
      for(i in m)
      {
        lf<-d[i]*en2(T[i],T[i+1])+r[i]*en2(T[i+1],Inf)
        lfem2 <- lfem2+lf
      }
      
      return(-(N*(log(theta[1])-theta[1]*log(theta[2]))+(theta[1]-1)*lfem1-lfem2))
    }
    
    m4<-suppressWarnings(nlm(semem_weibull, theta <- c(mu1[L],la[L]), hessian=F)$estimate)
    
    mu1<-c(mu1,m4[1])
    la<-c(la,m4[2])
    if((abs((mu1[L+1]-mu1[L]))<10^(-5))&(abs((la[L+1]-la[L]))<10^(-5))) break
  }
  b<-1
  },error=function(e){print("error")})
  
  if(b=="error"){
  semem.lambda<-c(semem.lambda,0)
  semem.mu<-c(semem.mu,0)}
  else{
    semem.lambda<-c(semem.lambda,la[length(la)])
    semem.mu<-c(semem.mu,mu1[length(mu1)])}

}



write.table(sem.lambda,"sem.lambda.txt",row.names=F,col.names=F,sep=",")
write.table(sem.mu,"sem.mu.txt",row.names=F,col.names=F,sep=",")
write.table(mle.lambda,"mle.lambda.txt",row.names=F,col.names=F,sep=",")
write.table(mle.mu,"mle.mu.txt",row.names=F,col.names=F,sep=",")
write.table(em.lambda,"em.lambda.txt",row.names=F,col.names=F,sep=",")
write.table(em.mu,"em.mu.txt",row.names=F,col.names=F,sep=",")
write.table(semem.lambda,"semem.lambda.txt",row.names=F,col.names=F,sep=",")
write.table(semem.mu,"semem.mu.txt",row.names=F,col.names=F,sep=",")
