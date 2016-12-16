########################################################################
## MCMC update
########################################################################

updateMCMC<-function(curr,prior,data,funcs=funcs){
  
  ## RJMCMC update
  
  u<-sample(1:3,size=1)
  if(curr$nbasis==0){
    u<-1 # birth for sure
  }
  if(curr$nbasis==prior$maxBasis){
    u<-sample(2:3,size=1) # no birth
  }
  if(u==1){ # birth
    curr<-funcs$birth(curr,prior,data)
  } else if(u==2){ # death
    curr<-funcs$death(curr,prior,data)
  } else{ # change
    curr<-funcs$change(curr,prior,data)
  }
  
  ## Gibbs updates
  
  # beta
  curr$beta<-curr$bhat/(1+curr$beta.prec)+curr$R.inv.t%*%rnorm(curr$nc)*sqrt(curr$s2/(1+curr$beta.prec)/data$temp.ladder[curr$temp.ind])
  
  # lambda
  lam.a<-prior$h1+curr$nbasis
  lam.b<-prior$h2+1
  curr$lam<-rgammaTemper(1,lam.a,lam.b,data$temp.ladder[curr$temp.ind])
  
  # s2
  qf2<-crossprod(curr$R%*%curr$beta)
  curr$s2.rate<-(data$ssy + (1+curr$beta.prec)*qf2 - 2*crossprod(curr$beta,curr$Xty[1:curr$nc]))/2
  s2.a<-prior$g1+(data$n+curr$nbasis+1)/2
  s2.b<-prior$g2+curr$s2.rate
  if(s2.b<=0){
    prior$g2<-prior$g2+1
    s2.b<-prior$g2+curr$s2.rate
    warning('Increased g2 for numerical stability')
  }
  curr$s2<-rigammaTemper(1,s2.a,s2.b,data$temp.ladder[curr$temp.ind])
  if(is.nan(curr$s2) | is.na(curr$s2)) # major variance inflation, get huge betas from curr$R.inv.t, everything becomes unstable
    browser()
  if(curr$s2==0 | curr$s2>1e10){ # tempering instability, this temperature too small
    curr$s2<-runif(1,0,1e6)
    prior$g2<-prior$g2+1
    warning('Small temperature too small...increased g2 for numerical stability')
    #browser()
  }
  
  # beta.prec
  beta.prec.a<-prior$a.beta.prec+(curr$nbasis+1)/2
  beta.prec.b<-prior$b.beta.prec+1/(2*curr$s2)*qf2
  curr$beta.prec<-rgammaTemper(1,beta.prec.a,beta.prec.b,data$temp.ladder[curr$temp.ind])
  
  ## save log posterior
  curr$lpost<-lp(curr,prior,data)
  
  return(curr)
}