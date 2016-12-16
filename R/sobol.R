############################################################
## get Sobol decomposition
############################################################

#' @title BASS Sensitivity Analysis
#'
#' @description Decomposes the variance of the BASS model into variance due to main effects, two way interactions, and so on, similar to the ANOVA decomposition for linear models.  Uses the Sobol' decomposition, which can be done analytically for MARS models.
#' @param bassMod a fitted model output from the \code{bass} function.
#' @param mcmc.use an integer vector indexing which MCMC iterations to use for sensitivity analysis.
#' @param func.var an integer indicating which functional variable to make sensitivity indices a function of.  Disregard if \code{bassMod} is non-functional or if scalar sensitivity indices are desired.
#' @param xx.func.var grid for functional variable specified by \code{func.var}.  Disregard if \code{func.var} is not specified.  If \code{func.var} is specified and \code{xx.func.var} not specified, the grid used to fit \code{bass} will be used.
#' @param verbose logical; should progress be displayed?
#' @details Performs analytical Sobol' decomposition for each MCMC iteration in mcmc.use (each corresponds to a MARS model), yeilding a posterior distribution of sensitivity indices.  Can obtain Sobol' indices as a function of one functional variable.
#' @return If non-functional (\code{func.var = NULL}), a list with two elements:
#'  \item{S}{a data frame of sensitivity indices with number of rows matching the length of \code{mcmc.use}.  The columns are named with a particular main effect or interaction.  The values are the proportion of variance in the model that is due to each main effect or interaction.}
#'  \item{T}{a data frame of total sensitivity indices with number of rows matching the length of \code{mcmc.use}.  The columns are named with a particular variable.}
#'  Otherwise, a list with four elements:
#'  \item{S}{an array with first dimension corresponding to MCMC samples (same length as \code{mcmc.use}), second dimension corresponding to different main effects and interactions (labeled in \code{names.ind}), and third dimension corresponding to the grid used for the functional variable.  The elements of the array are sensitivity indices.}
#'  \item{S.var}{same as \code{S}, but scaled in terms of total variance rather than percent of variance.}
#'  \item{names.ind}{a vector of names of the main effects and interactions used.}
#'  \item{xx}{the grid used for the functional variable.}
#'
#' @keywords Sobol decomposition
#' @seealso \link{bass} for model fitting and \link{predict.bass} for prediction.
#' @export
#' @examples
#' # See examples in bass documentation.
#'
sobol<-function(bassMod,mcmc.use=NULL,func.var=NULL,xx.func.var=NULL,verbose=TRUE){
  if(class(bassMod)!='bass')
    stop('First input needs to be a bass object')
  if(bassMod$p==1 & !bassMod$func)
    stop('Sobol only used for multiple input models')
  if(bassMod$p==1 & sum(func.var)>0)
    stop('Cannot decompose the variance in terms of only one variable')
  mcmc.use.poss<-1:((bassMod$nmcmc-bassMod$nburn)/bassMod$thin)
  if(any(!(mcmc.use%in%mcmc.use.poss))){
    mcmc.use<-mcmc.use.poss
    warning('disregarding mcmc.use because of bad values')
  }
  if(any(is.null(mcmc.use))){
    mcmc.use<-mcmc.use.poss
  }
  
  if(is.null(func.var)){
    func<-F
  } else{
    func<-T
    if(!bassMod$func){
      func<-F
      warning('disregarding func.var because bassMod parameter is not functional')
    }
  }

  if(!is.null(xx.func.var) & !func)
    warning('disregarding xx.func.var because bassMod parameter is not functional')
  
  if(func){
    if(!(func.var%in%(1:ncol(bassMod$xx.func))))
      stop('func.var in wrong range of values')
    if(is.null(xx.func.var)){
      xx.func.var<-bassMod$xx.func[,func.var,drop=F]
    } else{
      rr<-range(xx.func.var)
      if(rr[1]<bassMod$range.func[1,func.var] | rr[2]>bassMod$range.func[2,func.var])
        warning(paste('range of func.var in bass function (',bassMod$range.func[1,func.var],',',bassMod$range.func[2,func.var],') is smaller than range of xx.func.var (',rr[1],',',rr[2],'), indicating some extrapolation',sep=''))
      xx.func.var<-scale.range(xx.func.var,bassMod$range.func[,func.var])
    }
    return(sobol_des_func(bassMod=bassMod,mcmc.use=mcmc.use,verbose=verbose,func.var=func.var,xx.func.var=xx.func.var))
  } else{
    return(sobol_des(bassMod=bassMod,mcmc.use=mcmc.use,verbose=verbose)) # applies to both des & func as long as functional sobol indices are not desired
  }
}






## get sobol indices - no functional

sobol_des<-function(bassMod,mcmc.use,verbose){
  models<-bassMod$model.lookup[mcmc.use] # only do the heavy lifting once for each model
  uniq.models<-unique(models)
  nmodels<-length(uniq.models)
  maxInt.tot<-sum(bassMod$maxInt.des)+sum(bassMod$maxInt.cat)+sum(bassMod$maxInt.func)
  maxBasis<-max(bassMod$nbasis)
  q<-bassMod$degree
  pdes<-sum(bassMod$pdes)
  pcat<-sum(bassMod$pcat)
  pfunc<-sum(bassMod$pfunc)
  p<-pdes+pcat+pfunc


  ################################################
  # get combs & ll including functional variables
  allCombs<-getCombs(bassMod,uniq.models,nmodels,maxBasis,maxInt.tot)
  combs<-allCombs$combs
  names.ind<-allCombs$names.ind
  ll<-allCombs$ll
  num.ind<-allCombs$num.ind
  cs.num.ind<-allCombs$cs.num.ind
  ################################################
  sob<-array(0,dim=c(length(mcmc.use),sum(num.ind)))

  if(verbose)
    cat('Sobol Start',timestamp(prefix='#--',suffix='--#',quiet=T),'Models:',length(unique(models)),'\n')

  i<-1
  mod.count<-0
  for(mod in uniq.models){ #do this in parallel?
    mod.count<-mod.count+1
    mcmc.use.mod<-mcmc.use[models==mod] # which part of mcmc.use does this correspond to?
    mod.ind<-i:(i+length(mcmc.use.mod)-1)
    M<-bassMod$nbasis[mcmc.use.mod][1] # number of basis functions in this model
    if(M>0){
      # for each model, m, tl stores everything necessary to get sobol indices
      tl<-get_tl(bassMod,mcmc.use.mod,M,mod,p,q,cs.num.ind,combs)
      if(tl$cat)
        tl<-add_tlCat(tl,bassMod,mcmc.use.mod,mod)
      tl<-add_tl(tl,p)
      lens<-apply(tl$Kind,1,function(x) length(na.omit(x)))


      var.tot<-Vu(1:p,tl) # total variance
      vars.used<-unique(unlist(na.omit(c(tl$Kind)))) # which variables are used?
      vars.used<-sort(vars.used)

      tl$integrals<-matrix(0,nrow=length(mcmc.use.mod),ncol=max(cs.num.ind)) # where we store all the integrals (not normalized by subtraction) - matches dim of sob
      tl$integrals[,which(combs[[1]]%in%vars.used)]<-apply(t(vars.used),2,Vu,tl=tl)
      sob[mod.ind,1:cs.num.ind[1]]<-tl$integrals[,1:cs.num.ind[1]]

      if(max(lens)>1){ # if there are any interactions
        for(l in 2:max(lens)){ # must go in order for it to work (tl$integrals is made sequentially)
          int.l.ind<-(cs.num.ind[l-1]+1):cs.num.ind[l] # index for which elements columns of sob this corresponds to
          basis.int.l.use<-matrix(nrow=M,ncol=num.ind[l]) # mm[i,j]=1 if interaction combs[[l]][,j] occurs in ith basis function
          for(m in 1:M){
            #basis.int.l.use[m,]<-unlist(lapply(ll[[l]],function(el){prod(el%in%tl$Kind[m,])})) #all the variables in question must be in the basis
            basis.int.l.use[m,]<-apply(combs[[l]],2,function(el){prod(el%in%tl$Kind[m,])})
          }
          use<-colSums(basis.int.l.use)>0 # indicator for 
          tl$integrals[,int.l.ind[use]]<-apply(combs[[l]][,use,drop=F],2,Vu,tl=tl) # perform the necessary integration
          sob.l<-matrix(0,nrow=length(mcmc.use.mod),ncol=num.ind[l])
          for(l2 in (1:num.ind[l])[use]){ # only go through the interactions that are actually in Kind (but still allow for 2-way when actual is 3-way, etc.)
            sob.l[,l2]<-VuInt(combs[[l]][,l2],tl) # do the normalizing
          }
          sob[mod.ind,int.l.ind]<-sob.l
        }
      }
      sob[mod.ind,]<-sob[mod.ind,]/var.tot # sobol indices
      i<-i+length(mcmc.use.mod) # for index


      if(verbose & mod.count%%10==0)
        cat('Sobol',timestamp(prefix='#--',suffix='--#',quiet=T),'Model:',mod.count,'\n')


    }
  }

  sob<-as.data.frame(sob)
  names(sob)<-unlist(names.ind) # give labels

  
  if(verbose)
    cat('Total Sensitivity',timestamp(prefix='#--',suffix='--#',quiet=T),'\n')

  tot<-getTot(combs,sob,names.ind,p,maxInt.tot,allCombs$aa)
  
  
  sob.reorder<-NA
  sob.reorder[1:length(names.ind[[1]])]<-mixOrd(allCombs$dispNames[[1]])
  for(l in 2:length(names.ind)){
    sob.reorder[(cs.num.ind[l-1]+1):(cs.num.ind[l])]<-cs.num.ind[l-1]+mixOrd(allCombs$dispNames[[l]])
  }
  sob<-sob[,sob.reorder,drop=F]
  names(sob)<-unlist(allCombs$dispNames)[sob.reorder]
  tot<-tot[,sob.reorder[1:length(names.ind[[1]])],drop=F]
  names(tot)<-allCombs$dispNames[[1]][sob.reorder[1:length(names.ind[[1]])]]
  
  ret<-list(S=sob,T=tot,func=F)
  class(ret)<-'bassSob'
  
  #browser()
  
  return(ret)
}






## get sobol indices - functional

sobol_des_func<-function(bassMod,mcmc.use,verbose,func.var,xx.func.var){
  models<-bassMod$model.lookup[mcmc.use] # only do the heavy lifting once for each model
  uniq.models<-unique(models)
  nmodels<-length(uniq.models)
  maxInt.tot<-sum(bassMod$maxInt.des)+sum(bassMod$maxInt.cat)+sum(bassMod$maxInt.func)
  maxBasis<-max(bassMod$nbasis)
  q<-bassMod$degree
  pdes<-sum(bassMod$pdes)
  pcat<-sum(bassMod$pcat)
  pfunc<-sum(bassMod$pfunc)
  p<-pdes+pcat+pfunc
  
  ################################################
  # get combs & ll including functional variables
  allCombs<-getCombs(bassMod,uniq.models,nmodels,maxBasis,maxInt.tot,func.var)
  combs<-allCombs$combs
  names.ind<-allCombs$names.ind
  ll<-allCombs$ll
  num.ind<-allCombs$num.ind
  cs.num.ind<-allCombs$cs.num.ind
  ################################################
  sob<-sob2<-array(0,dim=c(length(mcmc.use),sum(num.ind),length(xx.func.var)))

  if(verbose)
    cat('Sobol Start',timestamp(prefix='#--',suffix='--#',quiet=T),'Models:',length(unique(models)),'\n')

  i<-1
  mod.count<-0
  for(mod in uniq.models){ #do this in parallel?
    mod.count<-mod.count+1
    mcmc.use.mod<-mcmc.use[models==mod] # which part of mcmc.use does this correspond to?
    mod.ind<-i:(i+length(mcmc.use.mod)-1)
    M<-bassMod$nbasis[mcmc.use.mod][1] # number of basis functions in this model
    if(M>0){
      tl<-get_tl(bassMod,mcmc.use.mod,M,mod,p,q,cs.num.ind,combs,func.var,xx.func.var)
      if(tl$cat)
        tl<-add_tlCat(tl,bassMod,mcmc.use.mod,mod)
      tl<-add_tl(tl,p)
      lens<-apply(tl$Kind,1,function(x) length(na.omit(x)))

      tl$tfunc.basis<-makeBasisMatrixVar(mod,M,vars=bassMod$vars.func,signs=bassMod$signs.func,knots.ind=bassMod$knotInd.func,q=bassMod$degree,xxt=t(tl$xx),n.int=bassMod$n.int.func,xx.train=bassMod$xx.func,var=func.var)[-1,,drop=F]

      var.tot<-Vu_des_func(1:p,tl) # total variance
      vars.used<-unique(unlist(na.omit(c(tl$Kind)))) # which variables are used?
      vars.used<-sort(vars.used)
      
      tl$integrals<-array(0,dim=c(length(mcmc.use.mod),max(cs.num.ind),length(xx.func.var))) # where we store all the integrals (not normalized by subtraction)
      jj=0
      for(pp in vars.used){
        jj=jj+1
        tl$integrals[,jj,]<-Vu_des_func(pp,tl)
      }
      
      sob[mod.ind,1:cs.num.ind[1],]<-tl$integrals[,1:cs.num.ind[1],]
      
      if(max(lens)>1){ # if there are any interactions
        for(l in 2:max(lens)){ # must go in order for it to work (tl$integrals is made sequentially)
          int.l.ind<-(cs.num.ind[l-1]+1):cs.num.ind[l]
          basis.int.l.use<-matrix(nrow=M,ncol=length(ll[[l]]))
          for(m in 1:M){
            basis.int.l.use[m,]<-unlist(lapply(ll[[l]],function(el){prod(el%in%tl$Kind[m,])})) # all the variables in question must be in the basis
          }
          use<-colSums(basis.int.l.use)>0
          for(pp in which(use)){
            tl$integrals[,int.l.ind[pp],]<-Vu_des_func(combs[[l]][,pp,drop=F],tl) # perform the necessary integration
          }
          sob.l<-array(0,dim=c(length(mcmc.use.mod),length(ll[[l]]),length(xx.func.var)))

          for(l2 in (1:length(ll[[l]]))[use]){ # only go through the interactions that are actually in Kind (but still allow for 2-way when actual is 3-way, etc.)

            sob.l[,l2,]<-VuInt_des_func(ll[[l]][[l2]],tl) # do the normalizing
          }
          sob[mod.ind,int.l.ind,]<-sob.l
        }
      }
      
      kk<-0
      for(ii in mod.ind){
        kk=kk+1
        sob2[ii,,]<-t(t(sob[ii,,])/var.tot[kk,]) # sobol indices
      }
      i<-i+length(mcmc.use.mod) # for index
      
      if(verbose & mod.count%%10==0)
        cat('Sobol',timestamp(prefix='#--',suffix='--#',quiet=T),'Model:',mod.count,'\n')
    }
  }
  
  # reorder for display
  sob.reorder<-NA
  sob.reorder[1:length(names.ind[[1]])]<-mixOrd(allCombs$dispNames[[1]])
  for(l in 2:length(names.ind)){
    sob.reorder[(cs.num.ind[l-1]+1):(cs.num.ind[l])]<-cs.num.ind[l-1]+mixOrd(allCombs$dispNames[[l]])
  }
  sob<-sob[,sob.reorder,,drop=F]
  sob2<-sob2[,sob.reorder,,drop=F]

  
  ret<-list(S=sob2,S.var=sob,names.ind=unlist(allCombs$dispNames)[sob.reorder],xx=unscale.range(tl$xx,bassMod$range.func),func=T)
  class(ret)<-'bassSob'
  return(ret)
}

# get total sensitivity indices
getTot<-function(combs,sob,names.ind,p,maxInt.tot,aa){ 
  vars.use<-unique(unlist(combs))
  puse<-length(vars.use)
  ncombs<-sapply(combs,ncol)
  tot<-sob[,1:ncombs[1]]
  ncombs[1]<-0
  int.use<-(2:maxInt.tot)[!is.na(aa[-c(1,maxInt.tot+1)])] # should get int.use somewhere else...
  if(maxInt.tot>1){
    for(pp in 1:puse){
      for(l in int.use){
        tot[,pp]<-tot[,pp]+rowSums(
          sob[,puse+ncombs[l-1]+which(apply(t(combs[[l]]),1,function(r){vars.use[pp]%in%r})),drop=F]
        )
      }
    }
  }
  return(tot)
}



########################################################################
## processing functions
########################################################################

## make for only one variable
makeBasisMatrixVar<-function(i,nbasis,vars,signs,knots.ind,q,xxt,n.int,xx.train,var){ 
  n<-ncol(xxt)
  tbasis.mat<-matrix(nrow=nbasis+1,ncol=n)
  tbasis.mat[1,]<-1
  if(nbasis>0){
    for(m in 1:nbasis){
      if(all(na.omit(vars[i,m,])!=var)){
        tbasis.mat[m+1,]<-1
      } else{
        use<-which(vars[i,m,]==var)
        knots<-xx.train[cbind(knots.ind[i,m,use],vars[i,m,use])]
        tbasis.mat[m+1,]<-makeBasis(signs[i,m,use],1,knots,xxt,q)
      }
    }
  }
  return(tbasis.mat)
}

## get all the variable combinations used in the models
getCombs<-function(bassMod,uniq.models,nmodels,maxBasis,maxInt.tot,func.var=NULL){
  des.labs<-which(bassMod$cx=='numeric')
  cat.labs<-which(bassMod$cx=='factor')
  func.labs<-letters[0:sum(bassMod$pfunc)]
  labs<-c(des.labs,func.labs,cat.labs) # this is the order things end up in
  vf<-bassMod$vars.func[uniq.models,,]
  sub<-0
  if(!is.null(func.var)){
    vf[vf==func.var]<-NA
    sub=-1
  }
  n.un<-array(c(as.integer(bassMod$vars.des[uniq.models,,]),as.integer(vf+sum(bassMod$pdes)),as.integer(bassMod$vars.cat[uniq.models,,]+sum(bassMod$pdes)+sum(bassMod$pfunc))),dim=c(nmodels,maxBasis,maxInt.tot))
  n.un<-apply(n.un,1:2,sort)
  n.un<-unique(c(n.un))
  n.un[sapply(n.un,length)==0]<-NULL
  temp<-list()
  for(ii in 1:length(n.un)){
    pp<-length(n.un[[ii]])
    if(pp==1){
      temp<-c(temp,n.un[[ii]])
    } else{
      temp<-c(temp,do.call(c,sapply(1:pp,function(x) combn(n.un[[ii]],x,simplify=F))))
    }
  }
  temp<-lapply(temp,as.integer)
  n.un<-unique(c(n.un,unique(temp)))
  ord<-order(sapply(n.un,length))
  n.un<-n.un[ord]
  a<-NA
  for(ii in 1:maxInt.tot)
    a[ii]<-which(sapply(n.un,length)==ii)[1]
  a[maxInt.tot+1]<-length(n.un)+1
  
  combs<-names.ind<-dispNames<-dispCombs<-ll<-mat<-list()
  aa<-a
  a<-na.omit(a)
  k<-0
  for(ii in (1:maxInt.tot)[!is.na(aa[-(maxInt.tot+1)])]){
    k<-k+1
    if(!is.na(a[ii])){
      mat<-do.call(rbind,n.un[a[k]:(a[k+1]-1)])
      mat<-mat[do.call(order, as.data.frame(mat)),,drop=F]
      combs[[ii]]<-t(mat)
      names.ind[[ii]]<-apply(combs[[ii]],2,paste,collapse='x') # labels for output
      dispCombs[[ii]]<-matrix(labs[combs[[ii]]],nrow = nrow(combs[[ii]]))
      dd<-dim(dispCombs[[ii]])
      dispCombs[[ii]]<-apply(dispCombs[[ii]],2,mixSort)
      dim(dispCombs[[ii]])<-dd
      dispNames[[ii]]<-apply(dispCombs[[ii]],2,paste,collapse='x')
      ll[[ii]]<-split(combs[[ii]],rep(1:ncol(combs[[ii]]),each=nrow(combs[[ii]])))
    }
  }

  num.ind<-sapply(combs,ncol)
  cs.num.ind<-cumsum(num.ind) # used for indexing
  return(list(combs=combs,names.ind=names.ind,ll=ll,num.ind=num.ind,cs.num.ind=cs.num.ind,aa=aa,dispCombs=dispCombs,dispNames=dispNames))
}

## process model information into a temporary list
get_tl<-function(bassMod,mcmc.use.mod,M,mod,p,q,cs.num.ind,combs,func.var=NULL,xx.func.var=NULL){
  a<-bassMod$beta[mcmc.use.mod,2:(M+1),drop=F] # basis coefficients excluding intercept
  vf<-bassMod$vars.func[mod,1:M,]
  pfunc<-sum(bassMod$pfunc)
  if(!is.null(func.var)){
    vf[vf==func.var]<-NA
    pfunc<-pfunc-1
  }
  
  if(bassMod$des){
    Kind<-cbind(bassMod$vars.des[mod,1:M,],vf+bassMod$pdes)
  } else if(bassMod$func){
    Kind<-matrix(vf)
  } else{
    Kind<-NA
  }
  
  if(M==1){
    Kind<-t(Kind)
  }
  p.df<-sum(bassMod$pdes)+sum(bassMod$pfunc)
  t<-s<-matrix(0,nrow=M,ncol=sum(bassMod$pdes)+sum(bassMod$pfunc))
  if(p.df>0){
    for(k in 1:M){ # these matrices mimic the output of earth
      if(bassMod$des){
        n.int.des<-bassMod$n.int.des[mod,k]
        knotInd.des<-bassMod$knotInd.des[mod,k,1:n.int.des]
        vars.des<-bassMod$vars.des[mod,k,1:n.int.des]
        t[k,vars.des]<-bassMod$xx.des[cbind(knotInd.des,vars.des)]
        s[k,vars.des]<-bassMod$signs.des[mod,k,1:n.int.des]
      }
      if(pfunc>0){
        if(bassMod$n.int.func[mod,k]>0){
          n.int.func<-bassMod$n.int.func[mod,k]
          knotInd.func<-bassMod$knotInd.func[mod,k,1:n.int.func]
          vars.func<-bassMod$vars.func[mod,k,1:n.int.func]
          t[k,vars.func+sum(bassMod$pdes)]<-bassMod$xx.func[cbind(knotInd.func,vars.func)]
          s[k,vars.func+sum(bassMod$pdes)]<-bassMod$signs.func[mod,k,1:n.int.func]
        }
      }
    }
  }
  if(!is.null(func.var)){
    s[,bassMod$pdes+func.var]<-t[,bassMod$pdes+func.var]<-0
  }
  tl<-list(s=s,t=t,q=q,a=a,M=M,Kind=Kind,cs.num.ind=cs.num.ind,combs=combs,xx=xx.func.var,pfunc=sum(bassMod$pfunc),cat=bassMod$cat,pdes=sum(bassMod$pdes)) #temporary list
  return(tl)
}

## process model information into a temporary list - categorical part
add_tlCat<-function(tl,bassMod,mcmc.use.mod,mod){
  tl$pcat<-bassMod$pcat
  tl$sub.cnt<-matrix(0,nrow=tl$M,ncol=tl$pcat)
  tl$sub<-list()
  for(mm in 1:tl$M){
    vars<-na.omit(bassMod$vars.cat[mod,mm,])
    tl$sub.cnt[mm,vars]<-bassMod$sub.size[mod,mm,bassMod$vars.cat[mod,mm,]%in%vars]/bassMod$nlevels[vars]
    tl$sub[[mm]]<-list()
    for(k in 1:tl$pcat){
      tl$sub[[mm]][[k]]<-NA
      if(k %in% vars)
        tl$sub[[mm]][[k]]<-bassMod$sub.list[[mod]][[mm]][[which(vars==k)]]
    }
  }
  p.df<-sum(bassMod$pdes)+sum(tl$pfunc)
  if(p.df>0)
    tl$Kind<-cbind(tl$Kind,bassMod$vars.cat[mod,1:tl$M,]+sum(bassMod$pdes)+sum(tl$pfunc))
  if(p.df==0)
    tl$Kind<-bassMod$vars.cat[mod,1:tl$M,]
  if(is.null(nrow(tl$Kind)))
    tl$Kind<-matrix(tl$Kind)
  tl$nlevels<-bassMod$nlevels
  return(tl)
}

## process model information into a temporary list - evaluate integrals (from Chen 2004, but vectorized as much as possible)
add_tl<-function(tl,p){
  p.df<-sum(tl$pdes)+sum(tl$pfunc)
  p.use<-p
  if(p.df==0){
    C1.all.cat<-tl$sub.cnt
    C1.all<-C1.all.cat
  } else{
    C1.all<-(1/(tl$q+1)*((tl$s+1)/2-tl$s*tl$t))*tl$s^2 # so I don't need C function anymore
    if(tl$cat){
      C1.all.cat<-tl$sub.cnt
      C1.all<-cbind(C1.all,C1.all.cat)
    }
  }
  
  C1.all2<-replace(C1.all,which(C1.all==0,arr.ind=T),1) # for products, make 0's 1's
  C1.all.prod<-apply(C1.all2,1,prod)
  tl$CC<-tcrossprod(C1.all.prod)
  C2.all<-C1.all.prod2<-both.ind<-array(0,dim=c(tl$M,tl$M,p.use))
  
  for(ii in 1:tl$M){
    for(jj in ii:tl$M){
      bb<-intersect(na.omit(tl$Kind[ii,]),na.omit(tl$Kind[jj,])) # variables that basis functions ii and jj have in common
      if(length(bb)>0){
        C1.all.prod2[ii,jj,]<-C1.all.prod2[jj,ii,]<-C1.all[ii,]*C1.all[jj,] # pairwise products of C1.all
        both.ind[ii,jj,bb]<-both.ind[jj,ii,bb]<-1
        bb.cat<-bb[bb>p.df]
        bb.des<-bb[bb<=p.df]
        temp<-rep(0,p.use)
        if(length(bb.des)>0){
          temp[1:p.df]<-apply(t(1:p.df),2,C2,m=ii,n=jj,tl=tl)
        }
        if(length(bb.cat)>0){
          temp[(p.df+1):p.use]<-apply(t(1:tl$pcat),2,C2Cat,m=ii,n=jj,tl=tl)
        }
        C2.all[ii,jj,]<-C2.all[jj,ii,]<-temp
      }
    }
  }
  tl$C1.all.prod3<-C1.all.prod2
  tl$C1.all.prod3[C1.all.prod2==0]<-1
  tl$C2.all2<-C2.all
  tl$C2.all2[!as.logical(both.ind)]<-1
  return(tl)
}

## sorting function with mixed numerical and character
mixSort<-function(x){
  ind<-is.na(suppressWarnings(as.numeric(x)))
  num<-which(!ind)
  char<-which(ind)
  return(c(sort(as.numeric(x[num])),x[char]))
}

## ordering function with mixed numerical and character
mixOrd<-function(x){
  ind<-is.na(suppressWarnings(as.numeric(x)))
  ord<-1:length(x)
  num<-which(!ind)
  char<-which(ind)
  return(c(ord[num][order(as.numeric(x[num]))],ord[char]))
}


########################################################################
## functions for Sobol decomposition - these all use scaling from const function
########################################################################
## refer to francom 2016 paper
pCoef<-function(i,q){ 
  factorial(q)^2*(-1)^i/(factorial(q-i)*factorial(q+1+i))
}
## integral from a to b of [(x-t1)(x-t2)]^q when q positive integer
intabq<-function(a,b,t1,t2,q){
  sum(pCoef(0:q,q)*(b-t1)^(q-0:q)*(b-t2)^(q+1+0:q)) - sum(pCoef(0:q,q)*(a-t1)^(q-0:q)*(a-t2)^(q+1+0:q))
}
## integral of two pieces of tensor that have same variable - deals with sign, truncation
C2<-function(k,m,n,tl){ 
  q<-tl$q
  t1<-tl$t[n,k]
  s1<-tl$s[n,k]
  t2<-tl$t[m,k]
  s2<-tl$s[m,k]
  cc<-const(signs=c(s1,s2),knots=c(t1,t2),degree=q)
  if((s1*s2)==0){
    return(0)
  }
  if(t2<t1){
    t1<-tl$t[m,k]
    s1<-tl$s[m,k]
    t2<-tl$t[n,k]
    s2<-tl$s[n,k]
  }
  if(m==n){ #t1=t2, s1=s2
    return(1/(2*q+1)*((s1+1)/2-s1*t1)^(2*q+1)/cc)
  } else{
    if(s1==1){
      if(s2==1){
        return(intabq(t2,1,t1,t2,q)/cc)
      } else{
        return(intabq(t1,t2,t1,t2,q)*(-1)^q/cc)
      }
    } else{
      if(s2==1){
        return(0)
      } else{
        return(intabq(0,t1,t1,t2,q)/cc)
      }
    }
  }
}

## same as C2, but categorical
C2Cat<-function(k,m,n,tl){ # k is variable (categorical), m & n are basis functions
  return(length(intersect(tl$sub[[m]][[k]],tl$sub[[n]][[k]]))/tl$nlevels[k])
}

## matrix used in sobol main effect variances - where most of the time is spent
VuMat<-function(u,tl){ 
  CCu<-apply(tl$C1.all.prod3[,,u,drop=F],1:2,prod) # TODO: utilize symmetry
  C2.temp<-apply(tl$C2.all2[,,u,drop=F],1:2,prod)
  mat<-tl$CC*(C2.temp/CCu-1)
  return(mat)
}

## sobol main effect variances
Vu<-function(u,tl){
  mat<-VuMat(u,tl)
  out<-apply(tl$a,1,function(x) t(x)%*%mat%*%x)
  if(any(out<0))
    browser()
  return(out)
}

## functional sobol main effect variances
Vu_des_func<-function(u,tl){ 
  mat<-VuMat(u,tl)
  nx<-length(tl$xx)
  nmodels<-length(tl$a[,1])
  out<-matrix(nrow=nmodels,ncol=nx)
  for(i in 1:nmodels){
    for(j in 1:nx){
      tt<-tl$a[i,]*tl$tfunc.basis[,j]
      out[i,j]<-t(tt)%*%mat%*%tt
    }
  }
  return(out)
}

## sobol interaction variances
VuInt<-function(u,tl){ 
  add<-0
  len<-length(u)
  for(l in 1:len){
    ind<-((sum(tl$cs.num.ind[l-1])+1):tl$cs.num.ind[l])[apply(tl$combs[[l]],2,function(x) all(x%in%u))] # this gets index for which combs are subsets of u, sum(cs.num.ind[l-1]) makes it 0 when it should be
    add<-add+(-1)^(len-l)*rowSums(tl$integrals[,ind,drop=F])
  }
  add[abs(add)<1e-15]<-0
  if(any(add<0))
    browser()
  return(add)
}

## sobol interaction variances - functional
VuInt_des_func<-function(u,tl){
  add<-0
  len<-length(u)
  #browser()
  for(l in 1:len){
    ind<-((sum(tl$cs.num.ind[l-1])+1):tl$cs.num.ind[l])[apply(tl$combs[[l]],2,function(x) all(x%in%u))] # this gets index for which combs are subsets of u, sum(cs.num.ind[l-1]) makes it 0 when it should be 
    add<-add+(-1)^(len-l)*apply(tl$integrals[,ind,,drop=F],c(1,3),sum)
  }
  return(add)
}

