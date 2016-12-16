########################################################################
## miscellaneous functions
########################################################################

## sample a tempered gamma
rgammaTemper<-function(n,shape,rate,temper){ 
  rgamma(n,temper*(shape-1)+1,temper*rate)
}
## sample a tempered IG
rigammaTemper<-function(n,shape,scale,temper){ 
  1/rgamma(n,temper*(shape+1)-1,rate=temper*scale)
}

## scale a vector to be between 0 and 1
scale.range<-function(x,r=NULL){ # x is a vector
  if(is.null(r))
    r<-range(x)
  (x-r[1])/(r[2]-r[1])
}
## rescale a vector between 0 and 1 to range r
unscale.range<-function(x,r){
  x*(r[2]-r[1])+r[1]
}

## get yhat under the different scenarios
getYhat_des<-function(curr,nb){
  curr$des.basis%*%curr$beta
}
getYhat_cat<-function(curr,nb){
  curr$cat.basis%*%curr$beta
}
getYhat_des_cat<-function(curr,nb){
  curr$dc.basis%*%curr$beta
}
getYhat_des_func<-function(curr,nb){
  tcrossprod(curr$des.basis%*%diag(c(curr$beta),nb+1),curr$func.basis)
}
getYhat_cat_func<-function(curr,nb){
  tcrossprod(curr$cat.basis%*%diag(c(curr$beta),nb+1),curr$func.basis)
}
getYhat_des_cat_func<-function(curr,nb){
  tcrossprod(curr$dc.basis%*%diag(c(curr$beta),nb+1),curr$func.basis)
}

## for checking inputs
posInt<-function(x){
  x==as.integer(x) & x>0
}
