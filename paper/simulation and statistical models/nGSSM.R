
#' Set the parameters of the nGSSM
#' 
#' @param t Maximumu life time of the pen
#' @param nf Number of weight observations
#' @param W System variance of the nGSSM
#' @param c0 Initial shape parameter of the prior distribution (Inv-Gamma distribution)
#' @param d0 Initial scale parameter of the prior distribution (Inv-Gamma distribution)
#' 
#' @return A list containing the recquired parameters. 
#' 
#'@author Reza Pourmoayed \email{rpourmoayed@econ.au.dk}
setmodel2<-function(t=12, nf, W, c0, d0){
  
  model<-list(t=t)  
  model$nf=35  # sampel size
  model$W=0    # system variance 
  model$c0<-(model$nf-1)/2  # initial shape parameter
  model$s0<-7.65 # initial sample variance in the pen at insertion time
  model$d0<-model$c0*model$s0 #initial scale parameter 
  model$m0= model$d0 / (model$c0 -1) # initial posterior mean 
  model$C0=  model$m0^2 / (model$c0 -2) # initial posterior variance
  g<-c()
  g[1]<-1
  for(j in 2:t){
          
     if(j==2 || j==3){
        g[j]=sqrt(1.8)
     }else{
        g[j]=j/(j-1)
     }  
  }  
  model$g<-g^2  
  return(model)
}

mod1<-setmodel2()

#' nGSSM filtering to estimate the weight variances during the growing period. 
#' 
#' @param mod1 set of parameters of the nGSSM
#' @param ob set of observations related to the sample variances of the weights during the growing period. 
#' @param W System variance of the nGSSM
#'
#' @return Updated information of posterior. 
#' 
#'@author Reza Pourmoayed \email{rpourmoayed@econ.au.dk}
DGLMfilter<-function(mod1,ob,W){
  mtt<-c()
  ctt<-c()
  Rtt<-c()
  error<-0
  
  for(i in 1:mod1$t){
    #Prior
    if(i==1){
      att<-(1) * mod1$m0
      Rtt[i]<-(1)^2 * mod1$C0 + W
    }else{
      att<- mod1$g[i]* mtt[i-1]
      Rtt[i]<-mod1$g[i]^2 * ctt[i-1] + W
    }
    ftt<-att
    error<-(ob[i]-ftt)^2+error
    qtt<-Rtt[i]
    alpha<-((ftt^3)/qtt) + ftt
    beta<-((ftt^2)/qtt) + 1
    alpha1<-alpha + ((mod1$nf-1)/2)*ob[i]
    beta1<-beta + ((mod1$nf-1)/2)
    ftt1<-(alpha1/(beta1))  
    qtt1<-(ftt1^2)/(beta1-1)
    mtt[i]<-ftt1
    ctt[i]<-qtt1  
  }
  
  dglm<-list()
  dglm$mtt<-mtt
  dglm$ctt<-ctt
  dglm$Rtt<-Rtt
  dglm$error<-error
  return(dglm)
}

   
