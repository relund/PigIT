
#DGLM mod1el for the variance companants 

setmodel2<-function(t=15,nf,W, c0, d0){
  
  model<-list(t=t)  
  model$nf=35  # sampel size
  model$W=0    # system variance 
  model$c0<-(model$nf-1)/2  # initial shape parameter 
  model$d0<-model$c0*7.65 #initial scale parameter 
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

   
