
#global parameters
setModel<-function(t=12,FF,GG,V,W,m0,C0,D){
   model<-list(t=t)
   model$FF<-matrix(data=c(1,0,0.044,1.549),ncol=2)
   model$GG<-matrix(data=c(1,0,1,1),nrow=2)
   model$V<-matrix(c(0.066,0.027,0.027,0.012),ncol=2)
   model$W<-matrix(data=c(0,0,0,0.12),ncol=2) # matrix(c(2.1,-0.124,-0.124,0.112),ncol=2)
   model$m0<-matrix(data=c(26.49,6),nrow=2)
   model$C0<-matrix(data=c(4.26,0.32,0.32,0.53),ncol=2)
   return(model)
}

mod<-setModel()


#DLM filtering 
DLMfilter<-function(mod,D,W){
   
   k1<-c()
   
   #The means of posterior for t =1 to t=12
   L1<-array(NA, dim=c(2,1,mod$t))
   
   
   #The variance matrices of posterior for t=1 to t=12
   L2<-array(NA, dim=c(2,2,mod$t))
   
   #The variance matrices fir the bivariate normal 
   L3<-array(NA, dim=c(2,2,mod$t))
   
   #The variance matrices of prior for t=1 to t=12
   Rt<-array(NA, dim=c(2,2,mod$t))
   
   for(i in 1:mod$t){
      #Prior
      if(i==1){
         at<-mod$GG %*% mod$m0
         Rt[,,i]<-mod$GG %*% mod$C0 %*% t(mod$GG) + W
      }
      else{
         at<-mod$GG %*% L1[,,i-1] 
         Rt[,,i]<-mod$GG %*% L2[,,i-1] %*% t(mod$GG) + W
      }
      # One step forcast
      k1[i]<-7*0.044/(at[1]^0.25)
      FF<-matrix(data=c(1,0,k1[i],1.549),ncol=2)
      ft<-t(FF) %*% at
      Qt<-t(FF) %*% Rt[,,i] %*% FF + mod$V
      
      #Posterior (we see y_t here)
      At<-Rt[,,i] %*% FF %*% solve(Qt)
      et<-D[,,i]-ft
      L1[,,i]<-at + At %*% et   
      L2[,,i]<-Rt[,,i] - At %*% Qt %*% t(At) 
      L3[,,i]<-Rt[,,i] - L2[,,i] 
   }
   dlm<-list()
   dlm$L1<-L1
   dlm$L2<-L2
   dlm$L3<-L3
   dlm$Rt<-Rt
   return(dlm)
}



