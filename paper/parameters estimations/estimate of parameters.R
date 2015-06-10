
# In this file we find the parameters of GSSM and nGSSM.  
# Prameters are observation variance of GSSM (V), system variance of GSSM(W), initial
# posterior parameters of GSSM (m0 and c0), and initial sample variance of weights in nGSSM (s0).
# To estimate m0, c0 and s0 we use the real data at the insertion time of piglets into the pen.
# To estimate V and W we use the the simulated data generated from the file "simulation.R".

#'@author Reza Pourmoayed \email{rpourmoayed@econ.au.dk}
#--------------------------------------------------------------------------------------------------

# Estimate of m0, c0 and s0: 
# The data of 7 pens in a standard danish herd are in the file "individualWeighings". 
load("individualWeighings")

# Using the file above, compute the mean and variance of weight in the pen level and store them in the array "weight": 
individualWeighings$Pig = as.integer(individualWeighings$GrisID)
individualWeighings = individualWeighings[,-1]

batches = unique(individualWeighings$OmsaetnDato)

disps = unique(individualWeighings$VentilNr)

weighttime = unique(individualWeighings$Tid)
weight<-array(NA, dim =c(length(batches),length(disps),length(weighttime), 2))

for (b in 1:length(batches)) {
  thisBatch = individualWeighings[individualWeighings$OmsaetnDato == batches[b], ]
  for (d in 1:length(disps)) {
    batchPen = thisBatch[thisBatch$VentilNr == disps[d], ]
    if (length(batchPen$Pig) > 0) {
      for(c in 1:length(weighttime)){
        batchPenwtime = batchPen[batchPen$Tid == weighttime[c], ]   
        if (length(batchPenwtime$Pig) > 0){
          
          weight[b,d,c,]<-c(mean(batchPenwtime$Vaegt),var(batchPenwtime$Vaegt))                  
        }
      }
    }
  }
}

# Extract the mean and variance of weight in the first and second week of growing period
# for all the pens and store them in the list "DL": 
DL<-list()

d<-weight[1,1,,]
DL[[1]]<-subset(d,(d[,1]>0 & d[,2]>0))
DL[[1]]<-rbind(DL[[1]][1,],DL[[1]][2,])

d<-weight[2,1,,]
DL[[2]]<-subset(d,(d[,1]>0 & d[,2]>0))
DL[[2]]<-rbind(DL[[2]][1,],DL[[2]][2,])

d<-weight[3,3,,]
DL[[3]]<-subset(d,(d[,1]>0 & d[,2]>0))
DL[[3]]<-rbind(DL[[3]][1,],DL[[3]][2,])

d<-weight[3,4,,]
DL[[4]]<-subset(d,(d[,1]>0 & d[,2]>0))
DL[[4]]<-rbind(DL[[4]][1,],DL[[4]][2,])

d<-weight[4,3,,]
DL[[5]]<-subset(d,(d[,1]>0 & d[,2]>0))
DL[[5]]<-rbind(DL[[5]][1,],DL[[5]][2,])

d<-weight[4,4,,]
DL[[6]]<-subset(d,(d[,1]>0 & d[,2]>0))
DL[[6]]<-rbind(DL[[6]][1,],DL[[6]][2,])

# Compute the average weekly gain of the first week for all the pens and save them in the list GL:  
GR<-function(weight){
  h<-c()  
  for(i in length(weight):2)
    h[i]<-weight[i]-weight[i-1] 
  h[1]<-h[2]
  return(h)
}

GL<-list()

for(i in 1:6){
  as.vector(DL[[i]][,1])
  GL[[i]]<-GR(as.vector(DL[[i]][,1]))   
}

# Estimate m0 and C0 of GSSM using the lists DL and GL :
m0_1 = mean(c(DL[[1]][1,1], DL[[2]][1,1],DL[[3]][1,1],DL[[4]][1,1],DL[[5]][1,1],DL[[6]][1,1] )) # mean of insertion weight
c0_1 = var(c(DL[[1]][1,1], DL[[2]][1,1],DL[[3]][1,1],DL[[4]][1,1],DL[[5]][1,1],DL[[6]][1,1] )) # variance of insertion weight
m0_2 = mean(c(GL[[1]][1], GL[[2]][1],GL[[3]][1],GL[[4]][1],GL[[5]][1],GL[[6]][1])) # mean of initial growth
c0_2 = var(c(GL[[1]][1], GL[[2]][1],GL[[3]][1],GL[[4]][1],GL[[5]][1],GL[[6]][1])) # variance of initial growth
c0_12 = cov(c(GL[[1]][1], GL[[2]][1],GL[[3]][1],GL[[4]][1],GL[[5]][1],GL[[6]][1]),c(DL[[1]][1,1], DL[[2]][1,1],DL[[3]][1,1],DL[[4]][1,1],DL[[5]][1,1],DL[[6]][1,1] )) # covariance between average insertion weight and initial growth

m0 = matrix(data = c(m0_1,m0_2), ncol=2)
c0 = matrix(data = c(c0_1,c0_12,c0_12,c0_2), ncol = 2)

# Estimate s0 of nGSSM using the avarage sample variance of insertion weighs in the list DL: 
s0 = mean (c(DL[[1]][1,2],DL[[2]][1,2],DL[[3]][1,2],DL[[4]][1,2],DL[[5]][1,2],DL[[6]][1,2]))

#--------------------------------------------------------------------------------------------------------------------------

# Estimation of V in the GSSM. 
# Matrix V includes the observational variance parameters of average weight(V_1), average feed intake(V_2) and
# covariance between average weight and feed intake(V_12) in the pen level. 

# Using the parameters in Kristensen et al. (2012) (doi:10.1016/j.livsci.2012.01.003), observation variance for
# weighting an individual pig is 1 and hence for the the observation variance of average weight in a pen with 15 pigs easily 
# we can conclude that V_1=1/15=0.066
V_1=1/15

# Using the parameters in Jørgensen (1993) (doi:10.1007/s10479-010-0688-z), observation variance 
# related to the daily feed intake of ane individula pig is (0.16)^2 and hence for the the observation variance of average
# weekly feed intake (for 7 days) in a pen with 15 pigs we easily can conclude that V_2 = ( (0.16^2)/15 )*7 = 0.012
V_2 = ( (0.16^2)/15 )*7 

# Using V_1 and V_2, we can compute V_12 by the correlation between average weight and feed intake in the pen level (r).
# To find the correlation r we use the simulated data of weight and feed intake generated from the file "simulation.R". 
# Generate simulated data:
#source("simulation.R")
numIteration<-10
MH<-list()
DL<-list()
corel<-c()
for(i in 1:numIteration){
  param<-setSimParam()
  dat<-simulatePen(param,1,T=11*7)
  iniWeight<-mean(dat$dtDailyPig[t==1]$OLW)
  iniFI<-sum(dat$dtDailyPig[t==1:6]$FI)/15
  DL[[i]]<-cbind(append(iniWeight,dat$dtWeekAve$aveOLWAll), append(iniFI,dat$dtWeekAve$aveFIAll))
  corel[i]<-cor(append(iniWeight,dat$dtWeekAve$aveOLWAll), append(iniFI,dat$dtWeekAve$aveFIAll))
}
r = mean(corel)
# The resulted value of r is 0.975 and now we can find the the value of V_12 with the statistical formula
# of correaltion: V_12=V_21=r*sqrt(V11)*sqrt(V22) = 0.027 
V_12=r*sqrt(V_1)*sqrt(V_2) 

# and hence matrix V is
V = matrix(c(V_1,V_12,V_12,V_2),ncol = 2)

#--------------------------------------------------------------------------------------------------------------------------

#Estimate of W in the GSSM
# To find the System variance of GSSM, we apply the EM algorithm (C. Dethlefsen. Space Time Problems and Applications. PhD thesis, Aalborg University,
# Denmark, 2001) using the simulated data. Simulate data is generated using the the file "simulation.R":
#source("simulation.R")
numIteration<-10
MH<-list()
DL<-list()
corel<c()
for(i in 1:numIteration){
  param<-setSimParam()
  dat<-simulatePen(param,1,T=11*7)
  iniWeight<-mean(dat$dtDailyPig[t==1]$OLW)
  iniFI<-sum(dat$dtDailyPig[t==1:6]$FI)/15
  DL[[i]]<-cbind(append(iniWeight,dat$dtWeekAve$aveOLWAll), append(iniFI,dat$dtWeekAve$aveFIAll))
}
for(j in 1:numIteration){
  D<-array(NA,dim=c(2,1,dim(DL[[j]])[1]))
  for(i in 1:dim(DL[[j]])[1]){
    D[1,1,i]<-DL[[j]][i,1]
    D[2,1,i]<-DL[[j]][i,2]
  }
  MH[[j]]<-D
}
# Next we define a function to set the real parameters of the GSSM except W:

#' Set the parameters of the GSSM. 
#' 
#' @param t Maximumu life time of the pen
#' @param FF design matrix of system equation.
#' @param GG design matrix of observation equation.
#' @param V Observation variance of the nGSSM
#' @param W System variance of the nGSSM (optional values)
#' @param m0 Initial mean of posterior distribution at the insertion time (t=0)
#' @param C0 Initial variance of posterior distribution at the insertion time (t=0)
#' 
#' @return A list containing the parameters of the GSSM. 
#' 
#'@author Reza Pourmoayed \email{rpourmoayed@econ.au.dk}
setModel<-function(t=12,FF,GG,V,W,m0,C0,D){
  model<-list(t=t)
  #model$FF<-matrix(data=c(1,0,0.044,1.549),ncol=2)
  model$GG<-matrix(data=c(1,0,1,1),nrow=2)
  model$V<-matrix(data=c(0.066,0.027,0.027,0.012),ncol=2)
  model$W<-matrix(data=c(0,0,0,0.1),ncol=2)
  model$m0<-matrix(data=c(26.49,5.8),nrow=2)
  model$C0<-matrix(data=c(4.2,0.32,0.32,0.53),ncol=2)
  k<-c(0.12900700, 0.12460995, 0.12071171,0.11713588,0.11390976,0.11098825,0.10836625,0.10598322,0.10382248,0.10185541,0.10007084, 0.09845869)
  FF<-array(NA, dim=c(2,2,t))
  for(z in 1:t){
     FF[1,1,z]<-1
     FF[1,2,z]<-k[z]
     FF[2,1,z]<-0
     FF[2,2,z]<-1.594        
  }
  model$FF<-FF  
  return(model)
}
mod<-setModel()

# EM algorithm needs to the filtering and smooting of GSSM. In below two functions are defined for filtering and smoothging:

#' GSSM filtering  
#' 
#' @param mod set of parameters needed in the GSSM model.
#' @param D Set of the weight and feed intake data generated from the simulation.
#' @param W System variance of the GSSM.
#' @param V Observation variance of the GSSM.
#' 
#' @return Updated information of filtering. 
#' 
#'@author Reza Pourmoayed \email{rpourmoayed@econ.au.dk} 
DLMfilter<-function(mod,D,W,V){
   
   k1<-c()
   
   #The means of posterior for t =1 to t=12
   L1<-array(NA, dim=c(2,1,mod$t))
   
   
   #The variance matrices of posterior for t=1 to t=12
   L2<-array(NA, dim=c(2,2,mod$t))
   
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
      #FF<-mod$FF[,,i] 
      FF<-matrix(data=c(1,0,k1[i],1.549),ncol=2)
      ft<-t(FF) %*% at
      Qt<-t(FF) %*% Rt[,,i] %*% FF + V
      
      #Posterior (we see y_t here)
      At<-Rt[,,i] %*% FF %*% solve(Qt)
      et<-D[,,i]-ft
      L1[,,i]<-at + At %*% et   
      L2[,,i]<-Rt[,,i] - At %*% Qt %*% t(At) 
      
   }
   dlm<-list()
   dlm$L1<-L1
   dlm$L2<-L2
   dlm$Rt<-Rt
   return(dlm)
}

#' GSSM smoothing  
#' 
#' @param mod set of parameters needed in the GSSM model.
#' @param fdlm1 output information of filtering the GSSM using the function DLMfilter.
#' @param W System variance of the GSSM.
#' 
#' @return Updated information of smoothing 
#' 
#'@author Reza Pourmoayed \email{rpourmoayed@econ.au.dk} 
Smoother<-function(mod,fdlm1,W){
   mts<-array(NA,dim=c(2,1,mod$t)) 
   Cts<-array(NA,dim=c(2,2,mod$t))
   Bt<-array(NA,dim=c(2,2,mod$t))
   
   mt<-fdlm1$L1
   Ct<-fdlm1$L2
   Rt<-fdlm1$Rt
      
   mts[,,mod$t]<-mt[,,mod$t]
   Cts[,,mod$t]<-Ct[,,mod$t]
   
   for(i in ((mod$t-1):1)){
      Bt[,,i]<-Ct[,,i] %*% t(mod$GG) %*% solve(Rt[,,i+1])
      mts[,,i]<-mt[,,i] + Bt[,,i] %*% (mts[,,i+1] - mod$GG %*%mt[,,i]) 
      Cts[,,i]<-Ct[,,i] + Bt[,,i] %*% (Cts[,,i+1] - Rt[,,i+1]) %*% t(Bt[,,i])
   }

   #for t=0
   Bt0<-mod$C0 %*% t(mod$GG) %*% solve(mod$GG %*% mod$C0 %*% t(mod$GG) + W)
   mts0<-mod$m0 + Bt0 %*% (mts[,,1] - mod$GG %*%mod$m0)
   Cts0<-mod$C0 + Bt0 %*% (Cts[,,1] - Rt[,,1]) %*% t(Bt0)

   smo<-list()
   smo$mts<-mts
   smo$Cts<-Cts
   smo$mts0<-mts0
   smo$Cts0<-Cts0
   
   return(smo)
}

# Finally the function of EM algorithm should be defined:

#' EM algorithm  
#' 
#' @param mod Set of parameters needed in the GSSM model.
#' @param fdlm1 Output information of filtering the GSSM using the function DLMfilter.
#' @param sdlm1 Output information of smoothing the GSSM using the function Smoother.
#' @param D Set of the weight and feed intake data generated from the simulation.
#' @param Wm Updated system variance of the GSSM.
#' 
#' @return Updated information of system variance W 
#' 
#'@author Reza Pourmoayed \email{rpourmoayed@econ.au.dk} 
EM<-function(mod,fdlm1,sdlm1,D,Wm){
   Bt<-array(NA,dim=c(2,2,mod$t))
   Lt<-array(NA,dim=c(2,2,mod$t))   
   W1<-array(0,dim=c(2,2))
   V1<-array(0,dim=c(2,2))
   result<-list()
   
   mt<-fdlm1$L1
   Ct<-fdlm1$L2
   mts<-sdlm1$mts
   Cts<-sdlm1$Cts
   mts0<-sdlm1$mts0
   Cts0<-sdlm1$Cts0
   
      
   for(i in 1:mod$t){      
      if(i==1){
         Bt0<-mod$C0 %*% t(mod$GG) %*% solve(mod$GG %*% mod$C0 %*% t(mod$GG) + Wm)
         Lt[,,i]<-Cts[,,i] + mod$GG %*% Cts0 %*% t(mod$GG) - Cts[,,i] %*% t(Bt0) - Bt0 %*% t(Cts[,,i])
         W1<-Lt[,,i]+(mts[,,i] - (mod$GG) %*% mts0) %*% t(mts[,,i] - (mod$GG) %*% mts0) + W1
        # V1<-t(mod$FF[,,i]) %*% Cts[,,i] %*% mod$FF[,,i] + (D[,,i]-t(mod$FF[,,i])%*%mts[,,i])%*%t((D[,,i]-t(mod$FF[,,i])%*%mts[,,i])) + V1
      }else{
         Bt[,,i-1]<-Ct[,,i-1] %*% t(mod$GG) %*% solve(mod$GG %*% Ct[,,i-1] %*% t(mod$GG) + Wm)
         Lt[,,i]<-Cts[,,i] + mod$GG %*% Cts[,,i-1] %*% t(mod$GG) - Cts[,,i] %*% t(Bt[,,i-1]) - Bt[,,i-1] %*% t(Cts[,,i])    
         W1<-Lt[,,i]+(mts[,,i] - (mod$GG) %*% mts[,,i-1]) %*% t(mts[,,i] - (mod$GG) %*% mts[,,i-1]) + W1
       #  V1<-t(mod$FF[,,i]) %*% Cts[,,i] %*% mod$FF[,,i] + (D[,,i]-t(mod$FF[,,i])%*%mts[,,i])%*%t((D[,,i]-t(mod$FF[,,i])%*%mts[,,i])) + V1
      }      
   }
   W1<-(W1+t(W1))/2
   #V1<-(V1+t(V1))/2
   W1<-W1/mod$t
   #V1<-V1/mod$t
   result[[1]]<-W1
   #result[[2]]<-V1
   return(result)
}
 
# No we implement the EM algorithm for z=1000 iterations : 
z<-1000
fdlm<-list()
sdlm<-list()
We<-array(NA,dim=c(2,2,z))
Ve<-array(NA,dim=c(2,2,z))
We1<-array(0,dim=c(2,2))
#Ve1<-array(0,dim=c(2,2))

for(h in 1:numIteration){
   for(u in 1:z){  
      if(u==1){      
         fdlm[[u]]<-DLMfilter(mod,MH[[h]],mod$W,mod$V)
         sdlm[[u]]<-Smoother(mod,fdlm[[u]],mod$W)
         We[,,u]<-EM(mod,fdlm[[u]],sdlm[[u]],MH[[h]],mod$W)[[1]]
         #Ve[,,u]<-EM(mod,fdlm[[u]],sdlm[[u]],MH[[h]],mod$W)[[2]]      
      }
      else{
         fdlm[[u]]<-DLMfilter(mod,MH[[h]],We[,,u-1],mod$V)
         sdlm[[u]]<-Smoother(mod,fdlm[[u]],We[,,u-1])
         We[,,u]<-EM(mod,fdlm[[u]],sdlm[[u]],MH[[h]],We[,,u-1])[[1]]
         #Ve[,,u]<-EM(mod,fdlm[[u]],sdlm[[u]],MH[[h]],We[,,u-1])[[2]]
      }  
   }
   We1<-We[,,z] + We1
#   Ve1<-Ve[,,z] + Ve1   
}

We1<-We1/numIteration
#Ve1<-Ve1/numIteration


