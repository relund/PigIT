#' Set the parameters used when build the HMDP model
#' 
#' @param tMax Maximum growth period in weeks.
#' @param tStartMarketing Week where start to consider marketing.
#' @param phases Maximum number of phases (possible changes in the ration).
#' @param pigs Number of pigs in a pen at insertion.
#' @param rations Number of possible rations.
#' @param rationCost A vector with ration costs (per feed unit) of length rations.
#' @param minPhaseT  A vector with minimum number of weeks used for each phase of feeding. (Used for three-level hmdp)
#' @param thresholds A vector of thresholds (in kg) used by the actions. The first threshold must be zero, i.e. the term. action. It must be increasing!
#' @param convRate Conversion rate between live weight and carcass weight  
#' @param convRateSd Standard deviation of the conversion rate of the live weight and the carcass weight 
#' @param priorGrowth Prior growth rate for each feed mix.
#' @param avgGRate Average growth rate in the herd.
#' @param avgLeanP Average leanness percent in the herd 
#' @param sDMeasure Std. dev. on measurements (vision).
#' @param CovWG Covariance between live weight and growth rate (we assume it is a constant value) 
#' @param avgInsWeight Avarage weight of the pigs at the insertion time into the pen  
#' @param iniFeedMix The first feed mix for starting the system.
#' @param disSD A vector with number of center points on each side of the center value and the length of each interval (standard deviation).
#' @param disWeight A vector with number of center points on each side of the center value and the length of each interval (weight)
#' @param disGrowth A vector with number of center points on each side of the center value and the length of each interval (growth)
#' @param centerWeight A vector with center values (prior mean weight at time t) (one for each stage 1,...,tMax).
#' @param centerSD A vector with center values (prior sd at time t) (one for each stage 1,...,tMax).
#' @param centerGrowth A vector with center values (prior growth) (one for each feed mix).
#' @param cullActions True is use cull actions instead of threshold actions.
#' @param samples Number of samples to use in the simulation.
#' @param check Check model e.g. do trans pr sum to one
#' 
#' @return A list containing all the parameters
#' @export 
#' @author Reza Pourmoayed \email{rpourmoayed@@econ.au.dk} and Lars Relund \email{lars@@relund.dk}
setParameters<-function(tMax=12,
                        tStartMarketing=9,
                        phases=4,
                        pigs=15,
                        rations=4,
                        iniFeedMix=2, 
                        rationCost=c(1.8,1.88,1.96), #c(0.5,0.5,0.5,0.5),#c(1.2,1.25,1.3,1.35), #c(1.65,1.7,1.75,1.8), #(1.2,1.25,1.3,1.35), #c(1.2,1.2,1.2,1.2), #c(1.85,1.9,1.95,2), #c(1.2,1.2,1.2,1.2), #c(1.2,1.7,2.5,2.8)
                        minPhaseT=c(1,4,7,10),
                        thresholds=c(0,seq(94,116,by=3),200),     # MUST be incresing!!
                        convRate=0.84,
                        convRateSd=1.4,
                        avgGRate = 6,
                        priorGrowth=seq(5.5, by=0.5, length=rations), #c(5.5, 6, 6.5, 7)
                        avgLeanP=61,
                        CovWG=8.66,
                        avgInsWeight=26,
                        sDMeasure=1,
                        disWeight=c(5,2), 
                        disSD=c(2,1),     
                        disGrowth=c(5,2), 
                        centerWeight=seq(26.4,130, by=7),#seq(avgInsWeight,avgInsWeight+(tMax-1)*avgGRate,by=avgGRate), 
                        centerSD=seq(4,11,by=0.5),#seq(4,by=0.5,length=tMax),
                        centerGrowth=priorGrowth, #c(5.2, 5.7, 6.2, 6.7, 6.8),
                        iniDist=c(26.49,4.26),
                        cullActions = TRUE,
                        samples = 1000,
                        check = TRUE
){
   model<-list(tMax=tMax)   
   model$tStartMarketing=tStartMarketing  
   model$phases<-phases
   model$pigs<-pigs  
   model$rations<-rations
   model$iniFeedMix<-iniFeedMix
   model$rationCost<-rationCost
   model$minPhaseT<-minPhaseT
   model$thresholds<-thresholds
   model$convRate<-convRate
   model$convRateSd<-convRateSd
   model$avgGRate<-avgGRate
   model$priorGrowth<-priorGrowth
   model$avgLeanP<-avgLeanP
   model$CovWG<-CovWG
   model$avgInsWeight<-avgInsWeight
   model$sDMeasure<-sDMeasure
#    model$disWeight<-disWeight
#    model$disSD<-disSD
#    model$disGrowth<-disGrowth        
#    model$centerWeights<-centerWeight
#    model$centerSD<-centerSD
#    model$centerGrowth<-centerGrowth
   model$iniDist<-iniDist                    
   model$cullActions <- cullActions
   model$samples <- samples
   model$check <- check
   
   require(discretizeGaussian)          # We start the discritization from here 
   obj<-Discretize()
   pointsWeight<-c()
   pointsSd<-c()
   wIntervals<-vector("list", tMax)    # matrices of the discritized intervals at the stages of the hmdp (weight) 
   sDIntervals<-vector("list", tMax)   # matrices of the discritized intervals at the stages of the hmdp (standard deviation)
   gIntervals<-vector("list", rations)
   for(t in 1:tMax){
#      if(t<9){
         pointsWeight<-PointsW(t,disWeight[1],disWeight[2],centerSD,centerWeight)
         pointsSd<-PointsV(t,disSD[1],disSD[2],centerSD)   
         wIntervals[[t]]<-as.matrix(obj$discretize1DVec(pointsWeight, inf=1000, asDF=F), ncol=3) 
         sDIntervals[[t]]<-as.matrix(obj$discretize1DVec(pointsSd, inf=100, mInf=0.01, asDF=F), ncol=3) 
#      }
#      if(t>8){
#         pointsWeight<-seq(70, 120, by=5)
#         pointsSd<-seq(5.5,12, by=1)   
#         wIntervals[[t]]<-as.matrix(obj$discretize1DVec(pointsWeight, inf=1000, asDF=F), ncol=3) 
#         sDIntervals[[t]]<-as.matrix(obj$discretize1DVec(pointsSd, inf=100, mInf=0.01, asDF=F), ncol=3)          
#      }
      #gIntervals[[t]]<-as.matrix(obj$discretize1DVec(centerPointsGrowth, inf=1000, asDF=F), ncol=3)
   }
   for(ration in 1:rations){
      #pointsGrowth<-PointsG(ration,disGrowth[1],disGrowth[2],centerGrowth)
      pointsGrowth<-c(4.4,5.8,6.3,6.8,8.2)
      gIntervals[[ration]]<-as.matrix(obj$discretize1DVec(pointsGrowth, inf=1000, asDF=F), ncol=3)
   }
   
   model$wIntervals<-wIntervals
   model$sDIntervals<-sDIntervals      
   model$gIntervals<-gIntervals
   model$rGIdx<-findInterval(rations,priorGrowth,gIntervals)
   #model$sW<-obj$discretize1DVec(centerPointsWeight, inf=1000, asDF=F)
   #model$sG<-obj$discretize1DVec(centerPointsGrowth, inf=1000, asDF=F)
   #model$sSd<-obj$discretize1DVec(centerPointsSd, inf=1000, mInf=0, asDF=F)
   
   return(model)
}

####################################################################

#' Find the intervals the priorGrowth belongs to and return the index (c++ style)
#' 
#' @param priorGrowth Prior growth estimates for each feed mix.
#' @param gIntervals Discretization list.
#' 
#' @author Lars Relund \email{lars@@relund.dk}
findInterval<-function(rations,priorGrowth, gIntervals) {
   res<-rep(NA,rations)
   for (i in 1:rations) {
      mat<-gIntervals[[i]]
      for (r in 1:nrow(mat)) {
         if (priorGrowth[i]>mat[r,2] & priorGrowth[i]<=mat[r,3]) {res[i]<-r-1; next}
      }
   }
   return(res)
}


####################################################################


#' Initialize the parameters for the DLM used for a specific ration.
#' 
#' @param tMax Maximum growth period in weeks.
#' @param k2 Value of \eqn{k_2} in the \eqn{F} matrix.
#' @param FF The \eqn{F} matrix.
#' @param GG The \eqn{G} matrix.
#' @param V The \eqn{V} matrix.
#' @param W The \eqn{W} matrix.
#' @param m0 Prior latent mean \eqn{m_0}.
#' @param C0 Prior latent covariance matrix \eqn{C_0}.
#' @return A list containing all the parameters
#' @export
#' @author Reza Pourmoayed \email{rpourmoayed@@econ.au.dk} and Lars Relund \email{lars@@relund.dk}
iniDLM<-function(tMax=15,
                   k2=1.549,
                   FF=matrix(data=c(1,0,0.044,k2),ncol=2), 
                   GG=matrix(c(1,0,1,1),nrow=2), 
                   V=matrix(c(0.066,0.027,0.027,0.012),ncol=2),#matrix(data=c(0.066,0.00013,0.00013,0.0019),ncol=2),#matrix(c(0.25,0.013,0.013,0.73),ncol=2), 
                   W=matrix(c(2.1,-0.124,-0.124,0.112),ncol=2),#matrix(c(0,0,0,0.112),ncol=2),#matrix(data=c(0,0,0,0.0066),ncol=2),#matrix(c(2.194,1.167,1.167,1.601),ncol=2), 
                   m0=matrix(c(26.49,6),nrow=2), 
                   C0=matrix(c(4.26,0.32,0.32,0.53),ncol=2)
){
   model<-list(t=t)
   model$k2<-k2
   model$FF<-t(FF)    # note FF in dlm package is actually F' in W&H book!
   model$GG<-GG
   model$V<-V
   model$W<-W
   model$m0<-m0
   model$C0<-C0
   return(model)
}


####################################################################


#' Estimate values of \eqn{k_{1,t}} in \eqn{F_t} of the DLM and calculate the
#' covariance matrices to be used when have to calculate the transition
#' probabilities.
#' 
#' Estimate using formula \eqn{k_{1,t} = k_1*a_t^(-0.25)} where \eqn{a_t} is the prior at time \eqn{t}. 
#' 
#' @param iniDLM The DLM created using \link{iniDLM}.
#' @param Y A matrix with 2 observation rows (weight and intake). One column for each week. Used to estimate
#'   \eqn{k_{1,t}}. Could be mean values observed in the herd.
#' @param k1 Value of \eqn{k_1} used in the estimation.
#' @return A list containing all the parameters. Element L contains the covariance matrix for transition probabilities. 
#'   Element JFF is a (2x2) matrix indicating where \eqn{F_t'} depends on \eqn{t} (i.e. a 1 is put where \eqn{k_t}) is
#'   included). Element X contains the \eqn{k_t} values. See the vignette (\code{vignette("dlm")}) in the \pkg{dlm}
#'   package for further info about the specification.
#' @export
#' @author Reza Pourmoayed \email{rpourmoayed@@econ.au.dk} and Lars Relund \email{lars@@relund.dk}
buildDLM<-function(iniDLM,Y,k1=0.044) {
   tMax<-dim(Y)[2]
   k<-rep(NA,tMax)
   m<-vector("list", tMax)    # means of posterior - preallocate an empty list of length tMax
   C<-vector("list", tMax)    # covariance matrices of posterior 
   L<-vector("list", tMax)    # covariance matrices of transition probabilities
   FF<-t(iniDLM$FF)           # F def used in W&H book
   for(t in 1:tMax){
      #Prior
      if(t==1){
         a<-iniDLM$GG %*% iniDLM$m0
         R<-iniDLM$GG %*% iniDLM$C0 %*% t(iniDLM$GG) + iniDLM$W
      }
      else{
         a<-iniDLM$GG %*% m[[t-1]] 
         R <-iniDLM$GG %*% C[[t-1]]  %*% t(iniDLM$GG) + iniDLM$W
      }
      # One step forcast
      k[t]<-7*k1/(a[1]^0.25)
      FF[1,2]<-k[t]
      f<-t(FF) %*% a
      Q<-t(FF) %*% R %*% FF + iniDLM$V
      #Posterior (we see Y_t here)
      A<-R %*% FF %*% solve(Q)
      m[[t]] <-a + A %*% (Y[,t]-f)   
      C[[t]] <-R  - A %*% Q %*% t(A) 
      #Compute the covariance matrix for transition probabilities
      L[[t]] <-R-C[[t]] 
   }
   #L[[1]]<-L[[2]]
   
   #iniDLM$m<-m
   iniDLM$C<-C
   #iniDLM$k<-k
   iniDLM$L<-L
   iniDLM$JFF<-matrix(c(0,1,0,0), nrow=2) # F' time dependent in entry (2,1) (dlm package specification way)
   iniDLM$X<-matrix(k,ncol=1) # value of k_t (dlm package specification way)
   return(iniDLM)
}


####################################################################


#' Initialize the parameters for the DGLM 
#' 
#' @param tMax Maximum growth period in weeks.
#' @param nf Sample size to estimate the sample variance in the pen (we suppose this value is constant and the sample is taken from the herd)
#' @param W System variance (we suppose this value is 0 in the DGLM)
#' @param alpha0 Initial shape parameter at time t_0
#' @param beta0 Initial rate parameter at time t_0
#' @param c0 coefficient of G_t in the dglm model (Based on the simulation test, this value is 1) 
#' 
#' @return A list containing all the parameters
#' @export
#' @author Reza Pourmoayed \email{rpourmoayed@@econ.au.dk} and Lars Relund \email{lars@@relund.dk}
iniDGLM<-function(tMax=15, nf=35, W=0, alpha0=130, beta0=17, c0=1) {
   model<-list(tMax=tMax)  
   model$nf=35
   model$W=0 #? I consider this number as zero  based on the Anders paper and talking with him. 
   model$alpha0=130
   model$beta0=17
   model$co=1
   return(model)
}


####################################################################

#' Descritize the state space of the HMDP based on the number of the stage: Subfunction : Find the centerpoints for the live weight
#' and stansard deviation  
#' 
#'@param Week The current stage (week) in the production system
#'@param length The length of the related intervals in the current styage
#'
#'@return A vector included the center points to discritize the state space in the current stage       
#'@author Reza Pourmoayed \email{rpourmoayed@@econ.au.dk}
PointsW<-function(week, numPoints, length, sdWeights, meanWeights){  
   centerPoints<-c()
   x<-2*numPoints+1; 
   for(i in 1:numPoints){
      centerPoints[i]= meanWeights[week] - length*( numPoints - (i-1) );
      centerPoints[x-(i-1)]= meanWeights[week] + length*( numPoints - (i-1) );
   }
   centerPoints[numPoints+1]= meanWeights[week];
   return(centerPoints)
}


PointsV<-function(week, numPoints, length, sdWeights){  
   centerPoints<-c()
   x<-2*numPoints+1; 
   for(i in 1:numPoints){
      centerPoints[i]= sdWeights[week] - length*( numPoints - (i-1) );
      centerPoints[x-(i-1)]= sdWeights[week] + length*( numPoints - (i-1) );
   }
   centerPoints[numPoints+1]= sdWeights[week];
   return(centerPoints)
}


PointsG<-function(ration, numPoints, length, growthWeights){  
   centerPoints<-c()
   x<-2*numPoints+1; 
   for(i in 1:numPoints){
      centerPoints[i]= growthWeights[ration] - length*( numPoints - (i-1) );
      centerPoints[x-(i-1)]= growthWeights[ration] + length*( numPoints - (i-1) );
   }
   centerPoints[numPoints+1]= growthWeights[ration];
   return(centerPoints)
}








