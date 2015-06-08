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
   
   # We start the discritization from here 
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








