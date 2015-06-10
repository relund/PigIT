# --------------------------------------------------------------------------------------------- #

#' Create a simulation of one pig on a daily basis.
#' The function uses the model described in the appendix of Jørgensen (1993) 
#' 'The Influence of Weighing Precision on Delivery Decisions in Slaughter Pig Production'
#'
#' @param pig The id to be used of the pig (integer). Usefull when merge multiple pigs into one data table.
#' @param measurementStd The std of the weight sensor (e.g vision).
#' @param T Last day to simulate
#' @param k4Mean Mean of growth rate parameter (Gompertz growth curve).
#' @param TLW Simulation start weight (day before start simulation t=0).
#' @param DFI Days from insertion into the pen (start simulation at t=DFI+1).
#' @param DGFactor Increase the daily gain (x2) with this factor. 
#'        Can be used to indicate genetics of higher growth.
#' @param culled Do we simulate a culled pig (true(1)/false(0)).
#' 
#' @return A \link(data.table) with columns \code{t} (days since insertion into the pen), 
#' \code{pig} (pig id), \code{FI} (feed intake in feed units), \code{DG} (daily gain), 
#' \code{TLW} (true live weight),
#' \code{OLW} (observed live weight), \code{SW} (slaughter weight),
#' \code{culled} (are we simulating a culled pig (=0, alive=1).
#' 
#' @author Lars Relund \email{lars@@relund.dk}
#' @seealso createDailySampleData.
#' @examples dat<-samplePig()
#' plotPigData(dat)
#' @export 
samplePig<-function(pig=1,measurementStd=2,T=7*11,k4Mean=0.0166, TLW=rnorm(1,30,measurementStd), DFI=0, DGFactor=1, culled=0) { # FI (x1), DG (x2), FC (x3), TLW (x4), OLW (x5), SW (x6)
  # parameters from Jørgensen93
  k1<-c(0.044,0.002)
  k2<-c(1.549,0.005)
  k3<-c(5.30,0.2667)
  k4=c(DGFactor*k4Mean,0.00005)
  k5<-c(0.763,0.01)
  mu<-c(k1[1],k2[1],k3[1],k4[1],k5[1])
  sigma<-diag(c(k1[2]^2,k2[2]^2,k3[2]^2,k4[2]^2,k5[2]^2))
  #sigma[3,4]<-sigma[4,3] <- -0.0431  # seems not to the pos-semi-def??
  eps1<-0.16
  eps2<-0.02
  eps5<-measurementStd
  eps6<-1.0
  # initialize
  dat<-data.table(t=(DFI+1):T)  
  T<-length(dat$t) # days to simulate
  dat[,`:=`(pig=as.integer(pig), FI=0, DG=0, TLW=0, OLW=0, SW=0)]
  k<-mvrnorm(1, mu, sigma)  # parameters for this pig
  #TLW <- max(0,rnorm(1,OLW,eps5))   # weight at simulation start t=DFI
  dat$FI[1] <- k[2]*k[4]*(k[3]-log(TLW))*TLW + k[1]*TLW^0.75 + rnorm(1,0,eps1)
  dat$DG[1] <- ( (dat$FI[1] - k[1]*TLW^0.75)/k[2] ) + rnorm(1,0,eps2)
  #dat$FC[1] <- dat$FI[1] + FC
  dat$TLW[1]<- max(0,dat$DG[1] + TLW)
  dat$OLW[1]<- dat$TLW[1]+rnorm(1,0,eps5)
  dat$SW[1] <- k[5]*dat$TLW[1] + rnorm(1,0,eps6)
  # simulate
  for (i in 2:T) {
    dat$FI[i] <- k[2]*k[4]*(k[3]-log(dat$TLW[i-1]))*dat$TLW[i-1] + k[1]*dat$TLW[i-1]^0.75 + rnorm(1,0,eps1)
    dat$DG[i] <- ( (dat$FI[i] - k[1]*dat$TLW[i-1]^0.75)/k[2] ) + rnorm(1,0,eps2)
    #dat$FC[i] <- dat$FI[i] + dat$FC[i-1]
    dat$TLW[i]<- max(0,dat$DG[i] + dat$TLW[i-1])
    dat$OLW[i]<- dat$TLW[i] + rnorm(1,0,eps5)
    dat$SW[i] <- k[5]*dat$TLW[i] + rnorm(1,0,eps6)
  }
  dat$culled<-culled
  return(dat)
}

# --------------------------------------------------------------------------------------------- #

#' Fit parameter k4 as a function of average daily growth (over a period of 84 days).
#' Used to find different k4's for different feed-mix.
#' Based on the appendix of Jørgensen (1993). 
#' 
#' @param ite Number of iterations (pigs used in the estimation per k4 value).
#' 
#' @return A data frame with k4 values given a set of daily gains 
#' @author Lars Relund \email{lars@@relund.dk}
#' @seealso samplePig.
#' @examples 
#' dat<-estimateK4()
#' plot(dat)
#' @export 
estimateK4<-function(ite=500) {
  k1<-c(0.044,0.002)
  k2<-c(1.549,0.005)
  k3<-c(5.30,0.2667)
  k5<-c(0.763,0.01)
  eps1<-0.16
  eps2<-0.02
  T<-7*12
  df<-data.frame(k4=NA, aveDG=NA)
  DG<-rep(NA,T)
  for (k4Mean in seq(0.005,0.02,by=0.0025) ) {
    k4=c(k4Mean,0.00005)
    mu<-c(k1[1],k2[1],k3[1],k4[1],k5[1])
    sigma<-diag(c(k1[2]^2,k2[2]^2,k3[2]^2,k4[2]^2,k5[2]^2))
    for (i in 1:ite) {
      k<-mvrnorm(1, mu, sigma)  # parameters for this pig
      TLW <- 30 
      for (i in 1:T) {
        FI<-k[2]*k[4]*(k[3]-log(TLW))*TLW + k[1]*TLW^0.75 + rnorm(1,0,eps1)
        DG[i] <- (FI - k[1]*TLW^0.75)/k[2]  + rnorm(1,0,eps2)
        TLW<- DG[1] + TLW
      }  
      df <- rbind(df,c(k4Mean,mean(DG)) )
    }
  }
  df<-df[-1,]
  lm.df<-lm(k4~aveDG+I(aveDG^2)+I(aveDG^3),data=df)
  #return(lm.df)
  aveDG<-seq(0.5,2,by=0.1)
  dat<-data.frame(aveDG, I(aveDG^2), I(aveDG^3))
  k4<-predict(lm.df,dat)
  return(data.frame(aveDG=dat$aveDG,k4))
}

# --------------------------------------------------------------------------------------------- #

#' Plot simulation data. 
#'
#' @param dat Simulation data obtained from \link{samplePig}.
#' @return NULL (a plot is output). 
#' @author Lars Relund \email{lars@@relund.dk}
#' @seealso createSampleData.
#' @examples plotPigData(dat)
#' @export 
plotPigData<-function(dat) {
  par(mfrow=c(2,2))
  lo <- loess(dat$OLW~dat$t)
  plot(dat$t,dat$OLW)
  lines(dat$t,predict(lo), col='red', lwd=1)
  lo <- loess(dat$TLW~dat$t)
  plot(dat$t,dat$TLW)
  lines(dat$t,predict(lo), col='red', lwd=1)
  lo <- loess(dat$DG~dat$t)
  plot(dat$t,dat$DG)
  lines(dat$t,predict(lo), col='red', lwd=1)
  lo <- loess(dat$FI~dat$t)
  plot(dat$t,dat$FI)
  lines(dat$t,predict(lo), col='red', lwd=1)
}
# dat1<-samplePig(TLW=80,T=110)
# plotPigData(dat1)

# dat<-samplePig(TLW=5,T=200)
# plotPigData(dat)

# --------------------------------------------------------------------------------------------- #

#' Create a data set simulating a pen from DFI and T days ahead. 
#' 
#' @param pen Pen id added to the output.
#' @param pigIds Id of the pigs simulated
#' @param TLW Start weight day before simulate (same length as pigIds).
#' @param culled Have the pig we simulate been culled (same length as pigIds).
#' @param measurementStd The std of the weight sensor (normally vision).
#' @param k4Mean Mean of growth rate parameter (Gompertz growth curve).
#' @param DGFactor Increase the daily gain (x2) with this factor.
#' @param T Number of days to simulate.
#' @param DFI Days from insertion into the pen (start simulation at t=DFI+1).
#' @param feedMix Feed-mix number used.
#' 
#' @return A data table containing a row for each pig on a daily basis with columns \code{t} (days
#'   since insertion into the pen), \code{pig} (pig id), \code{FI} (feed intake in feed units),
#'   \code{DG} (daily gain), \code{TLW} (true live weight), \code{OLW} (observed live weight),
#'   \code{SW} (slaughter weight), \code{culled} (are we simulating a culled pig (=0, alive=1),
#'   \code{week} (week number), \code{weekDay} (week day - monday = 1), \code{pen} (pen id),
#'   \code{feedMix} (feed-mix number used).
#' 
#' @author Lars Relund \email{lars@@relund.dk}
#' @examples dat<-createDailySampleData()
#' @export 
createDailySampleData<-function(pen,pigIds, TLW,
                                culled,
                                measurementStd,
                                k4Mean,
                                DGFactor,
                                T,
                                DFI, feedMix
) {
  ## create data for all pigs
  pigs<-length(pigIds)
  datPigs<-lapply(1:pigs,FUN=function(p) samplePig(pig=pigIds[p],measurementStd,k4Mean=k4Mean,DGFactor=DGFactor,TLW=TLW[p], culled=culled[p], T = T, DFI = DFI))
  l<-length(datPigs[[1]]$t) # timeperiods considered
  dtPigs<-rbindlist(datPigs)  # join in one DT
  dtPigs[,week:=(t-1)%/%7+1] 
  dtPigs[,weekDay:=(t-1)%%7+1]
  dtPigs[,pen:=pen]
  dtPigs[,feedMix:=feedMix]
  return(dtPigs)
}

# --------------------------------------------------------------------------------------------- #

#' Calc a weekly data set based on the daily samples. We assume that actions are
#' considered on Mondays (in the morning), i.e week w in the HMDP corresponds week w-1 and
#' weight data to data from sunday the week before!
#' 
#' @param measurementsPerDay Measurements per day from the weigthing sensor (normally vision).
#' @param measurementStd The std of the weight sensor (normally vision).
#' 
#' @return A data table with a row for each week with columns \code{aveOLW<X>} (ave observed
#'   liveweight at t), \code{sdOLW<X>} (sd observed live weight at t), \code{aveFI<X>} (ave feed
#'   intake per pig per week, monday-sunday, we consider sundays)
#' 
#' the number of pigs, the average feed intake,
#'   ave weight and sd at the given day in the week. Both data for the selected and uselected pen
#'   are given. Note the observations of ave weight and sd are based on pigs randomly selected and
#'   mesured (often more than once) and hence may vary. The true weight and sd are based on the true
#'   weights of each pig.
#' 
#' @author Lars Relund \email{lars@@relund.dk}
#' @examples dat<-calcWeekSampleData()
#' @export 
calcWeekSampleData<-function(dtPigs, measurementsPerDay, measurementStd) {
  setkey(dtPigs,t)
  ## create weight measurements data set for all pigs (unselected pen)
  pigsMeasured<-data.table() # pigs measured using vision each day
  for (tt in dtPigs[weekDay==7,unique(t)]) { # consider only data at days used in the HMDP (monday morning -> data at sundays)
    alive<-dtPigs[J(tt),pig]  # pigs current day
    if (length(alive)==1) tmp<-rep(alive,measurementsPerDay)
    else tmp<-sample(alive, measurementsPerDay, replace = TRUE)
    pigsMeasured<-rbindlist(list(pigsMeasured,data.table(t=tt,pig=tmp)))
  }  
  # now make sample
  keycols = c("t","pig")
  setkeyv(dtPigs,keycols)
  dtPigsSample<-dtPigs[pigsMeasured,]
  # recalculate OLW since want diff OLW obs for the same pig
  dtPigsSample$OLW<-dtPigsSample$TLW + rnorm(length(dtPigsSample$OLW),0,measurementStd)
  ## create weekly data for all the pigs (unselected pen) 
  tmp1<-dtPigsSample[,list(aveOLWAll=mean(OLW),sdOLWAll=sd(OLW)), by = t]
  # create feed data set
  dtFeed<-dtPigs[,list(pigs=length(unique(pig)), FIt=sum(FI)), by=t] # daily total feed intake (all pigs)
  dtFeed[,week:=(t-1)%/%7+1]
  dtFIPerWeekPerPig<-dtFeed[,list(pigs=mean(pigs), FI=sum(FIt)/mean(pigs)), by = week]
  tmp1[,aveFIAll:=dtFIPerWeekPerPig$FI]
  tmp1[,pigs:=dtFIPerWeekPerPig$pigs]
  tmp<-dtPigs[weekDay==7,list(aveTLWAll=mean(TLW),sdTLWAll=sd(TLW)), by = t]  # find true mean and sd
  setkey(tmp1,t)
  tmp1<-merge(tmp1,tmp)
  
  ## Now do the same for alive pigs (selected pen) 
  ## create weight measurements data set for alive pigs
  dtDaily<-dtPigs  # make a copy for later use
  dtPigs<-dtPigs[culled==0]
  pigsMeasured<-data.table() # pigs measured using vision each day   
  for (tt in dtPigs[weekDay==7,unique(t)]) { 
    alive<-dtPigs[J(tt),pig]  # pigs with culled = 0
    if (length(alive)==1) tmp<-rep(alive,measurementsPerDay)
    else tmp<-sample(alive, measurementsPerDay, replace = TRUE)
    pigsMeasured<-rbindlist(list(pigsMeasured,data.table(t=tt,pig=tmp)))
  }
  # now make sample
  keycols = c("t","pig")
  setkeyv(dtPigs,keycols)
  dtPigsSample<-dtPigs[pigsMeasured,]
  # recalculate OLW since want diff OLW obs for the same pig
  dtPigsSample$OLW<-dtPigsSample$TLW + rnorm(length(dtPigsSample$OLW),0,measurementStd)
  ## create weekly data for all the pigs (unselected pen) 
  setkey(dtPigsSample,weekDay)
  tmp2<-dtPigsSample[J(7)] # consider only data at days we consider in the HMDP (monday morning -> data at sundays)
  tmp2<-tmp2[,list(aveOLWAlive=mean(OLW), sdOLWAlive=sd(OLW)), by = t]
  ## create feed data set
  dtFeed<-dtPigs[,list(alive=length(unique(pig)),FIt=sum(FI)), by=t] # daily total feed intake (alive pigs)
  dtFeed[,week:=(t-1)%/%7+1]
  dtFIPerWeekPerPig<-dtFeed[,list(alive=mean(alive), FI=sum(FIt)/mean(alive)), by = week]
  tmp2[,aveFIAlive:=dtFIPerWeekPerPig$FI]
  tmp2[,alive:=dtFIPerWeekPerPig$alive]
  # add TLW alive
  tmp<-dtPigs[weekDay==7,list(aveTLWAlive=mean(TLW),sdTLWAlive=sd(TLW),feedMix=mean(feedMix) ), by = t]  # find true mean and sd
  setkey(tmp2,t)
  tmp2<-merge(tmp2,tmp)
  # now check if not have any pigs alive in some weeks (alive = 0)
  tt<-tmp1[t>tmp2[,max(t)],t]   # weeks with no pigs alive
  if (!length(tt)==0) {
    tmp<-data.table(t=tt, aveOLWAlive=0, sdOLWAlive=0, aveFIAlive=0, alive=0, aveTLWAlive=0, sdTLWAlive=0, feedMix=tmp[J(tmp[,max(t)]),feedMix])
    tmp2<-rbindlist(list(tmp2,tmp))
  }
  # add further columns
  tmp2[,pen:=dtPigs$pen[1]]
  tmp2[,week:=(t-1)%/%7+1]
  tmp2[,stage:=(t-1)%/%7+1+1]
  tmp2[,weekDay:=(t-1)%%7+1]
  tmp2[,culled:=tmp1$pigs-tmp2$alive]
  
  dtWeek<-merge(tmp1,tmp2) # merge data
  dtWeek<-calcFullObs(dtWeek)
  return(dtWeek)
  #setkey(dtPigs,TLW)  # set key for binary search
  # remove all observations with TLW above threshold
  #firstT<-dtPigs[TLW>threshold,list(t=min(t)),by=pig][order(pig)]  # first obs above threshold
  #dtPigs<-dtPigs[,.SD[t<firstT[pig,t]],by=pig] # save all obs below firstT for each pig
  #dtPigs<-subset(dtPigs, TLW<=threshold)
}

# --------------------------------------------------------------------------------------------- #

#' Cull pigs based on a threshold action at a given stage.
#' 
#' @param dTList A list similar to the one returned by \link{simulatePen}
#' @param stage Stage number.
#' @param th Threshold.
#' 
#' @return The list dTList modified so pigs have been culled and weekly data recalculated.
#' 
#' @author Lars Relund \email{lars@@relund.dk}
#' @examples dat<-thresholdAction(dat,5,90)
#' @export 
thresholdAction<-function(param, dTList, stage, th) { 
  tt<-(stage-1)*7 + 1    # time instance where apply action (cull monday)  
  x<-dTList$dtDailyPig[culled==0 & t==tt,]$OLW # Reza: Observed live weights
  x<-x[order(x,decreasing=T)[1:th]] # Reza: find the th highest value of OLW
  cull<-dTList$dtDailyPig[J(tt),][OLW==x]$pig # Reza: the id of th largest pigs     
  #cull<-dTList$dtDailyPig[J(tt),][OLW>th,pig]   # pigs to be culled
  dTList$dtDailyPig[t>=(tt) & pig %in% cull,culled:=1]   # cull pigs
  # recalculate weekly data
  datWeekly<-calcWeekSampleData(dTList$dtDailyPig, param$measurementsPerDay, param$measurementStd) 
  dTList$dtWeekAve<-datWeekly
  return(dTList)
}

# --------------------------------------------------------------------------------------------- #

#' Resimulate based on a new feed-mix action at a given stage.
#' 
#' @param dTList A list returned by \link{simulatePen}
#' @param stage Stage number.
#' @param feedMix New feed-mix.
#' 
#' @return The list dTList modified so pigs have been culled and weekly data recalculated.
#' 
#' @author Lars Relund \email{lars@@relund.dk}
#' @examples dat<-thresholdAction(dat,5,90)
#' @export 
changeFeedMixAction<-function(param, dTList, stage, feedMix) { 
  tt<-(stage-1)*7  # last time instance where use old feed-mix (sunday)
  T<-dTList$dtDailyPig[,max(t)]  # last time instance simulated
  datDaily<-dTList$dtDailyPig[t<=tt,] # all observations until now
  TLW<-datDaily[J(tt),TLW]
  culled<-datDaily[J(tt),culled]
  # resimulate
  tmp<-createDailySampleData(param$pen, param$pigIds, TLW, culled, param$measurementStd, param$k4Mean[feedMix], 
                             param$DGFactor, T = T, DFI = tt, feedMix = feedMix)
  datDaily<-rbindlist(list(datDaily,tmp) )
  datWeekly<-calcWeekSampleData(datDaily,  param$measurementsPerDay, param$measurementStd) 
  dTList$dtDailyPig<-datDaily
  dTList$dtWeekAve<-datWeekly
  return(dTList)
}

# --------------------------------------------------------------------------------------------- #

#' Update the sample data with "pseudo" observations (observations if all pigs in the pen 
#' were present) based on the formulars of the truncted normal distribution 
#' (Herd management book, Vol 1, p 128) 
#' 
#' @param dat Weekly sample data table.
#' 
#' @return Updated data table with additional observation columns (ending with Calc).
#' 
#' @author Lars Relund \email{lars@@relund.dk}
#' @export 
calcFullObs<-function(dat) {
  dat$sdOLWAllCalc<-dat$sdOLWAlive
  dat$aveOLWAllCalc<-dat$aveOLWAlive
  th<-qnorm(dat$alive/dat$pigs)
  cMean<-dnorm(th)/pnorm(th)
  cVar<-1-th*cMean-cMean^2
  idx<-which(th<Inf & th > -Inf)   # those weeks where pigs have been culled (not all culled)
  dat$sdOLWAllCalc[idx]<-dat$sdOLWAlive[idx]/sqrt(cVar[idx])
  dat$aveOLWAllCalc[idx]<-dat$aveOLWAlive[idx]+dat$sdOLWAllCalc[idx]*cMean[idx]     
  return(dat)
}

# --------------------------------------------------------------------------------------------- #

#' Set the simulation parameters for a pen.
#' 
#' @param pen Pen id added to the output.
#' @param pigIds Id of the pigs simulated
#' @param measurementsPerDay Measurements per day from the weigthing sensor (normally vision).
#' @param measurementStd The std of the weight sensor (normally vision).
#' @param k4Mean Mean of growth rate parameter (Gompertz growth curve). One value for each possible feed-mix.
#' @param DGFactor Increase the daily gain (x2) with this factor.
#' 
#' @return A list containing the parameters. 
#' 
#' @author Lars Relund \email{lars@@relund.dk}
#' @examples dat<-simulatePen()
#' @export 
setSimParam<-function(pen=1, pigIds=1:15, 
                      measurementsPerDay=30,  
                      measurementStd=1,
                      k4Mean=c(0.01012228,0.01150325,0.01387502,0.01487517),
                      DGFactor=1) {
  return(list(pen=pen, pigIds=pigIds, measurementsPerDay=measurementsPerDay, 
              measurementStd=measurementStd, DGFactor = DGFactor, k4Mean = k4Mean))
}

# --------------------------------------------------------------------------------------------- #

#' Create a data set simulating a pen from DFI and T days ahead. We assume that actions are
#' considered on Mondays (in the morning), i.e week w in the HMDP corresponds week w-1 and
#' weight data to data from sunday the week before!
#' 
#' @param param Pen simulation parameters.
#' @param feedMix Feed-mix used.
#' @param TLW Start weight day before simulate (same length as pigIds).
#' @param culled Have the pig we simulate been culled (same length as pigIds).
#' @param T Number of days to simulate.
#' @param DFI Days from insertion into the pen (start simulation at t=DFI+1)
#' 
#' @return A list with two data tables. Table \code{dtDailyPig} contains a row 
#' for each pig on a daily basis. Table \code{dtWeekAve}
#' contain a row for each week with the numbe of pigs, the average feed
#' intake, ave weight and sd at the given day in the week. 
#' 
#' @author Lars Relund \email{lars@@relund.dk}
#' @examples dat<-simulatePen()
#' @export 
simulatePen<-function(param, feedMix, 
                      TLW = rnorm(length(param$pigIds),30,param$measurementStd),
                      culled = rep(0,length(param$pigIds)),
                      T = 7*11, DFI = 0
) {
  datDaily<-createDailySampleData(param$pen, param$pigIds, TLW, culled, param$measurementStd, 
                                  param$k4Mean[feedMix], param$DGFactor, T = T, DFI = DFI, feedMix = feedMix)
  datWeekly<-calcWeekSampleData(datDaily, param$measurementsPerDay, param$measurementStd) 
  return(list(dtDailyPig=datDaily, dtWeekAve=datWeekly))
}

# --------------------------------------------------------------------------------------------- #

#' Add estimated posterior values to observed values of the GSSM and nGSSM.
#' 
#'@param dat Data table with observed values (columns aveOLW, sdOLW and aveFI).
#'@param m0 Initial mean in DLM. This value is equal to filtered mean at time "startT-1"
#'@param C0 Initial variance in DLM. This value is equal to filtered variance at time "startT-1"
#'@param var0 Initial variance in DGLM. This value is equal to the varince at time "startT-1" resulted from DGLM.
#'@param startT Starting time of using the DLM and the DGLM (starting time after having a new decision)

#'@return An extended data table with new columns eAveOLW, eSdOLW and eAve, and also the posterior variances to be used as the initial variances in the next run of SSMs 
#'
#'@author Reza Pourmoayed \email{rpourmoayed@@econ.au.dk}
#'@export 
estimatePosterior<-function(dat, m0, C0, var0, startT){  
  
  # compute the filtered data:
  #DLM model;
  # store the observations related to the DLM in the Array D   
  D<-array(NA,dim=c(2,1,(param$tMax-1)-startT+1))
  for(t in 1:((param$tMax-1)-startT+1) ){
    D[1,1,t]<-round(dat$aveOLWAll[t] , 3)
    D[2,1,t]<-round(dat$aveFIAll[t]  ,3)
  }
  #Use the filter source to filter the raw data 
  #source("GSSM.R")
  mod<-setGSSM()
  mod$m0<-m0
  mod$C0<-C0
  mod$t<-(param$tMax-1)-startT+1
  dlm<-DLMfilter(mod,D,mod$W)
  eAveOLW<-dlm$L1[1,,] # filtered live weight 
  eAveG<-dlm$L1[2,,] # filtered growth rate
  eVarOLW<-dlm$L2[1,1,] # Variance of estimation for OLW
  eVarG<-dlm$L2[2,2,] # Variance of estimation for G
  eCov<-dlm$L2[1,2,] # Covariance of estimation for OLW and G
  
  #DGLM model
  # store the observations related to the DGLM in the vector ob (variance components)
  #source("nGSSM.R")
  mod1<-setnGSSM()
  ob<-round(dat$sdOLWAll^2 , 3)
  mod1$m0=var0
  mod1$t<-(param$tMax-1)-startT+1
  dglm<-DGLMfilter(mod1,ob,mod1$W)
  eSdOLW<-sqrt(dglm$mtt)
  
  dat[, "eAveOLWAll"]<-round(eAveOLW,3)   
  dat[, "eAveGAll"]<-round(eAveG, 3)
  dat[, "eSdOLWAll"]<-round(eSdOLW, 3)
  dat[, "eVarOLWAll_cov"]<-round(eVarOLW, 3)
  dat[, "eVarGAll_cov"]<-round(eVarG, 3)
  dat[, "eCovAll_cov"]<-round(eCov, 3)
  
  return(dat)   
}

# --------------------------------------------------------------------------------------------- #

#' Add the optimal actions to a data set.
#' 
#' @param dat Data table with estimated values (columns eAveOLW, eSdOLW and eAveG).
#' @param dat Data table with observed values (columns aveOLW, sdOLW and aveFI).
#' @param policy Optimal policy acquired from the hmdp . We need to run HMDP and get the optimal policy
#' @param startT Starting time of using the DLM and the DGLM (starting time after having a new decision)
#' @param phase The phase number of feeding in the system at time "startT"
#' @param feedMix The feed-mix that is used in week "startT"
#' @param startFeed The time that the last feedMix has been started
#' 
#' @return An extended data table with new columns threshold (NA if not used) and feedMix (feed-mix used).
#' 
#' @author Reza Pourmoayed \email{rpourmoayed@@econ.au.dk} 
#' @export 
findActions<-function(pen, paramSim, param, policy, startT, phase, feedMix, startFeed) {
  # find the id of the states based on the true information 
  
  dat<-pen$dtWeekAve
  
  idxTrueWeight<-c()
  for(t in 1:((param$tMax-1)-startT+1) ){
    for(i in 1:dim(param$wIntervals[[startT+t]])[1]){
      if( (param$wIntervals[[startT+t]][i,2]<=dat$eAveOLWAll[t]) && (dat$eAveOLWAll[t]<param$wIntervals[[startT+t]][i,3]) )
        idxTrueWeight[t]<-i-1
    }
  }
  
  idxTrueGrowth<-c()
  for(t in 1:((param$tMax-1)-startT+1) ){
    for(i in 1:dim(param$gIntervals[[feedMix]])[1]){
      if( (param$gIntervals[[feedMix]][i,2]<=dat$eAveGAll[t]) && (dat$eAveGAll[t]<param$gIntervals[[feedMix]][i,3]) )
        idxTrueGrowth[t]<-i-1
    }
  }
  
  idxTrueSdWeight<-c()
  for(t in 1:((param$tMax-1)-startT+1) ){
    for(i in 1:dim(param$sDIntervals[[startT+t]])[1]){
      if( (param$sDIntervals[[startT+t]][i,2]<=dat$eSdOLWAll[t]) && (dat$eSdOLWAll[t]<param$sDIntervals[[startT+t]][i,3]) )
        idxTrueSdWeight[t]<-i-1
    }
  }
  
  #find the optimal decision per week based on the optimal policy acquired from the HMDP
  
  if(phase==1){
    ration=0 
  } else{
    ration=feedMix-1
  }
  opt<-c()
  optMChange<-c()
  optFeed<-c()
  optMarket<-c()
  for(t in 1:((param$tMax-1)-startT+1) ){
    
    if(dat$alive[t]==0) break
    
    opt[t]<-subset(policy, stateLabel==paste("(",dat$alive[t],",",idxTrueWeight[t],",",idxTrueGrowth[t],",",idxTrueSdWeight[t],",",startFeed+1,",",(phase-1),",",feedMix-1,",",startT+t,")",sep="") ) ["actionLabel"]  
    
    if( (opt[t]!="cont.") && (opt[t]!="term.") && (opt[t]!="change fm")  ){
      opt[t] = unlist(strsplit(paste(opt[t]), split='l ', fixed=TRUE))[2]
      optMarket[t]<-as.numeric(opt[t])
      pen<-thresholdAction(paramSim, dTList=pen, stage=startT+t, th=optMarket[t])
      dat<-pen$dtWeekAve
    }
    
    if(opt[t]=="term."){
      optMarket[t]<-dat$alive[t]
      pen<-thresholdAction(paramSim, dTList=pen, stage=startT+t, th=optMarket[t])
      dat<-pen$dtWeekAve      
    }
    
    if( (opt[t]!="change fm")  )
      optFeed[t]<-feedMix
    
    if(opt[t]=="change fm"){
      opt[t]<-subset(policy, stateLabel== paste("(",feedMix-1,",",startT+t,",",dat$alive[t],",",idxTrueWeight[t],",",idxTrueSdWeight[t],",",phase,")",sep="") ) ["actionLabel"] 
      opt[t] = unlist(strsplit(paste(opt[t]), split='= ', fixed=TRUE))[2]
      optFeed[t]<-as.numeric(opt[t]) + 1
      idGrowth<-param$rGIdx[optFeed[t]] #as.integer(dim(param$gIntervals[[optFeed[t]]])[1]/2)
      optMChange[t]<-subset(policy, stateLabel==paste("(",dat$alive[t],",",idxTrueWeight[t],",",idGrowth,",",idxTrueSdWeight[t],",",startT+t,",",phase,",",optFeed[t]-1,",",startT+t,")",sep="") ) ["actionLabel"] 
      if(optMChange[t]=="term."){
        optMarket[t]<-dat$alive[t]
        pen<-thresholdAction(paramSim, dTList=pen, stage=startT+t, th=optMarket[t])
        dat<-pen$dtWeekAve
      }
      if( (optMChange[t]!="cont.") && (optMChange[t]!="term.") ){
        optMarket[t] = unlist(strsplit(paste(optMChange[t]), split='l ', fixed=TRUE))[2]
        optMarket[t]<-as.numeric(optMarket[t])
        pen<-thresholdAction(paramSim, dTList=pen, stage=startT+t, th=optMarket[t])
        dat<-pen$dtWeekAve
      }        
    }      
  }   
  
  # Add 2 columns to the data frame pen$dtWeekAve  
  
  pen$dtWeekAve[, "OptFeedMix"]<-optFeed
  pen$dtWeekAve[, "optThreshold"]<-optMarket
  
  return(pen)   
}

# --------------------------------------------------------------------------------------------- #

#' Remove rows where re-simulation is needed. 
#' 
#' @param dat List with 2 data tables (see createSampleData) where the weekly data table 
#'        contains optimal actions (columns threshold and feedMix).
#' @param phase Number of phase feeding in the system
#' @param feedMix The feed-mix that is used in week "startT"
#' @param startFeed The time that the last feedMix has been started
#' 
#' @return A reduced data table.
#' @export 
cutOff<-function(pen, feedMix, m0, C0, var0, phase, startFeed, startT, PigsCull) {
  # some data to play with
  dat1<-pen$dtWeekAve
  dat1<-estimatePosterior(pen$dtWeekAve, m0, C0, var0, startT)
  dat2<-pen$dtDailyPig
  
  # find week where feedMix change
  mx<-feedMix; weekMx <- Inf;
  if(length(dat1$OptFeedMix)==1){
    weekMx <- Inf
  }else{
    for (i in 1:length(dat1$OptFeedMix)) {
      if (mx!=dat1$OptFeedMix[i]) {weekMx <- dat1$week[i]; mx1<-dat1$OptFeedMix[i]; break} #for the last week we had problem here and I change the method here
    }
  }   
  
  if (weekMx != Inf) { # change the feed-mix: do not remove pigs and re-simulate from weekMX + 1
    
    PigsCull<-dat2[t==7*weekMx+1,culled]
    
    dat2<-dat2[week<=weekMx]
    dat1<-dat1[week<=weekMx]
    
    
    m0<-matrix(data=c(dat1[week==weekMx]$eAveOLW,dat1[week==weekMx]$eAveG),nrow=2)
    C0<-matrix(data=c(dat1[week==weekMx]$eVarOLW,0.32,0.32,0.53),ncol=2)
    var0<-dat1[week==weekMx]$eSdOLW
    
    phase<-phase + 1
    startFeed<-weekMx
    feedMix<-mx1
  }
  dat1$eVarOLWAll_cov<-NULL
  dat1$eVarGAll_cov<-NULL
  dat1$eCovAll_cov<-NULL
  
  return(list( pen=list(dtDailyPig=dat2,dtWeekAve=dat1), feedMix=feedMix, m0=m0, C0=C0, var0=var0, phase=phase, startFeed=startFeed, PigsCull=PigsCull ))
}

# --------------------------------------------------------------------------------------------- #

#' Compute the growth rate parameter of the Gompertz function based on the current live weight and the daily gain. 
#' 
#' @param W weight at the current week in the production system (we suppose weight is the insertion weight in the system)
#' @param K Logarithm of the outgrowth weight in kg
#' @param G Daily gain in kg
#'
#' @return The growth rate parameter in the Gompertz function.
#' @author Reza Pourmoayed \email{rpourmoayed@@econ.au.dk}
#' @export 
GrowthParam<-function(W,G,K){
  return( G/(W*(K-log(W))) )
}
