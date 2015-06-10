## R script file for plotting the simulated data




# Plot the data: 

#Three simulated pen 

#Three simulated pen 

pen1<-pen1Weekly
pen2<-pen2Weekly
pen3<-pen3Weekly

insertRow<-function(pen, penDaily){
  
  for( i in 1:dim(pen)[1])
    pen$week[i]<-pen1$week[i]+1
  
  pen<-rbind(pen[1],pen)
  
  pen[1]$t<-as.integer(0)
  pen[1]$aveOLWAll<-mean(penDaily[t==1]$OLW)
  pen[1]$aveTLWAll<-mean(penDaily[t==1]$TLW)
  pen[1]$sdOLWAll<-sd(penDaily[t==1]$OLW)
  pen[1]$aveFIAll<-sum(penDaily[t==1]$FI)/15*7
  pen[1]$week<-1   
  pen[1]$stage<-1    
  pen[1]$eAveOLWAll<-26.49
  pen[1]$eAveGAll<-6
  pen[1]$eSdOLWAll<-2.850987
  pen[1]$eVarOLWAll<-4.26
  pen[1]$eVarGAll<-0.53
  pen[1]$eCovAll<-0.32
  return(pen)
}

pen1<-insertRow(pen1, pen1Daily)
pen2<-insertRow(pen2, pen2Daily)
pen3<-insertRow(pen3, pen3Daily)
#pen3<-pen3[-12,]

#' Find the times when we resimulate 
#' 
#' @param pen A data fram included our information
#' @return a vector with times when we resimulate (when a ration is changed or a threshould is selected)
changeTime<-function(pen){
  th<-c()
  fe<-c()
  change<-list()
  for(i in 1:dim(pen)[1]){
    if( !is.na(pen$optThreshold[i]) )
      th[i]<-i
  }
  for(j in 1:(dim(pen)[1]-1) ){
    if(pen$OptFeedMix[j]!=pen$OptFeedMix[j+1])      
      fe[j+1]<-j+1
  }
  
  th<- as.integer(th[!is.na(th)])
  if(length(fe)==0){
    fe<-0
  }else{
    fe<- as.integer(fe[!is.na(fe)])
  }   
  change[[1]]<-th
  change[[2]]<-fe
  #    change<-c(fe,th)
  #    change<- change[!is.na(change)]
  #    change<-sort(change)
  return(change)
}



dtPlot<-rbindlist(list(pen1,pen2,pen3))

# transform it to long format
library(reshape2)
library(ggplot2)
dat<-melt(data.frame(dtPlot),
          # ID variables - all the variables to keep but not split apart on
          id.vars=c("pen", "week","weekDay","t"),
          # The source columns
          measure.vars=c("aveOLWAll", "aveFIAll", "sdOLWAll", "eAveOLWAll", "eAveGAll", "eSdOLWAll" ),
          # Name of the destination column that will identify the original
          # column that the measurement came from
          variable.name="name",
          value.name="y"
)
dat$pen<-factor(dat$pen, labels=c("Pen 1 (low)","Pen 2 (normal)","Pen 3 (high)"))
library(plyr)
dat$name<-mapvalues(dat$name, from = c("aveOLWAll","aveFIAll","sdOLWAll","eAveOLWAll","eAveGAll","eSdOLWAll"), to = c("$ \\bar{w}_t $ \\small(kg)","$ \\bar{z}_t $ \\small(FEsv)","$ s^{2}_t $ \\small(kg)","$ \\hat{\\mu}^{w}_t $ \\small(kg)","$ \\hat{\\mu}^{g}_t $ \\small(kg)","$ \\hat{\\sigma}^{w}_t $ \\small(kg)") )




# vertical lines 

# tmp<-changeTime(pen1)
# vline.fm <- data.frame(w = tmp[[2]], pen = "Pen 1 (low)")
# vline.th <- data.frame(w = tmp[[1]], pen = "Pen 1 (low)")
# tmp<-changeTime(pen2)
# vline.fm <- rbind(vline.fm,data.frame(w = tmp[[2]], pen = "Pen 2 (average)"))
# vline.th <- rbind(vline.th,data.frame(w = tmp[[1]], pen = "Pen 2 (average)"))
# tmp<-changeTime(pen3)
# vline.fm <- rbind(vline.fm,data.frame(w = tmp[[2]], pen = "Pen 3 (high)"))
# vline.th <- rbind(vline.th,data.frame(w = tmp[[1]], pen = "Pen 3 (high)"))

#-------------------------

#Plot the data related to the simulation and the SSMs: 
# dat$linecol<-"black"
# dat$linecol[substr(dat$name,1,1)=="e"]<-"gray"



library(tikzDevice)
tikz("sim_plot1.tex", width = 10, height = 7, standAlone=T)

plot<-ggplot(data=dat,  aes(x=factor(week), y=y, group=name, shape=name, linetype=name )  ) + geom_line() +  geom_point() +  facet_grid(. ~ pen) + 
  xlab("week") + ylab(" ") 
#theme_bw() #+
#theme(panel.background = element_blank(), panel.grid.minor=element_blank(), panel.grid.major=element_blank())
g <- guide_legend("")
plot +    guides(shape = g, linetype=g) +
  geom_vline(data=subset(dat, pen=="Pen 1 (low)"), aes(xintercept=changeTime(pen1)[[1]][-length(changeTime(pen1)[[1]] ) ] ), linetype="solid", color="gray66") + # Added by Reza
  geom_vline(data=subset(dat, pen=="Pen 1 (low)"), aes(xintercept=changeTime(pen1)[[2]] ), linetype="solid", color="gray10") +
  
  geom_vline(data=subset(dat, pen=="Pen 2 (normal)"), aes(xintercept=changeTime(pen2)[[1]][-length(changeTime(pen2)[[1]] ) ] ), linetype="solid", color="gray66")+  # Added by Reza
  #geom_vline(data=subset(dat, pen=="Pen 2 (normal)"), aes(xintercept=changeTime(pen2)[[2]] ), linetype="solid", color="gray10")+ # 
  
  geom_vline(data=subset(dat, pen=="Pen 3 (high)"), aes(xintercept=changeTime(pen3)[[1]][-length(changeTime(pen3)[[1]] ) ] ), linetype="solid", color="gray66")+  # Added by Reza
  geom_vline(data=subset(dat, pen=="Pen 3 (high)"), aes(xintercept=changeTime(pen3)[[2]] ), linetype="solid", color="gray10")  

#geom_vline(aes(xintercept = w), data=vline.fm, color="gray") + 
#geom_vline(aes(xintercept = w), data=vline.th, color="gray", linetype="twodash") +  

dev.off()


#-------------------

# Plot optimal decisions

pdf("optimal.pdf", width = 15, height = 10)

library(tikzDevice)
tikz("optimal2.tex", width = 12, height = 8, standAlone=T)

par(mfrow=c(2,2))

titles<-c("Pen 1 (low)","Pen 2 (normal)","Pen 3 (high)")

for(hh in 1: 3){
  if(hh==1)
    pen<-pen1
  if(hh==2)
    pen<-pen2
  if(hh==3)
    pen<-pen3
  
  #Plot optimal decisions (Average growth rate)
  plot(c(0,13), c(1,4.7), yaxt="n", xlab='', ylab='', xaxt="n", bty='n', pch=NA)
  abline(v=1:12,lty= 2, col="gray")
  #axis(2, las=1, at=c(2,4), labels=c("Marketing\n decisions","Feeding \ndecisions"),lwd=15, line=-1, lty=0, cex.axis=0.8)
  axis(2, las=1, at=c(1.5, 2, 2.5), labels=c("Individual\n marketing", "Continue", "Terminating"), line=-3, cex.axis=1.36)
  axis(2, las=1, at=c(3.2,3.7,4.2,4.7), labels=c("Feed-mix 1","Feed-mix 2", "Feed-mix 3", "Feed-mix 4"), line=-3, cex.axis=1.36)
  axis(1, las=1, at=c(1,2,3,4,5,6,7,8,9,10,11,12), labels=c("1","2","3","4","5","6","7","8","9","10","11", "12"), cex.axis=1.15)
  mtext(text="Weeks after insertion into the pen", side = 1, line = 3, at=6.5, cex=1.1)
  title(main=titles[hh], cex.main=1.9)
  
  
  vecFM<-changeTime(pen)[[2]]
  vecTH<-changeTime(pen)[[1]]
  last<-vecTH[length(changeTime(pen)[[1]] ) ]
  
  vecFCor<-c(3.2,3.7,4.2,4.7)
  if( vecFM==0){
    datL3<-rep(vecFCor[3],last)
    lines(datL3, type="l", lwd=12)
  }else{
    for(j in 1:length(vecFM) ){
      if(j==1){
        datL<-rep(vecFCor[3],vecFM[1])
        lines(datL, type="l", lwd=12)
        counter<-vecFM[j]
      }
      datL1<-rep(NA,counter-1)
      if( j==(length(vecFM) ) ){
        dur<-last-vecFM[j]+1
        datL2<-rep(vecFCor[pen$feedMix[vecFM[j]]],dur)
      }else{
        dur<-vecFM[j+1]-vecFM[j]+1
        datL2<-rep(vecFCor[pen$feedMix[vecFM[j]]],dur)
      }
      datL<-c(datL1,datL2)
      lines(datL, type="l", lwd=12)
      counter<-counter+dur-1
    }
  }
  
  
  for(i in 1:length(vecTH)){  
    if(pen$threshold[vecTH[i]]==0){
      label=paste("term.")
      points(vecTH[i],2.5, type="p", lwd=4, pch=15) 
      text(vecTH[i],2.5, labels=label,pos=1, cex=1.2)
    }else{
      label=paste(pen$threshold[vecTH[i]])
      points(vecTH[i],1.5, type="p", lwd=6, pch=19) 
      text(vecTH[i],1.5, labels=label,pos=1, cex=1.5)
    }
  }
  
  for(k in 1:dim(pen)[1]){
    if( (is.na(pen$threshold[k])) && (k<9) && (k!=vecFM) ){
      label=paste("cont.")
      points(k,2, type="p", lwd=4, pch=17) #Individual marketing
      text(k,2, labels=label,pos=1, cex=1.2)
    }
  }
  
}


plot(c(0,13), c(1,4.7), yaxt="n", xlab='', ylab='', xaxt="n", bty='n', pch=NA)
legend(0,4.5,c("Feeding action - symbol in HMDP: $ a^f $"," Individual marketing action (kg) - symbol in HMDP: $ a^\\delta $", "Continuing action - symbol in HMDP: $ \\mathtt{\\small cont.} $", "Terminating action - symbol in HMDP: $ \\mathtt{\\small term.} $"), pch = c(NA,19,17,15),lty=c(1,NA,NA,NA), lwd=c(6,NA,NA,NA), cex=1.5,box.lty=0, y.intersp=2.5 )
dev.off()

