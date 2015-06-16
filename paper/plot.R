## R script file for plotting the simulated data
require(tools)
# Three simulated pens
pen1<-pen1Weekly
pen2<-pen2Weekly
pen3<-pen3Weekly

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
dat$name<-mapvalues(dat$name, from = c("aveOLWAll","aveFIAll","sdOLWAll","eAveOLWAll","eAveGAll","eSdOLWAll"), to = c("$ \\bar{w}_t $ \\small(kg)","$ \\bar{z}_t $ \\small(FEsv)","$ s_t $ \\small(kg)","$ \\hat{\\mu}_t $ \\small(kg)","$ \\hat{g_t} $ \\small(kg)","$ \\hat{\\sigma}_t $ \\small(kg)") )


# vertical lines 
tmp<-changeTime(pen1)
vline.fm <- data.frame(w = tmp[[2]], pen = "Pen 1 (low)")
vline.th <- data.frame(w = tmp[[1]], pen = "Pen 1 (low)")
tmp<-changeTime(pen2)
vline.fm <- rbind(vline.fm,data.frame(w = tmp[[2]], pen = "Pen 2 (normal)"))
vline.th <- rbind(vline.th,data.frame(w = tmp[[1]], pen = "Pen 2 (normal)"))
tmp<-changeTime(pen3)
vline.fm <- rbind(vline.fm,data.frame(w = tmp[[2]], pen = "Pen 3 (high)"))
vline.th <- rbind(vline.th,data.frame(w = tmp[[1]], pen = "Pen 3 (high)"))

# pigs in pen
datPigs<-dtPlot[,list(pen,week,alive)]
datPigs$pen<-factor(datPigs$pen, labels=c("Pen 1 (low)","Pen 2 (normal)","Pen 3 (high)"))
datPigs$y<-datPigs$alive
datPigs$name<-NA


#Plot the data related to the simulation and the SSMs: 
library(ggplot2)
library(grid)
library(tikzDevice)
tikz("sim_plot.tex", width = 10, height = 7, standAlone=T)
plot<-ggplot(data=dat, aes(x=factor(week), y=y, group=name, shape=name, linetype=name ) )+ 
  geom_line() + scale_y_continuous(breaks=seq(0,120,10)) +
  #geom_point() + 
  facet_grid(. ~ pen) + 
  xlab("week") + ylab(" ") 
g <- guide_legend("", override.aes = list(fill=NA))

plot + guides(shape = g, linetype=g) + 
  geom_histogram(stat="identity", data=datPigs, alpha = 1/4, colour=NA, width=0.25) + 
  geom_vline(aes(xintercept = w), data=vline.fm, color="gray") + 
  geom_vline(aes(xintercept = w), data=vline.th, color="gray", linetype="twodash") +
  geom_line() + 
  theme_bw() + 
  theme(legend.position="bottom", panel.background = element_blank(), 
        panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
        legend.key = element_rect(fill = NA, colour = NA), 
        legend.key.width = unit(1, "cm"), legend.text.align=0.5,
        strip.background=element_rect(fill = NA)) 
dev.off()

#------------------------------------------------------------------------------------------------
# Plot optimal decisions
#pdf("tmp.pdf", width = 12, height = 6)
#library(tikzDevice)
tikz("opt_plot.tex", width = 8.5, height = 3.5, standAlone =T)
titles<-c("Pen 1 (low)","Pen 2 (normal)","Pen 3 (high)")
#par(mfrow=c(1,3))
par(mar = c(4, 0, 3, 0))
layout(t(c(1,2,2,2,3,3,3,4,4,4)))
# make y axes
plot(c(0,1), c(1,4.4), yaxt="n", xlab='', ylab='', xaxt="n", bty='n', pch=NA)
axis(2, las=1, at=c(1.5, 2, 2.5), labels=c("cull", "continue", "terminate"), 
     line=-8, cex.axis=1.36)
axis(2, las=1, at=c(3.2,3.8,4.4), labels=c("feed-mix 1","feed-mix 2", "feed-mix 3"), 
     line=-8, cex.axis=1.36)
# plot the pens
for(hh in 1:3){
  if(hh==1)
    pen<-pen1
  if(hh==2)
    pen<-pen2
  if(hh==3)
    pen<-pen3
  
  plot(c(0,13), c(1,4.4), yaxt="n", xlab='', ylab='', xaxt="n", bty='n', pch=NA)
  abline(v=1:12,lty= 2, col="gray85")
  axis(1, las=1, at=c(1,2,3,4,5,6,7,8,9,10,11,12), labels=c("1","2","3","4","5","6","7","8","9","10","11","12"), 
       cex.axis=1)
  title(main=titles[hh], cex.main=1)
  
  vecFM<-changeTime(pen)[[2]]
  vecTH<-changeTime(pen)[[1]]
  last<-vecTH[length(changeTime(pen)[[1]] ) ]
  
  vecFCor<-c(3.2,3.8,4.4)
  if( vecFM[1]==0){
    datL3<-rep(vecFCor[1],last)
    lines(datL3, type="l", lwd=4)
  }else{
    for(j in 1:length(vecFM) ){
      if(j==1){
        datL<-rep(vecFCor[1],vecFM[1])
        lines(datL, type="l", lwd=4)
        counter<-vecFM[j]
      }
      datL1<-rep(NA,counter-1)
      if( j==(length(vecFM) ) ){
        dur<-last-vecFM[j]+1
        datL2<-rep(vecFCor[pen$OptFeedMix[vecFM[j]]],dur)
      }else{
        dur<-vecFM[j+1]-vecFM[j]+1
        datL2<-rep(vecFCor[pen$OptFeedMix[vecFM[j]]],dur)
      }
      datL<-c(datL1,datL2)
      lines(datL, type="l", lwd=4)
      counter<-counter+dur-1
    }
  }
  for(i in 1:length(vecTH)){  
    if(pen$optThreshold[vecTH[i]]==pen$alive[vecTH[i]]){
      label="" #paste("term.")
      points(vecTH[i],2.5, type="p", lwd=4, pch=15) 
      text(vecTH[i],2.5, labels=label,pos=1, cex=1)
    }else{
      label=paste(pen$optThreshold[vecTH[i]])
      points(vecTH[i],1.5, type="p", lwd=6, pch=19, cex=0.5) 
      text(vecTH[i],1.5, labels=label,pos=1, cex=1)
    }
  }
  for(k in 1:dim(pen)[1]){
    if( (is.na(pen$optThreshold[k]))  ){ #&& (k!=vecFM)
      label="" #paste("cont.")
      points(k,2, type="p", lwd=4, pch=17) 
      text(k,2, labels=label,pos=1, cex=1)
    }
  }
  if (hh==2) mtext(text="week", side = 1, line = 3, at=6.5, cex=1)
}

#legend(0,4.5,c("Feeding action - symbol in HMDP: $ a^f $"," Individual marketing action (kg) - symbol in HMDP: $ a^\\delta $", "Continuing action - symbol in HMDP: $ \\mathtt{\\small cont.} $", "Terminating action - symbol in HMDP: $ \\mathtt{\\small term.} $"), pch = c(NA,19,17,15),lty=c(1,NA,NA,NA), lwd=c(6,NA,NA,NA), cex=1.5,box.lty=0, y.intersp=2.5 )
dev.off()

#do.call(file.remove,list(list.files(pattern = ".tex")))
