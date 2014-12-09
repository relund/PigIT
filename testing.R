## Build the model (test)
library(PigIT)
require(Rcpp)   
require(RcppArmadillo)   
library(data.table)
set.seed(234567)

param<-setParameters(pigs=4,  # program crash if rations = 1!!
                     rations=3,
                     phases=3,
                     tMax=12,
                     iniFeedMix=2,
                     minPhaseT=c(1,4,7,10),
                     disWeight=c(5,3),
                     disSD=c(3,0.8),
                     disGrowth=c(2,0.1),
                     priorGrowth=c(5.2,5.7,6.2,6.7),
                     cullActions = T
)
# create DLMs
dlms<-list()
for (i in 1:param$rations){
   dlmMod<-iniDLM()
   Y<-matrix(c(28.91122,11.17050, 34.73722, 11.57909, 40.55647, 11.91425, 46.43588, 12.31422, 52.60592, 12.76750, 58.90283, 13.18425, 64.73567, 13.57245, 70.88710, 13.94777, 76.19884, 14.36087, 82.30708, 14.76617, 88.41322, 15.10727, 94.68169,  15.56838, 100.68169, 15.96838, 106.68169, 16.57838, 112.68169, 17.06838), nrow=2)
   dlms[[i]]<-buildDLM(dlmMod, Y)
}

# Get the DGLM parameters.
dglmParam<-iniDGLM()

prefix<-"test_"
BuildHMDP(prefix, param, dlms, dglmParam)


#------------------------------
#solve:
rm(mdp)
wLbl<-"Reward"
durLbl<-"Time"
mdp<-loadMDP(prefix, check = FALSE)
g<-policyIteAve(mdp,wLbl,durLbl)      # Finds the optimal policy using the average reward per week (g) criterion 



# After changes in setParam function (ok)
# Run policy iteration under average reward criterion using 
# reward 'Reward' over 'Time'. Iterations (g):
#    1 (1.3102) 2 (1.3102) finished.

# 1 - Results after
# Run policy iteration under average reward criterion using 
# reward 'Reward' over 'Time'. Iterations (g):
#    1 (1.4446) 2 (1.4446) finished.

# 0 - Results before changing + to - in SimulatePigs when calc lWM
# Run policy iteration under average reward criterion using 
# reward 'Reward' over 'Time'. Iterations (g):
#    1 (1.39684) 2 (1.39684) finished.




policy<-getPolicy(mdp, labels=TRUE)   # optimal policy for each sId
sIdx<-stateIdxDf(prefix)              # information about the states
policy<-merge(sIdx,policy)            # merge the two data frames
policyW<-getPolicyW(mdp, wLbl)        # the optimal rewards of the policy
policy<-merge(policy,policyW)         # add the rewards
#policy<-merge(policy,states,by.x="sId",by.y="id")   # add state variables


# transform label=(q,iSW,iSG,iSSd,RS, phase, iRation, t) into the (center)values. NOTE iRation is not the feed-mix but feedmix-1!!
label2Values<-function(label, param) {
   label<-gsub("[()]","",label)
   label<-strsplit(label, split=",")
   label<-lapply(label,as.numeric)
   values<-matrix(nrow=length(label),ncol=8)
   for (j in 1:length(label)) {
      if (length(label[[j]])<8) next
      values[j,1]<-label[[j]][1]    # q
      values[j,5]<-label[[j]][5]    # startWeek
      values[j,6]<-label[[j]][6]+1  # phase
      w<-label[[j]][8]  # week
      w1<-label[[j]][7] +1
      values[j,8]<-w
      values[j,2]<-param$wIntervals[[w]][label[[j]][2]+1,1]  # weight
      values[j,3]<-param$gIntervals[[w1]][label[[j]][3]+1,1]  # growth
      values[j,4]<-param$sDIntervals[[w]][label[[j]][4]+1,1] # sd
      values[j,7]<-label[[j]][7]+1 # feed-mix
   }   
   values<-as.data.table(values)
   setnames(values,c("q","sW","sG","sSd","startWeek", "phase", "ration", "week") )
   return(values)
}

#------------------------------
dat<-NULL
# Find the states with wrong optimal actions. 
# For instance a threshold 94 would never be optimal for pigs with low weight.
dat<-as.data.table( policy) # subset(policy, aLabel=="1")
values<-label2Values(dat$label,param)
values<-cbind(dat,values)
values[aLabel==2 & week>=9 & sW<=100 & sSd<=11.5 & ration==1 & sG>5.2,]  # & sG==5.2 &sSd<=7.5
values[phase==4 & week>=10 & sW==83 & ration==1,])
# sId = 126185 is choosen to examine
id<-   9836
v<-values[sId==id,]
policy<-as.data.table(policy)
policy[sId==id,w1]
# info about the actions (note string output)
a<-info(mdp,idS=id)
# extract action labels:
m<-regexpr('[\\S]*[)]$',a[[1]]$actions, perl=T)
labels<-regmatches(a[[1]]$actions, m)
# better info from the binary files (slow)
aInfo<-actionInfo(prefix)
aI<-aInfo[aInfo$sId==id,]
aI<-aI[,colSums(is.na(aI))<nrow(aI)]  # remove NA columns
prCols<-grep("pr",colnames(aI),value=T)
# some observations: 
# 1) Is some transPr = 0 in the hgf? We should not include transitions with pr = 0. If trans pr = 0
# it must be due to rounding. Is the trans pr in the HMDP the same as in the binray files?
getActionTransPr(mdp,id,0)
aI[1,prCols]
# NO small numbers have been rounded. Due to the way we load the multipliers in the hypergraph (as
# integers which are normalized) How big an effect do this have on calculations using binary files
# trans pr:
res1<-NULL
for (i in 0:1) {
   w<-aI$Reward[i+1]
   t<-aI$Time[i+1]
   idN<-getActionTransIdS(mdp, idS=id, idxA=i)
   prN<-aI[i+1,prCols]
   prN<-prN[!is.na(prN)]
   wN<-getPolicyW(mdp, wLbl, sId = idN)
   V<-w + sum(prN*wN[,2]) - g*t
   res1<-rbind(res1,c(V,w,V-w+g*t) )
}
res1<-cbind(as.data.frame(res1),labels)
colnames(res1) = c("V","r","Exp","label")
res1  # result using correct trans pr
 values[sId==id,]
 NEXT<-values[sId==idN,]
 NEXT$prN<-prN
NEXT<-NULL
# using HMDP trans pr
res1<-NULL
for (i in 0:2) {
   w<-aI$Reward[i+1]
   t<-aI$Time[i+1]
   idN<-getActionTransIdS(mdp, idS=id, idxA=i)
   prN<-getActionTransPr(mdp,id,i)
   wN<-getPolicyW(mdp, wLbl, sId = idN)
   V<-w + sum(prN*wN[,2]) - g*t
   res1<-rbind(res1,c(V,w,V-w+g*t) )
}
res1<-cbind(as.data.frame(res1),labels)
colnames(res1) = c("V","r","Exp","label")
res1  # result using HMDP trans pr
# calculate using policy ite
res<-NULL
for (i in 0:2) {
   resetActions(mdp)
   t<-aI$Time[i+1]
   w<-getActionW(mdp,idS = id,idxA = i)[2]
   fixAction(mdp, sId=id, iA=i)
   g<-policyIteAve(mdp,wLbl,durLbl)
   V<-getPolicyW(mdp, wLbl, sId = id)[1,2]
   res<-rbind(res,c(V,w,V-w+g*t) )
}
resetActions(mdp)
res<-cbind(as.data.frame(res),labels)
colnames(res) = c("V","r","Exp","label")
res
# OK - Seems to give the same results for all methods 2) For all actions we have the same              #i<-1
# transitions and trans pr. Why? Note since the SSMs model the unselected pen they are unaffected by
# the actions. That is only Pr(q_{n+1}|q_n) change. Let us calculate it for all the actions:
res2<-vector("list", 11)
for (i in 0:10) {
   res2[[i+1]]$label<-as.character(res$label[i+1])
   res2[[i+1]]$idN<-getActionTransIdS(mdp, idS=id, idxA=i)
   res2[[i+1]]$idNValues<-label2Values(policy[sId==res2[[i+1]]$idN]$label,param)
   res2[[i+1]]$pr<-getActionTransPr(mdp,id,i)
   if (i==0) { # cont.
      prN<-data.table(q=v$q,prN=1)
   } else {
      if (i==10) { # term.
         prN<-data.table(q=0,prN=1)
      } else {
         prN<-data.table(q=0:v$q,prN=NA)
         prS<-pnorm(eval(parse(text=labels[i+1])), v$sW, sqrt(v$sSd^2+1))
         for (m in 0:v$q) { # m = numb of pigs left at t+1  
            prN$prN[m+1]<-dbinom(m, size=param$pigs, prob=prS)
         }
      }
   }
   res2[[i+1]]$prN<-prN  
   res2[[i+1]]$prN$prN<-res2[[i+1]]$prN$prN/sum(res2[[i+1]]$prN$prN)
}
res2
# Let us have a look at the ones calculated in the package using the function:
# double TransPrN(int & t, int & iTH, int & iSWt, int & iSSdt, int & nt, int & n) {
#    double meanWeight = dis.sW[t-1](iSWt,0);
#    double sDWeight =  dis.sSd[t-1](iSSdt,0) +1; // 1 is observation error
#    double frac= (double) (nt)/(pigs);
#    double limit= R::qnorm(frac,meanWeight,sDWeight,1,0);
#    double probSuccess= 1 - cdflTail(tH[iTH],meanWeight,sDWeight,limit);
#    
#    return (R::dbinom((double)nt-n, (double)nt, probSuccess, 1/*use log transform*/) ); 
# }
library(truncnorm)
res3<-vector("list", 11)
for (i in 0:10) {
   res3[[i+1]]$label<-as.character(res$label[i+1])
   res3[[i+1]]$idN<-getActionTransIdS(mdp, idS=id, idxA=i)
   res3[[i+1]]$idNValues<-label2Values(policy[sId==res3[[i+1]]$idN]$label,param)
   res3[[i+1]]$pr<-getActionTransPr(mdp,id,i)
   if (i==0) { # cont.
      prN<-data.table(q=v$q,prN=1)
   } else {
      if (i==10) { # term.
         prN<-data.table(q=0,prN=1)
      } else {
         frac<-v$q/param$pigs
         limit<-qnorm(frac, v$sW, sqrt(v$sSd^2+1))
         prN<-data.table(q=0:v$q,prN=NA)
         prS<-1-ptruncnorm(eval(parse(text=labels[i+1])) , a=-Inf, b=limit, v$sW, sqrt(v$sSd^2+1)) 
         for (m in 0:v$q) { # m = numb of pigs left at t+1  
            prN$prN[m+1]<-dbinom(v$q-m, size=v$q, prob=prS)
         }
      }
   }
   res3[[i+1]]$prN<-prN  
   #res3[[i+1]]$prN$prN<-res2[[i+1]]$prN$prN/sum(res2[[i+1]]$prN$prN)
}
res3
# This is a totally different result!!!

KOMMET TIL



# we have transitions to the same states



# - For threshold > 94 the cost is equal. Why?
# - For all thresholds the transitions and trans pr are the same. Why?


# next step is to have a look at the code where the actions is generated, i.e. functions CalcTransPrTH
# We first take a look at the trans pr
aInfo<-actionInfo(prefix)
tmp<-aInfo[aInfo$sId==id,]
tmp
# Note that some trans pr is at e-13 (very low). Moreover, since the SSMs model the unselected pen they are unaffected by the actions. Let us have a look at the last pr Pr(q_{n+1}|q_n) which for continue is 1 for 
values[sId==id,] # we see that q_n = 3
for (i in 0:10) {
   getActionTransIdS(mdp, idS=id, idxA=i)
   w<-getActionW(mdp,idS = id,idxA = i)
   fixAction(mdp, sId=id, iA=i)
   policyIteAve(mdp,wLbl,durLbl)
   V<-getPolicyW(mdp, wLbl, sId = id)[1,2]
   res<-rbind(res,c(V,w[2],V-w[2]) )
}





# ... and threshold 94 would never be optimal for pigs with high weight, since the terminating action should be better.
values[values$phase==4 & values$week==10 & values$sW>=103 & values$ration==1,]
# when is terminate optimal
dat1<-subset(policy, aLabel=="terminate") 
values1<-label2Values(dat1$label,param)
values1[values1$phase==4 & values1$week==10 & values1$sW>=103 & values1$ration==1,]





# OLD STUFF






# We know the best live weight for slaughter, according to Anders's paper, is 109 and hence the threshold weights that are less
# than 109 are not reasonable. To find the states that the optimal action is a threshold less than 109 we follow this procedure:

subset(policy, n0==0 & s0==0 & a0==0 & n1==2  &  s1==2732 ) & a1==0 & n2==10-param$minPhaseT[n1+1] &  aLabel=="94") 


# The steps that we need to find the states with threshold action 94:

# Step 1: Use "subset" function to filter all the states with threshold 94 in the "policy" data frame:
# subset(policy, n0==0 & s0==0 & a0==0 & n1==?  &  s1==0 & a1==? & n2==t-param$minPhaseT[n1+1] &  aLabel=="94") 
# In the subset function: 
# n1 is the phase number (from 0 to 2); s1 is always 0 (index of dummy state in second level); a1 is the index of feed-mix (from 0 to 2)-exception: when n1 is 0 always a1 must be 0);
# n2 is the statge number in the last level of the HMDP (stage number: t-param$minPhaseT[n1+1] where t is the week number); aLabel shows the optimal action.




# Now we can do the above steps for threshold 94:

#Step1: we assume we are in week 10(t=10) of phase 0 (n1=0) with feed-mix 2 (a1=0) (Point: when we are in phase 0 (n1==0) then always a1 must be 0 because just we have one defined feed-Mix for phase 0) 
   
subset(policy, n0==0 & s0==0 & a0==0 & n1==0  &  s1==0 & a1==0 & n2==10-param$minPhaseT[n1+1] &  aLabel=="94") 

# Step 2: Now we select one of the labels. For instance label= (5,3,3,1,1,0,2,10).
# Step 2-1: This label shows that we have threshold 94 when number of pigs is 5 (n=5);
# index of weight, growth and sd are 3,3 and 1 respectively; starting time of using current feed-mix is week 1 (RS=1), phase number is 0 (phase=0); 
# current feed-mix is 2 (iRation=2); and week number is 10 (t=10).
# Step 2-2: Based on the index of weight, growth and sd (3,3,1), the related center points in week 10 are: 
# Weight: 98 kg in the list param$wInterval
# Growth: 6.8 kg in the list param$gInterval
# Sd: 9.5 in the list param$sDInterval


# Point 1: In the filtered labels (states) for threshold 94, when we have low velues for iSW it means that we have wrong decisions (threshold 94) with low weights and 
# when we have high velues for iSW in the related labels, it shows that we have wrong decisions (threshold 94) with high weights.   


# Point 2: If the model shows an optimal threshold less than 109, this shows that we have some errors according to Anders's comments.  


