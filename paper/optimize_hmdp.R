## R file for solving the HMDP 

library(hmdpFeedPigIT)
set.seed(234567)

# Set HMDP parameters: 
param<-setParameters(pigs=5, 
                     rations=3,
                     phases=4,
                     tMax=12,
                     tStartMarketing=9,
                     iniFeedMix=1,
                     minPhaseT=c(1,4,7,10),
                     disWeight=c(5,4),
                     disSD=c(3,0.8),
                     disGrowth=c(2,0.3),
                     priorGrowth=c(5.8,6.3,6.8),
                     cullActions = T
)

# Create GSSM
dlms<-list()
for (i in 1:param$rations){
  dlmMod<-iniDLM()
  Y<-matrix(c(28.91122,11.17050, 34.73722, 11.57909, 40.55647, 11.91425, 46.43588, 12.31422, 52.60592, 12.76750, 58.90283, 13.18425, 64.73567, 13.57245, 70.88710, 13.94777, 76.19884, 14.36087, 82.30708, 14.76617, 88.41322, 15.10727, 94.68169,  15.56838, 100.68169, 15.96838, 106.68169, 16.57838, 112.68169, 17.06838), nrow=2)
  dlms[[i]]<-buildDLM(dlmMod, Y)
}

# Get the nGSSM
dglmParam<-iniDGLM()

# Build the HMDP model
prefix<-"hmdp_"
BuildHMDP(prefix, param, dlms, dglmParam)

# Solve the HMDP model
wLbl<-"Reward"
durLbl<-"Time"
mdp<-loadMDP(prefix)
g<-policyIteAve(mdp,wLbl,durLbl)      # Finds the optimal policy using the average reward per week (g) criterion 
policy<-getPolicy(mdp) # optimal policy for each sId

do.call(file.remove,list(list.files(pattern = prefix)))