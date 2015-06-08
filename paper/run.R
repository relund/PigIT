## R script file for getting results used in the paper. Do the following
## 
## 1) Find the optimal policy 
## 2) Create data for the 3 pens (stored in 6 csv files)
## 3) Plot the results (stored in pdf files)

source("optimize_hmdp.R")   # find optimal policy
source("simulate.R")        # simulate the 3 pens based on the optimal policy
source("plot.R")            # plot the results


