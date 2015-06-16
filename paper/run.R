## R script file for getting results used in the paper. 
## Do the following
## 
## 1) Find the optimal policy 
## 2) Create data for the 3 pens (stored in 6 csv files)
## 3) Plot the results (stored in pdf files)

# remember to set the working dir til ./paper/
library(hmdpFeedPigIT)
useSimPaper <- TRUE   # use the already simulated data?
if (!useSimPaper) {
  message("Use a new simulation.")
  source("optimize_hmdp.R")   # find optimal policy
  source("simulate.R")        # simulate the 3 pens (use optimal policy to identify feed-mix and number of pigs)
} else {
  message("Use the data generated in the simulation_paper folder.")
  pen1Weekly <- read.csv2("simulation_paper/pen1Weekly.csv")
  pen1Daily <- read.csv2("simulation_paper/pen1Daily.csv")
  pen2Weekly <- read.csv2("simulation_paper/pen2Weekly.csv")
  pen2Daily <- read.csv2("simulation_paper/pen2Daily.csv")
  pen3Weekly <- read.csv2("simulation_paper/pen3Weekly.csv")
  pen3Daily <- read.csv2("simulation_paper/pen3Daily.csv")
}
source("plot.R", echo = TRUE)    # plot results as tex files
tools::texi2pdf(file = "sim_plot.tex", clean = T)
tools::texi2pdf(file = "opt_plot.tex", clean = T)

