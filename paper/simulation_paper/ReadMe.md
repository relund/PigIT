# Simulation results used in the paper

This folder contains simulated data of three pens with different genetics properties. Pen 2 has a normal growth rate and Pens 1 and 3 grow twenty percent slower and faster than Pen 2. Simulated data for each pen has been stored in the "penDaily" and "penWeekly" files.

## penDaily 

Includes daily information related to daily simulation of individual pigs in the pen. Observed data in this file is daily weight and feed intake of individual pigs. Main fields describing the data in this file are:  

  t: Day number
  pig: Id of an individual pig in the pen. 
  FI: Feed intake of an individual pig in the pen.
  OLW: Observed live weight of an individual pig in the pen.
  culled: A binary value showing that a specific  pig has been culled from the pen or not (resulted from the optimal policy of HMDP).
  week: Week number. 
  pen: Pen number. 
  feedMix: The feed-mix number that is currently used in the pen (resulted from the optimal policy of HMDP).

## penWeekly

Contains weekly information in the pen level. Observed weekly weigh and feed intake in the pen level were calculated using the file "penDaily" that are inputs to the GSSM and nGSSM. Estimated mean and standard deviation of weight and growth are the outputs of GSSM and nGSSM. Moreover the optimal decisions of marketing and feeding have been shown. Main fields describing the data in this file are:       

  t: Day number in the pen.
  aveOLWALL: The observed average live weight in the pen level (for an unselected pen).  
  sdOLWALL: The observed standard deviation of live weight in the pen level (for an unselected pen).
  aveFIAll: The average feed intake in the pen level (for an unselected pen).
  alive: Number of pigs that has been remained in the pen.
  pen: Pen number.
  week: Week number.
  weekDay: Number of days in a week (always equal to 7).
  stage: Stage number in the HMDP.
  eAveOLWALL: Estimated mean of live weight in the pen level (output of GSSM).  
  eAveGALL: Estimated mean of growth rate in the pen level (output of GSSM).
  eSdOLWALL: Estimated standard deviation of live weight in the pen level (output of nGSSM).
  optFeedMix: Optimal feeding decision resulted from the optimal policy of HMDP. 
  optThreshold: Optimal marketing decision resulted from the optimal policy of HMDP. 
