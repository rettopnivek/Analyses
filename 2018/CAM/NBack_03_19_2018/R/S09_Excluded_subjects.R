# Checks for subjects excluded from fNIRS analyses
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you 
# have any questions or comments
# Last updated 2018-05-25

# Table of contents
# 1) Initial setup
# 2) Script that extracted behavioral data to analyze
# 3) Comparison of behavioral/fNIRS subjects
# 4) 

###
### 1) Initial setup
###

# Save current working directory
orig_dir = getwd()

# Navigate to project directory
setwd( '..' )
proj_dir = getwd()

# Load in useful packages

# Collection of useful functions
# devtools::install_github("rettopnivek/utilityf")
library( utilityf )

# Package for easier manipulation of data frames
# install.packages( 'dplyr' )
library( dplyr )

# Package for estimating mixed effects models
# install.packages( 'lme4' )

# Load in data
setwd( 'Data' )
# Behavioral data
load( "NBack_3_19_2018.RData" )
# fNIRS data
setwd( 'Original_files' )
fd = read.csv( 
  file = "self-report_Mereyem.csv",
  header = T,
  stringsAsFactors = FALSE )
setwd( proj_dir )

# Load in useful functions
setwd( 'R' )
source( 'S02_Useful_functions.R' )
setwd( proj_dir )

###
### 2) Script that extracted behavioral data to analyze
###

# Exclude combined data
bd = dat %>% 
  filter( Task != 'Combined' )

# Drop rows with missing 
# self-report data for subject 
# FN_041 (66)
# Subject FN_041 did not have any 
# self report data on how high
# for the T1 drug condition
bd = bd %>% 
  filter( Self_report_on_high != '-' )

# Determine count data for misses/correct rejections
bd$Z = bd$Trials - bd$Counts
bd$Y = bd$Counts

# Break self-reported high into 
# discrete bins
SRH_range = seq( 0, 100, 20 )

bd$SRH_bins = 0
for ( i in 2:length( SRH_range ) ) {
  sel = bd$Self_report_on_high > SRH_range[i-1] & 
    bd$Self_report_on_high <= SRH_range[i]
  bd$SRH_bins[sel] = SRH_range[i]
}

###
### 3) Comparison of behavioral/fNIRS subjects
###

behav = unique( bd$ID )
fNIRS = unique( fd$Subject )
es = behav[ !(behav %in% fNIRS) ]

###
### 4) Details on excluded subjects
###

# Isolate excluded subjects
cd = bd[ bd$ID %in% es, ]
# Measure of accuracy
cd$Accuracy = ( cd$Y / cd$Trials ) * cd$Target + 
  ( cd$Z / cd$Trials ) * ( 1 - cd$Target )
check = cd %>% 
  group_by( ID, Condition ) %>% 
  summarize(
    Ac = mean( Accuracy ),
    fNIRS = mean( (R_DLPFC + L_DLPFC + MPFC + R_VLPFC + L_VLPFC)/5 ),
    SR = mean( Self_report_on_high )
  )

p = check$ID[ check$Condition == 'Placebo' ]
d = check$ID[ check$Condition == 'Drug' ]
# Placebo only
# print( p[ !(p %in% d) ] )
# 33, 34, 74, 95, and 96
# Drug only
# 58, 78, 86, 88, and 91
# print( d[ !(d %in% p) ] )
# Both
# print( p[ p %in% d ] )
# 41 and 69

# fNIRS subjects
cd = bd[ !(bd$ID %in% es), ]
# Measure of accuracy
cd$Accuracy = ( cd$Y / cd$Trials ) * cd$Target + 
  ( cd$Z / cd$Trials ) * ( 1 - cd$Target )
check2 = cd %>% 
  group_by( ID, Condition ) %>% 
  summarize(
    Ac = mean( Accuracy ),
    fNIRS = mean( (R_DLPFC + L_DLPFC + MPFC + R_VLPFC + L_VLPFC)/5 ),
    SR = mean( Self_report_on_high )
  )

p = check2$ID[ check2$Condition == 'Placebo' ]
d = check2$ID[ check2$Condition == 'Drug' ]
# Placebo only
# print( p[ !(p %in% d) ] )
# 33, 34, 74, 95, and 96
# Drug only
# 58, 78, 86, 88, and 91
# print( d[ !(d %in% p) ] )
# Both
# print( p[ p %in% d ] )
# 41 and 69

setwd( orig_dir )