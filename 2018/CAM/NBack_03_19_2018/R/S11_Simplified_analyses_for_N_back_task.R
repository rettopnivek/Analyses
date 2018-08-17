# Simplified analyses for the N-back task
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you 
# have any questions or comments
# Last updated 2018-08-06

# Table of contents
# 1) Initial setup
# 2) Analyses and figures

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

# Packages for easier manipulation of data frames
# install.packages( 'dplyr' )
library(dplyr)

# Package for Bayes factor analyses
# install.packages( 'BayesFactor' )
library(BayesFactor)

# Load in data
setwd( 'Data' )
load( "NBack_3_19_2018.RData" )
setwd( proj_dir )

# Custom functions
setwd( 'R' );
source( 'S02_Useful_functions.R' )
setwd( proj_dir )

# Create a new data frame that
# a) Excludes the combined trials
# b) Has separate variables for hits and false alarms
all_ptbp = dat %>% 
  filter( Task != 'Combined' & Response_type == 'Hits' ) %>% 
  arrange( Task, Condition, Timepoints, Subject ) %>% 
  mutate( H = Counts / Trials,
          fH = Counts, 
          Npos = Trials )
tmp = dat %>% 
  filter( Task != 'Combined' & Response_type == 'False_alarms' ) %>% 
  arrange( Task, Condition, Timepoints, Subject ) %>% 
  mutate( FA = Counts / Trials,
          fFA = Counts,
          Nneg = Trials )
all_ptbp$FA = tmp$FA
all_ptbp$fFA = tmp$fFA
all_ptbp$Nneg = tmp$Nneg

# Clean up workspace
rm( tmp )

# Convert the hit/false alarm rates into 
# estimates of d' and bias
all_ptbp = all_ptbp %>% 
  mutate( dp = sdt_calc_binary( fH, Npos, fFA, Nneg, 1, T ),
          crt = sdt_calc_binary( fH, Npos, fFA, Nneg, 1, F ) )

subj = unique( all_ptbp$ID )
dtba = data.frame(
  ID = rep( subj, 3 * 2 * 2 ),
  Time_point = rep( rep( 1:3, each = length(subj) ), 4 ),
  Task = rep( rep( c( '0-back', '2-back' ), each = length(subj) * 3 ), 2 ), 
  Condition = rep( c( 'Placebo', 'Drug' ), each = length( subj ) * 3 * 2 ),
  d_prime = NA,
  criterion = NA,
  Freq_hits = NA,
  Positive_trials = NA,
  Freq_false_alarms = NA,
  Negative_trials = NA,
  SRH = NA,
  Visit = NA, 
  stringsAsFactors = F
)

tsk = unique( all_ptbp$Task )
tm = unique( all_ptbp$Timepoints )
for ( n in 1:nrow( all_ptbp ) ) {
  # Determine row in main data frame
  if ( all_ptbp$Task[n] == tsk[1] ) tsk_cur = '0-back' else tsk_cur = '2-back'
  tm_cur = which( tm == all_ptbp$Timepoints[n] )
  sel = all_ptbp$ID[n] == dtba$ID & 
    all_ptbp$Condition[n] == dtba$Condition & 
    tsk_cur == dtba$Task & 
    tm_cur == dtba$Time_point
  
  # Fill in data
  dtba$d_prime[sel] = all_ptbp$dp[n]
  dtba$criterion[sel] = all_ptbp$crt[n]
  
  dtba$Freq_hits[sel] = all_ptbp$fH[n]
  dtba$Positive_trials[sel] = all_ptbp$Npos[n]
  dtba$Freq_false_alarms[sel] = all_ptbp$fFA[n]
  dtba$Negative_trials[sel] = all_ptbp$Nneg[n]
  
  dtba$SRH[sel] = all_ptbp$Self_report_on_high[n]
  
  dtba$Visit[sel] = all_ptbp$Visit[n]
}

# Standardize scores
sel = !is.na( dtba$SRH )
dtba$SRH[sel] = my_standardize( dtba$SRH[sel] )

# Clean up workspace
rm( subj, tsk_cur, tsk, tm, tm_cur, sel, n )

###
### 2) Analyses and figures
###

setwd( proj_dir )
setwd( 'Figures' )
pdf( file = 'Simple_analyses_of_N_back.pdf', width = 12 )
setwd( proj_dir )

# Compute difference scores between T1 and T2
T1 = dtba %>% 
  filter( Time_point == 1 & Task == '2-back' ) %>% 
  select( ID, Condition, d_prime, SRH )
T2 = dtba %>% 
  filter( Time_point == 2 & Task == '2-back' ) %>% 
  select( ID, Condition, d_prime, SRH )

by_time_point = data.frame(
  ID = T1$ID,
  Condition = T1$Condition,
  x = T2$d_prime - T1$d_prime,
  SRH = T2$SRH - T1$SRH,
  Visit = NA,
  stringsAsFactors = F
)
for ( i in 1:nrow( by_time_point ) ) {
  sel = by_time_point$ID[i] == dtba$ID & 
    by_time_point$Condition[i] == dtba$Condition
  by_time_point$Visit[i] = unique( dtba$Visit[sel] )
}
by_cond = data.frame(
  ID = by_time_point$ID[ by_time_point$Condition == 'Placebo' ],
  x = by_time_point$x[ by_time_point$Condition == 'Drug' ] - 
    by_time_point$x[ by_time_point$Condition == 'Placebo' ],
  SRH = by_time_point$SRH[ by_time_point$Condition == 'Drug' ] - 
    by_time_point$SRH[ by_time_point$Condition == 'Placebo' ],
  Visit = by_time_point$Visit[ by_time_point$Condition == 'Drug' ] - 
    by_time_point$Visit[ by_time_point$Condition == 'Placebo' ],
  stringsAsFactors = F
)
sel = !is.na( by_cond$x )
by_cond = by_cond[sel,]
by_cond$Visit = -by_cond$Visit

# Frequentist approach
tt = t.test( by_cond$x )
# Pretty results
paste( 't(', tt$parameter, ') = ', round( tt$statistic, 2 ),
       ', p = ', round( tt$p.value, 3 ), sep = '' )

# Plot results
# x11( width = 12 )
layout( cbind( 1, 2 ) )

dn = densityPoints( by_cond$x )
xl = lowerUpper( 1, dn$x )
yl = c( 0, lowerUpper( .5, dn$y )[2] )
blankPlot( xl, yl )

vertLines( 0, yl, lty = 2, lwd = 2 )
vertLines( mean( by_cond$x ), yl, col = 'black', lwd = 2 )

lines( dn$x, dn$y, type = 'b', lwd = 2, pch = 19 )

customAxes( xl, yl )
axis( 1, seq( xl[1], xl[2], 1 ),
      tick = F, line = -1.5, cex.axis = 1.25 )
mtext( "Difference scores for d'", side = 1, line = 2, cex = 1.25 )
mtext( "Density", side = 2, line = 2, cex = 1.25 )

legend( 'topright', 'Mean', lwd = 2, col = 'black', lty = 1,
        bty = 'n', cex = 1.25 )

# Controlling for self-reported high
lmf = lm( x ~ 1 + SRH, data = by_cond )
sm = summary( lmf )
# Pretty results
paste( 'Beta(Drug) = ', round( sm$coefficients[1,1], 2 ), 
       ', p = ', round( sm$coefficients[1,4], 3 ), sep = '' )
paste( 'Beta(Self-reported high) = ', round( sm$coefficients[2,1], 2 ), 
       ', p = ', round( sm$coefficients[2,4], 3 ), sep = '' )

xl = lowerUpper( 1, by_cond$SRH )
yl = c( -4, 4 )
blankPlot( xl, yl )

nd = by_cond
nd = nd %>% arrange( SRH )
nd$x = NA
pred = predict( lmf, newdata = nd, interval = 'prediction' )
polygon( c( nd$SRH, rev( nd$SRH ) ),
         c( pred[,2], rev( pred[,3] ) ),
         col = 'grey', border = NA )
horizLines( 0, xl, lty = 2, lwd = 2 )
lines( nd$SRH, pred[,1], lwd = 2 )

points( by_cond$SRH, by_cond$x, pch = 19 )

customAxes( xl, yl )
axis( 1, seq( xl[1], xl[2], 1 ),
      tick = F, line = -1.5, cex.axis = 1.25 )
mtext( 'Difference scores for self-reported high',
       side = 1, line = 2, cex = 1.25 )
axis( 2, seq( yl[1], yl[2], 1 ),
      tick = F, line = -1.5, cex.axis = 1.25 )
mtext( "Difference scores for d'",
       side = 2, line = 2, cex = 1.25 )

r = cor.test( by_cond$SRH, by_cond$x )
legend( 'topright', 
        paste( 'R = ', round( r$estimate, 2 ),
               ', p = ', round( r$p.value, 3 ), sep = '' ),
        bty = 'n', cex = 1.25 )

dev.off()

setwd( orig_dir )

