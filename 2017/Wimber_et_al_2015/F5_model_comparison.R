#---------------------#
# Model comparisons   #
# Kevin Potter        #
# Updated 01/16/2017  #
#---------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

# Indicate whether to save a pdf file
savePlot = T
if (savePlot) {
  setwd( 'Plots' )
  fName = 'Model_comparisons.pdf'
  pdf( fName, width = 12, height = 6 )
  rm( fName )
  setwd( orig_dir )
}

# Index
# Lookup - 01:  Load in useful packages and data

###
### Load in useful packages and data
###
# Lookup - 01

# For geting github packages
# install.packages(devtools)
# library(devtools)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

# Load in package for Bayesian estimation
library(rstan)
# For parallel processing
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Load in package for model comparison
library(loo)

# Load in data
setwd( 'Data' )
load( 'Original_all_data.RData' )
setwd( orig_dir )

# Extract observations for final recogniton memory 
# test and pre-process data

d = OriginalAllData[ OriginalAllData$Cond == 6, ]
# For easy manipulation
colnames( d ) = c('S','Tr','C','IN','Co','Ch','RT',
                  'Ac','IT','B','Bl','Cat','CR','fN')
d$I = createIncrement( d$IN )

# Use dummy coding for responses
d$Ch = d$Ch - 1; d$Co = d$Co - 1

# Define a variable denoting conditions

# Targets that underwent selective retrieval (T-SR)
d$Cnd = 1;
# Targets in the baseline condition (T-B)
d$Cnd[ d$IT == 1 & d$B == 1 ] = 2
# Competitors that underwent selective retrieval (C-SR)
d$Cnd[ d$IT == 2 & d$B == 0 ] = 3
# Competitors in the baseline condition (C-B)
d$Cnd[ d$IT == 2 & d$B == 1 ] = 4

# Set missing data to be incorrect
d$Ac[ is.na(d$RT) ] = 0
d$Ch[ is.na(d$RT) ] = 1 - d$Co[ is.na(d$RT) ]

###
### Load in posterior samples
###
# Lookup - 02

mNames = c(
  'SDT_model_M1_post.RData',
  'SDT_model_M2_post.RData',
  'SDT_model_M3_post.RData',
  'mSDT_model_no_item_post.RData' )

waic_all = c()
loo_all = c()

for (i in 1:4 ) {
  
  # Load in posterior
  fName = "C:/Users/Kevin/Documents/Posteriors from Stan"
  setwd(fName)
  setwd( "Wimber_et_al" )
  load( mNames[i] )
  setwd( orig_dir )
  
  waic_all = c( waic_all, list( waic( post$logLik ) ) )
  loo_all = c( loo_all, list( loo( post$logLik ) ) )
  
}
names( waic_all ) = c('M1','M2','M3','M4')

###
### Aikake weights
###

ll = c(
  # waic_all$M1$waic,
  waic_all$M2$waic,
  waic_all$M3$waic,
  waic_all$M4$waic
)
lbls = c(
  'Main effects only',
  'SR-C',
  'Mixture' )

ll = ll - min( ll )
w = exp( -.5*ll )
w = w/sum(w)

if (!savePlot) x11()
plot( c(1,length(w)), c(0,1), type = 'n', xlab = ' ',
      ylab = 'Akaike weights', xaxt = 'n', bty = 'n' )
abline( h = 0 )
axis( 1, 1:3, lbls, tick = F )
points( 1:length(w), w, pch = 19 )

###
### Pair-wise comparisons
###

m2_v_m4 = compare( waic_all$M2, waic_all$M4 )
m3_v_m2 = compare( waic_all$M3, waic_all$M2 )
m3_v_m4 = compare( waic_all$M3, waic_all$M4 )

ya = c(m2_v_m4[1],
       m3_v_m2[1],
       m3_v_m4[1])
u = c(m2_v_m4[2],
      m3_v_m2[2],
      m3_v_m4[2])
ui = cbind(
  ya - qnorm( .975 )*u,
  ya + qnorm( .975 )*u )

if (!savePlot) x11()
par( mar = c( 5, 4, 4, .5 ) )
plot( c(0,4), c(-100,20), type = 'n', bty = 'n',
      xlab = ' ', xaxt = 'n', ylab = 'WAIC difference' )
abline( h = 0 )
axis( 1, 1:3, lbls[c(2,2,1)], tick = F )
axis( 3, 1:3, lbls[c(3,1,3)], tick = F )

segments( 1:3, ui[,1], 1:3, ui[,2] )
points( 1:3, ya, pch = 19 )

if (savePlot) dev.off()