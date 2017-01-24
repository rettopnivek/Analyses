#---------------------------------#
# Data creation and preprocessing #
# Kevin Potter                    #
# Updated 01/18/2017              #
#---------------------------------#

# Clear workspace
rm( list = ls() )

# Save current directory
orig_dir = getwd()

# Load in useful packages

# For geting github packages
# install.packages(devtools)
# library(devtools)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

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

# Save results to file
setwd( 'Data' )
save( d, file = 'Recog_mem.RData' )
setwd( orig_dir )