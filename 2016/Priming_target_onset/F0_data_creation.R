#-----------------------------#
# Script to create R data set #
# Kevin Potter                #
# Updated 11/11/16            #
#-----------------------------#

# Clear workspace
rm(list = ls())

# Save current working directory
orig_dir = getwd()

###
### Load in useful packages
###

# For geting github packages
# install.packages(devtools)
# library(devtools)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

###
### Read in .csv file
###

# Change directory to where data is stored
setwd( 'Data' )

# Load in .csv file with data
rawDat = read.table( file = 'Offset.csv', sep = ',', header = F )

# Columns for rawDat
# 1: Subject index
# 2: Prime duration ( 0 = 50 ms, 1 = 2000 ms )
# 3: Which word was primed ( 0 = target, 1 = foil )
# 4: Offset duration ( 0 - 4: 0 = target 33 ms before foil; 
#    1 = target 17 ms before foil; 2 = same time; 
#    3 = target 17 ms after foil; 4 = target 33 ms after foil )
# 5: Correct answer ( 0 = left, 1 = right )
# 6: Response time (ms)
# 7: Choice ( 0 = left, 1 = right )

# Create meaningful column names
colnames( rawDat ) = c(
  'Subject',
  'PrimeDur',
  'PrimeType',
  'Offset',
  'Correct',
  'RT',
  'Choice'
)

###
### Add in additional variables
###

# Determine accuracy per trial ( 0 = error, 1 = correct )
rawDat$Accuracy = as.numeric( rawDat$Correct == rawDat$Choice )

# Include exact time of offsets in seconds
rawDat$OffsetTime = -33/1000
rawDat$OffsetTime[ rawDat$Offset == 1 ] = -17/1000
rawDat$OffsetTime[ rawDat$Offset == 2 ] = 0
rawDat$OffsetTime[ rawDat$Offset == 3 ] = 17/1000
rawDat$OffsetTime[ rawDat$Offset == 4 ] = 33/1000

# Convert response times to seconds
rawDat$RT = rawDat$RT/1000

# Change increment for offset
rawDat$Offset = rawDat$Offset + 1

# Create a set of covariates

# Covariate DurxType:
# Collapse over offset and type of correct response
CB = covCreate( cbind( rawDat$PrimeDur, rawDat$PrimeType ) );
tmp = aggregate( CB, list(rawDat$PrimeDur, rawDat$PrimeType), unique );
colnames( tmp ) = c('PD','PT','CB')
rawDat$DurxType = 1;
for (i in 2:nrow(tmp)) rawDat$DurxType[ CB == tmp$CB[i] ] = i;

# Prime duration | type   | DurxType
#          50 ms | Target | 1
#        2000 ms | Target | 2
#          50 ms | Foil   | 3
#        2000 ms | Foil   | 4

# Clean up workspace
rm( CB, tmp, i )

# Covariate CorxDurxType:
# Collapse over offset
CB = covCreate( cbind( rawDat$PrimeDur, rawDat$PrimeType, rawDat$Correct ) );
tmp = aggregate( CB, list(
  rawDat$Correct, rawDat$PrimeDur, rawDat$PrimeType ), unique );
colnames( tmp ) = c('Co','PD','PT','CB')
rawDat$CorxDurxType = 1;
for (i in 2:nrow(tmp)) rawDat$CorxDurxType[ CB == tmp$CB[i] ] = i;

# Correct side | Prime duration | type   | CorxDurxType
#         Left |          50 ms | Target | 1
#        Right |          50 ms | Target | 2
#         Left |        2000 ms | Target | 3
#        Right |        2000 ms | Target | 4
#         Left |          50 ms | Foil   | 5
#        Right |          50 ms | Foil   | 6
#         Left |        2000 ms | Foil   | 7
#        Right |        2000 ms | Foil   | 8

# Clean up workspace
rm( CB, tmp, i )

# Covariate Cond:
# All conditions
CB = covCreate( cbind( rawDat$PrimeDur, rawDat$PrimeType, 
                       rawDat$Correct, rawDat$Offset ) );
tmp = aggregate( CB, list(
  rawDat$Correct, rawDat$PrimeDur, rawDat$PrimeType, rawDat$Offset ), unique );
colnames( tmp ) = c('Co','PD','PT','O','CB')
rawDat$Cond = 1;
for (i in 2:nrow(tmp)) rawDat$Cond[ CB == tmp$CB[i] ] = i;

# Correct side | Prime duration | type   | Cond
#              |                |        | Offset 1  2  3  4  5
#--------------------------------------------------------------
#         Left |          50 ms | Target |        1  9 17 25 33
#        Right |          50 ms | Target |        2 10 18 26 34
#         Left |        2000 ms | Target |        3 11 19 27 35
#        Right |        2000 ms | Target |        4 12 20 28 36
#         Left |          50 ms | Foil   |        5 13 21 29 37
#        Right |          50 ms | Foil   |        6 14 22 30 38
#         Left |        2000 ms | Foil   |        7 15 23 31 39
#        Right |        2000 ms | Foil   |        8 16 24 32 40

# Clean up workspace
rm( CB, tmp, i )

# Determine total number of subjects
N = length( unique( rawDat$Subject ) )

# Save extracted data
fName = 'Priming_offset.RData'
save(rawDat,N,file=fName)