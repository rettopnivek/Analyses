#--------------------#
# Data creation      #
# Kevin Potter       #
# Updated 02/28/2017 #
#--------------------#

# Clear workspace
rm(list = ls())

# Save current working directory
orig_dir = getwd()

# Index
# Lookup - 01:  Load in useful packages
# Lookup - 02:  Read in data and create/adjust variables
# Lookup - 03:  Identify subjects who may have flipped responses
# Lookup - 04:  Data key
# Lookup - 05:  Extract hits/false alarms for subjects
# Lookup - 06:  Save results as an .RData file

###
### Load in useful packages
###
# Lookup - 01

# For geting github packages
# install.packages(devtools)
# library(devtools)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

###
### Read in data and create/adjust variables
###
# Lookup - 02

# Change directory to where data is stored
setwd( 'Data' )

# Read in .csv file
allDat = read.table( file = 'DRMdataforKevin4.csv',
                     sep = ',', header = T )

# Create variable storing original subject ID
allDat$ID = allDat$subject

# Adjust subject numbers to be sequential
allDat$subject = createIncrement( allDat$subject )

# Determine sample size
N = max( allDat$subject )

# Create numerical code for age group (dummy coding)
allDat$Age = 0
allDat$Age[ allDat$group == 'old' ] = 1

# Create numerical code for response category
allDat$Old = 0
allDat$Old[ allDat$old_new == 'old' ] = 1

# Create two numeric variables tracking conditions

# Cond14 allows unrelated lures to differ over 
# semantic/phonetic conditions to produce 14 conditions
# total:
tmp = aggregate( allDat$acc, 
                 list( allDat$class, 
                       allDat$list_type, 
                       allDat$list_length ), mean )
colnames( tmp ) = c( 'class', 'list_type', 'list_length', 'x' )

allDat$Cond14 = NA
for ( i in 1:14 ) {
  allDat$Cond14[ allDat$class == tmp[ i, 1 ] & 
                   allDat$list_type == tmp[ i, 2 ] & 
                   allDat$list_length == tmp[ i, 3 ] ] = i
}

# Cond13 collapses unrelated lures over semantic/phonetic
# conditions to produce 13 conditions total:
allDat$Cond13 = 1
for ( i in 1:12 ) {
  allDat$Cond13[ allDat$class == tmp[ i + 2, 1 ] & 
                 allDat$list_type == tmp[ i + 2, 2 ] & 
                 allDat$list_length == tmp[ i + 2, 3 ] ] = 1 + i
}

# Create a set of labels for the 14 conditions
no = nrow( allDat ) # Total number of observations
p1 = rep( 'URL', no ) # Abbreviate class type
p1[ allDat$class == 'critical' ] = 'CL'
p1[ allDat$class == 'related' ] = 'RL'
p1[ allDat$class == 'target' ] = 'T'
p2 = rep( 'P', no ) # Abbreviate list type
p2[ allDat$list_type == 'semantic' ] = 'S'
p3 = rep( '0', no ) # Abbreviate list length
p3[ allDat$list_length == 2 ] = '2'
p3[ allDat$list_length == 8 ] = '8'

# Collapse all vector of abbreviated strings
allDat$labels14 = paste( p1, '_', p2, '_', p3, sep = '' );

# To check if labels are correct
# aggregate( allDat$labels14, list( allDat$Cond14 ), unique )

# Create a set of labels for the 13 conditions
allDat$labels13 = allDat$labels14
allDat$labels13[ allDat$labels13 == 'URL_P_0' | 
                   allDat$labels13 == 'URL_S_0' ] = 'URL_0';

# To check if labels are correct
# aggregate( allDat$labels13, list( allDat$Cond13 ), unique )

###
### Identify subjects who may have flipped responses
###
# Lookup - 03

# Calculate average accuracy for target words
sel = allDat$class == 'target'
ac = aggregate( allDat$acc[sel], list( allDat$subject[sel] ), mean )
colnames( ac ) = c('S','A')
# Identify subjects who had accuracy significantly lower than chance
flipped = which( ac$A < .35 )
# Flip the responses of these subjects
for ( s in flipped ) {
  sel = allDat$subject == s
  allDat$acc[ sel ] = 
    1 - allDat$acc[ sel ]
  allDat$Old[ sel ] = 
    1 - allDat$Old[ sel ]  
  allDat$rating[ sel ] = 
    7 - allDat$rating[ sel ]  
  tmp = allDat$old_new[ sel ]
  allDat$old_new[ sel ][ tmp == 'new' ] = 'old'
  allDat$old_new[ sel ][ tmp == 'old' ] = 'new'
}

###
### Data key
###
# Lookup - 04

# subject:     Index for the subject number (recoded to be sequential)
# group:       Indicates whether participants are young adults or 
#              older
# hostname:    Name of computer on which experiment was conducted
# word:        The test word
# list_type:   The list type to which a word belonged (semantic versus
#              phonetic)
# list_length: The number of items on the original study list
# class:       The class of word (lure, related, critical, or target)
# rating:      The confidence rating response of the subject, where...
#              1 = Certain new,
#              2 = Probably new,
#              3 = Maybe new,
#              4 = Maybe old,
#              5 = Probably old,
#              6 = Certain old,
# old_new:     The binary classification of a subject's response as 
#              old or new
# acc:         Accuracy
# credit:      Indicates whether subject recieved class credit or 
#              not
# lab:         Indicates whether subject completed experiment in the 
#              lab or not
# ID:          Original subject ID (non-sequential)
# Age:         Dummy-coded variable for age (1 = Old)
# Old:         Dummy-coded variable for binary response (1 = Old) 
# Cond14:      Numerical index for all 14 conditions, where...
#              1  = Unrelated lures (Phonetic)
#              2  = Unrelated lures (Semantic)
#              3  = Critical lures  (Phonetic, 2 items)
#              4  = Related lures   (Phonetic, 2 items)
#              5  = Targets         (Phonetic, 2 items)
#              6  = Critical lures  (Semantic, 2 items)
#              7  = Related lures   (Semantic, 2 items)
#              8  = Targets         (Semantic, 2 items)
#              9  = Critical lures  (Phonetic, 8 items)
#              10 = Related lures   (Phonetic, 8 items)
#              11 = Targets         (Phonetic, 8 items)
#              12 = Critical lures  (Semantic, 8 items)
#              13 = Related lures   (Semantic, 8 items)
#              14 = Targets         (Semantic, 8 items)
# Cond13:      Numerical index for all 13 conditions, where...
#              1  = Unrelated lures 
#              2  = Critical lures  (Phonetic, 2 items)
#              3  = Related lures   (Phonetic, 2 items)
#              4  = Targets         (Phonetic, 2 items)
#              5  = Critical lures  (Semantic, 2 items)
#              6  = Related lures   (Semantic, 2 items)
#              7  = Targets         (Semantic, 2 items)
#              8  = Critical lures  (Phonetic, 8 items)
#              9  = Related lures   (Phonetic, 8 items)
#              10 = Targets         (Phonetic, 8 items)
#              11 = Critical lures  (Semantic, 8 items)
#              12 = Related lures   (Semantic, 8 items)
#              13 = Targets         (Semantic, 8 items)

###
### Extract hits/false alarms for subjects
###
# Lookup - 05

ratesDat = aggregate( allDat$Old,
                      list( allDat$Cond14, allDat$subject ),
                      length )
colnames( ratesDat ) = c( 'Cond14', 'S', 'N' )
ratesDat$Y = aggregate( allDat$Old, 
                      list( allDat$Cond14, allDat$subject ), 
                      sum )$x
ratesDat$P = aggregate( allDat$Old, 
                        list( allDat$Cond14, allDat$subject ), 
                        mean )$x
ratesDat$A = aggregate( allDat$Age, 
                        list( allDat$Cond14, allDat$subject ), 
                        unique )$x
ratesDat$Cond13 = aggregate( allDat$Cond13, 
                        list( allDat$Cond14, allDat$subject ), 
                        unique )$x
ratesDat$label14 = aggregate( allDat$labels14, 
                             list( allDat$Cond14, allDat$subject ), 
                             unique )$x
ratesDat$label13 = aggregate( allDat$labels13, 
                              list( allDat$Cond14, allDat$subject ), 
                              unique )$x

###
### Save results as an .RData file
###
# Lookup - 06

save( allDat, ratesDat, N, file = 'DRMdata.RData' )
setwd( orig_dir )



