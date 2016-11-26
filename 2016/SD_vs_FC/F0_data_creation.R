#-----------------------------#
# Script to create R data set #
# Kevin Potter                #
# Updated 11/23/16            #
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

# Index
# Lookup - 01:  Read in .csv file
# Lookup - 02:  Identify subjects to exclude based on Dave's cut-offs
# Lookup - 03:  Add in additional variables

###
### Read in .csv file
###
# Lookup - 01

# Change directory to where data is stored
setwd( 'Data' )

# Load in .csv file with data
rawDat = read.table( file = 'total_all.csv', sep = ',', header = F )

# Columns for rawDat
# 1: Subject index
# 2: Condition index, where...
#    (T = target word, F = foil word, one letter is 
#    a short duration prime (50 ms), two letters are 
#    a long duration prime (400 ms), and the final 
#    column gives the prime, target flash, and the 
#    on-screen alternatives:
#    2AFC task  (Two-alternative forced-choice)
#    0          LEFT  SHORT TAR   T  T T/F
#    1          RIGHT SHORT TAR   T  T F/T
#    2          LEFT  LONG  TAR   TT T T/F
#    3          RIGHT LONG  TAR   TT T F/T
#    4          LEFT  SHORT FOIL  F  T T/F
#    5          RIGHT SHORT FOIL  F  T F/T
#    6          LEFT  LONG  FOIL  FF T T/F
#    7          RIGHT LONG  FOIL  FF T F/T
#    SD task    (Same-different)
#    8          SAME  SHORT TAR   T  T T
#    9          DIFF  SHORT TAR   T  T F
#    10         SAME  LONG  TAR   TT T T
#    11         DIFF  LONG  TAR   TT T F
#    12         SAME  SHORT FOIL  F  T T
#    13         DIFF  SHORT FOIL  F  T F
#    14         SAME  LONG  FOIL  FF T T
#    15         DIFF  LONG  FOIL  FF T F
# 3: Accuracy ( 0 = incorrect, 1 = correct )
# 4: Response time (ms)

# Create meaningful column names
colnames(rawDat) = c('Subject','Condition','Accuracy','RT')

###
### Identify subjects to exclude based on criterions set by Dave
###
# Lookup - 02

# Define variable to identify which subjects have been excluded
# excl == 00000:  Subject was kept
# excl == 10000:  Subject was excluded due to computer error
# excl == 01000:  Subject was excluded due to missing second session
# excl == 00100:  Subject was excluded due to accuracy being too high 
#                 ( > 0.8 )
# excl == 00010:  Subject was excluded due to accuracy being too low 
#                 ( < 0.6 )
# excl == 00001:  Subject was excluded due to RT being too slow 
#                 ( mean RT > 1000 )

rawDat$excl = '00000'

rawDat$Subject[ rawDat$Subject==99 ] = 42 # Re-label last subject to be 42

# Identify subjects in which computer crashed and failed to record 
# RT properly
sel = which(rawDat$RT > 1000000)
excl = unique( rawDat$Subject[sel] )
for (n in 1:nrow(rawDat)) {
  for (ex in 1:length(excl)) {
    if ( rawDat$Subject[n] == excl[ex] ) {
      tmp = strsplit( rawDat$excl[n], split='' )[[1]]
      tmp[1] = '1';
      rawDat$excl[n] = paste( tmp, collapse='' )
    }
  }
}

# Identify subjects who did not return for a second session
cnt = rep(1,nrow(rawDat))
sel = aggregate(cnt,list(rawDat$Subject),sum)
excl = sel$Group.1[ sel$x < 1280 ]
for (n in 1:nrow(rawDat)) {
  for (ex in 1:length(excl)) {
    if ( rawDat$Subject[n] == excl[ex] ) {
      tmp = strsplit( rawDat$excl[n], split='' )[[1]]
      tmp[2] = '1';
      rawDat$excl[n] = paste( tmp, collapse='' )
    }
  }
}

# Remove the cases in which the computer failed to record RT
rawDat = rawDat[rawDat$RT < 1000000,]

# Identify subjects who had accuracy that was too high or too low
sel = aggregate(rawDat$Accuracy,list(rawDat$Subject),mean)
excl = sel$Group.1[ sel$x > 0.8 ]
for (n in 1:nrow(rawDat)) {
  for (ex in 1:length(excl)) {
    if ( rawDat$Subject[n] == excl[ex] ) {
      tmp = strsplit( rawDat$excl[n], split='' )[[1]]
      tmp[3] = '1';
      rawDat$excl[n] = paste( tmp, collapse='' )
    }
  }
}

excl = sel$Group.1[ sel$x < 0.6 ]
for (n in 1:nrow(rawDat)) {
  for (ex in 1:length(excl)) {
    if ( rawDat$Subject[n] == excl[ex] ) {
      tmp = strsplit( rawDat$excl[n], split='' )[[1]]
      tmp[4] = '1';
      rawDat$excl[n] = paste( tmp, collapse='' )
    }
  }
}

# Identify subjects who had an average RT that was too high
sel = aggregate(rawDat$RT,list(rawDat$Subject),mean)
excl = sel$Group.1[ sel$x > 1000 ]
for (n in 1:nrow(rawDat)) {
  for (ex in 1:length(excl)) {
    if ( rawDat$Subject[n] == excl[ex] ) {
      tmp = strsplit( rawDat$excl[n], split='' )[[1]]
      tmp[5] = '1';
      rawDat$excl[n] = paste( tmp, collapse='' )
    }
  }
}

# Clean up workspace
rm(cnt,sel,excl,ex,n,tmp)

###
### Add in additional variables
###
# Lookup - 03

# Define a dummy-coded variable indicating prime duration
# where 0 = short and 1 = long
rawDat$PrimeDur = 1;
ind = rawDat$Condition==0 | rawDat$Condition==1 | 
  rawDat$Condition==4 | rawDat$Condition==5 | 
  rawDat$Condition==8 | rawDat$Condition==9 |
  rawDat$Condition==12 | rawDat$Condition==13
rawDat$PrimeDur[ind] = 0

# Define a dummy-coded variable indicating which type
# of task was used where 0 = 2AFC and 1 = SD
rawDat$Task = 0
ind = rawDat$Condition==8 | rawDat$Condition==9 | 
  rawDat$Condition==10 | rawDat$Condition==11 | 
  rawDat$Condition==12 | rawDat$Condition==13 |
  rawDat$Condition==14 | rawDat$Condition==15
rawDat$Task[ind] = 1

# Define a dummy-coded variable indicating whether
# the correct choice was on the left or right for 
# the 2AFC task, and whether the answer was 'same' (0) 
# or 'different' (1) for the SD task
rawDat$Correct = 0
ind = rawDat$Condition==1 | rawDat$Condition==3 | 
  rawDat$Condition==5 | rawDat$Condition==7 | 
  rawDat$Condition==9 | rawDat$Condition==11 | 
  rawDat$Condition==13 | rawDat$Condition==15
rawDat$Correct[ind] = 1

# Define a dummy-coded variable indicating whether
# the target (1) or foil (0) was primed
rawDat$PrimeType = 0
ind = rawDat$Condition==0 | rawDat$Condition==1 | 
  rawDat$Condition==2 | rawDat$Condition==3 | 
  rawDat$Condition==8 | rawDat$Condition==9 | 
  rawDat$Condition==10 | rawDat$Condition==11
rawDat$PrimeType[ind] = 1

# Create a variable indicate whether a subject chose left/right or
# same/different (0 = left/same, 1 = right/different)
rawDat$Choice=rawDat$Correct
rawDat$Choice[rawDat$Accuracy==0] = 1-rawDat$Correct[rawDat$Accuracy==0]

# Shift the conditions so they start at 1
rawDat$Condition = rawDat$Condition + 1

# Rescale response times to be in seconds
rawDat$RT = rawDat$RT/1000

# Adjust subject indices
sel = rawDat$excl == '00000'
tmp = rawDat$Subject[ sel ]
tmp = createIncrement(tmp)
tmp2 = rawDat$Subject[ !sel ]
tmp2 = createIncrement(tmp2)
tmp2 = tmp2 + length( unique( tmp ) )
rawDat$Subject[ sel ] = tmp
rawDat$Subject[ !sel ] = tmp2

# Determine number of subjects
N = length( unique( rawDat$Subject ) )

# Clean up workspace
rm(ind, sel, tmp, tmp2)

"
Data key

42 subjects in rawDat
25 subjects in dat1
80 trials per condition
Trials with computer errors (RT == 86400000) were removed after noting 
  subject number
Prime durations were 50 ms (short) or 400 ms (long)
No counter-balancing, subjects pressed f = left/different or 
  j = right/same

Subject   - Indicates the subject to which a current trial belongs
Condition - Indicates which condition a trial belongs to, where...
            Condition
            2AFC task
            1          LEFT  SHORT TAR   T   T/F
            2          RIGHT SHORT TAR   T   F/T
            3          LEFT  LONG  TAR   TT  T/F
            4          RIGHT LONG  TAR   TT  F/T
            5          LEFT  SHORT FOIL  F   T/F
            6          RIGHT SHORT FOIL  F   F/T
            7          LEFT  LONG  FOIL  FF  T/F
            8          RIGHT LONG  FOIL  FF  F/T
            
            Same/different task
            9          SAME  SHORT TAR   T   T
            10         DIFF  SHORT TAR   T   F
            11         SAME  LONG  TAR   TT  T
            12         DIFF  LONG  TAR   TT  F
            13         SAME  SHORT FOIL  F   T
            14         DIFF  SHORT FOIL  F   F
            15         SAME  LONG  FOIL  FF  T
            16         DIFF  LONG  FOIL  FF  F
Accuracy  -   Indicates whether a subject made the correct choice
              (0 = inaccurate, 1 = accurate)
RT        -   How fast a subject responded (s)
excl      -   Indicates whether a subject was included (00000), or 
              excluded based on one of 5 reasons...
                10000 - Computer error in which RT was not measured 
                        correctly
                01000 - Subjects did not complete 2nd session
                00100 - Subjects were too accuracte (average > 0.8)
                00010 - Subjects were too inaccuracte (average < 0.6)
                00001 - Subjects were too slow (average RT > 1000)
              Note that some subjects had a combination of the above 
              (e.g. 10010 means that a subject had both the computer 
              error and was too inaccurate)
PrimeDur  -   A dummy-coded variable indicating whether prime duration
              was short (0) or long (1)
Task      -   A dummy-coded variable indicating whether the task at
              hand was 2AFC (0) or same-different (1)
Correct   -   A dummy-coded variable indicating whether the correct
              choice was shown on the left (0) or right (1) or was
              same (0) or different (1)
PrimeType -   A dummy-coded variable indicating whether the target (1)
              or foil (0) was primed
              1) Prime shown 2) Target/Foil shown 3) Choices to pick shown
              For primes: A = short duration; AA = long duration
              For target/foils: A = target primed; B = foil primed
              For choices: A/B = 2AFC (target on left);
              A = Same/different (target shown)
Choice    -   A variable indicating whether a subject chooses 
              left/right or same/different (0 = left/same, 1 = 
              right/different)
"

# Save data sets
save(rawDat,N,file='SD_v_FC.RData')