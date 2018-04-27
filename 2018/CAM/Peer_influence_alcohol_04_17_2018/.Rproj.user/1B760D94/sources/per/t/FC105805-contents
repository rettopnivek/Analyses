# Conversion from .csv to .RData
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directory if you 
# have any questions or comments
# Last updated 2018-04-25

# Table of contents
# 1) Initial setup
# 2) Randi's data
#   2.1) Demographic measures
#   2.2) Personality measure
#   2.3) Impulsivity measure
#   2.4) Drinking motive measures
#   2.5) Suggestibility measures
#   2.6) Alcohol expectancy measures
#   2.7) Complete data set
# 3) Compute summary scores
#   3.1) MISS
#   3.2) UPPS
#   3.3) DMQ
#   3.4) B-CEOA
#   3.5) Drinking behavior
# 4) Save results

###
### 1) Initial setup
###

# Save current working directory
orig_dir = getwd()

# Navigate to project directory
setwd( '..' ); proj_dir = getwd();

# Load in useful packages

# Collection of useful functions
# devtools::install_github("rettopnivek/utilityf")
library( utilityf )

# Packages for creating tidy data
# install.packages( 'tidyr' )
library(tidyr)
# install.packages( 'dplyr' )
library(dplyr)

###
### 2) Randi's data
###

# Navigate to data directory
setwd( 'Data' ); dat_dir = getwd();
setwd( 'Original_files' )

# Load in data set from Randi's study
Data_Randi = read.csv( file = "Revised Dataset for Kevin.csv",
                       header = TRUE, stringsAsFactors = FALSE )
Data_Randi$Subject = createIncrement( Data_Randi$id )
# Extract variable names
cn = colnames( Data_Randi )

# 2.1) 
### Demographic measures ###
Demographics = data.frame(
  ID = Data_Randi$id,
  Subject = Data_Randi$Subject,
  Acrostic = Data_Randi$acrostic,
  Age = Data_Randi$age_exact,
  Sex = 'Male',
  Total_days_alcohol_consumed_last_90_days = 
    Data_Randi$tlfb_etoh_12a_recalc,
  Total_number_of_drinks_last_90_days = 
    Data_Randi$tlfb_etoh_13a_recalc,
  stringsAsFactors = FALSE
)
Demographics$Sex[ Data_Randi$sex_conf == 2 ] = 'Female'

# 2.2)
### Personality measure ###

# Extract scores
tipi_sel = grep( 'tipi', cn )
scores = Data_Randi[ , tipi_sel ]

# Ten Item Personality Inventory (TIPI)
TIPI = list(
  # The personality trait to which subjects 
  # indicated their degree of match
  personality_traits_per_item = c(
    Q1  = "Extraverted, enthusiastic",
    Q2  = "Critical, quarrelsome (Recoded)",
    Q3  = "Dependable, self-disciplined",
    Q4  = "Anxious, easily upset (Recoded)",
    Q5  = "Open to new experiences, complex",
    Q6  = "Reserved, quiet (Recoded)",
    Q7  = "Sympathetic, warm",
    Q8  = "Disorganized, careless (Recoded)",
    Q9  = "Calm, emotionally stable",
    Q10 = "Conventional, uncreative (Recoded)" ),
  # The items whose likert scores were reversed
  flipped_items = c( 2, 4, 6, 8, 10 ), 
  # Item assignments for summed scores
  subscale_indices = list(
    extraversion = c( 1, 6 ),
    agreeableness = c( 2, 7 ),
    conscientiousness = c( 3, 8 ),
    emotional_stability = c( 4, 9 ),
    openness = c( 5, 10 )
  ), 
  # The response categories (likert scale)
  response_categories = c(
    'Disagree strongly',
    'Disagree moderately',
    'Disagree a little',
    'Neither disagree or agree',
    'Agree a little',
    'Agree moderately',
    'Agree strongly'
  ),
  # The observed data
  scores = NULL
)
TIPI$scores = scores[,1:10]
colnames( TIPI$scores ) = paste( 'TIPI_Q', 1:10, sep = '' )

# 2.3)
### Impulsivity measure ###

# Fix mislabeled item
upps_sel = grep( 'upps', cn )
sel = grep( 'upps_12_', cn )
cn[sel] = 'upps_12_rec'
colnames( Data_Randi ) = cn

# Extract scores
scores = Data_Randi[ , upps_sel ]

# Urgency, Premediation, Perseverance, Sensation seeking,
# and Positive urgency (UPPS-P) impulsivity inventory
UPPS = list(
  # The items whose likert scores were reversed
  flipped_items = c(
    2, 3, 5, 7, 8, 9, 10, 12, 13, 15, 17, 18, 
    20, 22, 23, 25, 26, 29, 30, 31, 34, 35, 36, 
    39, 40, 41, 44, 45, 46, 47, 49, 50, 51, 52, 
    54, 56, 57, 58, 59 ),
  # Item assignments for average scores
  subscale_indices = list(
    negative_urgency = c(
      2, 7, 12, 17, 22, 29, 34, 39, 44, 50, 53, 58 ),
    premediation = c(
      1, 6, 11, 16, 21, 28, 33, 38, 43, 48, 55 ),
    perseverance = c(
      4, 9, 14, 19, 24, 27, 32, 37, 42, 47 ),
    sensation_seeking = c(
      3, 8, 13, 18, 23, 26, 31, 36, 41, 46, 51, 56 ),
    positive_urgency = c(
      5, 10, 15, 20, 25, 30, 35, 40, 45, 49, 52, 54, 57, 59 )
  ),
  # The response categories (likert scale)
  response_categories = c(
    'Agree strongly',
    'Agree somewhat',
    'Disagree somewhat',
    'Disagree strongly'
  ),
  scores = NULL
)
item_labels = paste(
  'upps', 1:59, sep = '_' )
item_labels[ UPPS$flipped_items ] = paste(
  item_labels[ UPPS$flipped_items ],
  'rec', sep = '_' )

UPPS$scores = scores[,item_labels]
colnames( UPPS$scores ) = paste( 'UPPS_Q', 1:59, sep = '' )

# 2.4)
### Drinking motive measures ###

# Extract scores
dmq_sel = grep( 'dmq', cn )
scores = Data_Randi[ , dmq_sel ]

# Drinking Motives Questionnaire (DMQ)
DMQ = list(
  # Item assignments for summed scores
  subscale_indices = list(
    social_motives = c(3,5,11,14,16),
    coping = c(1,4,6,15,17),
    enhancement = c(7,9,10,13,18),
    conformity = c(2,8,12,19,20)
  ), 
  # The response categories (likert scale)
  response_categories = c(
    'Almost never / Never',
    'Some of the time',
    'Half of the time',
    'Most of the time',
    'Almost always / Always'
  ),
  # The observed data
  scores = NULL
)

DMQ$scores = scores[,1:20]
colnames( DMQ$scores ) = paste( 'DMQ_Q', 1:20, sep = '' )

# 2.5)
### Suggestibility measures ###

# Extract scores
miss_sel = grep( 'miss', cn )
scores = Data_Randi[ , miss_sel ]

# Multidimensional Iowa Suggestibility Test (MISS)
MISS = list(
  # The items whose likert scores negatively impact 
  # their construct
  flipped_items = c(
    34, 39 ),
  # Item assignments for summed scores
  subscale_indices = list(
    consumer_suggestibility = 
      c(2,10,14,20,24,32,45,51,57,63,70),
    persuadibility = 
      c(35,5,89,69,47,62,75,44,88,13,1,22,82,76), 
    physiological_suggestibility = 
      c(64,94,71,25,52,15,58,66,33,68,11,77),
    physiological_reactivity = 
      c(12,40,85,91,31,43,73,27,50,17,21,3,60),
    peer_conformity = 
      c(4,16,29,34,39,46,53,59,65,72,78,84,95,90),
    mental_control = 
      c(8,18,23,28,36,48,55,67,80,92,83,79,74,41,6),
    unpersuadibility = 
      c(7,9,19,37,38,42,49,87,81,86,93,26,56,30,54,61)
  ), 
  # The response categories (likert scale)
  response_categories = c(
    'Not at all / Slightly',
    'A little',
    'Somewhat',
    'Quite a bit',
    'A lot'
  ),
  # The observed data
  scores = NULL
)

MISS$scores = scores[,1:95]
colnames( MISS$scores ) = paste( 'MISS_Q', 1:95, sep = '' )

# 2.6)
### Alcohol expectancy measures ###

# Extract scores
bceoa_sel = grep( 'bceoa', cn )
scores = Data_Randi[ , bceoa_sel ]


# Brief Comprehensive Effects of Alcohol (B-CEOA)
BCEOA = list(
  # Item assignments for summed scores
  subscale_indices = list(
    liquid_courage = 
      c( 2, 5, 6, 8, 9, 13, 14 ),
    self_perception = 
      c( 10, 11, 12, 15 ),
    sexuality = 
      c( 1, 4 ),
    tension_reduction = 
      c( 2, 7 )
  ), 
  # The response categories (likert scale)
  response_categories = c(
    'Disagree',
    'Slightly disagree',
    'Slightly agree',
    'Agree'
  ),
  # The observed data
  scores = NULL
)

# 2.7) Complete data set
BCEOA$scores = scores[,1:15]
colnames( BCEOA$scores ) = paste( 'BCEOA_Q', 1:15, sep = '' )

# Complete data set
DR = list(
  Scores = 
    cbind(
      ID = Demographics$ID,
      Subject = Demographics$Subject,
      Age = Demographics$Age,
      Sex = Demographics$Sex,
      Total_days = 
        Demographics$Total_days_alcohol_consumed_last_90_days,
      Total_drinks = 
        Demographics$Total_number_of_drinks_last_90_days,
      TIPI$scores,
      UPPS$scores,
      DMQ$scores,
      MISS$scores,
      BCEOA$scores
    ),
  Demographics = Demographics,
  TIPI = TIPI,
  UPPS = UPPS,
  DMQ = DMQ,
  MISS = MISS,
  BCEOA = BCEOA
)


###
### 3) Compute summary scores
###

# 3.1) MISS

# Intialize data frame
SC = data.frame(
  ID = DR$Demographics$ID,
  Subject = DR$Demographics$Subject,
  Age = DR$Demographics$Age,
  Sex = DR$Demographics$Sex, 
  # Consumer suggestability
  MISS_CS = NA,
  # Persuadibility
  MISS_P = NA,
  # Physiological suggestibility
  MISS_PS = NA,
  # Physiological reactivity
  MISS_PR = NA,
  # Peer conformity
  MISS_PC = NA
)
N_trials = c(
  MISS_CS = NA,
  MISS_P = NA,
  MISS_PS = NA,
  MISS_PR = NA,
  MISS_PC = NA
)

tmp = DR$MISS$scores
# Flip certain scores for the MISS inventory
tmp[,DR$MISS$flipped_items] = 
  5 - tmp[,DR$MISS$flipped_items]
sel = colnames( SC )[ grep( 'MISS', colnames( SC ) ) ]
for ( i in 1:5 ) {
  SC[,sel[i]] = rowSums( tmp[,DR$MISS$subscale_indices[[i]]] )
  N_trials[sel[i]] = length( DR$MISS$response_categories ) * 
    length( DR$MISS$subscale_indices[[i]] )
}

# 3.2) UPPS

# Negative urgency
SC$UPPS_NU = NA
# Premediation
SC$UPPS_Pr = NA
# Perseverance
SC$UPPS_Pe = NA
# Sensation seeking
SC$UPPS_SS = NA
# Positive urgency
SC$UPPS_PU = NA

N_trials = c(
  N_trials, 
  UPPS_NU = NA,
  UPPS_Pr = NA,
  UPPS_Pe = NA,
  UPPS_SS = NA,
  UPPS_PU = NA
)
tmp = DR$UPPS$scores
sel = colnames( SC )[ grep( 'UPPS', colnames( SC ) ) ]
for ( i in 1:5 ) {
  SC[,sel[i]] = rowSums( tmp[,DR$UPPS$subscale_indices[[i]]] )
  N_trials[sel[i]] = length( DR$UPPS$response_categories ) * 
    length( DR$UPPS$subscale_indices[[i]] )
}

# 3.3) DMQ

# Social motives
SC$DMQ_SM = NA
# Coping
SC$DMQ_Cp = NA
# Enhancement
SC$DMQ_E = NA
# Conformity
SC$DMQ_Cn = NA

N_trials = c(
  N_trials,
  DMQ_SM = NA,
  DMQ_Cp = NA,
  DMQ_E = NA,
  DMQ_Cn = NA
)

tmp = DR$DMQ$scores
sel = colnames( SC )[ grep( 'DMQ', colnames( SC ) ) ]
for ( i in 1:4 ) {
  SC[,sel[i]] = rowSums( tmp[,DR$DMQ$subscale_indices[[i]]] )
  N_trials[sel[i]] = length( DR$DMQ$response_categories ) * 
    length( DR$DMQ$subscale_indices[[i]] )
}

# 3.4) B-CEOA

# Liquid courage
SC$BCEOA_LC = NA
# Self perception
SC$BCEOA_SP = NA
# Sexuality
SC$BCEOA_S = NA
# Tension reduction
SC$BCEOA_TR = NA

N_trials = c(
  N_trials,
  BCEOA_LC = NA,
  BCEOA_SP = NA,
  BCEOA_S = NA,
  BCEOA_TR = NA
)

tmp = DR$BCEOA$scores
sel = colnames( SC )[ grep( 'BCEOA', colnames( SC ) ) ]
for ( i in 1:4 ) {
  SC[,sel[i]] = rowSums( tmp[,DR$BCEOA$subscale_indices[[i]]] )
  N_trials[sel[i]] = length( DR$BCEOA$response_categories ) * 
    length( DR$BCEOA$subscale_indices[[i]] )
}

# 3.5) Drinking behavior

SC$Total_days = 
  DR$Demographics$Total_days_alcohol_consumed_last_90_days
SC$Total_drinks = 
  DR$Demographics$Total_number_of_drinks_last_90_days

###
### 4) Save results
###

setwd( proj_dir )
setwd( dat_dir )
save( DR, SC, file = 'Peer_influence_alcohol.RData' )
setwd( orig_dir )

