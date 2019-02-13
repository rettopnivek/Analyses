# Convert .csv to .RData
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you 
# have any questions or comments
# Last updated 2018-12-05

# Table of contents
# 1) Initial setup
# 2) Read in raw data
# 3) Convert to data frame for analyses
#   3.1) quick_proprogate
#   3.2) Extract variables of interest
#     3.2.1) Subject details
#     3.2.2) Measures of drug use
#     3.2.3) Primary variables of interest
#     3.2.4) Data issues
#     3.2.5) Create a subject index
# 4) Extract CUDIT and CWS scores
# 5) Save results

###
### 1) Initial setup
###

source( 'S01_Folder_paths.R' )

# Load in useful packages

# Collection of useful functions
my_package_load( 'utilityf', github = T )

# Packages for manipulating data frames
my_package_load( 'dplyr' )

###
### 2) Read in raw data
###

setwd( dat_dir )
setwd( 'Original_files' )

# Read in .csv file
all_files = dir()
fname = all_files[ grepl( 'Combined_data_', all_files ) ]
raw_dat = read.csv( 
  file = fname,
  header = T,
  stringsAsFactors = FALSE )
# Clean up workspace
rm( all_files, fname )

###
### 3) Convert to data frame for analyses
###

# 3.1)
quick_proprogate = function( vrb ) {
  # Purpose:
  # Function to proprogate singular values for a 
  # given subject over repeated measures.
  # Arguments:
  # vrb - Character string for variable of interest 
  #       in 'raw_dat'
  # Returns:
  # A vector with the repeated values for each subject.
  
  # Isolate subjects
  subj = unique( all_dat$ID )
  
  # Initialize output
  sel = all_dat$ID == subj[1]
  x = raw_dat[sel,vrb]
  no_na = !is.na( x )
  out = rep( x[no_na], nrow( all_dat ) )
  
  # Loop over remaining subjects
  for ( s in 2:length( subj ) ) {
    
    sel = all_dat$ID == subj[s]
    x = raw_dat[sel,vrb]
    no_na = !is.na( x )
    if ( any( no_na ) ) {
      out[sel] = x[no_na]
    } else {
      out[sel] = NA
    }
    
  }
  
  return( out )
}

# 3.2) Extract variables of interest

# Initialize data frame
all_dat = data.frame(
  ID = raw_dat$id,
  stringsAsFactors = F
)

# 3.2.1) Subject details

# Sex
all_dat$Sex = quick_proprogate( 'sex_conf' )
# Convert to effects coding
all_dat$Sex[ all_dat$Sex == 1 ] = -1
all_dat$Sex[ all_dat$Sex == 2 ] = 1
attributes( all_dat$Sex ) = list(
  Coding = data.frame(
    Sex = c( 'Male', 'Female' ),
    Value = c(-1,1),
    stringsAsFactors = F ) )

# Age
all_dat$Age = quick_proprogate( 'age_exact' )
attributes( all_dat$Age ) = list(
  Units = 'Years'
)

# Body mass index
all_dat$BMI = quick_proprogate( 'bmi' )
attributes( all_dat$BMI ) = list(
  Units = 'Kilograms/Meters^2'
)

# Race (White versus non-white)
all_dat$Race = quick_proprogate( 'race_revised' )
# Convert to effects coding
all_dat$Race[ all_dat$Race == 1 ] = -1
all_dat$Race[ all_dat$Race == 2 ] = 1
attributes( all_dat$Race ) = list(
  Coding = data.frame(
    Race = c( 'White', 'Non-white' ),
    Value = c(-1,1),
    stringsAsFactors = F ) )

# Height (in inches)
all_dat$Height = quick_proprogate( 'height' )
attributes( all_dat$Height ) = list(
  Units = 'Inches' )

# Weight (in pounds)
all_dat$Weight = quick_proprogate( 'weight' )
attributes( all_dat$Weight ) = list(
  Units = 'Pounds' )

# Race (Original coding scheme)
tmp_1 = quick_proprogate( 'race' )
tmp_2 = quick_proprogate( 'ethnicity' )
all_dat$Race_orig = 'W'
sel = tmp_1 == 2
all_dat$Race_orig[sel] = 'AA'
sel = tmp_1 == 3
all_dat$Race_orig[sel] = 'AS'
sel = tmp_1 == 4
all_dat$Race_orig[sel] = 'PI'
sel = tmp_1 == 5
all_dat$Race_orig[sel] = 'NA'
sel = tmp_1 == 6
all_dat$Race_orig[sel] = 'MOR'
sel = tmp_1 == 7
all_dat$Race_orig[sel] = 'OTH'

# Add ethnicity status
for ( i in 1:length( tmp_2 ) ) {
  
  if ( !is.na( tmp_2[i] ) & 
       tmp_2[i] == 2 ) {
    all_dat$Race_orig[i] = paste(
      all_dat$Race_orig[i],
      '(H)' )
  }
  
}
attributes( all_dat$Race_orig ) = list(
  Coding = data.frame(
    Race_orig = c( 
      'White',
      'Haitian, Black or African American',
      'Asian',
      'Hawaiian or other Pacific Islander',
      'American Indian/Alaska Native',
      'More than one race',
      'Other',
      'Hispanic' ),
    Value = c( 
      'W',
      'AA',
      'AS',
      'PI',
      'NA',
      'MOR',
      'OTH',
      'H' ),
    stringsAsFactors = F ) )
# Clean up workspace
rm( tmp_1, tmp_2, i )

# 3.2.2) Measures of drug use

# Years spent using cannabis
all_dat$Years_of_MJ_use = quick_proprogate( 'years_mj_use' )
attributes( all_dat$Years_of_MJ_use ) = list(
  Units = 'Number of years spent using MJ'
)

# Recency of cannabis use (0 - 2 days)
# Make sure subjects only have one initial value for 
# self-report on recency of use
# Loop over subjects
subj = unique( raw_dat$id )
for ( s in 1:length( subj ) ) {
  sel = raw_dat$id == subj[s]
  x = raw_dat$tlfb_mj_8[sel]
  x = unique( x[ !is.na( x ) ] )
  raw_dat$tlfb_mj_8[sel] = NA
  raw_dat$tlfb_mj_8[sel][1] = x
}
# Clean up workspace
rm( subj, s, sel, x )
all_dat$Recency_of_MJ_use = quick_proprogate( 'tlfb_mj_8' )
attributes( all_dat$Recency_of_MJ_use ) = list(
  Units = 'Number of days since last use of MJ'
)

# Level of cannabis use
all_dat$Level_of_MJ_use = quick_proprogate( 'tlfb_mj_14b' )
attributes( all_dat$Level_of_MJ_use ) = list(
  Units = 'Number of times MJ was used in last 30 days'
)

# Level of alcohol use
all_dat$Level_of_alcohol_use = quick_proprogate( 'tlfb_etoh_13b' )
attributes( all_dat$Level_of_alcohol_use ) = list(
  Units = 'Number of times alcohol was used in last 30 days'
)

# Nicotine use
all_dat$Nicotine_use = quick_proprogate( 'nic_user' )
# Convert to effects coding
all_dat$Nicotine_use[ all_dat$Nicotine_use == 0 ] = -1
all_dat$Nicotine_use[ all_dat$Nicotine_use == 1 ] = 1
attributes( all_dat$Nicotine_use ) = list(
  Coding = data.frame(
    Nicotine_use = c( 'No nicotine', 'Nicotine user' ),
    Value = c(-1,1),
    stringsAsFactors = F ) )

# Dipstick test
all_dat$Dipstick = raw_dat$thc_qual
attributes( all_dat$Dipstick ) = list(
  Coding = data.frame(
    Dipstick = c( 'Positive (>50)', 'Negative(<50)' ),
    Value = c(1,0),
    stringsAsFactors = F ) )

# Positive test based on federal guidelines
all_dat$Positive_test_fed = 0
attributes( all_dat$Positive_test_fed ) = list(
  Coding = data.frame(
    Dipstick = c( 'Positive (Dipstick and THCCOOH>15)', 
                  'Negative (No dipstick or THCCOOH<15)' ),
    Value = c(1,0),
    stringsAsFactors = F ) )
# Define results based on federal guidelines
# Positive tests need positive dipstick and 
# unadjusted THC above 15
sel = raw_dat$thc_qual == 1 & 
  as.numeric( raw_dat$ua_conf_result ) > 15
sel = sel & !is.na( sel )
all_dat$Positive_test_fed[sel] = 1
all_dat$Positive_test_fed[ is.na( raw_dat$thc_qual ) ] = NA

# Age when subject first used MJ
all_dat$Age_first_used_MJ = quick_proprogate( 'tlfb_mj_2' )
attributes( all_dat$Age_first_used_MJ ) = list(
  Units = 'Years'
)

# Urine pH
all_dat$Urine_pH = raw_dat$ua_thc_ph
# Urine specific gravity
all_dat$Urine_specific_gravity = raw_dat$ua_thc_grav

# THCCOOH no creatinine adjustment
all_dat$THCCOOH_no_CN = 
  as.numeric( raw_dat$ua_conf_result )

# Days spent using in past 30
all_dat$Days_used_past_30 = 
  quick_proprogate( 'tlfb_mj_14b' )
# Times used in past 30
all_dat$Times_used_past_30 = 
  quick_proprogate( 'tlfb_mj_15b' )
# Amount used in past 30
all_dat$Amount_used_past_30 = 
  quick_proprogate( 'tlfb_mj_16b' )
# Times per day smoked
all_dat$Times_smoked_per_day = 
  quick_proprogate( 'tlfb_mj_6' )

# 3.3.3) Primary variables of interest

# THCCOOH level (creatinine adjusted)
all_dat$THCCOOH = raw_dat$ua_thclcmsms_crtadj
attributes( all_dat$THCCOOH ) = list(
  Units = 'Nanograms per milliliter (ng/mL)'
)

# Time since baseline measurement

all_dat$Days_from_baseline = 
  raw_dat$days_since_baseline
attributes( all_dat$Days_from_baseline ) = list(
  Units = 'Days since first baseline measurement'
)

all_dat$Time = all_dat$Days_from_baseline + 
  all_dat$Recency_of_MJ_use
attributes( all_dat$Time ) = list(
  Units = 'Days since last estimated use of MJ'
)

# Visit number
all_dat$Visit_number = raw_dat$visit_number

# Withdrawal symptoms
all_dat$Withdrawal_intensity = raw_dat$cws_intensity
all_dat$Withdrawal_negative_impact = raw_dat$cws_negimpact

# 3.3.4) Data issues

# Define function that determines different types of 
# data issues for each row of the data set
check_for_data_issues = function( r ) {
  
  # Function to concatenate issues
  f = function(x,y) {
    if ( x == '0' ) {
      out = y
    } else {
      out = paste( x, y, sep = ';' )
    }
    return( out )
  }
  
  # Initialize output
  out = '0'
  
  # Flag observations with THC values over the readout ceiling 
  # of the instrument
  if ( !is.na( raw_dat$thc_500_flag[r] ) ) {
    if ( raw_dat$thc_500_flag[r] ) out = f(out,'1')
  }
  
  # Flag subjects who did not recently use THC
  if ( !is.na( raw_dat$thc_recency_flag[r] ) ) {
    if ( raw_dat$thc_recency_flag[r] ) out = f(out,'2')
  }
  
  # Check for missing data over variables
  
  # Dependent variables
  if ( is.na( all_dat$THCCOOH[r] ) ) out = f(out,'3')
  # if ( is.na( all_dat$Withdrawal_intensity[r] ) ) out = f(out,'4')
  # if ( is.na( all_dat$Withdrawal_negative_impact[r] ) ) out = f(out,'5')
  
  # Independent variables
  if ( is.na( all_dat$Time[r] ) ) out = f(out,'6')
  
  # Flag subject who did not appear to maintain abstinence
  if ( raw_dat$id[r] %in% 
       c( '10016',
          '060_COMM',
          '076_COMM',
          '10063' ) ) out = f(out,'7')
  
  # Flag specific time points in which subject did not 
  # maintain abstinence
  if ( raw_dat$id[r] == '003_COMM'
       & raw_dat$visit_number[r] > 5 ) out = f( out, '7' )
  if ( raw_dat$id[r] == '006_COMM'
       & raw_dat$visit_number[r] > 5 ) out = f( out, '7' )
  if ( raw_dat$id[r] == '084_COMM'
       & raw_dat$visit_number[r] > 3 ) out = f( out, '7' )
  if ( raw_dat$id[r] == '086_COMM'
       & raw_dat$visit_number[r] > 4 ) out = f( out, '7' )
  if ( raw_dat$id[r] == '10005'
       & raw_dat$visit_number[r] > 5 ) out = f( out, '7' )
  if ( raw_dat$id[r] == '10037'
       & raw_dat$visit_number[r] > 5 ) out = f( out, '7' )
  if ( raw_dat$id[r] == '20014'
       & raw_dat$visit_number[r] < 2 ) out = f( out, '7' )
  if ( raw_dat$id[r] == '20024'
       & raw_dat$visit_number[r] < 2 ) out = f( out, '7' )
  
  return( out )
}

# Apply function to all rows
all_dat$Data_issues = sapply( 1:nrow( all_dat ), 
                              check_for_data_issues )
attributes( all_dat$Data_issues ) = list(
  Coding = data.frame(
    Issue = c( 
      'None',
      'Could not read value over 500',
      'Did not use THC recently',
      'Missing data for THCCOOH',
      'Missing data for withdrawal intensity',
      'Missing data for withdrawal impact',
      'Missing data for days from baseline',
      'Did not maintain abstinence' ),
    Value = 0:7,
    stringsAsFactors = F ) )

# Track data issues by visit number and subject
di = all_dat %>% 
  group_by( Visit_number, ID ) %>% 
  select( Data_issues ) %>% 
  unique
# Remove cases with no issues
di = di[ di$Data_issues != '0', ]

# Track individual issues
di$Flag.No_reading_for_THC_over_500 = 0
di$Flag.No_recent_THC_use = 0
di$Missing.THCCOOH = 0
di$Missing.Withdrawal_intensity = 0
di$Missing.Withdrawal_impact = 0
di$Missing.Days_from_baseline = 0
di$Missing.Failed_abstinence = 0
for ( i in 1:nrow( di ) ) {
  sel = strsplit( di$Data_issues[i], split = ';' )[[1]]
  sel = as.numeric(sel)
  di[i,sel+3] = 1
}
# Arrange based on issue
di = di %>% arrange( Visit_number, Data_issues )
di$Data_issues = NULL

# Save as a table
setwd( dat_dir )
setwd( 'Original_files' )
write.table( di,
             sep = ',',
             row.names = F,
             quote = F,
             file = 'Data_issues.csv'
)
setwd( dat_dir )

# 3.3.5) Create a subject index

all_dat$Subjects = 0
sel1 = all_dat$Data_issues == '0'
sel2 = !sel1
subj = unique( all_dat$ID )
inc = 1
inc2 = length( unique( all_dat$ID[sel1] ) ) + 1
for ( s in 1:length( subj ) ) {
  
  sel = all_dat$ID == subj[s]
  if ( any( all_dat$ID[sel1] %in% subj ) ) {
    all_dat$Subjects[sel] = inc
    inc = inc + 1
  } else {
    all_dat$Subjects[sel] = inc2
    inc2 = inc2 + 1
  }
}

# Clean up workspace
rm( sel, s, sel1, sel2, subj, inc, inc2, di, i )

###
### 4) Extract CUDIT and CWS scores
###

# Extract the CUDIT scores
sel = paste( 'cudit', 1:8, sep = '_' )
init = matrix( NA, nrow( all_dat ), 8 )
colnames( init ) = paste( 'CUDIT', 1:8, sep = '.' )
all_dat = cbind( all_dat, init )

for ( i in 1:ncol( init ) ) {
  
  all_dat[[ colnames( init )[i] ]] = 
    quick_proprogate( sel[i] )
  
}

# Clean up workspace
rm( sel, init, i )

###
### 5) Save results
###

setwd( dat_dir )
save( raw_dat, all_dat, file = 'THC_decay.RData' )

# Clean up workspace
rm( check_for_data_issues, quick_proprogate )

setwd( R_dir )

