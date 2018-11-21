# Script to combine pilot/active studies
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you 
# have any questions or comments
# Last updated 2018-11-13

# Table of contents
# 1) Initial setup
#   1.1) update_cd
#   1.2) quick_proprogate
# 2) Initial creation of combined data set
#   2.1) Extract variable names
#   2.2) Add pilot data to combined data set
#   2.3) Assess differences between data sets
#   2.4) Add active data to combined data set
# 3) Add in missing variables and adjust values
#   3.1)  Adjust visit number
#   3.2)  Code race (white vs. non-white)
#   3.3)  Create variable for years of MJ use
#   3.4)  Days since baseline
#   3.5)  Days since last used nicotine
#   3.6)  Assorted logical variables
#   3.7)  Add more mis-matched variables
#   3.8)  Adjust units for years (pilot)
#   3.9) Non-critical variables
# 4) Adjustments following trimming of rows
#   4.1) Days since baseline
#   4.2) Adjust dipstick results for pilot group
#   4.3) Imput zero for THC values not detected
#   4.4) Misc. adjustments
#   4.5) Remove values for non-abstinent cases
#   4.6) Adjust values for late-start abstinence
#   4.7) Set censored values as missing
# 5) Save results

###
### 1) Initial setup
###

# Specify folder paths

# Current directory
pa_dir = getwd()
# Directory with data sets
setwd('..'); orig_dir = getwd()
# Overarching directory for data
setwd('..'); dat_dir = getwd()

setwd( pa_dir )

# Load in useful packages
library( utilityf )
library( dplyr )

# Read in data
all_files = dir()

# Remove any previous instances of combined data
fname = all_files[ grepl( 'Combined_data', all_files ) ]
if ( length( fname ) > 0 ) {
  file.remove( file = fname )
}

# Read in pilot study data
fname = all_files[ grepl( 'CognitionAndAdolesce_DATA', all_files ) ]
pilot_dat = read.csv( file = fname,
                      header = T,
                      stringsAsFactors = F
)

# Read in active study data
fname = all_files[ grepl( 'ARCHESK23ActiveTrial_DATA', all_files ) ]
active_dat = read.csv( file = fname,
                       header = T,
                       stringsAsFactors = F
)
# Remove subjects with no ID label
active_dat = 
  active_dat %>% 
  filter( id != "" )

# Read in original variable names from initial analysis
orig_vrb = read.csv( file = 'Original_variable_names.csv',
                     header = T,
                     stringsAsFactors = F
)
orig_vrb = orig_vrb$Variables
# Read in subject IDs
orig_id = read.csv( file = 'Original_subject_IDs.csv',
                    header = T,
                    stringsAsFactors = F
)
orig_id = orig_id$ID

# Read in data dictionaries for pilot/active data
dd_pilot = read.csv(
  file = "CognitionAndAdolescentHealthAR_DataDictionary_2018-11-02.csv",
  header = T,
  stringsAsFactors = F
)
dd_active = read.csv(
  file = "ARCHESK23ActiveTrialUpdated_DataDictionary_2018-11-02.csv",
  header = T,
  stringsAsFactors = F
)
# Rename variable field name
colnames( dd_active )[1] = 'Variable...Field.Name'

# Clean up workspace
rm( all_files, fname )

# 1.1) 
update_cd = function() {
  # Purpose:
  # Helper function to quickly update the 
  # short-hand variable containing 
  # column names for the combined dataset.
  # Returns:
  # A new vector of column names, 'cd'.
  
  cd <<- colnames( combined_dat )
  
}

# 1.2)
quick_proprogate = function( vrb, df ) {
  # Purpose:
  # Function to proprogate singular values for a 
  # given subject over repeated measures.
  # Arguments:
  # vrb - Character string for variable of interest 
  # df  - Data frame of interest
  # Returns:
  # A vector with the repeated values for each subject.
  
  # Isolate subjects
  subj = unique( df$id )
  
  # Initialize output
  out = rep( NA, nrow( df ) )
  
  # Loop over subjects
  for ( s in 1:length( subj ) ) {
    
    sel = df$id == subj[s]
    x = df[[ vrb ]][sel]
    no_na = !is.na( x )
    if ( any( no_na ) ) {
      out[sel] = x[no_na]
    }
    
  }
  
  return( out )
}

###
### 2) Initial creation of combined data set
###

# 2.1) Extract variable names

# Variables for pilot data set
pd = colnames( pilot_dat )
# Variables for active data set
ad = colnames( active_dat )

# Initial variables for combined data set
cd = pd[ pd %in% ad ]
cd = cd[ cd %in% orig_vrb ]

# Subjects to select from pilot data
subj = pilot_dat$id %in% orig_id

# 2.2) Add pilot data to combined data set

# Initialize data set starting with pilot data
combined_dat = pilot_dat[subj,cd]

# Use re-calculated variables 
# for days and times MJ used
combined_dat$tlfb_mj_14a = 
  pilot_dat$tlfb_mj_14a_recalc[subj]
combined_dat$tlfb_mj_14b = 
  pilot_dat$tlfb_mj_14b_recalc[subj]
combined_dat$tlfb_mj_15a = 
  pilot_dat$tlfb_mj_15a_recalc[subj]
combined_dat$tlfb_mj_15b = 
  pilot_dat$tlfb_mj_15b_recalc[subj]

# Use re-calculated variables 
# for MJ grams variable
combined_dat$tlfb_mj_16a = 
  pilot_dat$tlfb_mj_16a_recalc[subj]
combined_dat$tlfb_mj_16b = 
  pilot_dat$tlfb_mj_16b_recalc[subj]

# Define variable to track different data sets
active = rep( F, nrow( combined_dat ) )

# Update variable names in combined data set
cd = update_cd()

# 2.3) Assess differences between data sets

# Document similarities/differences between 
# data sets for TLFB
notes_on_tlfb_btw_datasets = data.frame(
  Variables = cd[ grep( 'tl', cd ) ],
  stringsAsFactors = F
)
# Add field labels for variables (Pilot)
sel = sapply( notes_on_tlfb_btw_datasets$Variables,
              function(x)
                which( dd_pilot$Variable...Field.Name == x ) )
notes_on_tlfb_btw_datasets$Pilot = 
  dd_pilot$Field.Label[sel]
# Add field labels for variables (Active)
sel = sapply( notes_on_tlfb_btw_datasets$Variables,
              function(x)
                which( dd_active$Variable...Field.Name == x ) )
notes_on_tlfb_btw_datasets$Active = 
  dd_active$Field.Label[sel]
# Track which variables are similar between data sets
notes_on_tlfb_btw_datasets$Same = FALSE
notes_on_tlfb_btw_datasets$Same[1:20] = T
# Indicate which variable matches in Active data set
notes_on_tlfb_btw_datasets$Active_that_matches = 
  notes_on_tlfb_btw_datasets$Variables
sel = !notes_on_tlfb_btw_datasets$Same
notes_on_tlfb_btw_datasets$Active_that_matches[sel] = 
  ""
# Match variables that differ
sel = notes_on_tlfb_btw_datasets$Variables == 'tlfb_nic_4'
notes_on_tlfb_btw_datasets$Active_that_matches[sel] = 'tlfb_nic_6'
sel = notes_on_tlfb_btw_datasets$Variables == 'tlfb_nic_5'
notes_on_tlfb_btw_datasets$Active_that_matches[sel] = 'tlfb_nic_7'
sel = notes_on_tlfb_btw_datasets$Variables == 'tlfb_nic_6'
notes_on_tlfb_btw_datasets$Active_that_matches[sel] = 'tlfb_nic_8'

# View( notes_on_tlfb_btw_datasets )

# 2.4) Add active data to combined data set

# Certain variables need to be proprogated over 
# all rows because they were measured during 
# screening
active_dat$sex_conf = 
  quick_proprogate( 'sex_conf', active_dat )
active_dat$age_exact = 
  quick_proprogate( 'age_exact', active_dat )
active_dat$bmi = 
  quick_proprogate( 'bmi', active_dat )

# Also adjust variables that are foundations 
# for computing new variables in combined 
# data set

# Switch 'race' to equal 'race_conf'
active_dat$race = active_dat$race_conf
active_dat$race = 
  quick_proprogate( 'race', active_dat )

# Set values for certain TLFB to 
# visit 2 in active data set
tlfb = ad[ grep( 'tl', ad ) ]
for ( vrb in tlfb ) {
  
  subj = unique( active_dat$id )
  sel_1 = active_dat$visit_number == 2
  v = active_dat[[ vrb ]][ sel_1 ]
  if ( all( is.na( v ) ) ) {
    # Loop over subjects
    for ( s in 1:length( subj ) ) {
      sel_1 = active_dat$id == subj[s] & 
        active_dat$visit_number == 1 & 
        !is.na( active_dat$visit_number )
      sel_2 = active_dat$id == subj[s] & 
        active_dat$visit_number == 2 & 
        !is.na( active_dat$visit_number )
      active_dat[[ vrb ]][ sel_2 ] = 
        active_dat[[ vrb ]][ sel_1 ]
    }
  }
  
}

# Subjects to select from active data
subj = active_dat$id %in% orig_id

# Concatenate with active data set
combined_dat = rbind(
  combined_dat,
  active_dat[subj,cd]
)
active = c( active, rep( T, sum( subj ) ) )
combined_dat$active = active

# Adjust mis-matched variables between 
# pilot and active studies
for ( i in which( !notes_on_tlfb_btw_datasets$Same  ) ) {
  
  vrb = notes_on_tlfb_btw_datasets$Variables[i]
  mtch = notes_on_tlfb_btw_datasets$Active_that_matches[i]
  if ( mtch != '' ) {
    combined_dat[[ vrb ]][ combined_dat$active ] = 
      active_dat[[ mtch ]][ subj ]
  }
}

# Update column names
update_cd()

# Clean up workspace
rm( subj, active, i, mtch, sel, vrb,
    sel_1, sel_2, s, v, tlfb )

###
### 3) Add in missing variables and adjust values
###

# 3.1) Adjust visit number

combined_dat$old_visit_number = combined_dat$visit_number
combined_dat$visit_number = NA

# Loop over visit numbers
for ( i in 1:7 ) {
  
  # Pilot study
  vrb = paste( 'v', i, sep = '' )
  sel = grepl( vrb, combined_dat$redcap_event_name ) & 
    !combined_dat$active
  combined_dat$visit_number[sel] = i
  # Active study
  vrb = paste( 'v', i+1, sep = '' )
  sel = grepl( vrb, combined_dat$redcap_event_name ) & 
    combined_dat$active
  combined_dat$visit_number[sel] = i
  
}

# 3.2) Code race (white vs. non-white)

# Initialize data
combined_dat$race_revised = rep( NA, nrow( combined_dat ) )
# Loop over subjects with data for race in the 
# pilot and active studies
tmp = combined_dat %>% 
  filter( !is.na( race ) ) %>% 
  group_by( id ) %>% 
  summarize( 
    R = unique( race )
  )
for ( i in 1:nrow( tmp ) ) {
  sel = combined_dat$id == tmp$id[i]
  r = tmp$R[i]
  if ( r == 1 ) {
    combined_dat$race_revised[sel] = 1
  } else {
    combined_dat$race_revised[sel] = 2
  }
}
# Some subjects had race information coded in 
# the varible 'school_race' for the pilot study
tmp = pilot_dat %>% 
  filter( !is.na( school_race ) & id %in% combined_dat$id ) %>% 
  group_by( id ) %>% 
  summarize( 
    R = unique( school_race )
  )
for ( i in 1:nrow( tmp ) ) {
  sel = combined_dat$id == tmp$id[i]
  r = tmp$R[i]
  if ( r == 1 ) {
    combined_dat$race_revised[sel] = 1
  } else {
    combined_dat$race_revised[sel] = 2
  }
}

# 3.3) Create variable for years of MJ use

# Pilot
combined_dat$years_mj_use = 
  (combined_dat$age_exact/365) - 
  combined_dat$tlfb_mj_2
sel = combined_dat$active
# Active
combined_dat$years_mj_use[sel] = 
  combined_dat$age_exact[sel] - 
  combined_dat$tlfb_mj_2[sel]

# 3.4) Recency of MJ use

# For pilot study, TLFB item 8 
# reports number of days since last use
rmju = combined_dat$tlfb_mj_8
# For active study, based on different 
# TLFB item (5) and depends potentially on 
# two time points (visits 1 and/or 2)
id_sel = active_dat$id %in% combined_dat$id
sel = combined_dat$active
rmju[sel] = active_dat$tlfb_mj_fu5[id_sel]
# Exclude values for cases past visit 1
rmju[sel & combined_dat$visit_number > 1 ] = NA
# Check cases in which subjects didn't use between 
# first and second visits
tmp = active_dat %>% 
  filter( id_sel ) %>% 
  group_by( id, redcap_event_name ) %>% 
  summarize(
    V = unique( tlfb_mj_fu_use ),
    D = date
  )
tmp$Diff_V1 = NA
# Loop over subject IDs
for ( s in unique( tmp$id ) ) {
  # Select subject
  sel = tmp$id == s
  # Isolate row for visit 1
  v1 = grepl( 'v1_', tmp$redcap_event_name[sel] )
  if ( !is.na( tmp$D[sel][v1] ) ) {
    # Compute difference
    tmp$Diff_V1[sel] = as.Date(tmp$D[sel]) - as.Date(tmp$D[sel][v1])
  } else {
    tmp$Diff_V1[sel] = NA
  }
}
tst = tmp %>% 
  filter( grepl( 'v2_', redcap_event_name ) )
tst = tst %>% 
  filter( V == 0 )
# Loop over selected subjects
for ( s in 1:nrow( tst ) ) {
  sel = combined_dat$id == tst$id[s]
  id_sel = active_dat$id == tst$id[s] #& 
    #active_dat$visit_number > 1 & 
    #active_dat$visit_number < 9
  rmju[sel] = 
    combined_dat$tlfb_mj_8[id_sel] + tst$Diff_V1[s]
}
# Recode TLFB item 8
combined_dat$tlfb_mj_8 = rmju

# 3.5) Days since last used nicotine

# For pilot study, TLFB item 8b 
# reports number of days since last use
combined_dat$tlfb_nic_8b = NA
sel = combined_dat$pilot
id_sel = pilot_dat$id %in% combined_dat$id
combined_dat$tlfb_nic_8b[sel] = pilot_dat$tlfb_nic_8b[id_sel]
# For active study, based on different 
# TLFB item (4) and depends potentially on 
# two time points (visits 1 and/or 2)
id_sel = active_dat$id %in% combined_dat$id
sel = combined_dat$active

combined_dat$tlfb_nic_8b[sel] = 
  active_dat$tlfb_nic_fu4[id_sel]
# Check cases in which subjects didn't use between 
# first and second visits
tmp = active_dat %>% 
  filter( id_sel ) %>% 
  group_by( id, redcap_event_name ) %>% 
  summarize(
    V = unique( tlfb_nic_fu_use ),
    D = date
  )
tmp$Diff_V1 = NA
# Loop over subject IDs
for ( s in unique( tmp$id ) ) {
  # Select subject
  sel = tmp$id == s
  # Isolate row for visit 1
  v1 = grepl( 'v1_', tmp$redcap_event_name[sel] )
  if ( !is.na( tmp$D[sel][v1] ) ) {
    # Compute difference
    tmp$Diff_V1[sel] = as.Date(tmp$D[sel]) - as.Date(tmp$D[sel][v1])
  } else {
    tmp$Diff_V1[sel] = NA
  }
}
tst = tmp %>% 
  filter( grepl( 'v2_', redcap_event_name ) )
tst = tst %>% 
  filter( V == 0 )
# Loop over selected subjects
for ( s in 1:nrow( tst ) ) {
  sel = combined_dat$id == tst$id[s]
  id_sel = active_dat$id == tst$id[s]
  combined_dat$tlfb_nic_8b[sel] = 
    active_dat$tlfb_nic_14[id_sel] + tst$Diff_V1[s]
}

# 3.6) Assorted logical variables

# Censored THC data
combined_dat$thc_500_flag = 
  combined_dat$ua_conf_result == '501'
# SET VALUES FLAGGED TO MISSING

# Flag for no recent use of THC
combined_dat$thc_recency_flag = 
  combined_dat$tlfb_mj_8 > 2

# Nicotine use
combined_dat$nic_user = NA
tmp = combined_dat %>% 
  filter( !active ) %>% 
  group_by( id ) %>% 
  summarize(
    U = unique( tlfb_nic_5 )[2]
  )
# Loop over subjects
for ( s in 1:nrow( tmp ) ) {
  sel = combined_dat$id == tmp$id[s] & 
    !is.na( combined_dat$visit_number )
  ind = min( which( sel ) )
  if ( !is.na( tmp$U[s] ) & 
       tmp$U[s] > 0 ) {
    combined_dat$nic_user[ind] = 1
  } else {
    combined_dat$nic_user[ind] = 0
  }
  
}
tmp = active_dat %>% 
  filter( id %in% combined_dat$id ) %>% 
  group_by( id ) %>% 
  summarize(
    U = unique( tlfb_nic_8 )[2]
  )
# Loop over subjects
for ( s in 1:nrow( tmp ) ) {
  sel = combined_dat$id == tmp$id[s] & 
    !is.na( combined_dat$visit_number )
  ind = min( which( sel ) )
  if ( !is.na( tmp$U[s] ) & 
       tmp$U[s] > 0 ) {
    combined_dat$nic_user[ind] = 1
  } else {
    combined_dat$nic_user[ind] = 0
  }
  
}

# 3.7) Add more mis-matched variables

combined_dat$tlfb_nic_7a = NA
sel = combined_dat$id %in% pilot_dat$id
id_sel = pilot_dat$id %in% combined_dat$id
combined_dat$tlfb_nic_7a[sel] = 
  pilot_dat$tlbf_nic_7a[id_sel]
sel = combined_dat$active
id_sel = active_dat$id %in% combined_dat$id
combined_dat$tlfb_nic_7a[sel] = 
  active_dat$tlbf_nic_9[id_sel]

combined_dat$tlfb_nic_10a = NA
sel = combined_dat$id %in% pilot_dat$id
id_sel = pilot_dat$id %in% combined_dat$id
combined_dat$tlfb_nic_10a[sel] = 
  pilot_dat$tlfb_nic_10a_recalc[id_sel]
sel = combined_dat$active
id_sel = active_dat$id %in% combined_dat$id
combined_dat$tlfb_nic_10a[sel] = 
  active_dat$tlfb_nic_10[id_sel]

combined_dat$tlfb_nic_10b = NA
sel = combined_dat$id %in% pilot_dat$id
id_sel = pilot_dat$id %in% combined_dat$id
combined_dat$tlfb_nic_10b[sel] = 
  pilot_dat$tlfb_nic_10b_recalc[id_sel]
sel = combined_dat$active
id_sel = active_dat$id %in% combined_dat$id
combined_dat$tlfb_nic_10b[sel] = 
  active_dat$tlfb_nic_11[id_sel]

# 3.8) Adjust units for years (pilot)

sel = !combined_dat$active
combined_dat$age_exact[ sel ] = 
  combined_dat$age_exact[ sel ]/365

# 3.9) Non-critical variables

combined_dat$X = NA
combined_dat$week = NA
combined_dat$height_ft = NA
combined_dat$height_in = NA
combined_dat$height_calc = NA
combined_dat$height = NA
combined_dat$bmiheightsquared = NA
combined_dat$urinalysis_comments = NA
combined_dat$thc_heavy_flag = NA
combined_dat$tlfb_nic_9_b = NA
combined_dat$weekly_use = NA

# Remove comments
combined_dat$ua_comment = " "

# Update names
cd = colnames( combined_dat )

# Remove unnecessary visits
combined_dat = 
  combined_dat %>% 
  filter( !is.na( visit_number ) )

# Remove extraneous rows
combined_dat = 
  combined_dat %>% 
  filter( id != "" )

# Clean up workspace
rm( i, tmp, tst, id_sel, ind, 
    r, s, sel, v1, vrb, rmju )

###
### 4) Adjustments following trimming of rows
###

# 4.1) Days since baseline

# Initialize variables
combined_dat$days_since_baseline = NA
subj = unique( combined_dat$id )
for ( s in 1:length( subj ) ) {
  
  # Compute difference between date for each visit  
  # and the date for visit 1
  sel = combined_dat$id == subj[s]
  v = combined_dat$date[sel & combined_dat$visit_number == 1]
  tmp = as.Date(combined_dat$date[sel]) - as.Date( v )
  combined_dat$days_since_baseline[sel] = 
    as.numeric( tmp )
  
}

# 4.2) Adjust dipstick results for pilot group

sel = !combined_dat$active
combined_dat$thc_qual[sel] = combined_dat$thc_qual[sel] - 1

# 4.3) Imput zero for THC values not detected

sel = combined_dat$ua_thc_flag == 1 & 
  !is.na( combined_dat$ua_thc_flag )
combined_dat$ua_thclcmsms_crtadj[sel] = 0

# 4.4) Misc. adjustments

# Set variables for sex, BMI, and age 
# to only first visit

subj = unique( combined_dat$id )
# Loop over variables
for ( vrb in c(
  'sex_conf', 
  'bmi', 
  'age_exact',
  'race_revised' ) ) {
  
  # Loop over subjects
  for ( s in 1:length( subj ) ) {
    sel = combined_dat$id == subj[s]
    x = combined_dat[[ vrb ]][sel]
    if ( sum( !is.na( x ) ) > 1 ) {
      combined_dat[[ vrb ]][sel] = NA
      combined_dat[[ vrb ]][sel][1] = unique( x[ !is.na(x) ] )
      
    }
  }
}

# 4.5) Remove values for non-abstinent cases

# Active study
# 003_COMM - Last visit w/ no issues - v6 (Redcap event name)
sel = combined_dat$id == '003_COMM' & 
  combined_dat$visit_number > 5
combined_dat$ua_thclcmsms_crtadj[sel] = NA
# Also add in missing value for recency
missing_val = 
  active_dat$tlfb_mj_8[ active_dat$id == '003_COMM' ][2]
combined_dat$tlfb_mj_8[ combined_dat$id == '003_COMM' ][1] = 
  missing_val

# 006_COMM - Last visit w/ no issues - v6 (Redcap event name)
sel = combined_dat$id == '006_COMM' & 
  combined_dat$visit_number > 5
combined_dat$ua_thclcmsms_crtadj[sel] = NA

# 065_COMM - Lost to follow-up, no data, should be fine

# 084_COMM - Last visit w/ no issues - v4 (Redcap event name)
sel = combined_dat$id == '084_COMM' & 
  combined_dat$visit_number > 3
combined_dat$ua_thclcmsms_crtadj[sel] = NA

# 086_COMM - Last visit w/ no issues - v5 (Redcap event name)
sel = combined_dat$id == '086_COMM' & 
  combined_dat$visit_number > 4
combined_dat$ua_thclcmsms_crtadj[sel] = NA

# Pilot study
# 10005 - Last visit w/ no issues - v5 (Redcap event name)
sel = combined_dat$id == '10005' & 
  combined_dat$visit_number > 5
combined_dat$ua_thclcmsms_crtadj[sel] = NA

# 10037 - Last visit w/ no issues - v5 (Redcap event name)
sel = combined_dat$id == '10037' & 
  combined_dat$visit_number > 5
combined_dat$ua_thclcmsms_crtadj[sel] = NA

# 4.6) Adjust values for late-start abstinence

# Manually change values for these subjects
# 20014
# tlbf_mj_8 == 1, visit_number = 2
# Change days_since_baseline relative to visit 2
sel = combined_dat$id == '20014' & 
  combined_dat$days_since_baseline == 0
combined_dat[ sel, grepl( 'ua_', cd ) ] = NA
combined_dat$days_since_baseline[sel] = NA
combined_dat$tlfb_mj_8[sel] = NA
combined_dat$thc_qual[sel] = NA
# Adjust days since baseline
sel = combined_dat$id == '20014' & 
  !is.na( combined_dat$days_since_baseline )
combined_dat$days_since_baseline[sel] = 
  combined_dat$days_since_baseline[sel] - 2
combined_dat$tlfb_mj_8[sel][1] = 1

# 20024
# tlbf_mj_8 == 1, visit_number = 2
# Change days_since_baseline relative to visit 2
sel = combined_dat$id == '20024' & 
  combined_dat$days_since_baseline == 0
combined_dat[ sel, grepl( 'ua_', cd ) ] = NA
combined_dat$days_since_baseline[sel] = NA
combined_dat$tlfb_mj_8[sel] = NA
combined_dat$thc_qual[sel] = NA
# Adjust days since baseline
sel = combined_dat$id == '20024' & 
  !is.na( combined_dat$days_since_baseline )
combined_dat$days_since_baseline[sel] = 
  combined_dat$days_since_baseline[sel] - 2
combined_dat$tlfb_mj_8[sel][1] = 1

# Clean up workspace
rm( s, sel, subj, vrb, x, missing_val )

# 4.7) Set censored values as missing
sel = combined_dat$thc_500_flag & !is.na( combined_dat$thc_500_flag )
# Set values to missing
combined_dat$thc_qual[sel] = NA
combined_dat$ua_conf_result[sel] = NA
combined_dat$ua_thclcmsms_crtadj[sel] = NA

###
### 5) Save results
###

combined_dat %>% 
  write.table(
    file = paste( 'Combined_data_',
                  Sys.Date(),
                  '.csv', 
                  sep = '' ),
    row.names = F,
    quote = F,
    sep = ',' )
# Save to upper directory as well
setwd( orig_dir )
# Remove any previous instances of combined data
all_files = dir()
fname = all_files[ grepl( 'Combined_data', all_files ) ]
if ( length( fname ) > 0 ) {
  file.remove( file = fname )
}
combined_dat %>% 
  write.table(
    file = paste( 'Combined_data_',
                  Sys.Date(),
                  '.csv', 
                  sep = '' ),
    row.names = F,
    quote = F,
    sep = ',' )
# Clean up workspace
rm( all_files, fname )

setwd( pa_dir )

