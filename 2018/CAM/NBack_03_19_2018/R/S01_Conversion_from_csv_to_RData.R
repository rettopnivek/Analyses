# Conversion of .csv to .RData
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you 
# have any questions or comments
# Last updated 2018-04-11

# Table of contents
# 1) Initial setup
# 2) Load in behavioral data
# 3) Modify variable names and create data key
# 4) Load in neural data and modify
# 5) Load in demographics and vitals
# 6) Convert raw data into a tidier form
# 7) Save as .RData file

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

# Packages for creating tidy data
# install.packages( 'tidyr' )
library(tidyr)
# install.packages( 'dplyr' )
library(dplyr)

###
### 2) Load in behavioral data
###

# Navigate to location of data
setwd( 'Data/Original_files' )

# Load in data
rawDat = read.csv( 
  file = 'Copy of Adrian_Redcap_Drug_Time.csv', 
  header = T,
  stringsAsFactors = F )
# Rename first column to remove 'i..'
colnames( rawDat )[1] = 'Drug_Or_Placebo'

# Subject FN_041 has missing data that proves 
# problematic for later analyses. Therefore 
# reorder 'Subject' variable so that FN_041
# equals 66, so for cases where FN_041 is 
# excluded from the data, the sequential
# order of 'Subject' is preserved with no skips

# Temporarily change FN_041 ID
rawDat$study_id[ rawDat$study_id == 'FN_041' ] = 'FN_141'

# Create a new variable labeling subjects from 1 to N
# (where N is the maximum number of subjects)
rawDat = rawDat %>% 
  mutate( Subject = utilityf::createIncrement( study_id ) )

# Convert back to correct ID for FN_041
rawDat$study_id[ rawDat$study_id == 'FN_141' ] = 'FN_041'

# Extract number of subjects
N = length( unique( rawDat$study_id ) )

###
### 3) Modify variable names and create data key
###

# Original labels
orig_var_names = colnames( rawDat )
# Include column numbers
names( orig_var_names ) = 
  as.character( 1:length( orig_var_names ) )

# Rename variables to make it easier 
# to use 'tidyr' functions
rawDat = rawDat %>% 
  rename(
    Subject = Subject,
    ID = study_id,
    Condition = Drug_Or_Placebo,
    Timepoints = Pre_Post_Post2,
    Self_report_on_high = SR.High,
    Hits.Combined = slnb3_tp,
    False_alarms.Combined = slnb3_fp,
    Correct_rejections.Combined = true_negative_responses, 
    Misses.Combined = false_negative_responses,
    Number_of_positive_trials.Combined = total_tp_poss,
    Median_RT_Correct.Combined = slnb3_rtc,
    Hits.NBack_0 = slnb3_tp2,
    False_alarms.NBack_0 = slnb3_fp2,
    Correct_rejections.NBack_0 = true_negative_responses_fo,
    Misses.NBack_0 = false_negatives_responses,
    Number_of_positive_trials.NBack_0 = total_tp_0back_poss,
    Median_RT_Correct.NBack_0 = slnb3_rtc2,
    Hits.NBack_2 = slnb3_tp3,
    False_alarms.NBack_2 = slnb3_fp3,
    Correct_rejections.NBack_2 = true_neg_responses,
    Misses.NBack_2 = false_negative_responses_f,
    Number_of_positive_trials.NBack_2 = total_tp_2back_poss,
    Median_RT_Correct.NBack_2 = slnb3_rtc3,
    Mean_for_median_RT_over_0_and_2_Back = slnb3_mrtc,
    Efficiency = slnb3_meff,
    Longest_run_of_non_responses = slnb3_lrnr,
    Longest_possible_run_of_non_responses = 
      the_longest_possible_run_o,
    Final_run_of_non_responses = slnb3_frnr,
    Possible_final_run_of_non_responses = 
      the_possible_final_run_of,
    Comments = nback_comments2,
    Completion = nback_abde_complete,
    nback = nback,
    Select_which_run_was_completed = select_which_run_was_compl___0
  )
# Create additional variables tracking number of trials 
# for false alarms
rawDat$Number_of_negative_trials.NBack_0 = 
  rawDat$False_alarms.NBack_0 + 
  rawDat$Correct_rejections.NBack_0
rawDat$Number_of_negative_trials.NBack_2 = 
  rawDat$False_alarms.NBack_2 + 
  rawDat$Correct_rejections.NBack_2
rawDat$Number_of_negative_trials.Combined = 
  rawDat$False_alarms.Combined + 
  rawDat$Correct_rejections.Combined

# Add variables to list of original names
nv = ncol( rawDat )
orig_var_names = c(
  orig_var_names,
  colnames( rawDat )[ (nv-2):nv ]
)

# Extract new labels
new_var_names = colnames( rawDat )

# Informative descriptions for variables
desc = c(
  # 1
  'Drug (THC) versus placebo condition',
  # 2
  'Measurement time point relative to drug administration',
  # 3
  'Subject ID number',
  # 4
  'Self report on how high, (0 is no high, 100 is max high)',
  # 5
  paste( 'Number of hits',
         '(correctly pressed for true positives for both tasks)' ),
  # 6
  'Total number of true positive trials for  both tasks',
  # 7
  paste( 'Number of false alarms',
         '(incorrectly pressed for true negative for  both tasks)' ),
  # 8
  paste( 'Number of correct rejections', 
         '(correct non-presses for true negatives for  both tasks)' ),
  # 9
  paste( 'Number of misses', 
         '(incorrect non-presses for', 
         'true positives for  both tasks)' ),
  # 10
  'Median response time for all correct responses (ms)',
  # 11
  paste( 'Number of hits', 
         '(correctly pressed for true positives for 0-back)' ),
  # 12
  'Total number of true positive trials for 0-back',
  # 13
  paste( 'Number of false alarms',
         '(incorrectly pressed for true negative for 0-back)' ),
  # 14
  paste( 'Number of correct rejections', 
         '(correct non-presses for true negatives for 0-back)' ),
  # 15
  paste( 'Number of misses', 
         '(incorrect non-presses for true positives for 0-back)' ),
  # 16
  'Median response time for correct 0-Back trials (ms)',
  # 17
  paste( 'Number of hits', 
         '(correctly pressed for true positives for 2-back)' ),
  # 18
  'Total number of true positive trials for 2-back',
  # 19
  paste( 'Number of false alarms',
         '(incorrectly pressed for true negative for 2-back)' ),
  # 20
  paste( 'Number of correct rejections', 
         '(correct non-presses for true negatives for 2-back)' ),
  # 21
  paste( 'Number of misses', 
         '(incorrect non-presses for true positives for 2-back)' ),
  # 22
  'Median response time for correct 2-Back trials (ms)',
  # 23
  paste( 'Mean of median RT for hits over the 0-Back', 
         'and 2-Back Trials (ms)' ),
  # 24
  paste( 'Efficiency for correct 0-Back and 2-Back trials', 
         '(efficiency =', 
         '(Hits for 3-back) /',
         'log( mean of median RT for 0 and 2-back) )' ),
  # 25
  'Longest run of non_responses ',
  # 26
  'The longest possible run of non_responses',
  # 27
  'Final run of non_responses',
  # 28
  'The possible final run of non_responses',
  # 29
  'Comments regarding the N-back task',
  # 30
  'Completion',
  # 31
  'Was the N-back test administered? (1 = yes)',
  # 32
  'Selection for which run was completed (Run 1)',
  # 33
  'Index for subjects',
  # 34
  'Total number of true negative trials for 0-back',
  # 35
  'Total number of true negative trials for 2-back',
  # 36
  'Total number of true negative trials for 3-back'
)

# Create a data key
data_key = data.frame(
  Original_labels = orig_var_names,
  New_labels = new_var_names,
  Description = desc
)

# Initialize data frame giving details 
# about which parts of the study subjects 
# completed
study_completion = data.frame(
  ID = '???', 
  Subject = rep( 1:N, each = 12 ), 
  Condition = '???', 
  Task = '???', 
  Timepoints = '???', 
  Completed = 0, 
  Positive_trials = 0, 
  Negative_trials = 0,
  Total_trials = 0,
  stringsAsFactors = FALSE
)

# Fill in conditions
tmp = rawDat %>% 
  group_by( Subject ) %>% 
  summarize( ID = unique( ID ) )
study_completion$ID = rep( tmp$ID, each = 12 )
study_completion$Condition = 
  rep( rep( c( 'Placebo', 'Drug' ), each = 6 ), N )
study_completion$Task = 
  rep( 
    rep( rep( c( '0-Back', '2-Back' ), each = 3 ), 2 ), N )
study_completion$Timepoints = 
  rep( 
    rep( c( 'T1_Pre_drug', 'T2_Post_drug', 'T3_Post_drug' ), 4 ),
    N )

# Loop over subjects
for ( s in 1:N ) {
  
  tmp = rawDat %>% 
    filter( Subject == s )
  
  Cnd = tmp$Condition
  TP = tmp$Timepoints
  Tsk0_pos = tmp$Number_of_positive_trials.NBack_0
  Tsk0_neg = tmp$Number_of_negative_trials.NBack_0
  Tsk2_pos = tmp$Number_of_positive_trials.NBack_2
  Tsk2_neg = tmp$Number_of_negative_trials.NBack_2
  
  sel = study_completion$Subject == s
  
  for ( i in 1:nrow( tmp ) ) {
    
    # Check whether a given set of conditions 
    # has been completed
    inc0 = 0
    inc2 = 0
    if ( Cnd[i] == 'Placebo' ) {
      
      if ( TP[i] == 'Pre' ) {
        if ( Tsk0_pos[i] > 0 ) inc0 = 1
        if ( Tsk2_pos[i] > 0 ) inc2 = 4
      }
      
      if ( TP[i] == 'Post' ) {
        if ( Tsk0_pos[i] > 0 ) inc0 = 2
        if ( Tsk2_pos[i] > 0 ) inc2 = 5
      }
      
      if ( TP[i] == 'Post2' ) {
        if ( Tsk0_pos[i] > 0 ) inc0 = 3
        if ( Tsk2_pos[i] > 0 ) inc2 = 6
      }
      
    }
    if ( Cnd[i] == 'Drug' ) {
      
      if ( TP[i] == 'Pre' ) {
        if ( Tsk0_pos[i] > 0 ) inc0 = 1 + 6
        if ( Tsk2_pos[i] > 0 ) inc2 = 4 + 6
      }
      
      if ( TP[i] == 'Post' ) {
        if ( Tsk0_pos[i] > 0 ) inc0 = 2 + 6
        if ( Tsk2_pos[i] > 0 ) inc2 = 5 + 6
      }
      
      if ( TP[i] == 'Post2' ) {
        if ( Tsk0_pos[i] > 0 ) inc0 = 3 + 6
        if ( Tsk2_pos[i] > 0 ) inc2 = 6 + 6
      }
      
    }
    
    # Fill in rows
    if ( inc0 != 0 ) {
      study_completion$Completed[sel][inc0] = 
        study_completion$Completed[sel][inc0] + 1
      study_completion$Positive_trials[sel][inc0] = 
        Tsk0_pos[i]
      study_completion$Negative_trials[sel][inc0] = 
        Tsk0_neg[i]
      study_completion$Total_trials[sel][inc0] = 
        Tsk0_pos[1] + Tsk0_neg[i]
    }
    if ( inc2 != 0 ) {
      study_completion$Completed[sel][inc2] = 
        study_completion$Completed[sel][inc2] + 1
      study_completion$Positive_trials[sel][inc2] = 
        Tsk2_pos[i]
      study_completion$Negative_trials[sel][inc2] = 
        Tsk2_neg[i]
      study_completion$Total_trials[sel][inc2] = 
        Tsk2_pos[1] + Tsk2_neg[i]
    }
    
  }
  
}

###
### 4) Load in neural data and modify
###

# Read in data with BOLD measurements
neurDat = read.csv( file = 'Neural_data.csv',
                    header = T, stringsAsFactors = FALSE )
# Note that not all subjects have both behavioral and 
# neural data

# The combined channels represent:
# Right dorsolateral prefrontal cortex (R. DLPFC)
#   -> Average of channels 10, 15, 17, and 18
# Left dorsolateral prefrontal cortex (L. DLPFC)
#  -> Average of channels 1, 2, 5, and 8
# Medial prefrontal cortex (MPFC)
#  -> Average of channels 7, 9, 12, and 14
# Right ventrolateral prefrontal cortex (R. VLPFC)
#  -> Average of channels 13, 16, 19, and 20
# Left ventrolateral prefrontal cortex (L. VLPFC)
#  -> Average of channels 3, 4, 6, and 11

# Recompute averages
roi = c( 'R_DLPFC', 'L_DLPFC', 'MPFC',
         'R_VLPFC', 'L_VLPFC' )
# Channels for each ROI
adj = min( grep( 'Channel', colnames( neurDat ) ) - 1 )
channel = list(
  c( 10, 15, 17, 18 ) + adj,
  c( 1, 2, 5, 8 ) + adj,
  c( 7, 9, 12, 14 ) + adj,
  c( 13, 16, 19, 20 ) + adj,
  c( 3, 4, 6, 11 ) + adj
)

# Initialize matrix
new_roi_avg = matrix( NA, nrow( neurDat ), 5 )
colnames( new_roi_avg ) = roi

# Average over channels for each ROI
for ( i in 1:5 ) {
  sel = channel[[i]]
  # Average over the set of 4 channels
  new_roi_avg[,roi[i]] = 
    rowMeans( neurDat[,sel], na.rm = T )
}

# Original values in .csv file and 
# recalculate values differ slightly 
# due to rounding error in former
b("
# Double check that averages for ROI are 
# correct

avg_correct = round( as.matrix( neurDat[,roi] ), 8 ) == 
  round( new_roi_avg, 8 )
i = which( rowSums( avg_correct ) != 5 )
")

for ( i in 1:5 ) 
  neurDat[,roi[i]] = new_roi_avg[,roi[i]]

# Add subject indices
neurDat$Subject = NA
for ( n in 1:nrow( neurDat ) ) {
  sel = rawDat$ID == neurDat$ID[n]
  s = unique( rawDat$Subject[sel] )
  neurDat$Subject[n] = s
}

###
### 5) Load in demographics and vitals
###

# Baseline pulse readings
baseline_pulse = read.csv( file = 'Adrian_Baseline_Vitals.csv',
                           header = FALSE, 
                           stringsAsFactors = FALSE )
colnames( baseline_pulse ) = c( 'ID', 'Pulse_string' )
# Convert into numeric values and code missing cases as NA
baseline_pulse$Pulse = NA
sel = baseline_pulse$Pulse_string != '-'
baseline_pulse$Pulse[sel] = 
  as.numeric( baseline_pulse$Pulse_string[sel] )

# Dosage information
dosage = read.csv( file = 'Adrian_Dosage.csv',
                   header = FALSE,
                   stringsAsFactors = FALSE )
colnames( dosage ) = c( 'ID', 'Dosage_string' )
# Correct first ID value
dosage$ID[1] = 'FN_030'
# Convert into numeric values and code missing cases as NA
dosage$Dosage = NA
sel = dosage$Dosage_string != '-'
dosage$Dosage[sel] = 
  as.numeric( dosage$Dosage_string[sel] )

# Extract demographics information
raw_demo = read.csv( file = 'Demographics_sheet_2.csv',
                         header = TRUE, 
                         stringsAsFactors = FALSE )
# Eliminate blank rows
raw_demo = raw_demo %>% 
  filter( study_id != '' )

demographics = raw_demo[,1:21]
demographics = rbind(
  demographics[1:63,],
  c( 'FN_095', rep( NA, 19 ) ),
  demographics[64:65,] )
colnames( demographics ) = c(
  'ID', 'Gender', 'Age', 
  'Body_mass_index', 'Years_of_education',
  'State_trace_anxiety_index', 
  'Cannabis_use_disorder_index_total_scores',
  'Alcohol_use_disorder_index_score',
  'Age_of_first_cannabis_use',
  'Cannabis_use_frequency',
  'Cannabis_use_per_day',
  'Joints_per_week',
  'Cannabis_quantity_per_time',
  'UA_conf_result',
  'Self_predicted_cognitive_and_behavioral_impact',
  'Self_predicted_relaxation_and_tension_reduction',
  'Self_predicted_social_and_sexual_facilitation',
  'Self_predicted_perceptual_and_cognitive_effect',
  'Self_predicted_global_negative_effects_score',
  'Self_predicted_craving_and_phaysical_effect',
  'Marijuana_Effect_Expectancy_Questionnaire_total_score'
)
demographics$Gender_old = demographics$Gender
demographics$Gender = 'Male'
demographics$Gender[ demographics$Gender_old == 0 ] = 'Female'
demographics$Gender[ is.na( demographics$Gender_old ) ] = NA
demographics$Gender_old = NULL

###
### 6) Convert raw data into a tidier form
###

# 'Tidy' data is arranged so that each column is a variable,
# each row an observation

# Create a long-form version of the data, 
# based on the different count data for 
# response types (Hits, etc...) and the 
# task (3-back, etc...)
dat = rawDat %>% 
  gather( key = Response_type,
          value = Counts, 
          factor_key = FALSE, 
          # 3-back
          Hits.Combined,
          False_alarms.Combined,
          # 0-back
          Hits.NBack_0,
          False_alarms.NBack_0,
          # 2-back
          Hits.NBack_2,
          False_alarms.NBack_2 )
# Split 'Response_type' to separately delineate 
# the response type and task
dat = dat %>%
  separate(Response_type, 
           into = c("Response_type", "Task"), sep = "\\.")

# Add in total number of trials associated with 
# hits and false alarms
dat$Trials = 0
for ( n in 1:nrow( dat ) ) {
  
  for ( tsk in unique( dat$Task ) ) {
    if ( dat$Task[n] == tsk & 
         dat$Response_type[n] == 'Hits' ) {
      sel = paste( 'Number_of_positive_trials',
                   tsk, sep = '.' )
      dat$Trials[n] = dat[n,sel]
    }
    if ( dat$Task[n] == tsk & 
         dat$Response_type[n] == 'False_alarms' ) {
      sel = paste( 'Number_of_negative_trials',
                   tsk, sep = '.' )
      dat$Trials[n] = dat[n,sel]
    }
  }
}

# Add in median RT for hits
dat$Median_RT = NA
for ( n in 1:nrow( dat ) ) {
  
  for ( tsk in unique( dat$Task ) ) {
    if ( dat$Task[n] == tsk & 
         dat$Response_type[n] == 'Hits' ) {
      sel = paste( 'Median_RT_Correct',
                   tsk, sep = '.' )
      dat$Median_RT[n] = dat[n,sel]
    }
  }
}

# Extract the variables of interest
dat = dat %>% 
  dplyr::select( Subject, 
          ID, 
          Condition, 
          Timepoints, 
          Task, 
          Response_type, 
          Counts, 
          Trials,
          Median_RT, 
          Self_report_on_high
  )

# Re-label time points
old_labels = unique( dat$Timepoints )
new_labels = old_labels
new_labels[ old_labels == 'Pre' ] = 'T1_Pre_drug'
new_labels[ old_labels == 'Post' ] = 'T2_Post_drug'
new_labels[ old_labels == 'Post2' ] = 'T3_Post_drug'
for ( i in 1:3) {
  sel = dat$Timepoints == old_labels[i]
  dat$Timepoints[sel] = new_labels[i]
}

# Add in numerical code for response type
dat$Target = 1
dat$Target[ dat$Response_type == 'False_alarms' ] = 0

# Add in neural data
dat$R_DLPFC = NA
dat$L_DLPFC = NA
dat$MPFC = NA
dat$R_VLPFC = NA
dat$L_VLPFC = NA

# Loop over neural data
for ( n in 1:nrow( neurDat ) ) {
  
  # Match conditions
  sel = dat$Subject == neurDat$Subject[n] & 
    dat$Timepoints == neurDat$Timepoints[n] & 
    dat$Condition == neurDat$Condition[n] & 
    dat$Task == neurDat$Task[n]
  
  dat$R_DLPFC[sel] = neurDat$R_DLPFC[n]
  dat$L_DLPFC[sel] = neurDat$L_DLPFC[n]
  dat$MPFC[sel] = neurDat$MPFC[n]
  dat$R_VLPFC[sel] = neurDat$R_VLPFC[n]
  dat$L_VLPFC[sel] = neurDat$L_VLPFC[n]
  
}

# Add in dosage information
dat$Dosage = NA
dat$Dosage[ dat$Condition == 'Placebo' ] = 0

# Match conditions and subjects
subj = sort( unique( dat$ID ) )
for ( s in 1:length( subj ) ) {
  
  sel = dat$Condition == 'Drug' & 
    dat$ID == subj[s]
  
  dat$Dosage[ sel ] = 
    dosage$Dosage[ dosage$ID == subj[s] ]
}
dat$Dosage[ dat$Timepoints == 'T1_Pre_drug' ] = 0

# Sort data
dat = dat %>% 
  arrange( Subject, 
           desc( Condition ), 
           Timepoints, 
           Task, 
           desc(Response_type) )

# Convert self-report values to be numeric
dat = dat %>% 
  mutate( Self_report_on_high = as.numeric( Self_report_on_high ) )

dat$Self_report_on_high[ dat$Timepoints == 'T1_Pre_drug' ] = 0

# Add in visit order
visit = read.csv( file = 'Visit_order.csv', 
                  header = T, stringsAsFactors = F )
# Initialize values
dat$Visit = 1
# Loop over subjects
for ( s in 1:nrow( visit ) ) {
  
  sel = dat$ID == visit$ID[s]
  if ( !is.na( visit$Visit.2[s] ) ) {
    dat$Visit[ sel & dat$Condition == visit$Visit.2[s] ] = 2
  }
  
}

# Incorporate visit info into neural 
# data data frame as well
neurDat$Visit = 1
for ( s in 1:nrow( visit ) ) {
  
  sel = neurDat$ID == visit$ID[s]
  if ( !is.na( visit$Visit.2[s] ) ) {
    neurDat$Visit[ sel & neurDat$Condition == visit$Visit.2[s] ] = 2
  }
  
}

###
### 7) Save as .RData file
###

setwd( '..' )
save( rawDat, neurDat, data_key, 
      N, dat, 
      study_completion, 
      baseline_pulse,
      dosage, 
      demographics, 
      file = 'NBack_3_19_2018.RData' )

# Return to original directory
setwd( orig_dir )

