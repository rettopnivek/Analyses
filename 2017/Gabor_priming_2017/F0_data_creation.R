#-----------------------------#
# Script to create R data set #
# Kevin Potter                #
# Updated 12/07/2017          #
#-----------------------------#

# Clear workspace
rm(list = ls())

# Save current working directory
orig_dir = getwd()

# Indicator for whether or not to save .RData files
saveFile = TRUE

# Index
# Lookup - 01:   Useful packages and functions
# Lookup - 01a:  extract_demographics
# Lookup - 01b:  extract_experiment_log
# Lookup - 01c:  extract_data
# Lookup - 01d:  replace_names
# Lookup - 01e:  create_condition_indices
# Lookup - 01f:  add_variables
# Lookup - 02:   Experiment 1
# Lookup - 03:   Experiment 2
# Lookup - 04:   Data key

###
### Useful packages and functions
###
# Lookup - 01

### Load in useful packages

# For geting github packages
# install.packages(devtools)
# library(devtools)

# Miscellanous functions for modeling/plotting
# install_github("rettopnivek/utilityf")
library(utilityf)

### Define useful functions

# Lookup - 01a
extract_demographics = function() {
  # Purpose: 
  # Extracts the demographic information from the 
  # relevant text files and creates a data frame 
  # with the results.
  # Returns: 
  # A data frame with the demographics information.
  
  # Extract file names
  all_files = dir()
  
  # Determine which files start with "Demographics"
  sel = grep( "Demographics", all_files )
  
  # Isolate files for demographics info
  demo_files = all_files[ sel ]
  
  # Determine sample size
  N = length( demo_files )
  
  # Initialize character matrix
  out = cbind( 
    ID = rep( "blank", N ), 
    Sex = rep( "blank", N ),
    Ethnicity = rep( "blank", N ),
    Race = rep( "blank", N ) )
  
  ID = numeric( N )
  
  for ( n in 1:N ) {
    
    # Subject ID
    tmp = strsplit( demo_files[n], '_' )[[1]][2]
    tmp = as.numeric( strsplit( tmp, '.txt' )[[1]] )
    ID[n] = tmp
    info = scan( file = demo_files[n], what = "character",
                 sep = "\n")
    out[n,2] = strsplit( info[1], split = ", " )[[1]][2]
    out[n,3] = strsplit( info[2], split = ", " )[[1]][2]
    out[n,4] = strsplit( info[3], split = ", " )[[1]][2]
  }
  
  # Convert to data-frame
  out = as.data.frame( out )
  # Ensure that ID variable is numeric
  out$ID = as.numeric( out$ID ); out$ID = ID
  
  # Sort data based on ID number
  out = out[ order( out$ID ), ]
  rownames( out ) = 1:nrow(out)
  
  return( out )
}

# Lookup - 01b
extract_experiment_log = function() {
  # Purpose: 
  # Extracts the experiment log and creates a 
  # data frame with the results, adding 
  # additional variables for subject ID numbers
  # and experiment duration.
  # Returns: 
  # A data frame with the experiment log.
  
  # Read in initial file
  experiment_log = read.csv( 'Experiment_log.csv' )
  
  # Add column with subject ID
  experiment_log$ID = NA
  
  # Loop through output files and extract ID numbers
  for ( i in 1:nrow( experiment_log ) ) {
    
    Subject_number = strsplit( as.character( 
      experiment_log$Output[i] ), split = '.csv' )[[1]][1]
    Subject_number = strsplit( Subject_number, split = '_' )[[1]][2]
    
    experiment_log$ID[i] = as.numeric( Subject_number )
  }
  
  # Remove backslashes from dates
  experiment_log$Start.time = 
    gsub( '\\\\', '_', as.character( experiment_log$Start.time ) )
  # Use forward slashes instead
  experiment_log$Start.time = 
    gsub( '_', '/', as.character( experiment_log$Start.time ) )
  
  # Remove backslashes from dates
  experiment_log$End.time = 
    gsub( '\\\\', '_', as.character( experiment_log$End.time ) )
  # Use forward slashes instead
  experiment_log$End.time = 
    gsub( '_', '/', as.character( experiment_log$End.time ) )
  
  # Function to extract starting and ending times 
  # from dates
  extract_time = function( string_time ) {
    out = numeric(3)
    tmp = strsplit( string_time, split = ':' )[[1]]
    out[1] = as.numeric( tmp[1] )
    tmp = strsplit( tmp[2], split = ' ' )[[1]]
    out[2] = as.numeric( tmp[1] );
    if ( tmp[2] == 'PM' ) out[3] = 0 else out[3] = 1
    
    return( out )
  }
  
  # Function to convert starting and ending times 
  # into duration in minutes
  convert_to_minutes = function( val ) {
    
    out = 0
    
    # Adjust for switch from morning to afternoon
    if ( val[3,2] != val[3,2] ) {
      val[1,] = c( 1, 2 );
    }
    
    h_change = diff( val[1,] ) # 0 if no change, one otherwise
    # Duration in minutes
    out = val[2,2] + 60 * h_change - val[2,1]
    names( out ) = ''
    
    return( out )
  }
  
  # Add column assessing duration of study
  experiment_log$Duration = NA
  
  for ( i in 1:nrow( experiment_log ) ) {
    
    string_time = character(2)
    string_time[1] = 
      strsplit( experiment_log$Start.time[i], split = '/' )[[1]][4]
    string_time[2] = 
      strsplit( experiment_log$End.time[i], split = '/' )[[1]][4]
    
    # Extract starting and ending times
    val = sapply( string_time, extract_time )
    
    experiment_log$Duration[i] = convert_to_minutes(val)
    
  }
  
  return( experiment_log )
}

# Lookup - 01c
extract_data = function() {
  # Purpose: 
  # Extracts data from .csv files and converts to 
  # data frame and integer for sample size.
  # Returns: 
  # A data frame with the raw data and an integer 
  # with the number of subjects.
  
  # Create empty variable to store subsequent data
  rawDat = c()
  
  # Extract filenames
  allFiles = dir()
  # Isolate subject files
  sel = grep( 'Subject', allFiles )
  datFiles = allFiles[ sel ]
  # Isolate .csv files
  sel = grep( '.csv', datFiles )
  datFiles = datFiles[ sel ]
  
  # Extract sample size
  N = length( datFiles )
  
  # Loop over subjects
  for (n in 1:N) {
    
    tmp = read.csv( file = datFiles[n], header = T )
    
    # If data is for experiment 1
    if ( !( 'Experiment' %in% colnames( tmp ) ) ) {
      sel = which( apply( tmp, 1, function(x) all( x[-1] == 0 ) ) )
      tmp = tmp[ -sel, ] # Trim out an empty row of zeros
    }
    
    # Collate data
    rawDat = rbind( rawDat, tmp )
    
  }
  
  # Clean up workspace
  rm( allFiles, sel, datFiles, n, tmp )
  
  return( list( N = N, rawDat = rawDat ) )
}

# Lookup - 01d
replace_names = function( rawDat, original_names, new_names ) {
  # Purpose: 
  # Generates a character vector with the new variable names 
  # for a data frame
  # Arguments: 
  # rawDat         - The data frame whose variable names should 
  #                  be adjusted
  # original_names - The original variable names to be replaced
  # new_names      - The new variable names replacing the originals
  #                  given above
  # Returns: 
  # A character vector with the adjusted variable names.
  
  if ( length( original_names ) != length( new_names ) ) 
    stop( 'Vectors for original and new names must match in length',
          call. = FALSE )
  
  variable_names = colnames( rawDat )
  # Loop over original names, replacing them with 
  # new names
  for ( i in 1:length( original_names ) ) {
    sel = variable_names == original_names[i]
    variable_names[sel] = new_names[i]
  }
  
  return( variable_names )
}

# Lookup - 01e
create_condition_indices = function( rawDat ) {
  # Purpose: 
  # Creates a set of indices giving the different 
  # level combinations over conditions of interest.
  # Arguments: 
  # rawDat - The updated data frame of observations 
  #          (i.e., after applying the 'add_variables'
  #           function).
  # Returns: 
  # A list consisting of two data frames:
  # 1) The set of indices and the associated condition levels;
  # 2) The original data frame, updated to include the indices.
  
  # Isolate trials for main study
  sel = rawDat$Block_type == 2
  
  trial_counter = rep( 1, sum( sel ) )
  
  # Full index of all conditions
  variables_of_interest = list(
    rawDat$Target_rotation_label[sel], 
    rawDat$Prime_rotation_label[sel], 
    rawDat$Onscreen_choices_label[sel],  
    rawDat$Correct_answer_label[sel], 
    rawDat$Prime_duration[sel], 
    rawDat$Prime_type_label[sel]
  )
  
  Full_index = aggregate(
    trial_counter, variables_of_interest, sum )
  colnames( Full_index ) = c(
    'Target',
    'Prime',
    'Onscreen_choice',
    'Correct',
    'Prime_duration',
    'Prime_type',
    'Trials'
  )
  
  Full_index$Full_index = 1:nrow( Full_index )
  
  Full_index$ChxDurxPri = NA
  for ( u in unique( Full_index$Onscreen_choice ) ) {
    sel = Full_index$Onscreen_choice == u
    Full_index$ChxDurxPri[sel] = 1:sum( sel )
  }
  
  Full_index$DurxPri = NA
  for ( u in unique( Full_index$Target ) ) {
    sel = Full_index$Target == u
    Full_index$DurxPri[sel] = 1:sum( sel )
  }
  
  rawDat$Full_index = NA
  rawDat$ChxDurxPri = NA
  rawDat$DurxPri = NA
  
  for ( i in 1:nrow( Full_index ) ) {
    
    sel = rawDat$Target_rotation_label == Full_index$Target[i] & 
      rawDat$Correct_answer_label == Full_index$Correct[i] & 
      rawDat$Prime_duration == Full_index$Prime_duration[i] & 
      rawDat$Prime_type_label == Full_index$Prime_type[i] & 
      rawDat$Block_type == 2
    
    rawDat$Full_index[sel] = Full_index$Full_index[i]
    rawDat$ChxDurxPri[sel] = Full_index$ChxDurxPri[i]
    rawDat$DurxPri[sel] = Full_index$DurxPri[i]
  }
  
  return( list( Full_index, rawDat ) )
}

# Lookup - 01f
add_variables = function( rawDat, Exp ) {
  # Purpose: 
  # Adds additional variables to the data frame of raw data.
  # Arguments: 
  # rawDat - The data frame of raw data generated from extract_data
  # Exp    - The experiment version, where...
  #            1 = The first experiment (A) with 5 prime durations
  #            2 = An additional pilot study with 3 prime durations
  #            3 = An additional pilot study that used same-different
  #                responding instead (Same subjects as 2)
  # Returns: 
  # A new data frame with additional variables.
  
  ### Adjustments to data from experiment 1 ###
  
  if ( Exp == 1 ) {
    
    # Adjust variable names
    original_names = c( 
      'PrimeType', 
      'PrimeDuration',
      'FoilContrast',
      'TargetContrast',
      'PP',
      'BlockType',
      'Target'
    )
    new_names = c( 
      'Prime_type', 
      'Prime_duration',
      'Foil_contrast',
      'Target_contrast',
      'Max_illuminance',
      'Block_type',
      'Target_rotation'
    )
    colnames( rawDat ) = replace_names( rawDat,
                                        original_names,
                                        new_names )
    
    # Adjust values of prime type to match 
    # those from experiment 2 and 3
    new_prime_type = rep( 0, nrow( rawDat ) )
    # Neither primed
    new_prime_type[ rawDat$Prime_type == 0 ] = 1
    # Foil primed
    new_prime_type[ rawDat$Prime_type == 2 ] = 2
    # Target primed
    new_prime_type[ rawDat$Prime_type == 1 ] = 3
    rawDat$Prime_type = new_prime_type
    
    # Add in missing variables that are included 
    # in experiments 2 and 3
    
    # Correct answer is based on rotation of target stripes
    rawDat$Correct_answer = rawDat$Target_rotation
    
    # Include ratio of target to foil contrasts
    rawDat$Target_to_foil_ratio = 
      rawDat$Target_contrast / rawDat$Foil_contrast
    
    # Include timing of the stimuli shown in the study
    frame_rate = 1000 / 60 # Duration in ms of a frame at 60 Hz
    rawDat$Fixation_time =  19 * frame_rate / 1000
    rawDat$Target_time = 5 * frame_rate / 1000
    rawDat$Mask_time = 0
    rawDat$Feedback_time = 24 * frame_rate / 1000
    # Compute more accurate values of the timing 
    # for the placeholder and prime stimuli
    nFrames = round( rawDat$Prime_duration / frame_rate )
    rawDat$Placeholder_time = ( 24 - nFrames ) * frame_rate / 1000
    rawDat$Placeholder_time[ rawDat$Placeholder_time < 0 ] = 0
    rawDat$Prime_time = nFrames * frame_rate / 1000
    
    # Indicate rotation of prime stimulus
    rawDat$Prime_rotation = rawDat$Target_rotation
    rawDat$Prime_rotation[ rawDat$Prime_type == 1 ] = 2
    sel = rawDat$Prime_type == 2
    rawDat$Prime_rotation[sel] = 1 - rawDat$Target_rotation[sel]
    
    # Provide information on choice display
    rawDat$Choice_rotation = 2
    rawDat$Choice_display = 3
    
    # Add in variable tracking blocks
    rawDat$Block = NA
    # Loop over subjects
    for ( s in unique( rawDat$Subject ) ) {
      
      sel = rawDat$Subject == s
      
      tmp = aggregate( rep( 1, sum(sel) ), 
                       list( rawDat$Block_type[sel] ), sum )
      block_inc = c(
        rep( 2, tmp$x[1] ),
        rep( 3, tmp$x[2] ),
        rep( 4, tmp$x[3] ),
        rep( 1:6 + 5, each = 75 ) )
      rawDat$Block[sel] = block_inc
      
    }
    
    # Add in meaningful labels for the numerically coded 
    # variables
    
    rawDat$Choice_label = 'Left'
    rawDat$Choice_label[ rawDat$Choice == 1 ] = 'Right'
    rawDat$Choice_label[ rawDat$Choice == 2 ] = 'Timeout'
    
    rawDat$Target_rotation_label = 'Left'
    sel = rawDat$Target_rotation == 1
    rawDat$Target_rotation_label[sel] = 'Right'
    
    rawDat$Correct_answer_label = rawDat$Target_rotation_label
    
    rawDat$Prime_rotation_label = 'Left'
    sel = rawDat$Prime_rotation == 1
    rawDat$Prime_rotation_label[sel] = 'Right'
    sel = rawDat$Prime_rotation == 2
    rawDat$Prime_rotation_label[sel] = 'Vertical'
    
    rawDat$Prime_type_label = 'Neither primed'
    sel = rawDat$Prime_type == 2
    rawDat$Prime_type_label[sel] = 'Foil primed'
    sel = rawDat$Prime_type == 3
    rawDat$Prime_type_label[sel] = 'Target primed'
    
    rawDat$Choice_rotation_label = 'None'
    
    rawDat$Onscreen_choices_label = 'None'
    
    rawDat$Block_type_label = 'Main study'
    sel = rawDat$Block_type == -1
    rawDat$Block_type_label[sel] = 'Exclusion practice'
    sel = rawDat$Block_type == 0
    rawDat$Block_type_label[sel] = 'Final training'
    sel = rawDat$Block_type == 1
    rawDat$Block_type_label[sel] = 'Staircase calibration'
    
  }
  
  ### Adjustments to data from experiment 2 ###
  
  if ( Exp == 2 ) {
    
    # Adjust variable names
    original_names = c( 
      'Target_to_foil_contrast',
      'Onscreen_choices_labels'
    )
    new_names = c( 
      'Target_to_foil_ratio',
      'Onscreen_choices_label'
    )
    colnames( rawDat ) = replace_names( rawDat,
                                        original_names,
                                        new_names )
    
    # Correct labels for choice rotation
    rawDat$Choice_rotation_label = 'None'
    
  }
  
  ### Adjustments to data from experiment 3 ###
  
  if ( Exp == 3 ) {
    
    # Adjust variable names
    original_names = c( 
      'Target_to_foil_contrast',
      'Onscreen_choices_labels'
    )
    new_names = c( 
      'Target_to_foil_ratio',
      'Onscreen_choices_label'
    )
    colnames( rawDat ) = replace_names( rawDat,
                                        original_names,
                                        new_names )
  }
  
  ### Adjustments to all experiments ###
  
  # Add variable indicating trial number
  trial_length = aggregate( rep( 1, nrow( rawDat ) ), 
                            list( rawDat$Subject ), length )
  colnames( trial_length ) = c( 'S', 'N' )
  rawDat$Trial = NA
  for ( s in unique( rawDat$Subject ) ) {
    sel1 = rawDat$Subject == s
    sel2 = trial_length$S == s
    rawDat$Trial[sel1] = 1:trial_length$N[sel2]
  }
  
  # Add labels for accuracy
  rawDat$Accuracy_label = 'Correct'
  sel = rawDat$Accuracy == 0
  rawDat$Accuracy_label[sel] = 'Incorrect'
  
  # Indicate experiment version
  rawDat$Experiment = NULL
  rawDat$Experiment = Exp
  
  # Store subject ID numbers
  rawDat$ID = rawDat$Subject
  # Re-label subject IDs to be in incremental, ascending order
  rawDat$Subject = createIncrement( rawDat$Subject )
  
  return( rawDat )
}

###
### Experiment 1 
###
# Lookup - 02

check_match_between_exp = list(
  Exp_1_variable_names = NULL,
  Exp_2_variable_names = NULL,
  Exp_3_variable_names = NULL
)

# Change working directory to folder location
setwd( 'Data/Raw_data' )

# Location of raw data files for 1st experiment
setwd( 'Experiment_1' )

# Extract demographic info
demographics = extract_demographics()

# Extract experiment log
experiment_log = extract_experiment_log()

# Read in data
tmp = extract_data()
N = tmp$N
rawDat = tmp$rawDat
# Add additional variables to rawDat file
rawDat = add_variables( rawDat, Exp = 1 )
rm( tmp )

# Add additional indices
tmp = create_condition_indices( rawDat )
Full_index = tmp[[1]]
rawDat = tmp[[2]]
rm( tmp )

# Indicate which subjects finished the study
demographics$Finished_study = 'No'
demographics$Finished_study[ 
  demographics$ID %in% sort( unique( rawDat$Subject ) ) ] = 'Yes'

# Indicate which subjects were correctly logged
demographics$In_exp_log = 'No'
demographics$In_exp_log[ 
  demographics$ID %in% experiment_log$ID ] = 'Yes'

check_match_between_exp$Exp_1_variable_names = 
  colnames( rawDat )

# Save data set
setwd( orig_dir ); setwd('Data')
if (saveFile) save (rawDat, N, demographics, experiment_log, 
                    file = 'Gabor_Exp_1_11-17.RData' )
setwd( orig_dir )

###
### Experiment 2
###
# Lookup - 03

# Change working directory to folder location
setwd( 'Data/Raw_data' )

# Location of raw data files for 1st experiment
setwd( 'Experiment_2' )

# Extract demographic info
demographics = extract_demographics()

# Extract experiment log
experiment_log = extract_experiment_log()

# Read in data
tmp = extract_data()
N = tmp$N
rawDat = tmp$rawDat
# Add additional variables to rawDat file
rawDat = add_variables( rawDat, Exp = 2 )
rm( tmp )

# Add additional indices
tmp = create_condition_indices( rawDat )
Full_index = tmp[[1]]
rawDat = tmp[[2]]
rm( tmp )

check_match_between_exp$Exp_2_variable_names = 
  colnames( rawDat )

# Save data set
setwd( orig_dir ); setwd('Data')
if (saveFile) save (rawDat, N, demographics, experiment_log, 
                    file = 'Gabor_Exp_2_11-17.RData' )
setwd( orig_dir )

###
### Data key
###
# Lookup - 04

# Check that all experiments have matching variable names
all_checks = c(
  Exp_1_and_2 = FALSE,
  Exp_2_and_1 = FALSE
)
all_checks[1] = 
  all( check_match_between_exp$Exp_1_variable_names %in% 
  check_match_between_exp$Exp_2_variable_names )
all_checks[2] = 
  all( check_match_between_exp$Exp_2_variable_names %in% 
       check_match_between_exp$Exp_1_variable_names )
print( all_checks )

# Data key
b('
Variable name          | Description
-----------------------|--------------------------------------
Subject                | The subject index
Trial                  | The trial index for each subject
Block                  | The block index for each subject
Block_type             | The type of block
Foil_contrast          | The Michelson contrast value for the
                       | foil stripes in the overlaid grid
Target_contrast        | The Michelson contrast value for the
                       | target stripes in the overlaid grid
Target_to_foil_ratio   | The ratio of the target and foil contrasts
Max_illuminance        | The max illuminance value in the overlaid grid
Angle                  | The absolute value for the angle of rotation 
                       | for the overlaid grid
Fixation_time          | The duration (in seconds) for the fixation dot
Placeholder_time       | The duration (in seconds) for the horizontal 
                       | stripes making up the placeholder stimulus
Prime_time             | The actual duration (in seconds) of the 
                       | prime stimulus
Target_time            | The duration (in seconds) of the overlaid grid
Mask_time              | The duration (in seconds) of the post-target 
                       | mask
Feedback_time          | The duration (in seconds) of the accuracy 
                       | feedback following a response
Prime_duration         | The duration of the prime stimulus, rounded 
                       | up to the nearest millisecond
Prime_rotation         | The overall angle of the stripes used for the 
                       | prime
Prime_type             | The type of prime, which could match the angle 
                       | of either the target or foil stripes in the grid, 
                       | or simply be vertical
Target_rotation        | The overall rotation of the darker stripes in 
                       | the overlaid grid (either Left or Right)
Correct_answer         | The correct answer for each trial
Choice_rotation        | The overall rotation for the onscreen choice 
                       | options (if shown)
Choice_display         | The type of choice options shown onscreen
RT                     | The response time (in seconds) - a trial would 
                       | time out after 5 seconds
Choice                 | The choice made by a subject using the "j" and 
                       | "k" keys
Accuracy               | Indicates whether the subject made the correct 
                       | choice
Choice_label           | Meaningful label for the choice values
Prime_rotation_label   | Meaningful label for the overall angle of the 
                       | prime stripes
Target_rotation_label  | Meaningful angle for the overall angle of the 
                       | target stripes in the grid
Choice_rotation_label  | Meaningful label for the (potential) angle of 
                       | the onscreen choice options
Prime_type_label       | Meaningful label for the type of prime
Correct_answer_label   | Meaningful label for the type of correct respnose
Onscreen_choices_label | Meaningful label for the onscreen choice options
Block_type_label       | Meaningful label for the type of block
Accuracy_label         | Meaningful label for the accuracy values
Experiment             | Indicates to which experiment type a set of data 
                       | belongs, where...
                       | 1) A study with 5 prime durations, 3 prime 
                       |    prime types (foil, target, and neither primed)
                       |    and a forced choice between whether the 
                       |    darker stripes were rotated to the left or right
                       | 2) A study with 3 prime durations, 2 prime 
                       |    prime types (foil and target primed)
                       |    and a forced choice between whether the 
                       |    darker stripes were rotated to the left or right
                       | 3) A study with 2 prime durations, 2 prime types
                       |    (foil and target primed), and a forced choice 
                       |    between whether a set of onscreen stripes were 
                       |    either the "same" or "different" compared to the 
                       |    darker stripes
ID                     | The original identifying number for a subject
Full_index             | A single index giving the differing combinations 
                       | of levels for the rotation of the target stripes,
                       | the type of correct answer, prime duration, and 
                       | prime type.
ChxDurxPri             | A single index giving the differing combinations 
                       | of levels for the type of correct answer, prime 
                       | duration, and prime type.
DurxPr   i             | A single index giving the differing combinations 
                       | of levels for the type of prime duration and prime 
                       | type.
')

# Return to original directory
setwd( orig_dir )