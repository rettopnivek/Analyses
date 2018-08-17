# Script to extract data from fNIRS .csv files
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directory if you 
# have any questions or comments
# Last updated 2018-04-04

# Table of contents
# 1) Function to extract fNIRS data
# 2) Placebo conditions
# 3) Drug conditions

###
### 1) Function to extract fNIRS data
###

read_fNIRS_csv = function( fNIRS, 
                           row_index,
                           condition,
                           timepoint ) {
  # Purpose:
  # Given a character vector containing 
  # the rows with fNIRS data, convert to 
  # a data frame.
  # Arguments:
  # fNIRS     - A character vector with the 
  #             rows from the .csv file with 
  #             the fNIRS data
  # row_index - The specific rows for the fNIRS 
  #             data
  # condition - The type of condition (Placebo 
  #             or drug)
  # timepoint - The timepoint label
  # Returns:
  # A data frame.
  
  # Number of observations
  nr = length( row_index )
  
  # Initialize output data frame
  df = data.frame(
    ID = 'NULL',
    Condition = rep( condition, nr ),
    Timepoints = timepoint,
    Task = 'NBack_2',
    stringsAsFactors = FALSE )
  
  tmp = matrix( NA, nr, 25 )
  colnames( tmp ) = c(
    paste( 'Channel', 1:20, sep = '_' ),
    'R_DLPFC',
    'L_DLPFC',
    'MPFC',
    'R_VLPFC',
    'L_VLPFC'
  )
  df = cbind( df, tmp )
  
  obs = fNIRS[ row_index ]
  
  for ( i in 1:nr) {
    
    val = strsplit( obs[i], split = ',' )[[1]]
    val = val[ nchar(val) != 0 ]
    
    # Subject ID
    df$ID[i] = paste( 'FN_0', val[1], sep = '' )
    
    # Read in fNIRS data
    if ( length( val ) == 21 ) {
      val = c( val, rep( 'NA', 5 ) )
    }
    df[i,5:ncol(df)] = as.numeric( val[-1] )
    
  }
  
  return( df )
}

###
### 2) Placebo conditions
###

# Scan in text from .csv file for fNIRS placebo data
placebo_fNIRS = scan( file = 'Placebo_fNIRS_03.29.18.csv',
                      what = 'character',
                      sep = '\n' )

# Determine row with header information
header_loc = grep( 'SUBJECT ID', placebo_fNIRS )

# Determine rows where the data for the three different 
# timepoints start
start_pos = c(
  grep( 'PRE', placebo_fNIRS ),
  grep( 'DIFF', placebo_fNIRS ),
  grep( 'POST', placebo_fNIRS )
)

# Determine rows without any actual data
gaps = grep( ",,,,,,,,,,,,,,,,,,,,,,,,,,", placebo_fNIRS )

# Based on the previous vectors, the rows containing the 
# fNIRS data are:
# T1 (Pre drug):    3 -  56
# T2 (Post drug):  61 - 114
# T3 (Post drug): 176 - 221
row_indices = list(
  pre = 3:56,
  post = 61:114,
  post2 = 176:221
)

# T1 pre drug (Placebo)
neurDat = read_fNIRS_csv( placebo_fNIRS, 
                          row_indices[[1]], 
                          'Placebo', 
                          'T1_Pre_drug' )
# T2 post drug (Placebo)
neurDat = rbind( neurDat,
                 read_fNIRS_csv( placebo_fNIRS, 
                          row_indices[[2]], 
                          'Placebo', 
                          'T2_Post_drug' ) )

# T2 post drug (Placebo)
neurDat = rbind( neurDat,
                 read_fNIRS_csv( placebo_fNIRS, 
                                 row_indices[[3]], 
                                 'Placebo', 
                                 'T3_Post_drug' ) )

###
### 3) Drug conditions
###

# Scan in text from .csv file for fNIRS drug data
drug_fNIRS = scan( file = 'Drug_fNIRS_03.29.18.csv',
                   what = 'character',
                   sep = '\n' )

# Determine row with header information
header_loc = grep( 'SUBJECT ID', drug_fNIRS )

# Determine rows where the data for the three different 
# timepoints start
start_pos = c(
  grep( 'PRE', drug_fNIRS ),
  grep( 'DIFF', drug_fNIRS ),
  grep( 'POST', drug_fNIRS )
)

# Determine rows without any actual data
gaps = grep( ",,,,,,,,,,,,,,,,,,,,,,,,,,", drug_fNIRS )

# Based on the previous vectors, the rows containing the 
# fNIRS data are:
# T1 (Pre drug):    3 -  56
# T2 (Post drug):  60 - 113
# T3 (Post drug): 176 - 222
row_indices = list(
  pre = 3:56,
  post = 60:113,
  post2 = 176:222
)

# T1 pre drug (Placebo)
neurDat = rbind( neurDat, 
                 read_fNIRS_csv( drug_fNIRS, 
                          row_indices[[1]], 
                          'Drug', 
                          'T1_Pre_drug' ) )
# T2 post drug (Placebo)
neurDat = rbind( neurDat,
                 read_fNIRS_csv( drug_fNIRS, 
                                 row_indices[[2]], 
                                 'Drug', 
                                 'T2_Post_drug' ) )

# T2 post drug (Placebo)
neurDat = rbind( neurDat,
                 read_fNIRS_csv( drug_fNIRS, 
                                 row_indices[[3]], 
                                 'Drug', 
                                 'T3_Post_drug' ) )

# Save as new .csv file
write.csv( neurDat, file = 'Neural_data.csv' )

