# Table creation
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you 
# have any questions or comments
# Last updated 2018-12-10

# Table of contents
# 1) Initial setup
# 2) Create Table 1
# 3) Save table 1
# 4) Create table 2
# 5) Save table 2

###
### 1) Initial setup
###

source( 'S01_Folder_paths.R' )

# Indicate code segment to run
run_code = c(
  F, # Table 1
  F  # Table 2
)

# Load in useful packages

# Collection of useful functions
my_package_load( 'utilityf', github = T )

# Package for working with data frames
my_package_load( 'dplyr' )

# Package for Bayesian estimation
Sys.setenv(USE_CXX14 = 1) # For Rstan to work
my_package_load( 'brms' )

# Package for creating Office documents
my_package_load( 'officer' )

# Package for creating nice tables
my_package_load( 'flextable' )

# Load in packages for exponential decay model
source( 'S03_Exponential_decay_functions.R' )

# Load in data
setwd( dat_dir )
load( 'THC_decay.RData' )
# Load in data with imputted values
load( 'Imputted_data.RData' )

# If present, load in 
# best-fitting model 
# results
if ( 'Best_fitting_model.RData' %in% dir() ) {
  load( 'Best_fitting_model.RData' )
}

###
### 2) Create Table 1
###

# Variable to order final table results
ord = numeric( nrow( dtbf ) )
for ( i in 1:length( LETTERS ) ) {
  sel = dtbf$ID_Alph == LETTERS[i]
  ord[sel] = i
  sel = dtbf$ID_Alph == 
    paste( rep( LETTERS[i], 2 ), collapse = '' )
  ord[sel] = i + length( LETTERS )
  sel = dtbf$ID_Alph == 
    paste( rep( LETTERS[i], 3 ), collapse = '' )
  ord[sel] = i + length( LETTERS )*2
}
dtbf$Order = ord

# Clean up workspace
rm( ord, i, sel )

# Initialize table
tbl_1 = dtbf %>% 
  group_by( Order ) %>% 
  summarize(
    id = unique( ID ),
    id2 = unique( ID_Alph )
  )
colnames( tbl_1 ) = c( 'Order', 'ID', 'ID_Alph' )

# Define function to quickly extract variable 
# values in correct order
quick_add = function( vrb, df, cur_id ) {
  
  out = rep( NA, length( cur_id ) )
  for ( i in 1:length( cur_id ) ) {
    
    sel = df$ID %in% cur_id[i]
    val = unique( df[[ vrb ]][sel] )
    out[i] = val
  }
  
  return( out )
}

# Age, years
tbl_1$Age = round( quick_add( 'Age', all_dat, tbl_1$ID ) )
# Sex
tbl_1$Sex = 'M'
tbl_1$Sex[ quick_add( 'Sex', all_dat, tbl_1$ID ) == 1 ] = 'F'
# Race (Hispanic Ethnicity)
tbl_1$Race = quick_add( 'Race_orig', all_dat, tbl_1$ID )
# Weight, lb
tbl_1$Weight = round( quick_add( 'Weight', all_dat, tbl_1$ID ) )
# Height, in
tbl_1$Height = round( quick_add( 'Height', all_dat, tbl_1$ID ) )
# BMI
tbl_1$BMI = round( quick_add( 'BMI', all_dat, tbl_1$ID ), 1 )
# Days used in past 30
tbl_1$Days_used = quick_add( 'Days_used_past_30', all_dat, tbl_1$ID )
# Times used in past 30
tbl_1$Time_used = round( 
  quick_add( 'Times_used_past_30', all_dat, tbl_1$ID ) )
# Amount used in past 30, grams
tbl_1$Amount_used = round( 
  quick_add( 'Amount_used_past_30', all_dat, tbl_1$ID ), 1 )
# Days since last used
tbl_1$Recency = quick_add( 'Recency_of_MJ_use', all_dat, tbl_1$ID )
# Age first used, years
tbl_1$First_used = quick_add( 'Age_first_used_MJ', all_dat, tbl_1$ID )
# Years used
tbl_1$Years_used = round( 
  quick_add( 'Years_of_MJ_use', all_dat, tbl_1$ID ), 1 )
# CUDIT-R
tbl_1$CUDIT = quick_add( 'CUDIT.SS', dtbf, tbl_1$ID )

###
### 3) Save table 1
###

if ( run_code[1] ) {
  
  # For easier formatting, convert columns 
  # to character strings
  ftbl = apply( tbl_1, 2, as.character )
  ftbl = as.data.frame( ftbl, stringsAsFactors = F )
  # Remove subject IDs and order
  ftbl = ftbl %>% 
    select( -Order, -ID )
  
  # unlist( strsplit( tbl_1$Race, split = ' (H)', fixed = T ) )
  tmp = tbl_1 %>% 
    select(
      -Order,
      -ID,
      -ID_Alph,
      -Sex,
      -Race
    )
  
  M = function(x) mean( x, na.rm = T )
  SD = function(x) sd( x, na.rm = T )
  MD = function(x) median( x, na.rm = T )
  IQR_25 = function(x) quantile( x, prob = .25, na.rm = T )
  IQR_75 = function(x) quantile( x, prob = .75, na.rm = T )
  
  sm = tmp %>% 
    summarize_all(
      funs(
        M,
        SD,
        MD,
        IQR_25,
        IQR_75
      ) )
  sm = matrix( sm, ncol(tmp), 5, byrow = F )
  colnames( sm ) = c( 'Mean', 'SD', 'Median', 'IQR_25', 'IQR_75' )
  rownames( sm ) = colnames( tmp )
  # Clean up workspace
  rm( M, SD, MD, IQR_25, IQR_75,
      tmp )
  sm = as.data.frame( sm )
  for ( i in 1:ncol(sm) ) sm[[ i ]] = unlist( sm[[ i ]] )
  
  # Percentages
  tmp = unlist( strsplit( tbl_1$Race, split = ' (H)', fixed = T ) )
  per = list(
    Sex = round( 
      catProp( tbl_1$Sex, unique( tbl_1$Sex ) )[1]*100, 1 ),
    Race = round( 
      catProp( tmp, unique( tmp ) )*100, 1 )
  )
  # Add percentage for Hispanic
  per$Race = c( per$Race,
                H = round( 100*sum( grepl( '(H)', tbl_1$Race ) )/
                  nrow( tbl_1 ) )
  )
  
  # Add footer with details about 
  # summary statistics
  tbl_footer = data.frame(
    ID_Alph = c( 'Mean', 'SD', 'Median', 'IQR', '%', 
                 rep( ' ', length( per$Race ) - 1 ) ),
    stringsAsFactors = F
  )
  m = matrix( ' ', nrow( tbl_footer ), ncol( ftbl ) - 1 )
  colnames( m ) = colnames( ftbl )[-1]
  tbl_footer = cbind( tbl_footer, m )
  for ( i in 1:ncol(tbl_footer) ) {
    tbl_footer[[i]] = as.character( tbl_footer[[i]] )
  }
  # Add values to footer
  
  # Means
  tbl_footer$Age[1] = 
    as.character( round( sm['Age',1], 1 ) )
  tbl_footer$Weight[1] = 
    as.character( round( sm['Weight',1], 1 ) )
  tbl_footer$Height[1] = 
    as.character( round( sm['Height',1], 1 ) )
  tbl_footer$BMI[1] = 
    as.character( round( sm['BMI',1], 1 ) )
  tbl_footer$Days_used[1] = 
    as.character( round( sm['Days_used',1], 1 ) )
  tbl_footer$Recency[1] = 
    as.character( round( sm['Recency',1], 1 ) )
  tbl_footer$First_used[1] = 
    as.character( round( sm['First_used',1], 1 ) )
  tbl_footer$Years_used[1] = 
    as.character( round( sm['Years_used',1], 1 ) )
  tbl_footer$CUDIT[1] = 
    as.character( round( sm['CUDIT',1], 1 ) )
  
  # Standard deviations
  tbl_footer$Age[2] = 
    as.character( round( sm['Age',2], 1 ) )
  tbl_footer$Weight[2] = 
    as.character( round( sm['Weight',2], 1 ) )
  tbl_footer$Height[2] = 
    as.character( round( sm['Height',2], 1 ) )
  tbl_footer$BMI[2] = 
    as.character( round( sm['BMI',2], 1 ) )
  tbl_footer$Days_used[2] = 
    as.character( round( sm['Days_used',2], 1 ) )
  tbl_footer$Recency[2] = 
    as.character( round( sm['Recency',2], 1 ) )
  tbl_footer$First_used[2] = 
    as.character( round( sm['First_used',2], 1 ) )
  tbl_footer$Years_used[2] = 
    as.character( round( sm['Years_used',2], 1 ) )
  tbl_footer$CUDIT[2] = 
    as.character( round( sm['CUDIT',2], 1 ) )
  
  # Medians
  tbl_footer$Time_used[3] = 
    as.character( round( sm['Time_used',3], 1 ) )
  tbl_footer$Amount_used[3] = 
    as.character( round( sm['Amount_used',3], 1 ) )
  
  # Inter-quartile ranges
  tbl_footer$Time_used[4] = 
    paste(
      '[', round( sm['Time_used',4], 1 ), 
      ', ',
      round( sm['Time_used',5], 1 ), ']', sep = '' )
  tbl_footer$Amount_used[4] = 
    paste(
      '[', round( sm['Amount_used',4], 1 ), 
      ', ',
      round( sm['Amount_used',5], 1 ), ']', sep = '' )
  
  # Percentages
  tbl_footer$Sex[5] = 
    paste( per$Sex, '% ', names( per$Sex ), sep = '' )
  tbl_footer$Race[5:nrow( tbl_footer )] = 
    paste( per$Race, '% ', names( per$Race ), sep = '' )
  
  # Add footer
  ftbl = rbind( ftbl, tbl_footer )
  
  # Create header for table
  tbl_header = data.frame(
    col_keys = colnames( ftbl ),
    colD = ' ',
    colC = ' ',
    colB = ' ',
    colA = ' ',
    stringsAsFactors = F
  )
  # Subject ID
  tbl_header$colC[1] = 'Study'
  tbl_header$colB[1] = 'participant'
  # Age
  tbl_header$colC[2] = 'Age,'
  tbl_header$colB[2] = 'years'
  # Sex
  tbl_header$colC[3] = 'Sex'
  # Race
  tbl_header$colC[4] = 'Race'
  tbl_header$colB[4] = '(Hispanic'
  tbl_header$colA[4] = 'ethnicity)'
  # Weight
  tbl_header$colC[5] = 'Weight,'
  tbl_header$colB[5] = 'lb'
  # Height
  tbl_header$colC[6] = 'Height,'
  tbl_header$colB[6] = 'in'
  # BMI
  tbl_header$colC[7] = 'BMI'
  # Days used
  tbl_header$colC[8] = 'Days'
  tbl_header$colB[8] = 'used in'
  tbl_header$colA[8] = 'past 30'
  # Times used
  tbl_header$colC[9] = 'Times used'
  tbl_header$colB[9] = 'in past 30'
  # Amount used
  tbl_header$colC[10] = 'Amount used'
  tbl_header$colB[10] = 'in past 30,'
  tbl_header$colA[10] = 'grams'
  # Recency
  tbl_header$colD[11] = 'Cannabis use'
  tbl_header$colC[11] = 'Days'
  tbl_header$colB[11] = 'since last'
  tbl_header$colA[11] = 'use'
  # First used
  tbl_header$colC[12] = 'Age first'
  tbl_header$colB[12] = 'used,'
  tbl_header$colA[12] = 'years'
  # Years used
  tbl_header$colC[13] = 'Years'
  tbl_header$colB[13] = 'used'
  # CUDIT
  tbl_header$colC[14] = 'CUDIT-R'
  
  # Create flextable object
  ft = regulartable(
    ftbl, 
    col_keys = tbl_header$col_keys
  )
  # Add header with detailed column labels
  ft = ft %>% 
    set_header_df( mapping = tbl_header,
                   key = 'col_keys'
    )
  # Center header
  ft = ft %>% 
    align( align = "center", part = "header" )
  # Center body
  ft = ft %>% 
    align( align = "center", part = "body" )
  # Change font size and type
  ft = ft %>% 
    fontsize( size = 9, part = "all" )
  ft = ft %>% 
    font( fontname = "Arial", part = "all" )
  # Remove grid lines
  ft = ft %>% 
    border_remove()
  # Add borders for header and bottom
  big_border = fp_border(color = "black", width = 2)
  std_border = fp_border(color = "black", width = 1)
  ft = ft %>%
    hline_bottom( border = big_border, part = 'body' )
  ft = ft %>%
    hline_top( border = std_border, part = 'body' )
  ft = ft %>%
    hline_top( border = big_border, part = 'header' )
  # Add border separating footer
  pos = grep( 'Mean', ftbl$ID_Alph )
  ft = ft %>% 
    hline( i = pos - 1, border = std_border, part = 'body' )
  # Add grey background for alternating subjects
  alt_subj = unique( ftbl$ID_Alph )
  # Remove last 9 entries
  alt_subj = 
    alt_subj[ -( grep( 'Mean', alt_subj ):length( alt_subj ) ) ]
  alt_subj = alt_subj[ seq( 1, length( alt_subj ), 2 ) ] 
  for ( i in 1:length( alt_subj ) ) {
    sel = which( ftbl$ID_Alph %in% alt_subj[i] )
    for ( j in 1:length( sel ) ) {
      ft = ft %>% 
        bg( i = sel[j], bg = colors()[238] )
    }
  }
  # Color footer
  for ( i in pos:nrow( ftbl ) ) {
    ft = ft %>% 
      bg( i = i, bg = 'grey' )
  }
  ft = ft %>% 
    autofit()
  
  # Add title and notes
  str1 = 'Table 1. Demographics and Cannabis Use Characteristics'
  str2 = paste(
    'Note. AA, African American; AS, Asian; BMI, Body mass ',
    'index; CUDIT-R, Cannabis Use Disorder Identification ',
    'Test-Revised; H, Hispanic; in, inches; IQR, interquartile ',
    'range; lb, pounds; MOR, more than one race; ND, no data; ',
    'OTH, Other race; W, White',
    sep = '' )
  
  # Create Word document
  setwd( proj_dir )
  setwd( 'Documents' )
  read_docx() %>% 
    body_add_par(str1, style = "Normal") %>% 
    body_add_flextable(ft) %>% 
    body_add_par(str2, style = "Normal") %>% 
    print( 'Table_1.docx' )
  setwd( R_dir )
  
}

###
### 4) Create table 2
###

# Variable to order final table results
ord = numeric( nrow( dtbf ) )
for ( i in 1:length( LETTERS ) ) {
  sel = dtbf$ID_Alph == LETTERS[i]
  ord[sel] = i
  sel = dtbf$ID_Alph == 
    paste( rep( LETTERS[i], 2 ), collapse = '' )
  ord[sel] = i + length( LETTERS )
  sel = dtbf$ID_Alph == 
    paste( rep( LETTERS[i], 3 ), collapse = '' )
  ord[sel] = i + length( LETTERS )*2
}
dtbf$Order = ord

# Clean up workspace
rm( ord, i, sel )

# Read in estimated starting levels
setwd( proj_dir )
setwd( 'Documents' )
est_sl = read.csv( file = 'Starting_levels.csv',
                   stringsAsFactors = F
)
setwd( R_dir )

# Define function to compute day of 
# last positive or day of first 
# negative
f = function( id, type = 1 ) {
  
  # Initialize output
  out = NA
  
  # Select participant
  sel = dtbf$ID == unique( id )
  
  # Day of last positive
  if ( type == 1 ) {
    if ( any( sel ) ) {
      pft = dtbf$Positive_test_fed[sel] == 1 & 
        !is.na( dtbf$Positive_test_fed[sel] )
      d = dtbf$Time[sel]
      if ( any( pft ) ) {
        out = max( d[pft] )
      }
    }
  }
  
  # Day of first negative
  if ( type == 2 ) {
    if ( any( sel ) ) {
      nft = dtbf$Positive_test_fed[sel] == 0 & 
        !is.na( dtbf$Positive_test_fed[sel] )
      d = dtbf$Time[sel]
      if ( any( nft ) ) {
        out = min( d[nft] )
      }
    }
  }
  
  return( out )
}

# Initialize table
tbl_2 = dtbf %>% 
  group_by( Order ) %>% 
  summarize(
    id = unique( ID ),
    id2 = unique( ID_Alph ),
    DoA = max( Time ),
    SC = length( ID ),
    DoLP = f( ID, 1 ),
    DoFN = f( ID, 2 ),
    Cest = est_sl$Starting_level[ est_sl$ID == unique( ID ) ],
    Cmo = max( THCCOOH_orig ),
    Tmo = Time[ which.max( THCCOOH_orig ) ],
    Cld = THCCOOH_orig[ which.max( Time[ THCCOOH_orig > 0 ] ) ],
    Tld = Time[ which.max( Time[ THCCOOH_orig > 0 ] ) ]
  )
colnames( tbl_2 )[1:3] = c( 'Order', 'ID', 'ID_Alph' )

###
### 5) Save table 2
###

if ( run_code[2] ) {
  
  # Round values as necessary
  tbl_2$Cest = 
    round( tbl_2$Cest, 1 )
  tbl_2$Cmo = 
    round( tbl_2$Cmo, 1 )
  tbl_2$Cld = 
    round( tbl_2$Cld, 1 )
  
  # For easier formatting, convert columns 
  # to character strings
  ftbl = apply( tbl_2, 2, as.character )
  ftbl = as.data.frame( ftbl, stringsAsFactors = F )
  # Remove subject IDs and order
  ftbl = ftbl %>% 
    select( -Order, -ID )
  
  # Remove subjects with readings > 500
  chk = dtbf %>% 
    group_by( ID ) %>% 
    summarize( 
      Alph = unique( ID_Alph ), 
      ND = any( grepl( '1', all_dat$Data_issues[ all_dat$ID %in% ID ] ) ),
      Nsp = length( ID ),
      Cmax = max( THCCOOH_orig ),
      Clast = THCCOOH_orig[ which( Time == max( Time ) ) ],
      Tlast = max( Time )
    )
  ftbl$Cmo[ ftbl$ID_Alph %in% chk$Alph[ chk$ND ] ] = '---f'
  ftbl$Tmo[ ftbl$ID_Alph %in% chk$Alph[ chk$ND ] ] = '---f'
  
  # Add labels for missing data
  ftbl$DoLP[ is.na( ftbl$DoLP ) ] = '---d'
  ftbl$DoFN[ is.na( ftbl$DoFN ) ] = '---e'
  
  tmp = tbl_2 %>% 
    select(
      -Order,
      -ID,
      -ID_Alph
    )
  
  M = function(x) mean( x, na.rm = T )
  SD = function(x) sd( x, na.rm = T )
  MD = function(x) median( x, na.rm = T )
  IQR_25 = function(x) quantile( x, prob = .25, na.rm = T )
  IQR_75 = function(x) quantile( x, prob = .75, na.rm = T )
  
  sm = tmp %>% 
    summarize_all(
      funs(
        M,
        SD,
        MD,
        IQR_25,
        IQR_75
      ) )
  sm = matrix( sm, ncol(tmp), 5, byrow = F )
  colnames( sm ) = c( 'Mean', 'SD', 'Median', 'IQR_25', 'IQR_75' )
  rownames( sm ) = colnames( tmp )
  # Clean up workspace
  rm( M, SD, MD, IQR_25, IQR_75,
      tmp )
  sm = as.data.frame( sm )
  for ( i in 1:ncol(sm) ) sm[[ i ]] = round( unlist( sm[[ i ]] ), 1 )
  
  # Add footer with details about 
  # summary statistics
  tbl_footer = data.frame(
    ID_Alph = c( 'Mean', 'SD', 'Median', 'IQR' ),
    stringsAsFactors = F
  )
  m = matrix( ' ', nrow( tbl_footer ), ncol( ftbl ) - 1 )
  colnames( m ) = colnames( ftbl )[-1]
  tbl_footer = cbind( tbl_footer, m )
  for ( i in 1:ncol(tbl_footer) ) {
    tbl_footer[[i]] = as.character( tbl_footer[[i]] )
  }
  
  # Add values to footer
  for ( i in 1:nrow( sm ) ) {
    tbl_footer[[i+1]][1] = as.character( sm[i,1] )
    tbl_footer[[i+1]][2] = as.character( sm[i,2] )
    if ( i > 2 ) {
      tbl_footer[[i+1]][3] = as.character(sm[i,3])
      tbl_footer[[i+1]][4] = 
        paste(
          '[', round( sm[i,4], 1 ), 
          ', ',
          round( sm[i,5], 1 ), ']', sep = '' )
    }
  }
  
  # Add footer
  ftbl = rbind( ftbl, tbl_footer )
  
  # Add gap between dipstick and CN-THCCOOH columns
  ftbl = cbind(
    ftbl[,1:5],
    Gap = " ",
    ftbl[,6:10]
  )
  
  # Define upper labels for header
  lbls = list(
    c(
      'Immunoasssay RDDT',
      'Screen (LOQ for ',
      'THCCOOH= 50 ng/mL)'
    ),
    c(
      'Urine Creatinine-Adjusted THCCOOH ',
      'LC/MS/MS (ng/mg; LOQ ',
      'for THCCOOH= 5 ng/mL)'
    )
  )
  
  # Create header for table
  tbl_header = data.frame(
    col_keys = colnames( ftbl ),
    colF = ' ',
    colE = ' ',
    colD = ' ',
    colC = ' ',
    colB = ' ',
    colA = ' ',
    stringsAsFactors = F
  )
  # ID_Alph
  tbl_header$colC[1] = 'Study'
  tbl_header$colB[1] = 'Participant'
  # DoA
  tbl_header$colC[2] = 'Days of'
  tbl_header$colB[2] = 'Abstinencea'
  # SC
  tbl_header$colC[3] = 'Specimens'
  tbl_header$colB[3] = 'Collectedb'
  # DoLP
  tbl_header$colC[4] = 'Day of Last'
  tbl_header$colB[4] = 'Positive'
  tbl_header$colA[4] = 'Specimenc'
  # DoFN
  tbl_header$colC[5] = 'Day of First'
  tbl_header$colB[5] = 'Negative'
  tbl_header$colA[5] = 'Specimenc '
  # Gap
  # Cest
  tbl_header$colC[7] = 'C-estimated'
  tbl_header$colB[7] = 'at last'
  tbl_header$colA[7] = 'exposure'
  # Cest - upper label
  tbl_header$colF[7] = lbls[[2]][1]
  tbl_header$colE[7] = lbls[[2]][2] 
  tbl_header$colD[7] = lbls[[2]][3]
  # Cmo
  tbl_header$colC[8] = 'C-max observed'
  # Cmo - upper label
  tbl_header$colF[8] = lbls[[2]][1]
  tbl_header$colE[8] = lbls[[2]][2] 
  tbl_header$colD[8] = lbls[[2]][3]
  # Tmo
  tbl_header$colC[9] = 'T-max observedc'
  # Tmo - upper label
  tbl_header$colF[9] = lbls[[2]][1]
  tbl_header$colE[9] = lbls[[2]][2] 
  tbl_header$colD[9] = lbls[[2]][3]
  # Cld
  tbl_header$colC[10] = 'C-last observed'
  # Cld - upper label
  tbl_header$colF[10] = lbls[[2]][1]
  tbl_header$colE[10] = lbls[[2]][2] 
  tbl_header$colD[10] = lbls[[2]][3]
  # Tld
  tbl_header$colC[11] = 'T-last observedc'
  # Tld - upper label
  tbl_header$colF[11] = lbls[[2]][1]
  tbl_header$colE[11] = lbls[[2]][2] 
  tbl_header$colD[11] = lbls[[2]][3]
  
  # Create flextable object
  ft = regulartable(
    ftbl, 
    col_keys = tbl_header$col_keys
  )
  # Add header with detailed column labels
  ft = ft %>% 
    set_header_df( mapping = tbl_header,
                   key = 'col_keys'
    )
  # Merge upper labels
  ft = ft %>% 
    merge_h( part = "header" )
  # ft = ft %>% 
  #   merge_at( i = 4:5, j = 1:3, part = "header" )
  # Center header
  ft = ft %>% 
    align( align = "center", part = "header" )
  # Center body
  ft = ft %>% 
    align( align = "center", part = "body" )
  # Change font size and type
  ft = ft %>% 
    fontsize( size = 9, part = "all" )
  ft = ft %>% 
    font( fontname = "Arial", part = "all" )
  # Remove grid lines
  ft = ft %>% 
    border_remove()
  # Add borders for header and bottom
  big_border = fp_border(color = "black", width = 2)
  std_border = fp_border(color = "black", width = 1)
  ft = ft %>%
    hline_bottom( border = big_border, part = 'body' )
  ft = ft %>%
    hline_top( border = std_border, part = 'body' )
  ft = ft %>%
    hline_top( border = big_border, part = 'header' )
  # Add border separating footer
  pos = grep( 'Mean', ftbl$ID_Alph )
  ft = ft %>% 
    hline( i = pos - 1, border = std_border, part = 'body' )
  # Add grey background for alternating subjects
  alt_subj = unique( ftbl$ID_Alph )
  # Remove last 9 entries
  alt_subj = 
    alt_subj[ -( grep( 'Mean', alt_subj ):length( alt_subj ) ) ]
  alt_subj = alt_subj[ seq( 1, length( alt_subj ), 2 ) ] 
  for ( i in 1:length( alt_subj ) ) {
    sel = which( ftbl$ID_Alph %in% alt_subj[i] )
    for ( j in 1:length( sel ) ) {
      ft = ft %>% 
        bg( i = sel[j], bg = colors()[238] )
    }
  }
  # Color footer
  for ( i in pos:nrow( ftbl ) ) {
    ft = ft %>% 
      bg( i = i, bg = 'grey' )
  }
  # Auto-adjust cell widths and heights
  ft = ft %>% 
    autofit()
  # Extract dimensions
  dm = dim_pretty( ft )
  # Manually adjust column width for gap
  ft = ft %>% 
    width( j = 6, width = .2 )
  # Manually adjust column widths for last 5 columns
  ft = ft %>% 
    width( j = 7:11, width = sum(dm$widths[-(1:6)])/5 )
  # Manually adjust column widths for columns 4 and 5
  ft = ft %>% 
    width( j = 4:5, width = sum(dm$widths[4:5])/2 )
  
  # Add title and notes
  str1 = 'Table 2. Summary of Urine THCCOOH Concentrations'
  str2 = paste(
    'a Total number of days without cannabis exposure,',
    'tabulated from estimated day of last cannabis use.' )
  str3 = paste(
    'b Includes only analyzable specimens. Samples that',
    'were low volume, were beyond the linear range of the',
    'assay and were not diluted, as well as those from',
    'participants without verifiable abstinence are not',
    'included in table. Detailed information on reasons',
    'for missing data per participant are provided in',
    'Supplemental Appendix 2.'
  )
  str4 = 'c Tabulated from estimated day of last cannabis exposure.'
  str5 = 'd No urines screened positive for cannabinoids at 50 ng/mL.'
  str6 = 'e No urines screened negative for cannabinoids at 5 0ng/mL.'
  str7 = paste(
    'f Not reported because concentrations of initial',
    'specimens collected could not be quantitated due',
    'to low volume or concentrations beyond the linear',
    'range of the assay and were not diluted. See',
    'Supplemental Appendix 2 for visit-by-visit',
    'descriptives of known specimen concentrations.'
  )
  
  # Create Word document
  setwd( proj_dir )
  setwd( 'Documents' )
  read_docx() %>% 
    body_add_par(str1, style = "Normal") %>% 
    body_add_flextable(ft) %>% 
    body_add_par(str2, style = "Normal") %>% 
    body_add_par(str3, style = "Normal") %>% 
    body_add_par(str4, style = "Normal") %>% 
    body_add_par(str5, style = "Normal") %>% 
    body_add_par(str6, style = "Normal") %>% 
    body_add_par(str7, style = "Normal") %>% 
    print( 'Table_2.docx' )
  setwd( R_dir )
  
}

setwd( R_dir )


