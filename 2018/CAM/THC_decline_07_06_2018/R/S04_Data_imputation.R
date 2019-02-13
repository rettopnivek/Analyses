# Data imputation
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you 
# have any questions or comments
# Last updated 2018-12-10

# Table of contents
# 1) Initial setup
# 2) Data issues
# 3) Plot of THCCOOH by day
# 4) Create supplemental Appendix 2
# 5) Data imputation based on tobit regression

###
### 1) Initial setup
###

source( 'S01_Folder_paths.R' )

# Indicate code segments to run
run_code = c(
  F,
  F,
  F
)

# Collection of useful functions
my_package_load( 'utilityf', github = T )

# Package for working with data frames
my_package_load( 'dplyr' )

# Package for creating nice tables
my_package_load( 'flextable' )

# Package for manipulating Word documents
my_package_load( 'officer' )

# Load in packages for exponential decay model
source( 'S03_Exponential_decay_functions.R' )

# Load in data
setwd( dat_dir )
load( 'THC_decay.RData' )

# Compute summed CUDIT score
all_dat$CUDIT.SS = 
  rowSums( all_dat[,grep('CUDIT',colnames(all_dat))] )
for ( s in 1:length( unique( all_dat$ID ) ) ) {
  sel = all_dat$ID == unique( all_dat$ID )[s]
  val =  all_dat$CUDIT.SS[sel]
  val = val[ !is.na( val ) ]
  if ( length( val ) > 0 )
    all_dat$CUDIT.SS[sel] = val
}

# Select data with no issues
cd = all_dat %>% 
  filter( Data_issues == '0' )

# Clean up workspace
rm( s, sel, val )

###
### 2) Data issues
###

# Summarize data issues
di = all_dat %>% 
  group_by( ID ) %>% 
  summarize(
    # Total number of observations
    No = length( ID[ !is.na( THCCOOH ) ] ), 
    # Number of missing values
    Nna = sum( is.na( THCCOOH ) ),
    # Number of undetectable values
    Nz = sum( log( THCCOOH ) == -Inf, na.rm = T ),
    # Estimated recency of use over 2 days
    Recency = unique( Recency_of_MJ_use > 2 ),
    # Data issues
    DI = unique( Data_issues )[ 
      which.max( nchar( unique( Data_issues ) ) ) ]
  )
# Subjects identified with having issues with 
# abstinence or sets of observations that don't 
# have sufficient information to capture decay
di$MI = F
di$MI[ di$ID %in% c( '10016', '060_COMM', '076_COMM', '10063' ) ] = T
# Compute number of usable observations
di$Nu = 
  di$No - di$Nz
# Label issues

# Temporarily assign value to ensure sorted as first
di$Issues = 'A'

sel = di$Nu < 2
di$Issues[sel] = 'Only one observation within detectable levels'

sel = di$Nna > 5
di$Issues[sel] = 'Six or more missing observations'

sel = di$Recency
di$Issues[sel] = 'Did not use THC recently'

sel = di$ID %in% c( '060_COMM', '10063' )
di$Issues[sel] = 'Did not maintain abstinence'

sel = di$ID %in% c( '10016', '076_COMM' )
di$Issues[sel] = 'Cannot model exponential decay'

# Sort data by issues
di = di %>% 
  arrange( Issues )
di$Issues[ di$Issues == 'A' ] = 'No issues'

plot_yes = T
if ( plot_yes ) {
  setwd( fig_dir )
  pdf( 'Individual_fits.pdf',
       width = 12, height = 6 )
  for ( s in 1:length(di$ID) ) {
    
    sel = all_dat$ID == di$ID[s]
    tmp = all_dat[sel,]
    tmp$Time = tmp$Time + tmp$Recency_of_MJ_use
    if ( di$Issues[s] %in% 
         c( 'No issues', 
            'Did not maintain abstinence' ) ) {
      val = T
    } else {
      val = F
    }
    tryCatch( {
        plot_decay( tmp, new = F, estimate = val )
      },
      error = function(e) {
        layout( cbind( 1, 2 ) ); blankPlot(); blankPlot();
      }
    )
    
    string = paste(
      di$ID[s],
      ' - ', di$Issues[s], sep = ''
    )
    mtext( string, side = 3, line = -2, outer = T )
    
  }
  dev.off()
  setwd( R_dir )
}

tmp = di %>% 
  select( ID, Issues )
# Flip order of data frame
tmp = tmp[ rev( 1:nrow( tmp ) ), ]

setwd( dat_dir )
setwd( 'Original_files' )
write.table( tmp,
             row.names = F,
             quote = F,
             sep = ',',
             file = 'Excluded_subjects.csv'
)
setwd( R_dir )

# Clean up workspace
rm( tmp, sel, val, s, string, plot_yes )

# Select subjects with no issues
dtbf = cd %>% 
  filter( cd$ID %in% di$ID[ di$Issues == 'No issues' ] )
# Create alphabet-coded ID values for selected subjects
dtbf$ID_Alph = NA
all_dat$ID_Alph = NA
tmp = c( LETTERS, 
         paste( LETTERS, LETTERS, sep = '' ),
         paste( LETTERS, LETTERS, LETTERS, sep = '' )
)
for ( s in 1:length( tmp ) ) {
  sel = dtbf$ID == unique( dtbf$ID )[s]
  dtbf$ID_Alph[sel] = tmp[s]
  # Add to full data set as well
  sel = all_dat$ID == unique( dtbf$ID )[s]
  all_dat$ID_Alph[sel] = tmp[s]
}

###
### 3) Plot of THCCOOH by day
###

if ( run_code[1] ) {
  
  plot_decay( cd )
  legend( 4, 0,
          c( 'Observed', 'Mean', 'Censored-corrected predicted' ),
          fill = c( 'grey', 'black', 'blue' ),
          bty = 'n', cex = 1.25 )
  title( paste( round( 100 * sum( cd$THCCOOH == 0 )/nrow(cd) ),
                '% total obs. equal 0', sep = '' ) )
  
}

###
### 4) Create supplemental Appendix 2
###

if ( run_code[2] ) {
  
  sel = !is.na( all_dat$ID_Alph )
  
  tbl = all_dat %>% 
    filter( sel ) %>% 
    select( ID,
            ID_Alph,
            Visit_number,
            Days_from_baseline,
            Time,
            Positive_test_fed,
            THCCOOH,
            Urine_pH,
            Urine_specific_gravity
    )
  
  M = function(x) mean( x, na.rm = T )
  SD = function(x) sd( x, na.rm = T )
  MD = function(x) median( x, na.rm = T )
  IQR_25 = function(x) quantile( x, prob = .25, na.rm = T )
  IQR_75 = function(x) quantile( x, prob = .75, na.rm = T )
  
  sm = tbl %>% 
    select( -ID,
            -ID_Alph, 
            -Visit_number,
            -Days_from_baseline,
            -Time ) %>% 
    summarize_all(
      funs(
        M,
        SD,
        MD,
        IQR_25,
        IQR_75
      ) )
  sm = matrix( sm, 4, 5 )
  colnames( sm ) = c( 'Mean', 'SD', 'Median', 'IQR_25', 'IQR_75' )
  rownames( sm ) = c(
    'Positive_test_fed',
    'THCCOOH',
    'Urine_pH',
    'Urine_specific_gravity'
  )
  # Clean up workspace
  rm( M, SD, MD, IQR_25, IQR_75 )
  
  # Convert results for dipstick to string
  tmp = rep( 'NEG', nrow( tbl ) )
  tmp[ tbl$Positive_test_fed == 1 ] = 'POS'
  tmp[ is.na( tbl$Positive_test_fed ) ] = NA
  tbl$Positive_test_fed = tmp
  
  # Loop over each column, convert to string 
  for ( i in 1:ncol( tbl ) ) {
    tbl[[i ]] = as.character( tbl[[ i ]] )
  }
  
  # Determine type of missing data
  tmp = all_dat %>% 
    filter( !is.na( ID_Alph ) ) %>% 
    select( Data_issues ) %>% 
    unlist()
  # Initialize list tracking details about 
  # missing data
  missing_data = list()
  
  # ---DC = Lost to follow-up/Discontinued due to non-compliance
  
  # Set code for missing data
  md_code = '---DC'
  
  sel = grepl( '7', tmp )
  sel[ tbl$ID == '10019' & tbl$Visit_number %in% 7 ] = T
  sel[ tbl$ID == '10056' & tbl$Visit_number %in% 3:7 ] = T
  sel[ tbl$ID == '065_COMM' & tbl$Visit_number %in% 6:7 ] = T
  sel[ tbl$ID %in% c( '20014', '20024' ) ] = F
  
  # Label all measurement columns
  tbl$Days_from_baseline[sel] = md_code
  tbl$Time[sel] = md_code
  tbl$Positive_test_fed[sel] = md_code
  tbl$THCCOOH[sel] = md_code
  tbl$Urine_pH[sel] = md_code
  tbl$Urine_specific_gravity[sel] = md_code
  
  # Number of observations/subjects
  missing_data$DC = c(
    Observations = sum(sel),
    Subjects = length( unique( tbl$ID[sel] ) )
  )
  
  # ---LV = Sample leaked and no assay due to low volume
  
  # Set code for missing data
  md_code = '---LV'
  
  sel = tbl$ID == '10037' & 
    tbl$Visit_number == 1
  
  # Label all columns for CN-THCCOOH
  tbl$Positive_test_fed[sel] = md_code
  tbl$THCCOOH[sel] = md_code
  tbl$Urine_pH[sel] = md_code
  tbl$Urine_specific_gravity[sel] = md_code
  
  # Number of observations/subjects
  missing_data$LV = c(
    Observations = 1,
    Subjects = 1
  )
  
  sel = tbl$ID == '003_COMM' & 
    tbl$Visit_number == 1
  
  # Label all columns for CN-THCCOOH
  tbl$Positive_test_fed[sel] = md_code
  tbl$THCCOOH[sel] = md_code
  tbl$Urine_pH[sel] = md_code
  tbl$Urine_specific_gravity[sel] = md_code
  
  # Number of observations/subjects
  missing_data$LV[1] = missing_data$LV[1] + 1
  missing_data$LV[2] = missing_data$LV[2] + 1
  
  # ---ND = Exceeded linear range and not diluted
  
  # Set code for missing data
  md_code = '---ND'
  
  sel = grepl( '1', tmp )
  
  # Label all columns for CN-THCCOOH
  tbl$Positive_test_fed[sel] = md_code
  tbl$THCCOOH[sel] = md_code
  tbl$Urine_pH[sel] = md_code
  tbl$Urine_specific_gravity[sel] = md_code
  
  # Number of observations/subjects
  missing_data$ND = c( 
    Observations = sum(sel),
    Subjects = length( unique( tbl$ID_Alph[sel] ) )
  )
  
  # ---NR = Not recorded
  
  # Set code for missing data
  md_code = '---NR'
  
  # First visit excluded because not abstinent, but 
  # abstinent by second visit
  sel = rep( F, nrow( tbl ) )
  sel[ tbl$ID == '20014' & tbl$Visit_number %in% 1 ] = T
  sel[ tbl$ID == '20024' & tbl$Visit_number %in% 1 ] = T
  
  # Label all measurement columns
  tbl$Days_from_baseline[sel] = md_code
  tbl$Time[sel] = md_code
  tbl$Positive_test_fed[sel] = md_code
  tbl$THCCOOH[sel] = md_code
  tbl$Urine_pH[sel] = md_code
  tbl$Urine_specific_gravity[sel] = md_code
  
  # Number of observations/subjects
  # Number of observations/subjects
  missing_data$NR = c( 
    Observations = sum(sel),
    Subjects = length( unique( tbl$ID_Alph[sel] ) )
  )
  
  # Re-label missing values for dipstick test
  sel = is.na( tbl$Positive_test_fed )
  tbl$Positive_test_fed[ sel ] = '---NR'
  
  # d (Below limit of detectability, imputted value instead)
  
  for ( i in 1:nrow( tbl ) ) {
    if ( tbl$THCCOOH[i] == '0' & !is.na( tbl$THCCOOH[i] ) ) {
      sel = dtbf$ID_Alph == tbl$ID_Alph[i] & 
        dtbf$Visit_number == as.numeric( tbl$Visit_number[i] )
      tbl$THCCOOH[i] = paste(
        round( dtbf$THCCOOH[sel], 1 ),
        '+d+', sep = '' )
    }
  }
  
  # Add footer with details about 
  # summary statistics
  tbl_footer = data.frame(
    ID = ' ', 
    ID_Alph = c( 'Mean', 'SD', 'Median', 'IQR', '%' ),
    Visit_number = ' ',
    Days_from_baseline = ' ',
    Time = ' ',
    Positive_test_fed = ' ',
    THCCOOH = ' ',
    Urine_pH = ' ',
    Urine_specific_gravity = ' ',
    stringsAsFactors = F
  )
  # Dipstick test
  tbl_footer$Positive_test_fed[5] = 
    paste( round( sm[1,]$Mean, 3 )*100, '% POS', sep = '' )
  # CN-THCCOOH
  tbl_footer$THCCOOH[1] = 
    as.character( round( sm[2,]$Mean, 1 ) )
  tbl_footer$THCCOOH[2] = 
    as.character( round( sm[2,]$SD, 1 ) )
  tbl_footer$THCCOOH[3] = 
    as.character( round( sm[2,]$Median, 1 ) )
  tbl_footer$THCCOOH[4] = 
    paste( '[',
           round( sm[2,]$IQR_25, 1 ), 
           ' - ',
           round( sm[2,]$IQR_75, 1 ),
           ']', sep = '' )
  # Urine pH
  tbl_footer$Urine_pH[1] = 
    as.character( round( sm[3,]$Mean, 1 ) )
  tbl_footer$Urine_pH[2] = 
    as.character( round( sm[3,]$SD, 1 ) )
  tbl_footer$Urine_pH[3] = 
    as.character( round( sm[3,]$Median, 1 ) )
  tbl_footer$Urine_pH[4] = 
    paste( '[',
           round( sm[3,]$IQR_25, 1 ), 
           ' - ',
           round( sm[3,]$IQR_75, 1 ),
           ']', sep = '' )
  # Urine specific gravity
  tbl_footer$Urine_specific_gravity[1] = 
    as.character( round( sm[4,]$Mean, 1 ) )
  tbl_footer$Urine_specific_gravity[2] = 
    as.character( round( sm[4,]$SD, 1 ) )
  
  # Add footer
  tbl = rbind( tbl, tbl_footer )
  
  # Remove column for ID
  tbl = tbl %>% 
    select( -ID )
  
  # Create header for table
  tbl_header = data.frame(
    col_keys = colnames( tbl ),
    colD = ' ',
    colC = ' ',
    colB = ' ',
    colA = ' ',
    stringsAsFactors = F
  )
  # Subject ID
  tbl_header$colD[1] = 'Study'
  tbl_header$colC[1] = 'Participant'
  # Visit number
  tbl_header$colD[2] = 'Visit'
  tbl_header$colC[2] = 'Number'
  # Days from baseline
  tbl_header$colD[3] = 'Study'
  tbl_header$colC[3] = 'Daya'
  # Days from baseline
  tbl_header$colD[4] = 'Days'
  tbl_header$colC[4] = 'Since Last'
  tbl_header$colB[4] = 'Cannabis'
  tbl_header$colA[4] = 'Exposureb'
  # Dipstick test
  tbl_header$colD[5] = 'Urine THCCOOH'
  tbl_header$colC[5] = 'RDDT'
  tbl_header$colB[5] = 'Immunoassay'
  tbl_header$colA[5] = '(LOQ = 50 ng/mL)'
  # CN-THCCOOH
  tbl_header$colD[6] = 'Urine Creatinine-Adjusted'
  tbl_header$colC[6] = 'THCCOOH LC/MS/MS'
  tbl_header$colB[6] = '(ng/mg; LOQ for'
  tbl_header$colA[6] = 'THCCOOH= 5 ng/mL)c'
  # Urine pH
  tbl_header$colD[7] = 'Urine pH'
  # Urine specific gravity
  tbl_header$colD[8] = 'Urine'
  tbl_header$colC[8] = 'specific'
  tbl_header$colB[8] = 'gravity'
  
  # Create flextable object
  ft = regulartable(
    tbl, 
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
  pos = grep( 'Mean', tbl$ID_Alph )
  ft = ft %>% 
    hline( i = pos - 1, border = std_border, part = 'body' )
  # Add grey background for alternating subjects
  alt_subj = unique( tbl$ID_Alph )
  # Remove last 5 entries
  alt_subj = 
    alt_subj[ -( grep( 'Mean', alt_subj ):length( alt_subj ) ) ]
  alt_subj = alt_subj[ seq( 1, length( alt_subj ), 2 ) ] 
  for ( i in 1:length( alt_subj ) ) {
    sel = which( tbl$ID_Alph %in% alt_subj[i] )
    for ( j in 1:length( sel ) ) {
      ft = ft %>% 
        bg( i = sel[j], bg = colors()[238] )
    }
  }
  # Color footer
  for ( i in pos:nrow( tbl ) ) {
    ft = ft %>% 
      bg( i = i, bg = 'grey' )
  }
  ft = ft %>% 
    autofit()
  
  # Title and notes for table
  str1 = paste(
    'Supplemental Appendix 2. Individual Level Descriptives',
    'of Urine Specimens for Participants with Verified Abstinence' )
  
  str2 = paste( 'a Tabulated from number of days since',
                'baseline visit (Visit 1).' )
  str3 = paste( 'b Tabulated from estimated day of last', 
                'cannabis exposure.' )
  str4 = paste( 'c Includes only analyzable specimens. ',
                'Samples that were low volume, were beyond ',
                'the linear range of the assay and were not ',
                'diluted, as well as those from participants ',
                'without verifiable abstinence are not included ',
                'in table. Missing data is specified as follows: ',
                '---DC, Lost to follow-up or discontinued due to ',
                'study non-compliance (data points: n = ',
                missing_data$DC[1],
                '; participants: n = ',
                missing_data$DC[2], '); ',
                '---LV, Sample leaked in transit and quantitative ',
                'assay could not be performed due to low volume ',
                '(data points: n = ',
                missing_data$LV[1],
                '; participants: n = ',
                missing_data$LV[2], '); ---ND, THCCOOH ',
                'concentration exceeded the linear range ',
                'of the assay and sample dilutions were not ',
                'performed (data points: n = ',
                missing_data$ND[1],
                '; participants; n = ',
                missing_data$ND[2],
                '); ---NR, Not recorded (data points: n = ',
                missing_data$NR[1],
                '; participants: n = ',
                missing_data$NR[2], ').', sep = '' )
  
  str5 = '+d+ THCCOOH < 5ng/mL. Value imputed from tobit regression. '
  
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
    print( 'Supplementary_Appendix_2.docx' )
  setwd( R_dir )
  
  # Save details for missing data
  setwd( dat_dir )
  save( missing_data, file = 'Summary_of_missing_data.RData' )
  setwd( R_dir )
  
}

###
### 5) Data imputation based on tobit regression
###

if ( run_code[3] ) {
  
  # Extract subject IDs
  subj = unique( dtbf$ID )
  Ns = length( subj )
  
  # Initialize output
  all_est = matrix( NA, Ns, 3 )
  all_output = c()
  for ( s in 1:Ns ) all_output = c( all_output, list(NULL) )
  
  # Measure run time
  tick()
  # Initialize progress bar
  pb = txtProgressBar( min = 0, max = Ns, style = 3 )
  # Loop over subjects
  for ( s in 1:Ns ) {
    
    sel = dtbf$ID == subj[s]
    est = estimate_ed_lm( dtbf$Time[sel],
                          dtbf$THCCOOH[sel] )
    chk = sum( dtbf$THCCOOH[sel] > 0 )
    if ( chk == 2 ) est$est[3] = 0.0
    all_est[s,] = est$est
    all_output[[s]] = list( est )
    
    # Update the progress bar
    setTxtProgressBar(pb,s)
  }
  close(pb)
  tock()
  print( run_time )
  
  # Data imputation for zeros
  dtbf$THCCOOH_orig = dtbf$THCCOOH
  for ( s in 1:Ns ) {
    sel = dtbf$ID == subj[s]
    is_zero = dtbf$THCCOOH[sel] == 0
    if ( any( is_zero ) ) {
      dtbf$THCCOOH[sel][is_zero] = 
        all_output[[s]][[1]][['pred']][is_zero]
    }
  }
  
  # Convert to data frame
  all_est = data.frame(
    ID = subj,
    Start_point = all_est[,1],
    Log_start_point = log( all_est[,1] ), 
    Elimination_rate = all_est[,2],
    Residual_SD = all_est[,3],
    stringsAsFactors = F
  )
  
  ### Set up primary variables
  
  dtbf$log_THCCOOH = 
    log( dtbf$THCCOOH )
  dtbf$Start_point = 1
  dtbf$Elimination_rate = -dtbf$Time
  
  ### Set up predictors for start point
  
  dtbf$Sex.SP = dtbf$Sex
  dtbf$Race.SP = dtbf$Race
  
  # Standardize predictors
  dtbf$zBMI.SP = 
    my_standardize( dtbf$BMI )
  dtbf$zYears_of_MJ_use.SP = 
    my_standardize( dtbf$Years_of_MJ_use )
  dtbf$zLevel_of_MJ_use.SP = 
    my_standardize( dtbf$Level_of_MJ_use )
  
  ### Set up predictors for elimination rate
  
  dtbf$Sex.ER = dtbf$Sex * dtbf$Time
  dtbf$Race.ER = dtbf$Race * dtbf$Time
  
  # Standardize predictors
  dtbf$zBMI.ER = 
    dtbf$zBMI * dtbf$Time
  dtbf$zYears_of_MJ_use.ER = 
    dtbf$zYears_of_MJ_use * dtbf$Time
  dtbf$zLevel_of_MJ_use.ER = 
    dtbf$zLevel_of_MJ_use * dtbf$Time
  
  ### Save results
  
  setwd( dat_dir )
  save( dtbf, all_est, all_output, di, 
        file = 'Imputted_data.RData' )
  setwd( R_dir )
  
  # Clean up workspace
  rm( cd, est, pb, 
      chk, is_zero, Ns, run_code,
      s, sel, subj, tmp, run_time )
  
}

setwd( R_dir )

