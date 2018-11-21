# Metrics for data checking
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you 
# have any questions or comments
# Last updated 2018-11-07

# Table of contents
# 1) Initial setup
# 2) Define functions
#   2.1) ss_tbl_quant
# 3) Table of summary statistics for most important variables
# 4) Histograms for variables

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
library( officer )
library( flextable )

# Read in data
all_files = dir()

fname = all_files[ grepl( 'Combined_data', all_files ) ]
raw_dat = read.csv( file = fname,
                    header = T,
                    stringsAsFactors = F
)

# Create variable tracking number of 
# non-zero and non-NA values for 
# primary CN-THCCOOH values
raw_dat$Check_nz_nna = NA
tmp = raw_dat %>% 
  group_by( id ) %>% 
  summarize(
    V = sum( is.na( ua_thclcmsms_crtadj ) | 
           ua_thclcmsms_crtadj == 0 )
  )
for ( s in 1:nrow( tmp ) ) {
  sel = raw_dat$id == tmp$id[s] & 
    raw_dat$visit_number == 1
  raw_dat$Check_nz_nna[sel] = tmp$V[s]
}

# Clean up workspace
rm( all_files, fname )

###
### 2) Define functions
###

# 2.1) 
ss_tbl_quant = function( vrb, digits = 1 ) {
  # Purpose:
  # Computes summary statistics for a quantitative 
  # variable.
  # Arguments:
  # vrb    - The variable to examine
  # digits - The number of digits to round to
  # Returns:
  # A vector with the desired statistics.
  
  # Extract variable
  raw_dat$vrb = as.numeric( raw_dat[[ vrb ]] )
  
  # Function to quickly round a variable to a 
  # desired number of digits
  qr = function( x ) {
    return( round( x, digits ) )
  }
  
  # Determine number of observations per subject
  counts = raw_dat %>% 
    group_by( id ) %>% 
    summarize(
      No = length( vrb ),
      N_na = sum( is.na( vrb ) )
    )
  
  # Compute summary statistics
  out = raw_dat %>% 
    summarize(
      N = length( unique( id ) ),
      N_with_NA = length( unique( id[ is.na( vrb ) ] ) ),
      Mean = qr( mean( vrb, na.rm = T ) ),
      SD = qr( sd( vrb, na.rm = T ) ),
      Min = qr( min( vrb, na.rm = T ) ),
      Max = qr( max( vrb, na.rm = T ) )
    )
  # Include summary statistics for number of 
  # observations and NA values per subject
  out = c( out,
           Avg_N_obs = qr( mean( counts$No ) ),
           Min_N_obs = min( counts$No ),
           Max_N_obs = max( counts$No ),
          Avg_NA_obs = qr( mean( counts$N_na ) ),
           Min_NA_obs = min( counts$N_na ),
           Max_NA_obs = max( counts$N_na )
  )
  
  return( out )
}

###
### 3) Table of summary statistics for most important variables
###

variables_to_check = 
  data.frame(
    label = c(
      'Sex',
      'Age',
      'BMI',
      'Race',
      'Years of MJ use',
      'Recency of MJ use',
      'Level of MJ use',
      'Level of alcohol use',
      'Nicotine users',
      'Dipstick test',
      'Age first used MJ',
      'CN-THCCOOH',
      'Days since baseline',
      'Visit number',
      'Flag for censored THC',
      'Flag for recency',
      'THC of NA or 0'
    ),
    variable = c(
      'sex_conf',
      'age_exact',
      'bmi',
      'race_revised',
      'years_mj_use',
      'tlfb_mj_8',
      'tlfb_mj_14b',
      'tlfb_etoh_13b',
      'nic_user',
      'thc_qual',
      'tlfb_mj_2',
      'ua_thclcmsms_crtadj',
      'days_since_baseline',
      'visit_number',
      'thc_500_flag',
      'thc_recency_flag',
      'Check_nz_nna'
    ), 
    stringsAsFactors = F
  )

# Compute summary statistics over specified 
# variables
check = t( sapply( variables_to_check$variable,
        ss_tbl_quant ) )
check = as.data.frame( check, stringsAsFactors = F )
check$Variable = variables_to_check$label
check$Column_name = variables_to_check$variable
ck = colnames( check )
ck = c( ck[13:14], ck[1:12] )
check = check[,ck]
rownames( check ) = 1:nrow( check )

# Create a flextable object
ft = check %>% 
  regulartable( col_keys = ck )

# Specify interpretable header labels
ft_header = data.frame(
  col_keys = ck,
  stringsAsFactors = F
)
ft_header$ColA = c(
  'Variable',
  'Raw column name',
  'N',
  'N w/ NA',
  'Mean',
  'SD',
  'Min',
  'Max',
  'Avg. # of obs.',
  'Min # of obs.',
  'Max # of obs.',
  'Avg. # of NA',
  'Min # of NA',
  'Max # of NA'
)
ft = ft %>% 
  set_header_df( mapping = ft_header, 
                 key = "col_keys" )

# Create initial borders
ft = ft %>% 
  theme_box()

# Create alternating gray lines
ft = ft %>% 
  theme_zebra( odd_header = "transparent", 
               even_header = "transparent" )
# Add upper/lower borders
ft = ft %>% 
  hline_bottom( border = 
                  fp_border( width = .75, 
                             color = "black"), 
                part = "header" )
ft = ft %>% 
  hline_top( border = fp_border( width = 2, 
                             color = "black"), 
         part = "header" )
ft = ft %>% 
  hline_bottom( border = fp_border( width = 2, 
                                 color = "black"), 
             part = "body" )
# Reduce font size
ft = ft %>% 
  fontsize( part = "header", size = 8)
ft = ft %>% 
  fontsize( part = "body", size = 8)
# Auto-adjust column widths
ft = ft %>% 
  autofit()

# View table before saving
# ft

# Save table
read_docx() %>% 
  body_add_flextable(ft) %>% 
  print(target = "Summary_statistics.docx")

###
### 4) Histograms for variables
###

# Read in subject IDs to be analyzed
s_ID = read.csv( file = "Subjects_to_analyze.csv",
                 header = T,
                 stringsAsFactors = F
)
s_ID = s_ID$V1

# Generate a PDF file
pdf( 'Plots_of_raw_data.pdf', width = 12 )

# Loop over variables
for ( v in 1:nrow( variables_to_check ) ) {
  
  # x11( width = 12 )
  
  # Initialize data frame with data to be plotted
  dtbp = raw_dat
  dtbp$V = raw_dat[[ variables_to_check$variable[v] ]]
  dtbp$Analyzed = !(dtbp$id %in% s_ID)
  
  # Arrange data
  dtbp = dtbp %>% 
    arrange( Analyzed, id, visit_number )
  # Track individual observations
  dtbp$Obs = 1:nrow( dtbp )
  
  # Set color scheme to group subjects into sets of 4
  clr = rep( 'black', nrow( dtbp ) )
  cur_id = unique( dtbp$id )
  all_clr = c(
    'black',
    colors()[179],
    colors()[205],
    colors()[231]
  )
  for( k in 1:4 ) {
    col_id = cur_id[ seq( k, length(cur_id), 4 ) ]
    clr[ dtbp$id %in% col_id ] = all_clr[k]
  }
  
  # Separate active from pilot groups
  pts = rep( 21, nrow( dtbp ) )
  pts[ dtbp$active ] = 24
  
  # Indicate which subjects will actually 
  # be analyzed
  sel = dtbp$id %in% s_ID & 
    pts == 21
  pts[sel] = 19
  sel = dtbp$id %in% s_ID & 
    pts == 24
  pts[sel] = 17
  
  # Generate plot
  plot( dtbp$V, 
        # Aesthetics
        pch = pts, col = clr, 
        # Axis labels and structure
        xlab = 'Subject',
        ylab = variables_to_check$label[v],
        bty = 'l', xaxt = 'n' )
  
  # Internal plotting dimensions
  ipr = par("usr")
  xl = ipr[1:2]
  yl = ipr[3:4]
  
  # Add vertical lines to group subjects into 
  # sets of 4
  dtbp$ID = NA
  inc = 1
  subj = unique( dtbp$id )
  for ( s in 1:length( subj ) ) {
    sel = dtbp$id == subj[s]
    dtbp$ID[sel] = inc
    inc = inc + 1
  }
  
  wm = sapply( dtbp$Obs, function (x) {
    cur = dtbp$Obs == x
    cur_id = dtbp$ID[ cur ]
    sel = dtbp$ID == cur_id
    wm = which( !is.na( dtbp$V[sel] ) )
    if ( length( wm ) > 0 ) {
      if ( x == dtbp$Obs[sel][ min( wm ) ] ) {
        return( T )
      } else {
        return( F )
      }
    } else {
      return( F )
    }
    
  } )
  sel = !is.na( dtbp$V ) & 
    (dtbp$ID %% 4) == 0 & 
    wm
  xa = dtbp$Obs[sel]
  vertLines( c( -2.5, xa+3.5 ), yl, 
             col = 'grey' )
  
  # x-axis
  axis( 1, xa, dtbp$ID[sel],
        tick = F, line = -.5 )
  
  # Add useful indicators
  if ( variables_to_check$variable[v] == 'tlfb_mj_8' ) {
    horizLines( 2, xl, lty = 2 )
  }
  if ( variables_to_check$variable[v] == 'tlfb_mj_14b' ) {
    horizLines( 4, xl, lty = 2 )
  }
  if (variables_to_check$variable[v] == 'Check_nz_nna') {
    horizLines( 6, xl, lty = 2 )
  }
  
  # Add legend
  legend(
    xl[1] + diff(xl)*.05,
    yl[2] + diff(yl)*.1,
    c(
      'Active',
      'Excluded',
      'Pilot',
      'Excluded'
    ),
    pch = c( 17, 24, 19, 21 ),
    bty = 'n',
    horiz = T,
    xpd = T
  )
  
}
dev.off()

setwd( pa_dir )

