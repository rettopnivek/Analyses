#--------------------------#
# Extract/Save drift rates # 
# Kevin Potter             #
# Updated 06/24/2017       #
#--------------------------#

# Get current working directory
cur_dir = getwd()

# Change directory to where model results are saved
setwd( 'Estimation_results' )

# Load in relevant model results
models = c(
  'DRM_SD_Orig_results.RData',
  'WP_SD_Comp_v1_results.RData',
  'WP_SD_Comp_v2_results.RData'
)

for ( i in 1:length( models ) ) {
  
  load( models[i] )
  
  # Extract drift rates
  chk = grep( 'DRM_SD_Orig', models[i] ) == 1 
  if ( length( chk ) > 0) {
    sel = grep( 'xi', colnames( raw_prm ) )
    xi_DRM_SD_Orig = exp( raw_prm[,sel[1:8]] )
    save( xi_DRM_SD_Orig, file = 'DRM_SD_Orig_xi.RData' )
  }
  
  # Extract drift rates
  chk = grep( 'WP_SD_Comp_v1', models[i] ) == 1 
  if ( length( chk ) > 0) {
    sel = grep( 'xi', colnames( raw_prm ) )
    xi_WP_SD_Comp_v1 = raw_prm[,sel[1:8]]
    save( xi_WP_SD_Comp_v1, file = 'WP_SD_Comp_v1_xi.RData' )
  }
  
  # Extract drift rates
  chk = grep( 'WP_SD_Comp_v2', models[i] ) == 1 
  if ( length( chk ) > 0) {
    sel = grep( 'xi', colnames( raw_prm ) )
    xi_WP_SD_Comp_v2 = raw_prm[,sel[1:8]]
    save( xi_WP_SD_Comp_v2, file = 'WP_SD_Comp_v2_xi.RData' )
  }
  
}

setwd( cur_dir )