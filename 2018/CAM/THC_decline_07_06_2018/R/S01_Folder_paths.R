# Folder paths
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you 
# have any questions or comments
# Last updated 2018-07-07

# Table of contents
# 1) Assignment of directory pathways to variables
# 2) Convenience function for package loading

###
### 1) Assignment of directory pathways to variables
###

# R scripts directory
R_dir = getwd()

# Project directory
setwd( '..' )
proj_dir = getwd()

# Data directory
setwd( 'Data' )
dat_dir = getwd()
setwd( proj_dir )

# Figure directory
setwd( 'Figures' )
fig_dir = getwd()
setwd( proj_dir )

###
### 2) Convenience function for package loading
###

my_package_load = function( pck, github = F, repos = NULL ) {
  # Purpose:
  # A convenience function to load packages if 
  # they are installed, and to install and then 
  # load them otherwise.
  # Arguments:
  # pck    - A character string giving the package name
  # github - Logical; if true, will if necessary install 
  #          packages from the author's personal Github 
  #          account
  # repos  - An optional character string to indicate 
  #          the repository from which to install the 
  #          package
  
  # Install packages from personal Github account
  if ( github ) {
    
    # Check if package is already installed
    if ( !( pck %in% rownames(installed.packages()) ) ) {
      
      # Check if devtools is installed
      if (!require("devtools")) {
        
        install.packages( "devtools" )
        
      }
      
      # Install package using devtools
      devtools::install_github(
        paste( "rettopnivek", pck, sep="/" ) )
      
    } else library( pck, character.only = T )
    
  } else {
    
    # Check if package is installed
    if ( !require( pck, character.only = T ) ) {
      
      install.packages( pck, repos = repos )
      
    }
    
  }
  
}

setwd( R_dir )
