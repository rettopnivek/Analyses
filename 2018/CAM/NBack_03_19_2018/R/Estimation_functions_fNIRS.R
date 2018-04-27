# Estimation functions (fNIRS data)
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directory if you 
# have any questions or comments
# Last updated 2018-04-12

# Table of contents
# 1) Initial setup
# 2) Functions


###
### 1) Initial setup
###

setwd( proj_dir )
setwd( 'Stan_scripts' )
stan_dir = getwd()

# Compile stan script for estimating SDT 
# model with mixed effects
if ( !exists('fNIRS_fit') ) {
  message( 'Estimation script' )
  comp_time = Sys.time()
  fNIRS_fit = stan_model(stanc_ret = "fNIRS_fit.stan", 
                         stanc_builder(".stan"))
  comp_time = Sys.time() - comp_time
  message(
    paste( 'Compilation time:', round(comp_time,2), 
           attributes(comp_time)$units ) )
  rm( comp_time )
}

setwd( proj_dir )

###
### 2) Functions
###

