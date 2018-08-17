# Manuscript figures for N-back tasks
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you 
# have any questions or comments
# Last updated 2018-05-25

# Table of contents
# 1) Initial setup
#   1.1) lmer_standardized_coef
# 2) Diagnostic functions
#   2.1) convergenceExtract
#   2.2) findDec
#   2.3) quick_CI_plot
#   2.4) extract_model_results
# 3) Figures for convergence diagnostics
# 4) Model comparisons
# 5) Posterior retrodictive checks (Model 3)
#   5.1) quick_pred
# 6) Marginal posteriors (d')
# 7) Marginal posteriors (response bias)
# 8) Correlation between Self-report and THC administration
# 9) Extraction of baseline vitals
# 10) DEQ correlations
# 11) Table of results for SDT model

# Save current working directory
orig_dir = getwd()

# Navigate to project directory
setwd( '..' )
proj_dir = getwd()

# Indicate which code segments to run
run_code = c(
  F, # Convergence diagnostics
  F, # Model comparisons
  F, # Post. pred. check (Model 3)
  F, # Marginal posteriors (d')
  F, # Marginal posteriors (response bias)
  F, # Correlation between Self-report and THC administration
  F, # Extraction of baseline vitals
  F, # DEQ correlations
  F  # Table of results for SDT model
)

# Indicate whether plots should be saved as PDFs
savePlot = F

# Load in useful packages

# Collection of useful functions
# devtools::install_github("rettopnivek/utilityf")
library( utilityf )

# Package for easier manipulation of data frames
# install.packages( 'dplyr' )
library( dplyr )

# Package for estimating mixed effects models
# install.packages( 'lme4' )

# 1.1)
lmer_standardized_coef = function( fit, intercept = 1 ) {
  # Purpose:
  # Computes the standardized fixed effects for a 
  # lmer fit object.
  # Arguments:
  # fit       - A lmer fit object
  # intercept - The location(s) of the intercept(s)
  # Reference:
  # https://stats.stackexchange.com/questions/123366/
  #   lmer-standardized-regression-coefficients
  # Returns:
  # The standardized fixed effects.
  
  # Extract fixed effect estimates
  b = lme4::fixef(fit)[-intercept]
  
  # Compute standard deviations for predictors
  X = matrix( lme4::getME( mod, "X" )[,-intercept], 
              ncol = length(b) )
  sd_x = apply( X, 2, sd )
  # Compute standard deviations for dependent variable
  sd_y = sd( lme4::getME( fit, "y" ) )
  
  # Standardize coefficients
  out = b * sd_x/sd_y
  
  return( out )
}

# Package for estimating mixed effects models (Bayesian)
# install.packages( 'rstanarm' )
library( rstanarm )
# For parallel processing
options(mc.cores = parallel::detectCores())

# Package for Bayes factors for Stan models
# install.packages( 'bridgesampling' )
library( bridgesampling )

# Load in data
setwd( 'Data' )
load( "NBack_3_19_2018.RData" )
setwd( proj_dir )

# Load in useful functions
setwd( 'R' )
source( 'S02_Useful_functions.R' )
setwd( proj_dir )

# Exclude combined data
cd = dat %>% 
  filter( Task != 'Combined' )

# Drop rows with missing 
# self-report data for subject 
# FN_041 (66)
# Subject FN_041 did not have any 
# self report data on how high
# for the T1 drug condition
cd = cd %>% 
  filter( Self_report_on_high != '-' )

# Determine count data for misses/correct rejections
cd$Z = cd$Trials - cd$Counts
cd$Y = cd$Counts

# Break self-reported high into 
SRH_range = seq( 0, 100, 20 )

cd$SRH_bins = 0
for ( i in 2:length( SRH_range ) ) {
  sel = cd$Self_report_on_high > SRH_range[i-1] & 
    cd$Self_report_on_high <= SRH_range[i]
  cd$SRH_bins[sel] = SRH_range[i]
}

# Create a new data frame that
# b) Has separate variables for hits and false alarms
all_ptbp = cd %>% 
  filter( Response_type == 'Hits' ) %>% 
  arrange( Task, Condition, Timepoints, Subject ) %>% 
  mutate( H = Counts / Trials,
          fH = Counts, 
          Npos = Trials )
tmp = cd %>% 
  filter( Response_type == 'False_alarms' ) %>% 
  arrange( Task, Condition, Timepoints, Subject ) %>% 
  mutate( FA = Counts / Trials,
          fFA = Counts,
          Nneg = Trials )
all_ptbp$FA = tmp$FA
all_ptbp$fFA = tmp$fFA
all_ptbp$Nneg = tmp$Nneg

# Add in correct rejections and 
# misses
all_ptbp$CR = 1 - all_ptbp$FA
all_ptbp$fCR = all_ptbp$Nneg - all_ptbp$fFA
all_ptbp$M = 1 - all_ptbp$H
all_ptbp$fM = all_ptbp$Npos - all_ptbp$fH

# Clean up workspace
rm( tmp )

# Convert the hit/false alarm rates into 
# estimates of d' and bias
all_ptbp = all_ptbp %>% 
  mutate( dp = sdt_calc_binary( fH, Npos, fFA, Nneg, 1, T ),
          crt = sdt_calc_binary( fH, Npos, fFA, Nneg, 1, F ) )


# Create variables for constructing 
# design matrices for d' and bias
init_X = create_design_mat( cd )

# Path to folder with model results
setwd( 'Data/Posterior_estimates' )
model_dir = getwd()
# Filenames for model results
fname_models = paste( 'M', 1:3, 
         '_Nback_posteriors.RData', sep = '' )
fname_models = c( fname_models,
                  'High_M4_Nback_posteriors.RData',
                  'Not_high_M4_Nback_posteriors.RData' )
setwd( proj_dir )

# Path to directory for figures
setwd( 'Figures' )
plot_dir = getwd()
setwd( proj_dir )

# Initialize variable to store model results
mdl = NA

###
### 2) Additional functions
###

# 2.1)
convergenceExtract = function( fit, parName = NULL ) {
  # Purpose:
  # Extract convergence diagnostics from a Stan fit object.
  # Arguments:
  # fit     - A Stan fit object
  # parName - An optional string giving the final parameter label 
  #           of the subset of the output to include
  # Notes:
  # Extracting the convergence statistics can be slow, especially 
  # when a large number of parameters were stored.
  # Returns:
  # A list with the Gelman-Rubin convergence statistic, the 
  # effective number of samples, and the total number of samples.
  
  Rhat = summary(fit)[,"Rhat"]
  n_eff = summary(fit)[,"n_eff"]
  totSampSize = nrow( as.matrix( fit ) )
  
  return( list( Rhat = Rhat, n_eff = n_eff, 
                totSampSize = totSampSize ) )
}

# 2.2)
findDec = function( x, spacing = 10 ) {
  # Purpose:
  # Determines the rounded leading digit and the 
  # number of trailing zeros for a number or the 
  # number of decimal places.
  # Arguments:
  # x       - A vector of values
  # spacing - The value whose exponent should be increased
  # Returns:
  # A vector giving the leading digit, the number of 
  # trailing/leading zeros, the same but in scientific 
  # notation, and 1 if it's trailing zeros, -1 if it's 
  # decimal places.
  
  mx = max( x )
  rnd = mx
  
  if ( round( mx ) > 1 ) {
    
    inc = 0;
    
    while ( rnd > 1 ) {
      inc = inc + 1;
      rnd = round( mx/( spacing^inc ) )
    }
    
    v = round( mx/spacing^(inc-1) )
    f = spacing^(inc-1)
    out = c( v,f,inc-1,1)
  }
  
  if ( round( mx ) == 1 ) {
    
    out = c( 1, 1, 1, 0 )
    
  }
  
  if ( round( mx ) == 0 ) {
    
    inc = 0;
    
    while ( rnd < 1 ) {
      inc = inc + 1;
      rnd = round( mx*( spacing^inc ) )
    }
    
    v = round( mx*spacing^(inc) )
    f = spacing^(inc)
    out = c( v,f,inc,-1)
    
  }
  
  return( out )
}

# 2.3)
quick_CI_plot = function( pos, 
                          pred,
                          lbls = paste( 'CI',
                                        c(2.5,16,84,97.5), 
                                        sep = '_' ),
                          pt = 21, 
                          lt = .05, 
                          lnSz = 2, 
                          clr = 'grey' ) {
  # Purpose:
  # A convenience function to quickly add 
  # inner and outer credible intervals to an 
  # already existing plot.
  # Arguments:
  # pos    - The x-axis position for the intervals
  # pred   - A labeled matrix with the lower and upper 
  #          boundes for the intervals
  # lbls   - A character vector giving the labels 
  #          for the outer and inner intervals
  # lt     - The width of the error bars
  # lnSz   - The width of the lines
  # clr    - The color of the intervals
  
  segments( pos, pred[,lbls[1]],
            pos, pred[,lbls[4]], lwd = lnSz, col = clr )
  errorBars( pos, rbind( pred[,lbls[2]], pred[,lbls[3]] ),
             length = lt, lwd = lnSz, col = clr )
  points( pos, pred[,'Mean'], pch = pt, bg = clr, col = clr )
  
}

# 2.4)
extract_model_results = function( model_num, extract = T ) {
  # Purpose:
  # A convenience function to load in the model estimation 
  # results for a specified model.
  # Arguments:
  # model_num - The model number (1 - 4)
  # extract   - Logical; if true, loads in the results, 
  #             otherwise deletes the existing results
  # Returns:
  # Either 1) loads in the list of results and priors and 
  # updates the variable 'mdl' (which should already 
  # exist in the workspace), or 2) removes the list 
  # of results and priors and resets 'mdl' to equal NA.
  
  if ( extract ) {
    
    # Load in model results
    setwd( model_dir )
    load( fname_models[ model_num ] )
    setwd( proj_dir )
    
    if ( exists( 'm1' ) ) mdl <<- m1
    if ( exists( 'm2' ) ) mdl <<- m2
    if ( exists( 'm3' ) ) mdl <<- m3
    if ( exists( 'm4' ) ) mdl <<- m4
    
  } else {
    
    # Clean up workspace
    mdl <<- NA
    objs = ls(pos = ".GlobalEnv")
    obj_to_remove = c(
      paste( 'm', model_num, sep = '' ),
      paste( 'm', model_num, '_priors', sep = '' )
    )
    rm( list = objs[objs %in% obj_to_remove ], 
        pos = ".GlobalEnv" )
  }
  
}

###
### 3) Figures for convergence diagnostics
###

if ( run_code[1] ) {
  
  if ( savePlot ) {
    
    setwd( plot_dir )
    pdf( 'Convergence_check.pdf', width = 12 )
    setwd( proj_dir )
    
    # Create a table of contents for figures
    # Initialize table of contents
    tbl_cnt = c(
      'Model 1 (Null)',
      'Model 2 (Intermediary)',
      'Model 3 (Full)',
      'Model 4 (Self-report' )
    
    # Generate table of contents
    tableContents( tbl_cnt )
    
  }
  
  # Loop over models
  for ( i in 1:4 ) {
    
    extract_model_results( i )
    
    conv = convergenceExtract( mdl$est, 'Tau' )
    
    # Subject-level parameters
    sp = grep( 'Subject', names( conv$Rhat ) )
    
    # Plotting characteristics
    lnSz = 2
    txtSz = 1.5
    axSz = 1.25
    axPos = -1.5
    
    if ( !savePlot ) x11( width = 12 )
    layout( cbind( 1, 2 ) )
    
    ### Gelman-Rubin
    
    xl = c( 0, length( conv$Rhat ) + 1 )
    yl = lowerUpper( .05, conv$Rhat )
    if ( yl[2] < 1.1 ) yl[2] = 1.11
    
    # Create a blank plot
    blankPlot( xl, yl )
    
    # Grid lines
    
    # Group-level means
    vertLines( min( sp ) - .5, yl, col = 'grey', lwd = 2 )
    # Subject-level parameters
    vertLines( max( sp ) + .5, yl, col = 'grey', lwd = 2 )
    
    # Poor convergence
    horizLines( 1, xl, col = 'grey', lwd = 2 )
    # Cut-off for poor convergence
    horizLines( 1.1, xl, lty = 2, lwd = 2 )
    text( length( conv$Rhat )/2, 
          1.1, 'Poor convergence', cex = axSz, pos = 3 )
    
    # Gelman-Rubin statistics
    points( 1:length( conv$Rhat ), conv$Rhat, pch = 19 )
    
    # Axes and labels
    customAxes( xl, yl,
                label = c( ' ', 
                           'Gelman-Rubin Statistic' ),
                lnSz = lnSz, lbSz = txtSz )
    axis( 2, c( yl[1], 1, 1.1, yl[2] ),
          tick = F, line = axPos, cex.axis = axSz )
    
    title( 'Convergence diagnostics',
           cex = txtSz )
    
    ### Effective sample size
    
    xl = c( 0, length( conv$n_eff ) + 1 )
    yl = c(0,conv$totSampSize)
    
    # Create a blank plot
    blankPlot( xl, yl )
    
    # Grid lines
    
    lsp = 1000
    horizLines( seq( yl[1], yl[2], lsp ),
                xl, col = 'grey', lwd = lnSz )
    
    # Group-level means
    vertLines( min( sp ) - .5, yl, col = 'grey', lwd = 2 )
    # Subject-level parameters
    vertLines( max( sp ) + .5, yl, col = 'grey', lwd = 2 )
    
    # Effective sample-size
    points( 1:length( conv$n_eff ), conv$n_eff, pch = 19 )
    
    # Axes and labels
    customAxes( xl, yl,
                label = c( ' ', 
                           'Effective sample size' ),
                lnSz = lnSz, lbSz = txtSz )
    axis( 2, seq( 0, yl[2], lsp ),
          tick = F, line = axPos, cex.axis = axSz )
    
    title( 'MCMC samples',
           cex = txtSz )
    
    mtext( 'Parameters', side = 1, cex = txtSz,
           line = -2, outer = T )
    
    mtext( paste( 'Model', i ), side = 3, cex = txtSz,
           line = -2, outer = T )
    
    extract_model_results( i, extract = F )
    
  }
  
  if ( savePlot ) dev.off()
}

###
### 4) Model comparisons
###

if ( run_code[2] ) {
    
    # Initialize list for 
    # leave-one_out cross validation
    # measure
    loo_cv = list(
      M1 = c(),
      M2 = c(),
      M3 = c(),
      M4 = c()
    )
    bs = list(
      M1 = c(),
      M2 = c(),
      M3 = c(),
      M4 = c()
    )
    
    # Loop over models
    for ( i in 1:4 ) {
      
      # Load in model results
      extract_model_results( i )

      loo_cv[[i]] = mdl$LOOCV
      bs[[i]] = mdl$bs
      
      # Clean up workspace
      extract_model_results( i, extract = F )
      
    }
    
    # Data frame with model comparison results
    mc = 
      data.frame(
        Model = 1:4,
        LOO_IC = c( loo_cv[[1]]$estimates[3,1],
                    loo_cv[[2]]$estimates[3,1],
                    loo_cv[[3]]$estimates[3,1],
                    loo_cv[[4]]$estimates[3,1] ),
        LML = c( bs[[1]]$logml,
                 bs[[2]]$logml,
                 bs[[3]]$logml,
                 bs[[4]]$logml )
      )
    # Compute Akaike weights
    mc$AW = icWeights( mc$LOO_IC )
    # Compute stacking weights
    mc$SW = as.numeric( loo_model_weights( loo_cv ) )
    # Compute Bayes factor
    ref = which.min( mc$LML )
    mc$BF = sapply(
      bs, function(x) bf( x, bs[[ref]] )$bf )
    # Convert the log of the marginal probabilities 
    # into relative probabilities
    mc$PML = mc$LML - min( mc$LML )
    mc$PML = exp( mc$PML )
    mc$PML = mc$PML / sum( mc$PML )
    
    # Create a similar table for only for the 
    # first three models
    mc2 = mc[1:3,]
    # Compute Akaike weights
    mc2$AW = icWeights( mc2$LOO_IC )
    # Compute stacking weights
    mc2$SW = as.numeric( loo_model_weights( loo_cv[1:3] ) )
    # Compute Bayes factor
    ref = which.min( mc$LML )
    
    mc2$BF = sapply(
      bs[1:3], function(x) bf( x, bs[[ref]] )$bf )
    # Convert the log of the marginal probabilities 
    # into relative probabilities
    mc2$PML = mc2$LML - min( mc2$LML )
    mc2$PML = exp( mc2$PML )
    mc2$PML = mc2$PML / sum( mc2$PML )
    
}

###
### 5) Posterior retrodictive checks (Model 3)
###

# 5.1)
quick_pred = function( ppc ) {
  # Purpose:
  # A convenience function to quickly generate 
  # model predictions for the overall hit and 
  # correct rejection rates.
  # Arguments:
  # ppc      - A matrix with the simulated observations 
  #            over the posterior samples
  
  h = cd$Response_type == 'Hits'
  fn = function(x) by( h * (1 - x/cd$Trials ) + 
                         (!h)*x/cd$Trials, 
                       list( 
                         cd$Timepoints,
                         cd$Condition,
                         cd$Task,
                         cd$Response_type ), mean )[1:24]
  ppc_1 = apply( apply( ppc, 1, fn ), 1, 
                 quick_desc )
  
  out = aggregate( h * (1 - cd$Z/cd$Trials ) + 
                     (!h)*cd$Z/cd$Trials,
                   list( cd$Timepoints,
                         cd$Condition,
                         cd$Task,
                         cd$Response_type ),
                   mean )
  colnames( out ) = c(
    'Timepoints', 'Condition', 
    'Task', 'Response_type', 'Observed' )
  out = cbind( out, t( ppc_1 ) )
  sel = out$Response_type == 'Hits'
  out$Response_type[!sel] = 'Correct rejections'
  
  return( out )
}

if ( run_code[3] ) {
  
  # Indicate which model to plot
  mn = 3
  
  # PDF for posterior predictive checks
  if ( savePlot ) {
    
    setwd( plot_dir )
    pdf( 'Observed_versus_predicted.pdf', width = 12 )
    setwd( proj_dir )
    
  }
  
  # Extract results
  extract_model_results( mn )
  
  # Extract posterior predictive distributions
  pred = quick_pred( mdl$post_pc )
  
  # Indicate whether new plotting window should 
  # be generated or if figure should be saved as PDF
  if (!savePlot) {
    x11(width = 12 )
  }
  
  # Plotting characteristics
  lnSz = 2
  lnSz2 = 3
  txtSz = 1.5
  ptSz = 1.25
  axPos = -1.75
  pts = rep( 24, 24 )
  clr = rep( rep( c( 'black', 'white' ), each = 3 ), 4 )
  
  # Extract elements to plot
  mean_val = numeric(24)
  error_bars = matrix( NA, 4, 24 )
  lbls = matrix( ' ', 4, 2 )
  cnd = rep( ' ', 24 )
  rpst = rep( ' ', 24 )
  # Hits
  vh = c( 4:6, 1:3, 10:12, 7:9 ) + 12
  vco = vh - 12
  for ( i in 1:4 ) {
    # Hits
    v1 = vh[ 1:3 + 3*(i-1) ]
    v3 = 1:3 + 6*( i - 1 )
    # Correct omissions
    v2 = vco[ 1:3 + 3*(i-1) ]
    v4 = 4:6 + 6*( i - 1 )
    # Hits
    mean_val[v3] = pred$Observed[v1]
    error_bars[1,v3] = pred$`Q2.5%`[v1]
    error_bars[2,v3] = pred$`Q97.5%`[v1]
    error_bars[3,v3] = pred$`Q16%`[v1]
    error_bars[4,v3] = pred$`Q84%`[v1]
    # Correct omissions
    mean_val[v4] = pred$Observed[v2]
    error_bars[1,v4] = pred$`Q2.5%`[v2]
    error_bars[2,v4] = pred$`Q97.5%`[v2]
    error_bars[3,v4] = pred$`Q16%`[v2]
    error_bars[4,v4] = pred$`Q84%`[v2]
    # Labels
    cnd[v3] = unique( pred$Condition[v1] )
    rpst[v3] = unique( pred$Response_type[v1] )
    cnd[v4] = unique( pred$Condition[v2] )
    rpst[v4] = unique( pred$Response_type[v2] )
    lbls[i,2] = unique( pred$Task[v1] )
  }
  lbls[,1] = cnd[ 3 + 6*(0:3) ]
  lbls[ lbls[,2] == 'NBack_0', 2 ] = '0-back'
  lbls[ lbls[,2] == 'NBack_2', 2 ] = '2-back'
  pts[ rpst == 'Hits' ] = 22
  
  # Create a blank plot
  xl = c( .5, 24.5 ); yl = c( .7, 1 )
  blankPlot( xl, yl )
  
  # Add grid lines
  horizLines( seq( .7, 1, .05 ), xl, col = 'grey80', 
              lwd = lnSz )
  customAxes( xl, yl )
  # Separate task and conditions
  vertLines( c( 6.5, 12.5, 18.5 ), yl, lwd = lnSz )
  
  # Add axis labels
  axis( 3, 24/4 * 1:4 - 2.5, 
        lbls[,1],
        tick = F, line = axPos, cex.axis = txtSz )
  axis( 1, 24/2 * 1:2 - 5.5, 
        unique( lbls[,2] ),
        tick = F, line = axPos, cex.axis = txtSz )
  axis( 2, seq( .7, 1, .1 ),
        seq( .7, 1, .1 ) * 100, 
        tick = F, line = axPos, cex.axis = txtSz )
  mtext( 'Average percentage', side = 2, line = 1.5, cex = txtSz )
  
  # Add plotting elements
  segments( 1:24, error_bars[1,],
            1:24, error_bars[2,], 
            lwd = lnSz, col = 'grey50' )
  errorBars( 1:24, error_bars[3:4,], 
             lwd = lnSz, length = .05, col = 'grey50' )
  
  # Lines by condition
  for ( i in 1:8 ) {
    v1 = 1:3 + 3*(i-1)
    lines( v1, mean_val[v1], lwd = lnSz )
  }
  
  # Add average percentages
  points( 1:24, mean_val, pch = pts, bg = clr )
  
  # Add labels for timepoints
  text( 1:3, rep( .95, 3 ), c( expression(T[1]),
                               expression(T[2]),
                               expression(T[3]) ),
        pos = 1, cex = txtSz )
  
  legend( -.5, yl[1] - diff(yl)*.1, 
          c( 'Correct Omissions', 'Hits' ),
          pch = rev( unique( pts ) ), 
          pt.bg = c( 'white', 'black' ),
          bty = 'n', horiz = T, xpd = T, cex = txtSz )
  
  legend( 7.5, yl[1] - diff(yl)*.1, 
          c( 'T1: Pre', 'T2: Peak', 'T3: Post' ),
          bty = 'n', horiz = T, xpd = T, cex = txtSz )
  
  legend( 7.5, yl[2] + diff(yl)*.2,
          c( 'Observed', 'Predicted' ),
          fill = c( 'black', 'grey50' ), bty = 'n',
          horiz = T, xpd = T, cex = txtSz )
  
  # Clean up workspace
  extract_model_results( mn, extract = F )
  
  if ( savePlot ) dev.off()
  
}

###
### 6) Marginal posteriors (d')
###

if ( run_code[4] ) {
  
  # PDF for d' effect sizes
  if ( savePlot ) {
    
    setwd( plot_dir )
    pdf( 'Effect_sizes_d_prime.pdf', width = 12 )
    setwd( proj_dir )
    
  }
  
  # Model number
  m = 3
  
  # Extract results
  extract_model_results( m )
  
  # 95% credible intervals
  ol = posterior_interval( mdl$est,
                            prob = .95 )
  # 68% credible intervals
  il = posterior_interval( mdl$est,
                           prob = .68 )
  # Posterior means
  mn = rstanarm::fixef( mdl$est )
  
  # Isolate d' estimates
  sel = 1:nrow( ol )
  sel = sel[ ( sel %in% grep( 'dp', rownames( ol ) ) ) & 
               !( sel %in% grep( 'Subject', rownames( ol ) ) ) ]
  # Marginal posteriors to plot
  ptbp = data.frame(
    Variable = rownames( ol )[sel],
    stringsAsFactors = FALSE )
  ptbp = cbind( ptbp, ol[sel,], il[sel,] )
  rownames( ptbp ) = rep( NULL, nrow( ptbp ) )
  ptbp$Mean = mn[ grep( 'dp', names( mn ) ) ]
  # Change column names for easier extracting
  colnames( ptbp ) = c( 'Variables',
                        'CI_2.5', 'CI_97.5',
                        'CI_16', 'CI_84',
                        'Mean' )
  
  # Remove baseline estimates
  ptbp = ptbp[ -grep( 'Baseline', ptbp$Variable ), ]
  
  # Indicate whether new plotting window should 
  # be generated or if figure should be saved as PDF
  if (!savePlot) {
    x11(width = 12 )
  } else {
    
  }
  par( mar = c( 7, 4, 3, 2 ) )
  
  # Plotting characteristics
  lnSz = 2
  txtSz = 1.5
  ptSz = 1.25
  axPos = -1.75
  
  yl = lowerUpper( .5, as.matrix( ptbp[,-1] ) )
  xl = c( .5, .5 + nrow( ptbp ) )
  blankPlot( xl, yl )
  horizLines( seq( yl[1], yl[2], .5 ), 
              xl, lwd = lnSz, col = 'grey80' )
  horizLines( 0, xl, lwd = lnSz )
  customAxes( xl, yl, lnSz = lnSz )
  
  axis( 2, seq( yl[1], yl[2], .5 ), 
        tick = F, line = -1.75, cex.axis = txtSz )
  mtext( "Estimated effect size (d')", 
         side = 2, line = 1.5, cex = txtSz )
  
  # Identify significant effects
  ssp = ptbp$CI_2.5 > 0 | ptbp$CI_97.5 < 0
  clr = rep( 'black', nrow( ptbp ) )
  clr[ ssp ] = 'blue'
  
  quick_CI_plot( 1:nrow( ptbp ),
                 ptbp, clr = clr )
  
  # Labels for regression coefficients
  lbl = c(
    'Pre (Placebo, 0-back)',
    'Peak (Placebo, 0-back)',
    'Post (Placebo, 0-back)',
    'Peak (Drug, 0-back)',
    'Post (Drug, 0-back)',
    'Pre (Placebo, 2-back)',
    'Peak (Placebo, 2-back)',
    'Post (Placebo, 2-back)',
    'Peak (Drug, 2-back)',
    'Post (Drug, 2-back)',
    'Visit order effect',
    'Intoxication effect' )
  
  xa = 1:nrow( ptbp )
  axis( 1, at = xa, tick = F, labels = FALSE )
  text( x = xa, 
        y = par()$usr[3]+0.0*(par()$usr[4]-par()$usr[3]),
        labels = lbl, srt = 45, adj = 1, xpd = TRUE )
  
  dl = max( cumsum( sapply( lbl, 
                            function(x) 
                              length( grep( '0-back', x ) ) ) ) )
  dl = c( dl,
          which( ptbp$Variables == 'dp_Visit_2' ) - 1 )
  vertLines( dl + .5, yl, lwd = lnSz, col = c( 'grey', 'black' ) )
  
  if ( savePlot ) dev.off()
  
}

###
### 7) Marginal posteriors (response bias)
###

if ( run_code[5] ) {
  
  # PDF for d' effect sizes
  if ( savePlot ) {
    
    setwd( plot_dir )
    pdf( 'Effect_sizes_bias.pdf', width = 12 )
    setwd( proj_dir )
    
  }
  
  # Model number
  m = 3
  
  # Extract results
  extract_model_results( m )
  
  # 95% credible intervals
  ol = posterior_interval( mdl$est,
                           prob = .95 )
  # 68% credible intervals
  il = posterior_interval( mdl$est,
                           prob = .68 )
  # Posterior means
  mn = rstanarm::fixef( mdl$est )
  
  # Isolate d' estimates
  sel = 1:nrow( ol )
  sel = sel[ ( sel %in% grep( 'crt', rownames( ol ) ) ) & 
               !( sel %in% grep( 'Subject', rownames( ol ) ) ) ]
  # Marginal posteriors to plot
  ptbp = data.frame(
    Variable = rownames( ol )[sel],
    stringsAsFactors = FALSE )
  ptbp = cbind( ptbp, ol[sel,], il[sel,] )
  rownames( ptbp ) = rep( NULL, nrow( ptbp ) )
  ptbp$Mean = mn[ grep( 'crt', names( mn ) ) ]
  # Change column names for easier extracting
  colnames( ptbp ) = c( 'Variables',
                        'CI_2.5', 'CI_97.5',
                        'CI_16', 'CI_84',
                        'Mean' )
  
  # Indicate whether new plotting window should 
  # be generated or if figure should be saved as PDF
  if (!savePlot) {
    x11(width = 12 )
  }
  
  par( mar = c( 7, 4, 3, 2 ) )
  
  # Plotting characteristics
  lnSz = 2
  txtSz = 1.5
  ptSz = 1.25
  axPos = -1.75
  
  yl = lowerUpper( .5, as.matrix( ptbp[,-1] ) )
  xl = c( .5, .5 + nrow( ptbp ) )
  blankPlot( xl, yl )
  horizLines( seq( yl[1], yl[2], .25 ), 
              xl, lwd = lnSz, col = 'grey80' )
  horizLines( 0, xl, lwd = lnSz )
  customAxes( xl, yl, lnSz = lnSz )
  
  axis( 2, seq( yl[1], yl[2], .25 ), 
        tick = F, line = -1.75, cex.axis = txtSz )
  mtext( "Estimated effect size (d')", 
         side = 2, line = 1.5, cex = txtSz )
  
  # Identify significant effects
  ssp = ptbp$CI_2.5 > 0 | ptbp$CI_97.5 < 0
  clr = rep( 'black', nrow( ptbp ) )
  clr[ ssp ] = 'blue'
  
  quick_CI_plot( 1:nrow( ptbp ),
                 ptbp, clr = clr )
  
  # Labels for regression coefficients
  lbl = c(
    'Bias (0-back)',
    'Pre (Placebo, 0-back)',
    'Peak (Placebo, 0-back)',
    'Post (Placebo, 0-back)',
    'Peak (Drug, 0-back)',
    'Post (Drug, 0-back)',
    'Bias (2-back)',
    'Pre (Placebo, 2-back)',
    'Peak (Placebo, 2-back)',
    'Post (Placebo, 2-back)',
    'Peak (Drug, 2-back)',
    'Post (Drug, 2-back)',
    'Visit order effect',
    'Intoxication effect' )
  
  xa = 1:nrow( ptbp )
  axis( 1, at = xa, tick = F, labels = FALSE )
  text( x = xa, 
        y = par()$usr[3]+0.0*(par()$usr[4]-par()$usr[3]),
        labels = lbl, srt = 45, adj = 1, xpd = TRUE )
  
  dl = max( cumsum( sapply( lbl, 
                            function(x) 
                              length( grep( '0-back', x ) ) ) ) )
  dl = c( dl,
          which( ptbp$Variables == 'crt_Visit_2' ) - 1 )
  vertLines( dl + .5, yl, lwd = lnSz, col = c( 'grey', 'black' ) )
  
  if ( savePlot ) dev.off()
  
}

###
### 8) Correlation between Self-report and THC administration
###

if ( run_code[6] ) {
  
  # Package for flexible Bayesian modeling
  # install.packages( 'brms' )
  library( brms )
  
  # Create variables for constructing 
  # design matrices
  init_X = create_design_mat( cd )
  
  # Data frame 
  dtbf = data.frame(
    Y = cd$pSRH,
    Drug = init_X$Drug,
    S = cd$Subject,
    stringsAsFactors = F
  )
  
  # Self-report measures are on the interval [0,1]
  # Beta regression needs values in the interval (0,1)
  # Use a transformation to slightly shift extreme 
  # values.
  # Reference:
  # Smithson, M., & Verkuilen, J. (2006). A better 
  # lemon squeezer? Maximum-likelihood regression 
  # with beta-distributed dependent variables.
  # Psychological Methods, 11, 54 - 71.
  # http://dx.doi.org/10.1037/1082-989X.11.1.54.
  trn = function( x ) {
    
    n = length(x)
    s = .5
    out = ( x * ( n - 1 ) + s )/n
    
    return( out )
  }
  
  # Apply transformation to data
  dtbf$Yt = trn( dtbf$Y )
  
  # Fit a simple multilevel model 
  # using the 'brms' package
  setwd( proj_dir )
  setwd( 'Data/Posterior_estimates' )
  
  fit_model = F
  
  if ( fit_model ) {
    
    # Fit a multilevel beta regression 
    # model that estimates:
    # 1) The correlation between self-report
    #    and THC administration
    # 2) A random intercept and slope for 
    #    each subject
    SRH_drug_cor = brm( 
      Yt ~ Drug + (Drug|S),
      data = dtbf,
      family = Beta(),
      warmup = 1000,
      iter = 3500,
      cores = 4,
      chains = 4 )
    
    # Posterior predictive check
    chk = brms::posterior_predict( SRH_drug_cor )
    f = function( x ) {
      out = as.numeric( by( x, list( dtbf$Drug ), mean ) )
      out
    }
    ppc = apply( chk, 1, f )
    f = function( x ) {
      out = unlist( by( x, list( dtbf$Drug ), 
                            quantile, prob = c( .025, .975 ) ) )
      out
    }
    ppc2 = apply( chk, 1, f )
    
    # Save results
    save( SRH_drug_cor, dtbf, ppc, ppc2, chk, 
          file = 'SRH_drug_correlation.RData' )
    setwd( proj_dir )
    
  } else {
    load( 'SRH_drug_correlation.RData' )
    setwd( proj_dir )
  }
  
  # Plot the results
  val = by( dtbf$Yt, list( dtbf$Drug ), table )
  xa = as.numeric( names( val[[1]] ) )
  ya = .2 * ( as.numeric( val[[1]] ) ) / 
    sum( as.numeric( val[[1]] ) )
  cnd = length( xa )
  xa = c( xa, as.numeric( names( val[[2]] ) ) )
  ya = c( ya, .2 * ( as.numeric( val[[2]] ) ) / 
    sum( as.numeric( val[[2]] ) ) )
  cnd = c( rep( 0, cnd ), rep( 1, length(xa) - cnd ) )
  xa1 = rep( 0, length( ya ) )
  xa1[ cnd == 1 ] = .8
  xa2 = ya + .8 * cnd
  
  # Indicate whether new plotting window should 
  # be generated or if figure should be saved as PDF
  if (!savePlot) {
    x11()
  }
  
  # Plotting characteristics
  ptSz = 1.25
  lnSz = 2
  axPos = -1
  axSz = 1.5
  txtSz = 1.5
  
  # Create a blank plot
  xl = c( -.025, 1.025 )
  yl = c( -.025, 1.025 )
  blankPlot( xl, yl )
  
  # 95% prediction interval
  prd = apply( ppc2, 1, mean )
  polygon( c( .2, .2, .8, .8 ),
           c( prd[1], prd[2], prd[4], prd[3] ),
           col = 'grey80', border = NA )
  
  # Guidelines
  horizLines( seq( 0, 1, .2 ), c(0,.2), col = 'grey', lwd = lnSz )
  horizLines( seq( 0, 1, .2 ), c(.8,1), col = 'grey', lwd = lnSz )
  
  # Add observations
  segments( xa1, xa, xa2, xa, lwd = 4 )
  
  # Add labels and axes
  customAxes( xl, yl, lnSz = lnSz )
  axis( 1, c( 0, 1 ), c( 'Placebo', 'Drug' ),
        tick = F, line = axPos, cex.axis = axSz )
  axis( 2, seq( 0, 1, .2 ), seq( 0, 1, .2 )*100,
        tick = F, line = axPos, cex.axis = axSz )
  mtext( 'Self-reported intoxication',
         side = 2, cex = txtSz, line = 2.5 )
  
  # Extract estimates
  est = brms::fixef( SRH_drug_cor )
  # Posterior p-value
  pst = as.matrix( SRH_drug_cor )
  pv = sum( pst[,'b_Drug'] < 0 )
  
  # Best-fitting line
  ya = apply( ppc, 1, mean )
  segments( .2, ya[1], .8, ya[2], lwd = lnSz, lty = 2 )
  
  title(
      substitute(
        paste( beta, " = ", v, ", p < 0.0001", sep = ""),
        list(v=round(est[2,1],2))
      ), cex = txtSz )
  
}

###
### 9) Extraction of baseline vitals
###

if ( run_code[7] ) {
  
  # Load in data with all vitals
  setwd( 'Data/Original_files' )
  av = read.csv( 
    file = "Vitals_All.csv",
    header = T,
    stringsAsFactors = FALSE )
  setwd( proj_dir )
  
  # Data frame with elements to plot
  dtbp = data.frame(
    ID = av$study_id,
    Subject = NA,
    stringsAsFactors = FALSE
  )
  # Match subject numbers
  subject_ID = dat %>% 
    group_by( Subject ) %>% 
    summarize( ID = unique( ID ) )
  for ( i in 1:nrow( subject_ID ) ) {
    if ( any( dtbp$ID == subject_ID$ID[i] ) ) {
      sel = dtbp$ID == subject_ID$ID[i]
      dtbp$Subject[sel] = subject_ID$Subject[i]
    }
  }
  # Session number
  dtbp$Visit = sapply( av$redcap_event_name, 
                       function(x) 
                         as.numeric( 
                           strsplit( x, split = '_' )[[1]][2] ) )
  
  # Time (in minutes) since initial measurement
  
  # Function to convert date into number of minutes
  convert_to_minutes = function( tm ) {
    
    tm = strsplit( tm, split = ' ' )[[1]][2]
    tm = strsplit( tm, split = ':' )[[1]]
    out = as.numeric( tm[1] ) * 60 + as.numeric( tm[2] )
    
    return( out )
  }
  dtbp$Time = as.numeric( 
    sapply( av$reading_time, convert_to_minutes ) )
  
  # Determine minimum time
  min_time = dtbp %>% 
    group_by( Subject, Visit ) %>% 
    summarize( M = min( Time ) )
  # Adjust times by minimum
  for ( i in 1:nrow( min_time ) ) {
    if ( any( dtbp$Subject == min_time$Subject[i] &
              dtbp$Visit == min_time$Visit[i] ) ) {
      sel = dtbp$Subject == min_time$Subject[i] &
        dtbp$Visit == min_time$Visit[i]
      dtbp$Time[sel] = dtbp$Time[sel] - min_time$M[i]
    }
  }
  
  # Discrete measurement points
  # (For aligning observations over subjects)
  dtbp$Measurement = 1
  val = paste( 'q', 1:16, sep = '' )
  for ( i in 1:length( val ) ) {
    sel = grep( val[i], av$redcap_event_name )
    dtbp$Measurement[sel] = i + 2
  }
  sel = grep( 'baseline_2', av$redcap_event_name )
  dtbp$Measurement[sel] = 2
  
  # Measurements for vitals
  
  # Number of heart beats per minute
  dtbp$Pulse = av$pulse
  # Blood pressure when heart is contracting
  dtbp$Systolic = av$bp_systolic
  # Blood pressure when heart is between beats
  dtbp$Diastolic = av$bp_diastolic
  
  # Self-report
  dtbp$DEQ_1 = av$deq_1
  dtbp$DEQ_2 = av$deq_2
  
  # Indicator for whether measurements had 
  # to be repeated
  dtbp$Repeat = 0
  # The number of attempts for a measurement
  dtbp$Attempt = 1
  
  # Extract repeated measurements
  sel = av$were_the_vitals_repeated_w == 1 & 
    !is.na( av$were_the_vitals_repeated_w )
  tmp = dtbp[sel,]
  tmp$Attempt = 2
  tmp$Pulse = av$pulse2[sel]
  tmp$Systolic = av$bp_systolic2[sel]
  tmp$Diastolic = av$bp_diastolic2[sel]
  dtbp$Repeat[sel] = 1
  dtbp = rbind( dtbp, tmp )
  
  # Match visits with drug/placebo condition
  condition_by_visit = cd %>% 
    group_by( Subject, Visit, Condition ) %>% 
    summarize( ID = unique( ID ) )
  dtbp$Condition = 'Placebo'
  for ( i in 1:nrow( condition_by_visit ) ) {
    if ( any( dtbp$Subject == condition_by_visit$Subject[i] &
              dtbp$Visit == condition_by_visit$Visit[i] ) ) {
      sel = dtbp$Subject == condition_by_visit$Subject[i] &
        dtbp$Visit == condition_by_visit$Visit[i]
      dtbp$Condition[sel] = condition_by_visit$Condition[i]
    }
  }
  
  # Sort data
  dtbp = dtbp %>% 
    arrange( Subject, Visit, Measurement )
  
  # Exclude NA responses
  dtbf = dtbp[ !is.na( dtbp$DEQ_2 ), ]
  dtbf$Y = ( dtbf$DEQ_2 - mean( dtbf$DEQ_2 ) ) / sd( dtbf$DEQ_2 )
  dtbf$X = ( dtbf$Pulse - mean( dtbf$Pulse ) ) / sd( dtbf$Pulse )
  # Apply simple mixed effects model to compute correlation
  res = stan_lmer( Y ~ -1 + X + (1|Subject),
                   dat = dtbf,
                   iter = 3500,
                   warmup = 1000,
                   chains = 4,
                   cores = 4 )
  # Results
  output = list(
    R = round( fixef( res ), 2 ),
    SD = round( sd( as.matrix( res )[,1] ), 2 ), 
    UI = round( posterior_interval( res, pars = 'X' ), 2 ),
    p_value = sum( as.matrix( res )[,1] < 0 )/10000
  )
  
  dtbf = cd[,c('Subject','Condition','Timepoints')]
  dtbf$Y = cd$Self_report_on_high
  dtbf = dtbf[ dtbf$Timepoints != 'T1_Pre_drug', ]
  dtbf$Y = ( dtbf$Y - mean( dtbf$Y ) ) / sd( dtbf$Y )
  dtbf$DvP = 0; dtbf$DvP[ dtbf$Condition == 'Drug' ] = 1
  dtbf$T3 = 0; dtbf$T3[ dtbf$Timepoints == 'T3_Post_drug' ] = 1
  dtbf$I = dtbf$DvP * dtbf$T3
  # Apply simple mixed effects model for self-report on high
  res2 = stan_lmer( Y ~ 1 + DvP + T3 + I + (1|Subject),
                    dat = dtbf,
                    iter = 3500,
                    warmup = 1000,
                    chains = 4,
                    cores = 4 )
  output = list(
    B = round( fixef( res2 ), 2 )[-1],
    SD = round( sd( as.matrix( res2 )[,2] ), 2 ), 
    UI = round( posterior_interval( res2, pars = c( 'DvP', 'T3', 'I') ), 2 ),
    p_value = sum( as.matrix( res2 )[,2] < 0 )/10000
  )
  
  # Function for descriptive statistics
  f = function(x) {
    is_na = is.na( x )
    N = length( !is_na )
    out = c(
      mean( x, na.rm = T ),
      sd( x, na.rm = T )
    )
    out[3] = out[2] / sqrt( N )
    names( out ) = c( 'M', 'SD', 'SEM' )
    
    return( out )
  }
  
  # Group means
  plt = aggregate( 
    dtbp$DEQ_2, list( dtbp$Measurement, dtbp$Condition ), f )
  colnames( plt ) = c( 'TP', 'Cnd', 'DS' )
  
  # Heart rate differences between neural data groups
  
  groups = list(
    U50 = as.numeric( 
      strsplit( paste( '30 32 36 39 42 43 46 48 49 51 52 53 56',
             '57 59 60 61 66 70 75 81 83 84 87 90 93 94' ), split = ' ' )[[1]] ),
    L50 = as.numeric(
      strsplit( paste( '31 35 37 40 45 47 50 54 55 62 63 64 65', 
                       '67 68 71 72 73 76 77 79 80 82 85 89 92 97' ), 
                split = ' ')[[1]] )
  )
  
  b("
  groups = list(
    U50 = as.numeric( 
      strsplit( paste( '39 42 43 46 48 49 51 52 53 56 57 59 60 61',
                       '66 70 75 81 83 84 87 90 93 94' ), split = ' ' )[[1]] ),
    L50 = as.numeric(
      strsplit( paste( '40 45 47 50 54 55 62 63 65 67 68 72 73 76', 
                       '77 79 80 82 85 89 92 97' ), 
                split = ' ')[[1]] )
  )
  ")
  
}

###
### 10) DEQ correlations
###

if ( run_code[8] ) {
  
  setwd( 'Data/Original_files' )
  av = read.csv( 
    file = "Copy of Adrian_fNIRS_Vitals_DEQ_Kevin.csv",
    header = T,
    stringsAsFactors = FALSE )
  colnames( av )[1] = 'study_id'
  setwd( proj_dir )
  
  # Data frame with elements to plot
  dtbp = data.frame(
    ID = av$study_id,
    Subject = NA,
    stringsAsFactors = FALSE
  )
  # Match subject numbers
  subject_ID = dat %>% 
    group_by( Subject ) %>% 
    summarize( ID = unique( ID ) )
  for ( i in 1:nrow( subject_ID ) ) {
    if ( any( dtbp$ID == subject_ID$ID[i] ) ) {
      sel = dtbp$ID == subject_ID$ID[i]
      dtbp$Subject[sel] = subject_ID$Subject[i]
    }
  }
  # Session number
  dtbp$Visit = sapply( av$redcap_event_name, 
                       function(x) 
                         as.numeric( 
                           strsplit( x, split = '_' )[[1]][2] ) )
  
  # Add variable for condition
  dtbp$Condition = 'Placebo'
  dtbp$Condition[ av$Drug == 1 ] = 'Drug'
  
  # Measurements for vitals
  
  # Number of heart beats per minute
  dtbp$Pulse = av$pulse
  # Blood pressure when heart is contracting
  dtbp$Systolic = av$bp_systolic
  # Blood pressure when heart is between beats
  dtbp$Diastolic = av$bp_diastolic
  
  # Self-report
  dtbp$DEQ_1 = av$deq_1
  dtbp$DEQ_2 = av$deq_2
  dtbp$DEQ_3 = av$deq_3
  dtbp$DEQ_4 = av$deq_4
  dtbp$DEQ_5 = av$deq_5
  
  # Indicator for whether measurements had 
  # to be repeated
  dtbp$Repeat = 0
  # The number of attempts for a measurement
  dtbp$Attempt = 1
  
  # Extract repeated measurements
  sel = av$were_the_vitals_repeated_w == 1 & 
    !is.na( av$were_the_vitals_repeated_w )
  tmp = dtbp[sel,]
  tmp$Attempt = 2
  tmp$Pulse = av$pulse2[sel]
  tmp$Systolic = av$bp_systolic2[sel]
  tmp$Diastolic = av$bp_diastolic2[sel]
  dtbp$Repeat[sel] = 1
  dtbp = rbind( dtbp, tmp )
  
  M = dtbp[ dtbp$Repeat == 0, grep( 'DEQ', colnames( dtbp ) ) ]
  
  if( savePlot ) {
    setwd( 'Figures' )
    pdf( 'DEQ_correlations.pdf' )
    setwd( proj_dir )
  }
  
  if ( !savePlot ) x11()
  lyt = matrix( 0, 5, 5 )
  diag( lyt ) = 1:5
  lyt[2,1] = 6
  lyt[3,1] = 7; lyt[3,2] = 10;
  lyt[4,1] = 8; lyt[4,2] = 11; lyt[4,3] = 13;
  lyt[5,1] = 9; lyt[5,2] = 12; lyt[5,3] = 14; lyt[5,4] = 15;
  lyt[ upper.tri(lyt) ] = 16:25
  layout( lyt )
  
  mn = c( 3, 3, 3, 3 )
  for ( i in 1:5 ) {
    par( mar = mn )
    hist( M[,i], xlab = ' ', ylab = ' ',
          bty = 'l', main = ' ',
          col = 'grey', border = 'white',
          yaxt = 'n', 
          xaxt = 'n' )
    axis( 1, c(0,50,100) )
    if ( i == 1 ) 
      mtext( 'DEQ (1)', side = 3, line = 1.5 )
    if ( i == 5 ) 
      mtext( 'DEQ (5)', side = 4, line = 1.5 )
  }
  
  R = numeric( 10 )
  vn = rbind(
    # x, y
    c( 'DEQ_1', 'DEQ_2' ),
    c( 'DEQ_1', 'DEQ_3' ),
    c( 'DEQ_1', 'DEQ_4' ),
    c( 'DEQ_1', 'DEQ_5' ),
    c( 'DEQ_2', 'DEQ_3' ),
    c( 'DEQ_2', 'DEQ_4' ),
    c( 'DEQ_2', 'DEQ_5' ),
    c( 'DEQ_3', 'DEQ_4' ),
    c( 'DEQ_3', 'DEQ_5' ),
    c( 'DEQ_4', 'DEQ_5' ) )
  for ( i in 1:nrow( vn ) ) {
    par( mar = mn )
    plot( M[,vn[i,1]], M[,vn[i,2]], pch = 19,
          xlab = ' ', ylab = ' ', bty = 'l',
          xaxt = 'n', yaxt = 'n',
          ylim = c(0,100), xlim = c(0,100) )
    axis( 1, c(0,50,100) )
    axis( 2, c(0,50,100) )
    frm = paste( vn[i,1], '~', vn[i,2],
                 '+ (1|Subject)' )
    mdl = lme4::lmer( as.formula( frm ), data = dtbp )
    b = lme4::fixef( mdl )
    abline( a = b[1], b = b[2], col = 'grey', lwd = 2, lty = 2 )
    R[i] = b[2]
  }
  
  pos = c( 1, 2, 5, 3, 6, 8, 4, 7, 9, 10 )
  inc = 1
  for ( i in 1:10 ) {
    par( mar = mn )
    blankPlot()
    # customAxes( c(0,1), c(0,1), pos = 1:4 )
    legend( 'center', paste( 'R =', round( R[ pos[i] ], 2 ) ),
            bty = 'n', cex = 1.25,
            xpd = T )
    if ( pos[i] %in% 1:4 ) {
      mtext( paste( 'DEQ (', pos[i]+1, ')', sep = '' ),
             side = 3, line = 1.5 )
    }
    
    if ( pos[i] %in% c( 4, 7, 9, 10 ) ) {
      mtext( paste( 'DEQ (', inc, ')', sep = '' ),
             side = 4, line = 1.5 )
      inc = inc + 1
    }
    
  }
  
  cor_mat = matrix( ' ', 5, 5 )
  rownames(cor_mat) = paste( 'DEQ (', 1:5, ')', sep = '' )
  colnames(cor_mat) = paste( 'DEQ (', 1:5, ')', sep = '' )
  diag( cor_mat ) = as.character( 
    round( sqrt( diag( var( M ) ) ) ) )
  sel = rbind(
    c( 1, 2 ),
    c( 1, 3 ),
    c( 1, 4 ),
    c( 1, 5 ),
    c( 2, 3 ),
    c( 2, 4 ),
    c( 2, 5 ),
    c( 3, 4 ),
    c( 3, 5 ),
    c( 4, 5 ) )
  for ( i in 1:10 ) {
    cor_mat[ sel[i,1], sel[i,2] ] = as.character(
      round( R[i], 2 ) )
  }
  
  if ( savePlot ) dev.off()
}

###
### 11) Table of results for SDT model
###

if ( run_code[9] ) {
  
  
  # Extract average visit order and SR rating per effect
  Xc = cbind(
    Visit_2 = init_X$O, 
    Self_report_high = init_X$SR
  )
  check = create_design_mat( cd, Xc )
  
  # 11.1)
  quick_results = function( mdl ) {
    # Purpose:
    # Computes peak - pre and post - peak 
    # comparisons for models 3, 4, and 5
    # Arguments:
    # mdl - Output from the 'estimate_sdt' function
    # Returns:
    # A table with the task, condition, comparison, 
    # effect size, posterior standard deviation, 
    # and posterior p-value.
    
    post = as.matrix( mdl$est )
    prm = colnames( post )
    sel = prm[ grep( 'dp', prm ) ]
    sel = sel[ -grep( 'Subject', sel ) ]
    
    # Function to correctly extract a desired condition 
    # from the parameeter names
    mc = function( Task, Condition, Timepoint ) {
      
      tmp = strsplit( sel, split = '_' )
      M = sapply( tmp, function(x) c( x[2], x[3], x[4] ) )
      M[ is.na( M ) ] = ''
      
      out = 
        apply( M, 2, function(x) 
          sum( x == c( Condition, Task, Timepoint ) ) == 3 )
      return( which( out == T ) )
    }
    
    # Function to compute the one-sided posterior p-value
    fp = function( x ) {
      
      if ( mean(x) > 0 ) {
        out = sum( x < 0 )/length( x )
      } else {
        out = sum( x > 0 )/length( x )
      }
      return( out )
    }
    
    # Extract the estimated average d' for 
    # each condition
    M = matrix( NA, nrow( post ), 12 )
    
    vrb = cbind(
      rep( c( '0','2'), each = 6 ),
      rep( rep( c( 'Placebo', 'Drug' ), each = 3 ), 2 ),
      rep( c( 'T1', 'T2', 'T3' ), 4 )
    )
    colnames( M ) = paste(
      paste( vrb[,1], 'back', sep = '-' ),
      vrb[,2],
      vrb[,3], sep = '_' )
    
    # 0-back
    # Determine baseline column
    bs = mc('0','Baseline','')
    M[,'0-back_Drug_T1'] = post[,sel[bs]]
    # Loop over remaining conditions
    for ( i in c(1:3,5:6) ) {
      nc = mc(vrb[i,1],vrb[i,2], vrb[i,3] )
      lbl = paste( vrb[i,1], '-back_', vrb[i,2], '_', vrb[i,3], sep = '' )
      M[,lbl] = rowSums( post[ , sel[ c( bs, nc ) ] ] )
    }
    
    # 2-back
    # Determine baseline column
    bs = mc('2','Baseline','')
    M[,'2-back_Drug_T1'] = post[,sel[bs]]
    # Loop over remaining conditions
    for ( i in ( c(1:3,5:6) + 6 ) ) {
      nc = mc(vrb[i,1],vrb[i,2], vrb[i,3] )
      lbl = paste( vrb[i,1], '-back_', vrb[i,2], '_', vrb[i,3], sep = '' )
      M[,lbl] = rowSums( post[ , sel[ c( bs, nc ) ] ] )
    }
    
    # Adjust posteriors by effect of visit order
    M2 = sapply( 1:12, function(x) {
      B = post[,'dp_Visit_2'];
      out = B * check$X[x,1];
      return( out )
    } )
    M = M + M2
    
    # If self-reported high is included
    if ( 'dp_Self_report_high' %in% sel ) {
      
      M2 = sapply( 1:12, function(x) {
        B = post[,'dp_Self_report_high'];
        out = B * check$X[x,2];
        return( out )
      } )
      M = M + M2
    }
    
    
    # Compute peak - pre and post - peak differences
    out = data.frame(
      Task = rep( c('0-back','2-back'), each = 4 ),
      Condition = rep( rep( c( 'Placebo', 'Drug' ), each = 2 ), 2 ),
      Comparison = rep( c( 'Peak - Pre', 'Post - Peak' ), 4 ),
      Effect_size = NA,
      SD = NA,
      p_value = NA,
      stringsAsFactors = F
    )
    fm = function(x,v1,v2) round( mean(x[,v2] - x[,v1]), 2 )
    out$Effect_size = apply( cbind( c( 1, 2, 4, 5, 7, 8, 10, 11 ),
                                    c( 2, 3, 5, 6, 8, 9, 11, 12 ) ), 1, 
                             function(x) fm( M, x[1], x[2] ) )
    fm = function(x,v1,v2) round( sd(x[,v2] - x[,v1]), 2 )
    out$SD = apply( cbind( c( 1, 2, 4, 5, 7, 8, 10, 11 ),
                           c( 2, 3, 5, 6, 8, 9, 11, 12 ) ), 1, 
                    function(x) fm( M, x[1], x[2] ) )
    fm = function(x,v1,v2) round( fp(x[,v2] - x[,v1]), 4 )
    out$p_value = apply( cbind( c( 1, 2, 4, 5, 7, 8, 10, 11 ),
                                c( 2, 3, 5, 6, 8, 9, 11, 12 ) ), 1, 
                         function(x) fm( M, x[1], x[2] ) )
    
    return( out )
  }
  
  for ( m in 3:5 ) {
    
    setwd( proj_dir )
    
    # Extract results
    extract_model_results( m )
    tbl = quick_results( mdl )
    # print( tbl )
    setwd( proj_dir )
    setwd( 'Documents' )
    write.table( tbl, 
                 row.names = F, quote = F, sep = ',',
                 file = paste( 'Results_M', m, '.csv', sep = '' ) )
    setwd( proj_dir )
    
  }
  
  
  
}


setwd( orig_dir )