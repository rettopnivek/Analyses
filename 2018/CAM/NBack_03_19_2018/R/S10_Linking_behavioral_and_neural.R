# Linking the behavioral and neural data
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you 
# have any questions or comments
# Last updated 2018-06-27

# Table of contents
# 1) Initial setup
# 2) Data processing
# 3) Estimation functions
#   3.1) estimate_sdt
#   3.2) quick_pred_plot
# 4) Model estimation
# 5) Scatter plots of neural data against behavioral
#   5.1) neural_scatter

###
### 1) Initial setup
###

# Save current working directory
orig_dir = getwd()

# Navigate to project directory
setwd( '..' )
proj_dir = getwd()

# Indicate which code segments to run
run_code = c(
  T, # Model estimation
  F  # Scatter plots of neural data against behavioral
)

# Indicate type of model estimation
est_type = c(
  T, # lme4 check
  F, # Prior predictive check
  F, # Estimation
  F, # Bridge sampling for Bayes factor
  F, # LOO-CV
  F  # Posterior retrodictive check
)

# Load in useful packages

# Collection of useful functions
# devtools::install_github("rettopnivek/utilityf")
library( utilityf )

# Package for easier manipulation of data frames
# install.packages( 'dplyr' )
library( dplyr )

# Package for estimating mixed effects models
# install.packages( 'lme4' )

# Package for estimating mixed effects models (Bayesian)
# install.packages( 'rstanarm' )
library( rstanarm )
# For parallel processing
options(mc.cores = parallel::detectCores())

# Package for computing Bayes factors using bridge sampling
# install.packages( 'bridge_sampling' )
library( bridgesampling )

# Load in data
setwd( 'Data' )
load( "NBack_3_19_2018.RData" )
setwd( proj_dir )

# Load in useful functions
setwd( 'R' )
source( 'S02_Useful_functions.R' )
setwd( proj_dir )

###
### 2) Data processing
###

# Only analyze 2-back data
cd = dat %>% 
  filter( Task == 'NBack_2' )

# Only analyze subjects included in Meryem's analyses
neural_subj = 
  strsplit( paste( '30 32 36 39 42 43 46 48 49 51 52 53 56',
                   '57 59 60 61 66 70 75 81 83 84 87 90 93 94',
                   '31 35 37 40 45 47 50 54 55 62 63 64 65', 
                   '67 68 71 72 73 76 77 79 80 82 85 89 92 97' ),
              split = ' ' )[[1]]
neural_subj = paste(
  'FN_0', neural_subj, sep = '' )
cd = cd[ cd$ID %in% neural_subj, ]

# Drop rows with missing 
# self-report data for subject 
# FN_041 (66)
# Subject FN_041 did not have any 
# self report data on how high
# for the T1 drug condition
cd = cd %>% 
  filter( Self_report_on_high != '-' )

# Exclude any missing data from fNIRS data
ndc = c( 'R_DLPFC', 'L_DLPFC', 'MPFC', 'R_VLPFC', 'L_VLPFC' )
is_na = apply( cd[,ndc], 1, function(x) any( is.na(x) ) )
cd = cd[ !is_na, ]

# Create new incremental subject labels
cd$Subject = createIncrement( cd$ID )

# Determine count data for misses/correct rejections
cd$Z = cd$Trials - cd$Counts
cd$Y = cd$Counts

# Initialize change scores
cd$R_DLPFC_D = NA
cd$L_DLPFC_D = NA
cd$MPFC_D = NA
cd$R_VLPFC_D = NA
cd$L_VLPFC_D = NA

# Compute change scores for fNIRS data
change_score_function = function( vrb ) {
  
  # Initialize a wide format data frame
  wf = data.frame(
    Subject = sort( unique( cd$Subject ) ),
    stringsAsFactors = F
  )
  # Initialize columns
  wf$ID = 'NA'
  wf$T1_P = NA
  wf$T2_P = NA
  wf$T3_P = NA
  wf$T1_D = NA
  wf$T2_D = NA
  wf$T3_D = NA
  
  # Specify conditions
  tm = c( 'T1_Pre_drug', 'T2_Post_drug', 'T3_Post_drug' )
  cnd = c( 'Placebo', 'Drug' )
  # Loop over subjects
  for ( s in 1:nrow( wf ) ) {
    
    # Track ID numbers
    subj = cd$Subject == wf$Subject[s]
    wf$ID[s] = unique( cd$ID[ subj ] )
    
    # Loop over conditions
    inc = 1
    for ( i in 1:2 ) {
      
      for ( j in 1:3 ) {
        sel = cd$Condition == cnd[i] & 
          cd$Timepoints == tm[j] & subj & 
          cd$Response_type == 'Hits'
        
        if ( any(sel) ) 
          wf[s,inc+2] = cd[sel,vrb]
        inc = inc + 1
      }
      
    }
    
  }
  
  # Loop over conditions
  for ( i in 1:2 ) {
    
    for ( j in 2:3 ) {
      
      sel1 = cd$Condition == cnd[i] & 
        cd$Timepoints == tm[1] & 
        cd$Response_type == 'Hits'
      
      sel2 = cd$Condition == cnd[i] & 
        cd$Timepoints == tm[j] & 
        cd$Response_type == 'Hits'
      
      cd[ sel1, paste( vrb, 'D', sep = '_' ) ] <<- 0
      
      adj = 'P'; if ( cnd[i] == 'Drug' ) adj = 'D'
      if ( j == 2 ) adj2 = 'T2' else adj2 = 'T3'
      
      df_sc = 
        wf[, paste( adj2, adj, sep = '_' ) ] - 
        wf[, paste( 'T1', adj, sep = '_' ) ]
      
      cd[ sel2, paste( vrb, 'D', sep = '_' ) ] <<- df_sc[ cd$Subject[ sel2 ] ]
      
      sel2 = cd$Condition == cnd[i] & 
        cd$Timepoints == tm[j] & 
        cd$Response_type == 'False_alarms'
      
      cd[ sel2, paste( vrb, 'D', sep = '_' ) ] <<- df_sc[ cd$Subject[ sel2 ] ]
      
      sel1 = cd$Condition == cnd[i] & 
        cd$Timepoints == tm[1] & 
        cd$Response_type == 'False_alarms'
      
      cd[ sel1, paste( vrb, 'D', sep = '_' ) ] <<- 0
      
    }
    
  }
  
}

# COmpute difference scores
for ( i in ndc ) change_score_function( i )

# Standardize neural measures
my_standardize = function( M ) {
  v = as.vector( M )
  z = ( v - mean(v ) ) / sd( v )
  return( matrix( z, nrow( M ), ncol( M ) ) )
}
sel = paste( ndc, 'D', sep = '_' )
cd[,sel] = my_standardize( as.matrix( cd[,sel] ) )

# Create variables for constructing 
# design matrices for d' and bias
init_X = create_design_mat( cd )
init_X[,ndc] = cd[,paste(ndc,'D',sep='_')]

###
### 3)
###

# 3.1)
estimate_sdt = function( dtbf, 
                         glmer_formula, 
                         priors, 
                         est_type,
                         seed,
                         niter = 3500,
                         warmup = 1000 ) {
  # Purpose:
  # Estimates an equal-variance SDT model 
  # using multi-level probit regression.
  # Arguments:
  # dtbf          - A data frame with the 
  #                 observed data and predictors
  # glmer_formula - The formula object to pass 
  #                 to the 'glmer' function
  # priors        - A list with a matrix of the 
  #                 means and standard deviations 
  #                 for the normal priors on the 
  #                 group-level regression coefficients,
  #                 and the value of the prior on the 
  #                 trace (sum of the diagonals) for the 
  #                 covariance matrix governing 
  #                 subject-level parameters
  # est_type      - A logical vector, indicating whether 
  #                 1) to estimate the model using the 
  #                    'lme4' package (maximum likelihood)
  #                 2) to estimate the prior predictive 
  #                    check
  #                 3) to estimate the model using the 
  #                    'rstanarm' package (Bayesian)
  #                 4) to compute the log marginal posterior 
  #                    using bridge sampling
  #                 5) to compute the estimate of 
  #                    the Leave-One-Out 
  #                    Cross-Validation metric
  #                 6) to estimate the posterior 
  #                    retrodictive check
  # seed          - The random seed to use for the 
  #                 'stan_glmer' function (for 
  #                 reproducibility)
  # niter         - The number of iterations to run 
  #                 per chain
  # warmup        - The number of warm-up iterations 
  #                 to use for the NUTS sampler
  # Returns:
  # A list with the output of either the 'glmer' function 
  # from the 'lme4' package (out$lme4), the results of the prior 
  # predictive check (out$prior_pc), the output of the 
  # 'stan_glmer' function from the 'rstanarm' package 
  # (out$est), and the results of the posterior predictive 
  # check.
  
  # Initialize output
  out = list( run_time = Sys.time() )
  
  # Check for model formula
  if ( est_type[1] ) {
    
    est = lme4::glmer( glmer_formula, 
                       family = binomial(link = "probit"),
                       data = dtbf,
                       control = 
                         lme4::glmerControl( optimizer="bobyqa", 
                                             optCtrl = 
                                               list( maxfun=2e5 ) )
    )
    
    out$lme4 = est
    
  }
  
  # Prior predictive check
  if ( est_type[2] ) {
    
    est = rstanarm::stan_glmer( glmer_formula,
                                family = 
                                  binomial(link = "probit"),
                                prior = 
                                  normal( priors[[1]][,1],
                                          priors[[1]][,2] ),
                                prior_covariance = 
                                  decov( 1, 1, 1, priors[[2]] ), 
                                data = dtbf,
                                iter = 5000, 
                                chains = 4,
                                cores = 4, 
                                prior_PD = T )
    
    out$prior_pc = posterior_predict(est)
    
  }
  
  # Bayesian estimation
  if ( est_type[3] ) {
    
    est = rstanarm::stan_glmer( glmer_formula,
                                family = 
                                  binomial(link = "probit"),
                                prior = 
                                  normal( priors[[1]][,1],
                                          priors[[1]][,2] ),
                                prior_covariance = 
                                  decov( 1, 1, 1, priors[[2]] ), 
                                data = dtbf,
                                warm = warmup,
                                iter = niter, 
                                chains = 4,
                                cores = 4, 
                                adapt_delta = .9999,
                                seed = seed, 
                                diagnostic_file = 
                                  file.path(tempdir(), "df.csv"))
    
    out$est = est
    out$glmer_formula = glmer_formula
    
  }
  
  # Bridge sampling
  if ( all( est_type[c(3,4)] ) ) {
    
    out$bs = bridge_sampler( est )
    
  }
  
  # LOO-CV
  if ( all( est_type[c(3,5)] ) ) {
    
    out$LOOCV = loo( est )
    
  }
  
  # Posterior predictive check
  if ( all( est_type[c(3,6)] ) ) {
    
    post_pc = posterior_predict( est )
    out$post_pc = post_pc
    
  }
  
  # Compute run time
  out$run_time = Sys.time() - out$run_time
  
  return( out )
}


# 3.2)
quick_pred_plot = function( ppc,
                            plot_yes = T,
                            new_plot = T ) {
  # Purpose:
  # A convenience function to quickly plot the 
  # model predictions for the overall hit and 
  # correct rejection rates.
  # Arguments:
  # ppc      - A matrix with the simulated observations 
  #            over the posterior samples
  # plot_yes - Logical; if true, plots the results
  # new_plot - Logical; if true, generates a new plotting 
  #            window
  
  h = dtbf$Response_type == 'Hits'
  fn = function(x) by( h * (1 - x/dtbf$Trials ) + 
                         (!h)*x/dtbf$Trials, 
                       list( 
                         dtbf$Timepoints,
                         dtbf$Condition,
                         dtbf$Task,
                         dtbf$Response_type ), mean )[1:24]
  ppc_1 = apply( apply( ppc, 1, fn ), 1, 
                 quick_desc )
  
  out = aggregate( h * (1 - dtbf$Z/dtbf$Trials ) + 
                     (!h)*dtbf$Z/dtbf$Trials,
                   list( dtbf$Timepoints,
                         dtbf$Condition,
                         dtbf$Task,
                         dtbf$Response_type ),
                   mean )
  colnames( out ) = c(
    'Timepoints', 'Condition', 
    'Task', 'Response_type', 'Observed' )
  
  out = cbind( out, t( ppc_1 ) )
  
  if ( plot_yes ) {
    
    if ( new_plot ) x11()
    xl = c( 0, 25 ); yl = c( 0, 1 )
    blankPlot( xl, yl )
    customAxes( xl, yl )
    
    horizLines( .5, xl, lwd = 2, lty = 2 )
    vertLines( seq( 3, 21, 3 ) + .5, yl, lwd = 2 )
    
    segments( 1:24, out$`Q2.5%`,
              1:24, out$`Q97.5%`,
              lwd = 2, col = 'grey' )
    errorBars( 1:24, rbind( out$`Q16%`, out$`Q84%` ),
               length = .05, col = 'grey' )
    points( 1:24, out$Mean, pch = 19, col = 'grey' )
    points( 1:24, out$Mode, pch = 17, col = 'grey' )
    
    points( 1:24, out$Observed, pch = 19 )
    
  }
  
  return( out )
}

###
### 4) Model estimation
###

if ( run_code[1] ) {
  
  # Navigate to folder to save posterior estimates
  setwd( proj_dir )
  setwd( 'Data' )
  setwd( 'Posterior_estimates' )
  
  ROI_labels = c( 'R_DLPFC',
                  'L_DLPFC',
                  'MPFC',
                  'R_VLPFC',
                  'L_VLPFC' )
  
  out = list(
    R_DLPFC = NULL,
    L_DLPFC = NULL,
    MPFC = NULL,
    R_VLPFC = NULL,
    L_VLPFC = NULL
  )
  
  for ( i in 1:length( ROI_labels ) ) {
    
    # Criterion
    # Specify separate dummy-coded variables 
    # for each condition and timepoint, 
    # along with coefficients for an 
    # order effect and self-reported high
    Xc = cbind(
      Baseline_2 = init_X$NBack_2,
      Placebo_2_T1 = init_X$T1_2_P,
      Placebo_2_T2 =  init_X$T2_2_P,
      Placebo_2_T3 = init_X$T3_2_P,
      Drug_2_T2 =  init_X$T2_2_D,
      Drug_2_T3 = init_X$T3_2_D,
      Visit_2 = init_X$O, 
      Self_report_high = init_X$SR,
      ROI = init_X[,ROI_labels[i]]
    )
    
    check = create_design_mat( cd, Xc )
    
    # d'
    # For d' variables, weight dummy coded 
    # values by .5 (False alarms) or -.5 (Hits)
    Xd = Xc; cn = colnames( Xd )
    sel = c( grep( 'Visit', cn ),
             grep( 'Self_report', cn ) )
    Xd[,-sel] = Xd[,-sel] * .5
    sel = cd$Response_type == 'Hits'
    Xd[sel] = -1 * Xd[sel]
    
    colnames( Xc ) = paste( 'crt', colnames( Xc ), sep = '_' )
    colnames( Xd ) = paste( 'dp', colnames( Xd ), sep = '_' )
    
    # Define design matrix
    DM = cbind( Xc, Xd )
    
    # Define random effects
    re_pos = c( 1 )
    re_pos = c( re_pos, re_pos + ncol(Xc) )
    
    # Create data to be fitted
    dtbf = cd %>% 
      select( Task, Condition, Timepoints, Response_type, 
              Subject, Trials, Z, Y )
    dtbf = cbind( dtbf, DM )
    # Make sure there are no missing values
    is_na = apply( dtbf, 1, function(x) any( is.na(x) ) )
    if ( sum( is_na ) > 0 ) {
      # Pairwise deletion of missing data
      dtbf = dtbf[ !is_na, ]
    }
    
    m1_priors = list(
      rbind(
        # Criterion
        c( 0.37, 0.15 ), # Baseline (2) *R - 1*
        c( 0.00, 0.15 ), # Placebo (2 - T1 )
        c( 0.00, 0.15 ), # Placebo (2 - T2 )
        c( 0.00, 0.15 ), # Placebo (2 - T3 )
        c( 0.00, 0.15 ), # Drug (2 - T2)
        c( 0.00, 0.15 ), # Drug (2 - T3)
        c( 0.00, 0.15 ), # Visit 2
        c( 0.00, 0.15 ), # Self-reported high
        c( 0.00, 0.15 ), # fNIRS z-scores
        # d'
        c( 2.9, 0.2 ), # Baseline (2) *R - 7*
        c( 0.0, 0.2 ), # Placebo (2 - T1 )
        c( 0.0, 0.2 ), # Placebo (2 - T2 )
        c( 0.0, 0.2 ), # Placebo (2 - T3 )
        c( 0.0, 0.2 ), # Drug (2 - T2)
        c( 0.0, 0.2 ), # Drug (2 - T3)
        c( 0.0, 0.2 ), # Visit 2
        c( 0.0, 0.2 ), # Self-reported high
        c( 0.0, 0.2 ) # fNIRS z-scores
      ),
      # Prior on the trace (sum of the diagonal)
      # for the subject level variabilities
      prior_tau = sqrt( .2*4 ) )
    
    glmer_formula = paste(
      # Dependent variable
      'cbind( Z, Y ) ~ -1 + ',
      # Group-level parameters
      paste( colnames( DM ), collapse = ' + ' ),
      # Subject-level parameters
      ' + ( -1 + ',
      paste( colnames( DM )[re_pos], collapse = ' + ' ),
      ' |Subject)',
      sep = '' )
    glmer_formula = as.formula( glmer_formula )
    
    m1 = estimate_sdt( dtbf, glmer_formula,
                       m3_priors, est_type )
    
    out[[i]] = m1$lme4
    
  }
  
  setwd( proj_dir )
  
}

###
### 5) Scatter plots of neural data against behavioral
###

# 5.1)
neural_scatter = function( vrb, lbl, new = T ) {
  # Purpose:
  # Generates a set of scatter plots comparing the 
  # standardized scores for the BOLD response of 
  # the 5 ROIs and a specified variable.
  # Arguments:
  # vrb - The y-axis variable name in the data frame
  # lbl - The label to use for the y-axis variable
  # new - Logical; if true, a new plot is generated
  
  for ( nd in paste( ndc, 'D', sep = '_' ) ) {
    
    ndl = strsplit( nd, split = '_' )[[1]]
    ndl = paste( ndl[1], '. ', ndl[2], sep = '' )
    
    if ( new ) x11( width = 12 )
    layout( cbind( 1, 2 ) )
    
    tm = c( 'T2_Post_drug', 'T3_Post_drug' )
    
    yz = my_standardize( cbind( cd[ cd$Response_type == 'Hits', vrb ] ) )
    yz = rep( as.vector( yz ), each = 2 )
    for ( i in 1:2 ) {
      
      # Isolate conditions of interest
      sel = cd$Timepoints == tm[i] & 
        cd$Response_type == 'Hits'
      
      # Separate drug and placebo conditions
      clr = rep( 'black', sum( sel ) )
      clr[ cd$Condition[sel] == 'Placebo' ] = 'grey60'
      
      # Create a blank plot
      xl = c( -5.5, 5.5 )
      yl = c( -5.5, 5.5 )
      blankPlot( xl, yl )
      
      # Compute correlation
      dtbf = data.frame(
        x = cd[sel,nd],
        y = yz[sel] )
      lmf = lm( y ~ -1 + x, data = dtbf )
      
      pid = data.frame( x = seq(-5.5,5.5,length=100),
                        y = 0 )
      pid$y = predict( lmf, newdata = pid, interval = 'predict' )
      polygon( c( pid$x, rev( pid$x ) ),
               c( pid$y[,2], rev( pid$y[,3] ) ),
               col = 'grey80', border = NA )
      
      
      segments( xl[1], coef( lmf )*xl[1],
                xl[2], coef( lmf )*xl[2],
                lwd = 2 )
      
      # Plot observations
      points( cd[sel,nd], yz[sel], pch = 19, col = clr, cex = ptSz )
      
      # Add axes and labels
      customAxes( xl, yl, lnSz = lnSz )
      axis( 1, seq( -5, 5, 1 ),
            tick = F, line = axPos, cex.axis = axSz )
      mtext( paste( 'z-scores for ROI (', ndl, ')', sep = '' ),
             side = 1, line = 2, cex = axSz )
      axis( 2, seq( -5, 5, 1 ),
            tick = F, line = axPos, cex.axis = axSz )
      mtext( paste( "z-scores for", lbl, "values" ),
             side = 2, line = 2, cex = axSz )
      
      if ( tm[i] == 'T2_Post_drug' ) {
        title( 'Post drug' )
      } else {
        title( 'Post drug (2)' )
      }
      
      sm = summary( lmf )
      legend( 'topright', 
              paste( 'R = ', 
                     round( sm$coefficients[1,1], 2 ),
                     ', p = ', 
                     round( sm$coefficients[1,4], 3 ),
                     sep = '' ), 
              cex = axSz, bty = 'n' )
      
      if ( i == 1 ) {
        legend( -5, 6, c( 'Placebo', 'Drug' ),
                fill = c( 'grey60', 'black' ),
                cex = axSz, bty = 'n' )
      }
      
    }
    
  }
  
}


if ( run_code[2] ) {
  
  # Compute descriptive d' and bias values
  sel = cd$Response_type == 'Hits'
  fH = cd$Counts[sel]
  Npos = cd$Trials[sel]
  sel = cd$Response_type == 'False_alarms'
  fFA = cd$Counts[sel]
  Nneg = cd$Trials[sel]
  dp = sdt_calc_binary( fH, Npos, fFA, Nneg, 1, T )
  crt = sdt_calc_binary( fH, Npos, fFA, Nneg, 1, F )
  sel = cd$Response_type == 'Hits'
  cd$dp[sel] = dp; cd$crt[sel] = crt
  sel = cd$Response_type == 'False_alarms'
  cd$dp[sel] = dp; cd$crt[sel] = crt
  
  # Plot characteristics
  lnSz = 2
  axPos = -1.5
  axSz = 1
  ptSz = 1.5
  
  setwd( proj_dir )
  setwd( 'Figures' )
  
  pdf( file = 'Scatter_plots_for_ROIs_and_dprime.pdf',
       width = 12 )
  neural_scatter( 'dp', "d'", new = F )
  dev.off()
  
  pdf( file = 'Scatter_plots_for_ROIs_and_self_report.pdf',
       width = 12 )
  neural_scatter( 'Self_report_on_high', "self-reported high", new = F )
  dev.off()
  
  setwd( proj_dir )
  
}

setwd( orig_dir )

