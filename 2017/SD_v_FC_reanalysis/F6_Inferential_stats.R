#-----------------------------#
# Inferential statistics with # 
# observed response times     #
# Kevin Potter                #
# Updated 08/21/2017          #
#-----------------------------#

# Initialize script
source('F3_starting_script.R')

# Load in package for Bayes factor tests
# install.packages( BayesFactor )
library( BayesFactor )

# Indicate which code segments to run
runCode = c( F, F, T, T, F )

aqs = function( string, obj ) {
  # Purpose: 
  # A function to allow for quick selection of a subset 
  # of a data frame.
  # Arguments: 
  # string - A string of shorthand conditions and the desired 
  #          subset, separated by commas (e.g., 'Co_L,PT_FP')
  # obj    - A data frame with 5 grouping columns
  # Returns: 
  # A logical vector.
  
  # Convert labels of current object to shorthand values
  SH = matrix( ' ', nrow( obj ), 5 )
  colnames( SH ) = c( 'PD', 'PT', 'Co', 'Ta', 'Ac' )
  
  SH[,1] = 'S'; SH[,2] = 'FP'; SH[,3] = 'L'; 
  SH[,4] = 'FC'; SH[,5] = '1'
  SH[ obj$PD == .4, 1 ] = 'L';
  SH[ obj$PT == 'Target-primed', 2 ] = 'TP';
  SH[ obj$Co == 'Right', 3 ] = 'R';
  SH[ obj$Co == 'Same', 3 ] = 'S';
  SH[ obj$Co == 'Different', 3 ] = 'D';
  SH[ obj$Ta == 'Same-different', 4 ] = 'SD';
  SH[ obj$Ac == 0, 5 ] = '0';
  
  vrb = strsplit( string, split = ',' )[[1]]
  
  K = length( vrb )
  which_true = matrix( NA, nrow(obj), K )
  for ( k in 1:K ) {
    
    val = strsplit( vrb[k], '_' )[[1]]
    which_true[,k] = SH[,val[1]] == val[2]
  }
  
  out = as.vector( apply( which_true, 1, all ) )
  
  return( out )
}

###
### Tests on proportion correct
###
# Lookup - 01

if ( runCode[1] ) {
  
  # Compute frequency for accuracy over conditions
  ac = aggregate( cbind( P = d$Ac, N = rep(1,nrow(d) ) ), 
                  list( d$PDL, d$PTL, d$CoL, d$TaL, d$S ), 
                  sum )
  colnames( ac ) = c( 'PD', 'PT', 'Co', 'Ta', 'S', 'F', 'N' )
  # Compute proportion correct
  ac$P = ac$F/ac$N
  # Apply correction for 100% accuracy
  ac$Pc = ( ac$F + .5 )/( ac$N + 1 )
  # Convert to unbounded measure via logit transform
  ac$lP = logit( ac$Pc )
  # Convert covaraties to factors
  ac$S = as.factor( ac$S )
  ac$Co = as.factor( ac$Co )
  ac$PD = as.factor( ac$PD )
  ac$PT = as.factor( ac$PT )
  ac$Ta = as.factor( ac$Ta )
  
  # Question:
  # Is there a difference in accuracy performance between 
  # tasks?
  
  # Fit repeated-measures ANOVA to logit transformed accuracy
  # testing all possible combinations of main effects and 
  # interactions for Task x Prime duration x Prime type
  # against the null model of a random effect for subjects.
  bf_task = anovaBF( lP ~ PD*PT*Ta + S, data = ac, 
                     whichRandom="S" )
  
  # Plot Bayes factors
  x11( width = 12 ); plot( bf_task )
  
  # A model with main effects for prime duration and type and their
  # interaction is preferred over all other fixed effects.
  
  # The next most likely model includes a main effect for task.
  # The ratio of the Bayes factors is:
  bf_task[4] / bf_task[9]
  
  # Question:
  # Are there differences in accuracy performance based on the 
  # type of correct answer (Left versus right)?
  
  # Fit repeated-measures ANOVA to logit transformed accuracy
  # for the forced-choice task only testing all possible 
  # combinations of main effects and interactions for 
  # Correct x Prime duration x Prime type against the null 
  # model of a random effect for subjects.
  ac_fc = ac[ ac$Ta == 'Forced-choice', ]
  bf_fc = anovaBF( lP ~ PD*PT*Co + S, data = ac_fc, 
                     whichRandom="S" )
  
  # Plot Bayes factors
  x11( width = 12 ); plot( bf_fc )
  
  # The most likely model only includes main effects and an 
  # interaction for prime duration and type. The next most 
  # likely model includes a main effect for position correct.
  # The ratio of the Bayes factors is:
  bf_fc[4] / bf_fc[9]
  
  # Question:
  # Are there differences in accuracy performance based on the 
  # type of correct answer (Same versus different)?
  
  # Fit repeated-measures ANOVA to logit transformed accuracy
  # for the same-different task only testing all possible 
  # combinations of main effects and interactions for 
  # Correct x Prime duration x Prime type against the null 
  # model of a random effect for subjects.
  ac_sd = ac[ ac$Ta == 'Same-different', ]
  bf_sd = anovaBF( lP ~ PD*PT*Co + S, data = ac_sd, 
                   whichRandom="S" )
  
  # Plot Bayes factors
  x11( width = 12 ); plot( bf_sd )
  
  # The most likely model includes main effects and an 
  # interaction for prime duration and type and a main 
  # effect for position correct. We can compare this 
  # against a model with no effect of position correct.
  # The ratio of the Bayes factors is:
  bf_sd[9] / bf_sd[4]
  
}

###
### Tests on median RTs
###
# Lookup - 02

if ( runCode[2] ) {
  
  # Compute median RT over relevant conditions
  md_rt = aggregate( d$RT, list( d$PDL, d$PTL, d$CoL, d$TaL,  
                                 d$Ch, d$Ac, d$S ), median )
  colnames( md_rt ) = c( 'PD', 'PT', 'Co', 'Ta', 'Ch', 
                         'Ac', 'S', 'MRT' )
  md_rt$lMRT = log( md_rt$MRT )
  for ( i in 1:7 ) md_rt[,i] = as.factor( md_rt[,i] )
  
  # Question:
  # Is there a difference between median response times 
  # based on whether the subject picked the left or 
  # right choice for the 2AFC task?
  
  ### For correct choices ###
  
  # Fit repeated-measures ANOVA to log of the median response 
  # times for the forced-choice task only testing all possible 
  # combinations of main effects and interactions for 
  # Choice x Prime duration x Prime type against the null 
  # model of a random effect for subjects.
  sel = md_rt$Ta == 'Forced-choice' & md_rt$Ac == 1
  dtbf = md_rt[ sel, ]
  bf_fc = anovaBF( lMRT ~ PD*PT*Ch + S, data = dtbf,  
                   whichRandom = "S" )
  
  # Plot Bayes factors
  x11( width = 12 ); plot( bf_fc )
  
  # The most likely model only includes main effects and an 
  # interaction for prime duration and type. The next most 
  # likely model includes a main effect for type of choice.
  # The ratio of the Bayes factors is:
  bf_fc[4] / bf_fc[9]
  
  ### For incorrect choices ###
  
  # Fit repeated-measures ANOVA to log of the median response 
  # times for the forced-choice task only testing all possible 
  # combinations of main effects and interactions for 
  # Choice x Prime duration x Prime type against the null 
  # model of a random effect for subjects.
  sel = md_rt$Ta == 'Forced-choice' & md_rt$Ac == 0
  dtbf = md_rt[ sel, ]
  bf_fc = anovaBF( lMRT ~ PD*PT*Ch + S, data = dtbf,  
                   whichRandom = "S" )
  
  # Plot Bayes factors
  x11( width = 12 ); plot( bf_fc )
  
  # The most likely model only includes a main effects of 
  # prime type. We can compare that against the most likely 
  # model that includes an effect fo choice: a main effect 
  # of prime type and choice.
  # The ratio of the Bayes factors is:
  bf_fc[2] / bf_fc[7]
  
  # Question:
  # Are response times for the same-different task slower 
  # than those from the 2AFC task?
  
  # Compute median RT over relevant conditions
  md_rt = aggregate( d$RT, list( d$TaL, d$Ac, d$S ), median )
  colnames( md_rt ) = c( 'Ta', 'Ac', 'S', 'MRT' )
  md_rt$lMRT = log( md_rt$MRT )
  for ( i in 1:3 ) md_rt[,i] = as.factor( md_rt[,i] )
  
  # Fit repeated-measures ANOVA to log of the median response 
  # times testing all possible combinations of main effects 
  # and interactions for Task x Accuracy against the null 
  # model of a random effect for subjects.
  bf_task = anovaBF( lMRT ~ Ta*Ac + S, data = md_rt, 
                   whichRandom = "S" )
  
  # Plot Bayes factors
  x11( width = 12 ); plot( bf_task )
  
}

###
### Tests on MLE for thresholds
###
# Lookup - 03

if ( runCode[3] ) {
  
  # Change directory to where model results are saved
  setwd( orig_dir )
  setwd( 'RT_and_choice_MLE' )
  setwd( 'Estimation_results' )
  
  # Load in parameter estimate
  load( 'DRM_SD_M7_PM_results.RData' )
  setwd( orig_dir )
  
  # Extract thresholds
  kappa = data.frame( S = rep(1:N,8) )
  kappa$kappa = NA
  for ( i in 1:8 ) kappa$kappa[ 1:N + N * (i-1) ] = all_prm[,i]
  # Add covariates
  kappa$PDL = rep( rep( c(.05,.4),4 ), each = N )
  kappa$PM = rep( rep( 
    c('Non-Match','Match','Match','Non-Match'), each = 2 ), each = N )
  kappa$R = rep( 
    rep( c('Same','Different'), each = 4 ), each = N )
  kappa$Cnd = rep( 1:8, each = N )
  dtbf = kappa
  
  # Plot predictions
  x11();
  tmp = aggregate( dtbf$kappa, list( dtbf$Cnd ), mean )
  xl = c( .5, 8.5 ); yl = c( .5, 2 ); blankPlot( xl, yl )
  customAxes( xl, yl );
  pts = c( 21, 21, 24, 24, 24, 24, 21, 21 )
  clr = rep( c( 'white', 'black' ), each = 4 )
  points( tmp[,1], tmp[,2], pch = pts, bg = clr )
  
  # Axes
  axis( 1, 1:8, rep( c(.05, .4 ), 4 ), 
        tick = F, cex = 1, line = -1.5 )
  axis( 1, seq(1,7,2)+.5, c( 'Non-match', 'Match', 'Match', 'Non-match' ), 
        tick = F, cex = 1, line = -.5 )
  axis( 1, seq(2,6,4)+.5, c('Same','Different'), 
        tick = F, cex = 1, line = .5 )
  axis( 2, seq(.5,2,.5), tick = F, line = -1.5, cex = 1.2 )
  mtext( 'Mean threshold', side = 2, line = 2, cex = 1.5 )
  
  # Create specific contrasts
  dtbf$SDB = -1 # Same
  dtbf$SDB[ dtbf$Cnd > 4 ] = 1 # Different
  dtbf$SDLD = 0 # No discounting for short durations
  # Biased against test word that matches with prime
  dtbf$SDLD[ dtbf$Cnd %in% c(4,6) ] = -1
  # Biased toward test word that does not match with prime
  dtbf$SDLD[ dtbf$Cnd %in% c(2,8) ] = 1
  dtbf$S = as.factor( dtbf$S )
  # Interaction
  dtbf$SDBxSDLD = dtbf$SDB * dtbf$SDLD
  
  ### Test three models ###
  
  # Null model
  # Random effect of subject
  Null = lmBF( kappa ~ S, whichRandom = 'S', data = dtbf )
  
  # Generate model predictions
  post = posterior( Null, iterations = 10000 )
  prm = colMeans( post )['mu']
  X = cbind( rep( 1, 8 ) )
  pred = X %*% cbind( prm )
  lines( 1:8, pred, col = 'purple' )
  
  # Model 1
  # Same-different bias
  # (1-4) vs. (5:8)
  M1 = lmBF( kappa ~ SDB + S, whichRandom = 'S', data = dtbf )
  
  # Generate model predictions
  post = posterior( M1, iterations = 10000 )
  prm = colMeans( post )[c('mu','SDB')]
  X = cbind( 1, aggregate( dtbf[,c('SDB')], list( dtbf$Cnd ), unique )[,-1] )
  pred = X %*% cbind( prm )
  lines( 1:8, pred, col = 'blue' )
  
  # Model 2
  # Same different bias
  # (1-4) vs. (5:8)
  # Strategic discounting with long durations
  # (2,8) vs. (4,6)
  M2 = lmBF( kappa ~ SDB + SDLD + SDBxSDLD + S, whichRandom = 'S', data = dtbf )
  
  # Generate model predictions
  post = posterior( M2, iterations = 10000 )
  prm = colMeans( post )[ c('mu','SDB','SDLD','SDBxSDLD') ]
  prm = as.vector( prm )
  X = cbind( 1, aggregate( dtbf[,c('SDB','SDLD')], list( dtbf$Cnd ), unique )[,-1] )
  X = cbind( X, X[,2] * X[,3] )
  X = as.matrix( X )
  colnames( X ) = c( 'I', 'SDB', 'SDLD', 'SDBxSDLD' )
  pred = X %*% cbind( prm )
  lines( 1:8, pred, col = 'orange' )
  
  # Credible intervals around regression coefficients
  round( apply( post[,c('SDB','SDLD','SDBxSDLD')], 2, 
                quantile, prob = c(.025,.5,.975) ), 3 )
  
  # Legend
  legend( 'topright', 
          c( 'Null model', 'Same-different bias', 'Strategic discounting' ),
          fill = c( 'purple', 'blue', 'orange' ), bty = 'n' )
  
  # Model comparison
  print( M1/M2 )
  # Evidence is equivocal - neither for or against.
  
  tmp = strsplit( capture.output( print( M2/M1 ) )[3], split = ':' )[[1]][2]
  tmp = as.numeric( strsplit( tmp, split = ' ' )[[1]][2] )
  text( 6, .6, paste( 'Bayes factor ratio: ', round( tmp, 2 ) ) )
  
  
}

###
### Tests on MLE for drift rates
###
# Lookup - 04

if ( runCode[4] ) {
  
  # Change directory to where model results are saved
  setwd( orig_dir )
  setwd( 'RT_and_choice_MLE' )
  setwd( 'Estimation_results' )
  
  # Load in parameter estimate
  load( 'DRM_SD_M7_PM_results.RData' )
  setwd( orig_dir )
  
  # Extract drift rates
  tmp = all_prm[,grep('xi',colnames(all_prm))]
  xi = data.frame( S = rep(1:N,16) )
  xi$xi = NA
  for ( i in 1:16 ) xi$xi[ 1:N + N * (i-1) ] = tmp[,i]
  # Add covariates
  xi$Co = rep( 
    rep( c('Same','Different', 'Different', 'Same'), each = 4 ), each = N )
  xi$PDL = rep( rep( c(.05,.4), 8 ), each = N )
  xi$PTL = rep( rep( rep( c('Foil','Target'), each = 2 ), 4 ), each = N )
  xi$R = rep( rep( c('Same','Different'), each = 8 ), each = N )
  xi$Cnd = rep( c(1:8,5:8+8,1:4+8), each = N )
  dtbf = xi
  
  # Plot predictions
  x11();
  tmp = aggregate( dtbf$xi, list( dtbf$Cnd ), mean )
  xl = c( .5, 8.5 ); yl = c( 0, 5 ); blankPlot( xl, yl )
  segments( c( 2.5, 4.5, 6.5 ), rep( yl[1], 3 ),
            c( 2.5, 4.5, 6.5 ), rep( yl[2], 3 ),
            col = 'grey', lwd = 2 )
  customAxes( xl, yl );
  
  # Add points
  pts = rep( c(22,21), each = 8 )
  clr = rep( c( 'white', 'black' ), each = 8 )
  points( 1:8 - .1, tmp$x[1:8], pch = pts[1:8], bg = clr[1:8] )
  points( 1:8 + .1, tmp$x[9:16], pch = pts[9:16], bg = clr[9:16] )
  
  # Axes
  axis( 1, 1:8, rep( c(.05, .4 ), 4 ), 
        tick = F, cex = 1, line = -1.5 )
  axis( 1, seq(1,7,2)+.5, rep( c('Foil', 'Target' ), 2 ), 
        tick = F, cex = 1, line = -.5 )
  axis( 1, seq(2,6,4)+.5, c('Same','Different'), 
        tick = F, cex = 1, line = .5 )
  axis( 2, seq(0,5,1), tick = F, line = -1.5, cex = 1.2 )
  mtext( 'Mean drift rate', side = 2, line = 2, cex = 1.5 )
  
  # Difference between primed and unprimed long duration 'same'
  # drift rates
  x = xi$xi[ xi$PDL == .4 & xi$PTL == 'Target' & xi$R == 'Same' ] - 
    xi$xi[ xi$PDL == .4 & xi$PTL == 'Foil' & xi$R == 'Same' ]
  # Positve indicates primed drift rates are higher, whereas 
  # negative drift rates indicate unprimed drifts are higher
  bf = ttestBF(x,nullInterval = c(0,Inf))
  print( 1/bf )
  
  # Predict 'different' drift rates from 'same' drift rates
  ord = order( xi$Cnd )
  dtbf = xi[ ord, ]
  tmp = dtbf$xi[ dtbf$R == 'Same' ]
  dtbf = dtbf[ dtbf$R == 'Different', ]
  dtbf$IV = tmp
  dtbf$IV2 = -1; dtbf$IV2[ dtbf$Co == 'Different' ] = 1
  dtbf$IV1x2 = dtbf$IV * dtbf$IV2
  dtbf$S = as.factor( dtbf$S )
  
  x11();
  xl = c(.5,8.5); yl = c(0,5); blankPlot( xl, yl )
  tmp = aggregate( dtbf[,c('xi','IV','IV2','IV1x2')], list( dtbf$Cnd ), mean )
  colnames( tmp ) = c('Cnd','xi','IV','IV2','IV1x2')
  points( 1:8, tmp$xi, pch = 19 )
  points( 1:8, tmp$IV, pch = 22, bg = 'white' )
  
  # Model 1
  M1 = lm( xi ~ 1 + IV, data = dtbf )
  prm = coef( M1 )
  X = cbind( 1, tmp$IV )
  pred = X %*% cbind( prm )
  lines( 1:8, pred, col = 'blue' )
  
}
