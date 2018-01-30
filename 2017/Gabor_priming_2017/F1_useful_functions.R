#--------------------#
# Useful functions   #
# Kevin Potter       #
# Updated 12/17/2017 #
#--------------------#

# Index
# Lookup - 01:  horizLines
# Lookup - 02:  tableContents
# Lookup - 03:  Colorblind palette
# Lookup - 04:  quickBar
# Lookup - 05:  quickCoverage
# Lookup - 06:  simplify_variable_names
# Lookup - 07:  vertLines
# Lookup - 08:  quickCAF
# Lookup - 09:  convergenceExtract
# Lookup - 10:  plotConvergence
# Lookup - 11:  createConditionIndex
# Lookup - 12:  
# Lookup - 13:  

# Lookup - 01
horizLines = function( yval, xl, ... ) {
  # Purpose: 
  # Draws a set of horizontal lines onto 
  # an existing plot.
  # Arguments: 
  # yval - A vector of y-axis values at which to 
  #        draw the lines.
  # xl   - A vector giving the lower and upper bounds 
  #        of the x-axis.
  
  l = length( yval )
  
  segments( rep( xl[1], l ), yval,
            rep( xl[2], l ), yval, ... )
  
}

# Lookup - 02
tableContents = function( string_vector, txtSz = 1.5, 
                          shift = NULL, ... ) {
  # Purpose: 
  # Creates a table of contents at the start of a 
  # PDF file with figures.
  # Arguments: 
  # string_vector - A vector of each line for the table 
  #                 table of contents
  # txtSz         - Text size
  # shift         - Optional shift parameter to page numbers
  # ...           - Additional parameters for the legend
  
  blankPlot()
  
  nLines = length( string_vector )
  
  if ( is.null( shift ) ) shift = rep( 0, nLines )
  
  shift = shift + 1;
  
  all_text = paste( 'Page ', 1:nLines + shift, '.......  ', string_vector,
                    sep = '' )
  
  title( 'Table of contents', cex = txtSz * 1.2 )
  legend( 0, 1, all_text, bty = 'n', cex = txtSz, ... )
  
}

# Lookup - 03
# Colorblind palette
cbbPalette = c( "#000000", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )

# Lookup - 04
quickBar = function( pos, w, h, ... ) {
  # Purpose: 
  # Generates a single bar for a custom bar graph.
  # Arguments: 
  # pos - The center of the x-axis position for the bar
  # w   - The width of the bar
  # h   - The height of the bar (from 0)
  # ... - Additional parameters for the 'polygon' function
  
  x = pos + c( -w/2, -w/2, w/2, w/2 )
  y = c( 0, h, h, 0 )
  polygon( x, y, ... )
  
}

# Lookup - 05
quickCoverage = function( coverage_width ) {
  # Purpose: 
  # Generates the lower and upper boundaries for an 
  # uncertainty interval given a coverage width.
  # Arguments: 
  # coverage_width - The coverage width for the uncertainty interval
  # Returns
  # A vector with the lower and upper boundaries.
  
  interval = numeric(2)
  interval[1] = (1 - coverage_width)/2
  interval[2] = interval[1] + coverage_width
  
  return( interval )
}

# Lookup - 06
simplify_variable_names = function( df ) {
  # Purpose: 
  # Replaces the variable names in the data frames 
  # with the experiment data with shorter, easier 
  # to type names.
  # Arguments: 
  # df
  # Returns: 
  # A character vector with the updated variable names. 
  
  original_names = c(
    "Subject", "Trial", "Block", "Block_type",
    "Foil_contrast", "Target_contrast", 
    "Target_to_foil_ratio", "Max_illuminance", 
    "Angle", "Fixation_time", "Placeholder_time", 
    "Prime_time", "Target_time", "Mask_time",
    "Feedback_time", "Prime_duration", "Prime_rotation", 
    "Prime_type",  "Target_rotation", "Correct_answer", 
    "Choice_rotation", "Choice_display", "RT", "Choice", 
    "Accuracy", "Choice_label", "Prime_rotation_label", 
    "Target_rotation_label", "Choice_rotation_label", 
    "Prime_type_label", "Correct_answer_label", 
    "Onscreen_choices_label", "Block_type_label", 
    "Accuracy_label", "Experiment", "ID", 
    "Full_index", "ChxDurxPri", "DurxPri" )
  
  new_names = c(
    "S", "Tr", "B", "BT",
    "FC", "TC", 
    "ToR", "MI", 
    "A", "FiTi", "PlTi", 
    "PTi", "TTi", "MTi",
    "FeTi", "PD", "PR", 
    "PT",  "TR", "Co", 
    "ChR", "ChD", "RT", "Ch", 
    "Ac", "ChL", "PRL", 
    "TRL", "ChRL", 
    "PTL", "CoL", 
    "OCL", "BTL", 
    "AcL", "Exp", "ID", 
    "FI", "CDP", "DP" )
  
  variable_names = colnames( df )
  
  # Make sure the variable names match the original names
  if ( !all( original_names %in% variable_names ) ) 
    stop( 'Mismatch between data frame column names and internal labels',
          call. = FALSE )
  
  # Loop over original names, replacing them with 
  # new names
  for ( i in 1:length( original_names ) ) {
    sel = variable_names == original_names[i]
    variable_names[sel] = new_names[i]
  }
  
  return( variable_names )
  
}

# Lookup - 07
vertLines = function( xval, yl, ... ) {
  # Purpose: 
  # Draws a set of vertical lines onto 
  # an existing plot.
  # Arguments: 
  # xval - A vector of x-axis values at which to 
  #        draw the lines.
  # yl   - A vector giving the lower and upper bounds 
  #        of the y-axis.
  
  l = length( xval )
  
  segments( xval, rep( yl[1], l ),
            xval, rep( yl[2], l ), ... )
  
}

# Lookup - 08
quickCAF = function( s, alpha = pnorm(2) - pnorm(-2),
                     axSz = 1.2, ptSz = 1.5, lnSz = 2 ) {
  # Purpose: 
  # Given a subject index, plots an estimate of the conditional 
  # accuracy function for a subject, based on a thin plate spline 
  # regression model.
  # Arguments: 
  # s     - The subject index
  # alpha - The coverage width for the prediction interval
  # axSz  - The size of the axis labels
  # ptSz  - The size of the figure points
  # lnSz  - The width of the figure lines
  
  # Consider main trials only
  sel = rd$BT == 2 & rd$S == s
  dtbf = rd[sel,]
  
  # Sort by RT
  ord = order( dtbf$RT )
  dtbf = dtbf[ ord, ]
  
  # Fit a generalized additive model predicting accuracy 
  # from RT using thin plate splines for smoothing
  # (Since only one regressor, fitting a basic spline model)
  gam_fit = mgcv::gam( Ac ~ s( RT ), data = dtbf, 
                       family = binomial(link = 'logit' ) )
  
  # Extract parameters
  beta = coef( gam_fit )
  Vb = vcov( gam_fit )
  
  # Predicted accuracy based on RT
  pred = predict( gam_fit, type = 'lpmatrix' )
  # Apply smoothing function
  s = pred %*% beta
  # Convert to predicted accuracy
  y = logistic( s )
  
  # Prediction interval (This can be slow)
  
  # Simulate replicate beta vectors from posterior
  Cv = chol( Vb )
  nRep = 10000
  nb = length( beta )
  br = t( Cv ) %*% matrix( rnorm( nRep * nb ), nb, nRep ) + beta
  
  s_sim = pred %*% br
  y_sim = apply( s_sim, 1, logistic )
  
  # Compute prediction interval via quantiles
  interval = quickCoverage( alpha )
  pi = apply( y_sim, 2, quantile, prob = interval )
  
  # Create a blank plot
  xl = lowerUpper( .2, dtbf$RT ); yl = c( 0, 1 )
  blankPlot( xl, yl )
  
  # Prediction interval
  xa = c( dtbf$RT, rev( dtbf$RT ) )
  ya = c( pi[1,], rev( pi[2,] ) )
  polygon( xa, ya, border = NA, col = 'grey80' )
  
  # Add predictions for best-fitting line
  lines( dtbf$RT, y, lwd = lnSz )
  
  # Add empirical estimate of CAF
  bins = quantile( dtbf$RT, seq( 0, 1, .1 ) )
  p = sapply( bins, function(x) sum( dtbf$Ac[ dtbf$RT < x ] ) )
  n = sapply( bins, function(x) sum( dtbf$RT < x ) )
  p[1] = 0; caf = diff( p )/diff( n )
  points( bins[-1], caf, pch = 19, cex = ptSz )
  
  # Add grid lines and axes
  horizLines( .5, xl, lwd = lnSz, lty = 2 )
  vertLines( .2, yl, lwd = lnSz, lty = 2 )
  customAxes( xl, yl, lnSz = lnSz )
  axis( 1, round( seq( xl[1], xl[2], length = 5 ), 3 ),
        tick = F, line = -.8, cex.axis = axSz )
  axis( 2, seq( 0, 1, .25 ), 
        tick = F, line = -.8, cex.axis = axSz )
  
}

# Lookup - 09
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
  
  Rhat = summary(fit)$summary[,"Rhat"]
  n_eff = summary(fit)$summary[,"n_eff"]
  totSampSize = length(extract(fit, pars = "lp__")[[1]])
  # We're only interested in a subset of parameters
  if ( length( parName ) == 0 ) 
    parName = names( Rhat )[ 
      which( names( Rhat ) == "logLik[1]" ) - 1 ]
  sel = 1:max( grep( parName, names(Rhat) ) )
  Rhat = Rhat[sel]; n_eff = n_eff[sel];
  
  return( list( Rhat = Rhat, n_eff = n_eff, totSampSize = totSampSize ) )
}

# Lookup - 10
plotConvergence = function( conv, savePlot, parName = NULL ) {
  # Purpose:
  # Generates a plot of the Gelman-Rubin statistics 
  # and the effective number of samples for the 
  # marginal posterior samples of the parameters.
  # Arguments:
  # fit      - A stan fit object
  # savePlot - A logical value, which when false generates a 
  #            new plotting window
  # parName  - An optional string giving the final parameter label 
  #            of the subset of the output to include
  
  # Remove NaN for correlation matrix
  sel = which( is.na( conv$Rhat ) )
  if ( length( sel ) == 0 ) sel = -(1:length(conv$Rhat))
  conv$Rhat = conv$Rhat[ -sel ]
  conv$n_eff = conv$n_eff[ -sel ]
  
  if (!savePlot) x11(width=12)
  layout( cbind(1,2) )
  
  # Plot a histogram of the Gelman-Rubin statistics for the 
  # marginal posterior samples of the parameters
  
  tmp = hist( conv$Rhat, plot = F )
  scl = findDec( tmp$density )
  if ( scl[4] == 1 ) scl = scl[1]*(scl[2]/10) else 
    scl = scl[1]/(scl[2])
  
  yl = lowerUpper( scl, tmp$density )
  yl[1] = 0
  
  xl = lowerUpper( .1, conv$Rhat )
  xl[2] = max( xl[2], 1.12 )
  xl[1] = min( xl[1], .98 )
  
  plot( xl, yl, type = 'n', cex.axis = 1.5, cex.lab = 1.5,
        xlab = expression( hat(R) ), ylab = 'Density',
        bty = 'l', main = 'Gelman-Rubin Statistic' )
  
  segments( tmp$mids, rep(0,length(tmp$mids)),
            tmp$mids, tmp$density, lwd = 3,
            col = 'grey' )
  abline( v = 1.1, lty = 2, lwd = 2 )
  
  # Plot a histogram of the effective sample size for the 
  # set of parameters
  
  tmp = hist( conv$n_eff, plot = F )
  scl = findDec( tmp$density )
  if ( scl[4] == 1 ) scl = scl[1]*(scl[2]/10) else 
    scl = scl[1]/(scl[2])
  
  yl = lowerUpper( scl, tmp$density )
  yl[1] = 0
  
  xl=c(0,conv$totSampSize)
  plot( xl, yl, type = 'n', cex.axis = 1.5, cex.lab = 1.5,
        xlab = expression(N[sample]), ylab = 'Density',
        bty = 'l', main = 'Effective sample size' )
  
  segments( tmp$mids, rep(0,length(tmp$mids)),
            tmp$mids, tmp$density, lwd = 3,
            col = 'grey' )
  
}

# Lookup - 11
createConditionIndex = function( varNames, dat ) {
  # Purpose: 
  # Creates an index of the unique level combinations 
  # given a set of variables in a data frame.
  # Arguments: 
  # varNames - A character vector of variable names
  # dat      - A data frame
  # Returns: 
  # A vector of the indices.
  
  # Create list of input variables
  lst = lapply( varNames, function(x) dat[[x]] )
  names( lst ) = varNames
  
  # Create index of conditions
  cnd = aggregate( rep( 1, nrow(dat) ), lst, sum )
  cnd = cnd[ , -ncol(cnd) ]
  colnames( cnd ) = varNames
  cnd$Cnd = 1:nrow( cnd )
  
  # Initialize output
  out = rep( NA, nrow( dat ) )
  
  for ( i in 1:nrow(cnd) ) {
    
    sel = matrix( FALSE, nrow( dat ), length( varNames ) )
    for ( v in 1:length( varNames ) ) {
      sel[,v] = dat[[varNames[v]]] == cnd[i,varNames[v]]
    }
    out[rowSums(sel)==length(varNames)] = cnd$Cnd[i]
  }
  
  return( out )
}

# Lookup - 12
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

# Lookup - 13
rotatedAxesLabels = function( pos, lab, adj = 0.1, angle = 45,
                              cex = 1, col = NULL, 
                              font = NULL, vfont = NULL, ... ) {
  
  axis( 1, at = pos, labels = FALSE, ... )
  
  text( x = pos, y = par()$usr[3]-adj*(par()$usr[4]-par()$usr[3]),
       labels = lab, srt = angle, adj=1, xpd = TRUE,
       cex = cex, col = col, font = font, vfont = vfont )
  
}
