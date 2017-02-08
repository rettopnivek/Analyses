#--------------------#
# Useful functions   #
# Kevin Potter       #
# Updated 01/19/2017 #
#--------------------#

# Index
# Lookup - 01:  draw_figure
# Lookup - 02:  convergence_extract
# Lookup - 03:  plot_conv
# Lookup - 04:  find_dec
# Lookup - 05:  SDT_prob
# Lookup - 06:  mSDT_prob
# Lookup - 07:  RvSDT_lpmf
# Lookup - 08:  p_f

# Lookup - 01
draw_figure = function( prp, lnSz = 1, lnTp = 1, ptSz = 1, 
                        pts = c(19,22,19,22),
                        clr = c('black','white'),
                        xa = NULL ) {
  # Purpose:
  # Adds a set of lines giving the proportion of times the 
  # response on the right was picked.
  # Arguments:
  # prp  - A set of 8 proportions, sorted by the correct position 
  #        of the target, the type of image (target versus competitor),
  #        and the condition (selective retrieval versus baseline)
  # lnSz - The width of the line
  # lnTp - The type of line
  # ptSz - The size of the point
  # pts  - A vector of 4 point types
  # clr  - The foreground and background color for the points
  # xa   - An optional list giving the x-axis positions for each 
  #        condition
  
  if ( length( xa ) == 0 ) 
    xa = list(
      c(.9,2.1),
      c(1.1,1.9),
      c(.9,2.1),
      c(1.1,1.9) )
  
  for ( i in 1:4 ) {
    
    sel = 1:2 + 2*(i-1)
    lines( xa[[i]], prp[sel], lwd = lnSz, lty = lnTp, col = clr[1] )
    points( xa[[i]], prp[sel], pch = pts[i], col = clr[1],
            bg = clr[2], cex = ptSz )
    
  }
  
}

# Lookup - 02
convergence_extract = function( fit, par_name = NULL ) {
  # Purpose:
  # Extract convergence diagnostics from a Stan fit object.
  # Arguments:
  # fit      - A Stan fit object
  # par_name - An optional string giving the final parameter label 
  #            of the subset of the output to include
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
  if ( length( par_name ) == 0 ) 
    par_name = names( Rhat )[ 
      which( names( Rhat ) == "logLik[1]" ) - 1 ]
  sel = 1:which( names(Rhat) == par_name )
  Rhat = Rhat[sel]; n_eff = n_eff[sel];
  
  return( list( Rhat = Rhat, n_eff = n_eff, totSampSize = totSampSize ) )
}

# Lookup - 03
plot_conv = function( fit, savePlot, parName = NULL ) {
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
  
  conv = convergence_extract(fit,parName)
  
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
  scl = find_dec( tmp$density )
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
  scl = find_dec( tmp$density )
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

# Lookup - 04
find_dec = function( x, spacing = 10 ) {
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

# Lookup - 05
SDT_prob = function( dp, crt, Co ) {
  # Purpose:
  # A function to calculate the probability of picking 
  # an alternative on the right in a two-alternative 
  # forced-choice task based on an equal-variance SDT model.
  # Arguments:
  # dp  - The d' value(s), the separation between means for the 
  #       strength distributions for the left and right 
  #       choices
  # crt - The criterion value(s), denoting the bias towards 
  #       left (positive values) or right (negative values)
  # Co  - The position(s) of the correct alternative, where 
  #       0 = left, and 1 = right.
  # Returns: 
  # The probability (or if vectors are inputs, a vector of 
  # probabilities) for picking the alternative on the right.
  
  # Calculate the probability of picking right given the 
  # correct answer is on the right.
  P_RR = 1.0 - pnorm( crt, dp/2.0 )
  # Calculate the probability of picking right given the 
  # correct answer is on the left.
  P_RL = 1.0 - pnorm( crt, -dp/2.0 )
  
  # Select the appropriate probability for the given trial
  theta = Co*P_RR + (1-Co)*P_RL
  
  return( theta )
}

# Lookup - 06
mSDT_prob = function( dp, crt, lmb, Co ) {
  # Purpose:
  # A function to calculate the probability of picking 
  # an alternative on the right in a two-alternative 
  # forced-choice task based on an mixture SDT model.
  # Arguments:
  # dp  - The d' value(s), the separation between means for the 
  #       strength distributions for the left and right 
  #       choices.
  # crt - The criterion value(s), denoting the bias towards 
  #       left (positive values) or right (negative values)
  # lmb - The mixture probability (represents properly attending 
  #       to the stimulus)
  # Co  - The position(s) of the correct alternative, where 
  #       0 = left, and 1 = right
  # Returns: 
  # The probability (or if vectors are inputs, a vector of 
  # probabilities) for picking the alternative on the right.
  
  # Calculate the probability of picking right given the 
  # correct answer is on the right.
  P_RR = lmb*( 1.0 - pnorm( crt, dp/2.0, 1.0 ) ) + 
    ( 1.0 - lmb )*( 1.0 - pnorm( crt, 0.0, 1.0 ) )
  
  # Calculate the probability of picking right given the 
  # correct answer is on the left.
  P_RL = lmb*( 1.0 - pnorm( crt, -dp/2.0, 1.0 ) ) + 
    ( 1.0 - lmb )*( 1.0 - pnorm( crt, 0.0, 1.0 ) )
  
  # Select the appropriate probability for the given trial
  theta = Co*P_RR + (1-Co)*P_RL
  
  return( theta )
}

# Lookup - 07
RvSDT_lpmf2 = function( y, theta, crt, dp, Co ) {
  # Purpose:
  # Calculates a mixture between a recall process and a SDT 
  # comparison process, where the recall process operates ony 
  # for the left choice.
  # Arguments:
  # y     - An observed response (0 = left, 1 = right)
  # theta - The probability of recall
  # crt   - The criterion for the SDT process
  # dp    - The d' value
  # Co    - The position of the correct response
  # Returns:
  # The log-density.
  
  # Define function to take the log of summed exponentiated vectors
  log_sum_exp = function(x) {
    out = log( exp(x[1]) + exp(x[2]) )
    return( out )
  }
  
  # Define probability of recall
  PR = (1.0 - Co) * theta;
  
  # Define log-density for recall process
  lcd = log( PR ) + dbinom( y, 1, 1.0, log = T );
  
  # Define log-density for comparison process
  sdt = Co * log( ( 1 - pnorm( crt, dp/2.0, 1.0 ) ) ) + 
    (1.0 - Co) * log( ( 1 - pnorm( crt, -dp/2.0, 1.0 ) ) );
  lcd = c( lcd, log( 1.0 - PR ) + sdt );
  
  # Return log-density for mixture
  out = log_sum_exp( lcd );
  
  return( out );
}

# Lookup - 08
p_f = function(x) {
  # Purpose;
  # A function to calculate choice probabilities.
  # Arguments:
  # x - A vector of binary choices (0 or 1)
  # Returns:
  # A vector of two values, the proportions for 
  # picking 0 and 1 respectively.
  
  out = c(0,0); names(out) = c("0","1")
  p = table(x);
  out[ names(p) ] = as.vector( p )
  return( out/sum(out) )
}
