#--------------------#
# Useful functions   #
# Kevin Potter       #
# Updated 04/30/2017 #
#--------------------#

# Index
# Lookup - 01:  design
# Lookup - 02:  create_grid
# Lookup - 03:  meta_subject_create
# Lookup - 04:  plot_jcdf_defaults
# Lookup - 05:  plot_jcdf

# Lookup - 01
design = function( Prime_time, display, cnd_color, 
                   yax = c(.3,.35), 
                   clr = c( 'orange', 'blue', 'red', 'green' ),
                   shft = c(.15,0),
                   txtSz = 1.5 ) {
  # Purpose:
  # A function to plot the design of the experiment based on 
  # a particular prime duration.
  # Arguments:
  # Prime_time
  # display   - 
  # cnd_color - 
  # yax       - 
  # clr       - 
  # shft      - 
  # txtSz     - 
  
  # Determine stimulus display times
  Total_time = 2900
  Target_display = 1900
  Fixation_time = Target_display - Prime_time
  
  # Labels for each stimulus display
  string = c( 'Fixation', 'Prime', 
              'Target until response', 'Feedback' )
  if ( Prime_time == 0 ) string[2] = ""
  
  # Boundaries for start/stop points of each display
  lb = c( 0, 
          Fixation_time/Total_time,
          Target_display/Total_time,
          (Target_display+500)/Total_time )
  ub = c( Fixation_time/Total_time,
          (Fixation_time+Prime_time)/Total_time,
          (Target_display+500)/Total_time,
          1 )
  
  for ( i in 1:4 ) {
    polygon( c( lb[i], lb[i], ub[i], ub[i] ),
             yax[ c(1,2,2,1) ], border = NA,
             col = clr[i] )
    text( lb[i] + (ub[i]-lb[i])/2, yax[1], 
          string[i], cex = txtSz, 
          pos = 1 )
    if ( i == 1 )
      text( ub[i], yax[2],
            paste( Fixation_time, 'ms' ),
            pos = 3, cex = .9 )
    if ( i == 2 & Prime_time > 0 )
      text( ub[i], yax[2],
            paste( Prime_time, 'ms' ),
            pos = 3, cex = .9 )
    if ( i == 4 )
      text( ub[i], yax[2],
            '500 ms',
            pos = 3, cex = .9 )
  }
  
  inc = yax[2] + shft[2]
  for ( ind in 1:length(display) ) {
    
    for ( i in c(1,3) ) {
      text( lb[i] + (ub[i]-lb[i])/2, inc, 
            display[[ind]][i], col = cnd_clr[ind],
            pos = 3, cex = txtSz )
    }
    
    if ( Prime_time != 0 )
      text( lb[2] + (ub[2]-lb[2])/2, inc, 
            display[[ind]][2], col = cnd_clr[ind],
            pos = 3, cex = txtSz )
    
    inc = inc + shft[1]
    
  }
  
}

# Lookup - 02
create_grid = function( xl, yl, xinc, yinc, ... ) {
  # Purpose:
  # A function to add a grid to an already existing plot.
  # Arguments:
  # xl   - The lower and upper limits for the x-axis
  # yl   - The lower and upper limites for the y-axis
  # xinc - The increments at which to draw lines for the x-axis
  # yinc - The increments at which to draw lines for the y-axis
  # ...  - Additional plotting parameters for the 'segments'
  #        function
  
  xa = seq( xl[1], xl[2], xinc )
  segments( xa, rep( yl[1], length(xa) ), 
            xa, rep( yl[2], length(xa) ), ... )
  xa = seq( yl[1], yl[2], yinc )
  segments( rep( xl[1], length(xa) ), xa,
            rep( xl[2], length(xa) ), xa, ... )
}

# Lookup - 03
meta_subject_create = function( d, Nt = NULL ) {
  # Purpose:
  # Creates a meta-subject by averaging the joint CDF functions 
  # over subjects for each condition and simulating data via 
  # an approximation to the inverse CDF.
  # Arguments:
  # d - A data-frame the appropriate variables
  # Returns:
  # A data-frame.
  
  # Define a function to determine either the earliest position 
  # in a vector larger than a value, or the latest position in 
  # a vector smaller than a value
  minMax = function( comp, x, Direction ) {
    
    L = length( x )
    
    if ( Direction == "larger" ) {
      chk = x > comp
      if ( sum( chk ) == 0 )
        out = L else
          out = min( which( chk == 1 ) )
    }
    if ( Direction == "smaller") {
      chk = x < comp
      if ( sum( chk ) == 0 )
        out = 1 else
          out = max( which( chk == 1 ) )
    }
    
    return( out )
  }
  
  # Define a function for linear interpolation given 
  # a pair of x and y-axis values
  linInterp = function( val, y0, y1, x0, x1, x = T ) {
    
    b1 = ( y1 - y0 ) / ( x1 - x0 ); # Slope
    b0 = y1 - b1*x1; # Intercept
    
    if ( x ) {
      num = val - b0;
      new = ( num )/b1;
    } else {
      new = b0 + b1 * val
    }
    
    return( new );
  }
  
  # Define a function to estimate a response time based on 
  # a given probability and an empirical CDF curve
  est_rt = function( q, cdf ) {
    
    y = cdf$y;
    x = cdf$x
    dff = 1
    inc = 1
    
    lb = minMax( q, y, "smaller" )
    ub = minMax( q, y, "larger" )
    
    xN = linInterp( q, y[lb], y[ub], x[lb], x[ub] )
    
    return( xN )
  }
  
  # If no trial numbers are provided, determine average 
  # number of trials per condition
  if ( is.null( Nt ) ) {
    Nt = aggregate( rep( 1, nrow( d ) ), 
                    list( d$S, d$Co, d$Cnd ), sum )
    Nt = aggregate( Nt$x, list( Nt[,2], Nt[,3] ), mean )
    colnames( Nt ) = c( 'Co', 'Cnd', 'N' )
    Nt$N = round( Nt$N )
  }
  
  # Initialize output
  out = matrix( NA, sum( Nt$N ), 4 )
  colnames( out ) = c( 'RT', 'Ac', 'Co', 'Cnd' )
  out = as.data.frame( out )
  for ( n in 1:nrow( Nt ) ) {
    
    if ( n == 1 ) {
      init = 1
      fin = Nt$N[n]
    } else {
      init = 1 + sum( Nt$N[ 1:(n-1) ] )
      fin = Nt$N[n] + init - 1
    }
    
    out$Co[ init:fin ] = Nt$Co[n];
    out$Cnd[ init:fin ] = Nt$Cnd[n];
  }
  
  # Loop over conditions
  for ( i in 1:nrow( Nt ) ) {
    
    # Extract relevant data
    sel = d$Cnd == Nt$Cnd[i] & d$Co == Nt$Co[i]
    
    rt = d$RT[ sel ];
    ac = d$Ac[ sel ];
    grp = d$S[ sel ];
    
    # Determine number of trials to simulate
    N_cur = Nt$N[i]
    
    # Determine number of actual trials
    N_obs = c( sum( ac == 1 ), sum( ac == 0 ) )
    
    # Extract overall accuracy
    acc = c( mean( ac == 1 ), mean( ac == 0 ) )
    
    # Determine CDF curves
    if ( acc[1] > 1/N_obs[1] ) {
      out1 = cdf_curve( rt, ac, sel = 1, grp = grp,
                        opt = list( draw = F, out = T ) )
      cdf1 = list( x = c( out1$pv$x, max( rt[ ac == 1 ] ) ),
                   y = c( out1$pv$y/acc[1], 1 ) )
    } else {
      if ( acc[1] == 0 ) cdf1 = 
          list( x = NULL, y = NULL )
      if ( acc[1] > 0 ) cdf1 = 
          list( x = mean( rt[ ac == 1 ] ), y = 1 )
    }
    
    if ( acc[2] > 1/N_obs[2] ) {
    out0 = cdf_curve( rt, ac, sel = 0, grp = grp,
                      opt = list( draw = F, out = T ) )
    cdf0 = list( x = c( out0$pv$x, max( rt[ ac == 0 ] ) ),
                 y = c( out0$pv$y/acc[2], 1 ) )
    } else {
      if ( acc[2] == 0 ) cdf0 = 
          list( x = NULL, y = NULL )
      if ( acc[2] > 0 ) cdf0 = 
          list( x = mean( rt[ ac == 0 ] ), y = 1 )
    }
    
    # Simulate data
    
    # Accuracy
    ac_sim = as.numeric( runif( N_cur ) < acc[1] )
    rt_sim = numeric( N_cur )
    
    u = runif( N_cur )
    
    if ( mean( ac_sim == 1 ) > 0 ) {
      rt_sim[ ac_sim == 1 ] = sapply( u[ ac_sim == 1 ], 
                                      est_rt, cdf = cdf1 )
      rt_sim[ is.na( rt_sim ) ] = min( rt[ ac == 1 ] )
    }
    if ( mean( ac_sim == 0 ) > 0 ) {
      rt_sim[ ac_sim == 0 ] = sapply( u[ ac_sim == 0 ], 
                                      est_rt, cdf = cdf1 )
      rt_sim[ is.na( rt_sim ) ] = min( rt[ ac == 0 ] )
    }
    
    sel = out$Co == Nt$Co[i] & out$Cnd == Nt$Cnd[i]
    out$RT[ sel ] = rt_sim
    out$Ac[ sel ] = ac_sim
    
  }
  
  # Add additional labels to simulated data
  
  # Label for flanker
  out$FL = 'None'
  sel = out$Cnd == 2 | out$Cnd == 8 | 
    out$Cnd == 14
  out$FL[ sel ] = 'Identical'
  
  sel = out$Cnd == 3 | out$Cnd == 9 | 
    out$Cnd == 15
  out$FL[ sel ] = 'Same'
  
  sel = out$Cnd == 4 | out$Cnd == 10 | 
    out$Cnd == 16
  out$FL[ sel ] = 'Incompatible'
  
  # Label for prime type
  out$PTL = 'None'
  inc = 0
  tmp = c( 'Identical', 'Same', 'Incompatible' )
  for (i in 2:4) {
    sel = out$Cnd == (5 + inc) | 
      out$Cnd == (8 + inc) | 
      out$Cnd == (11 + inc) | 
      out$Cnd == (14 + inc)
    out$PTL[ sel ] = tmp[ i - 1 ]
    inc = inc + 1
  }
  
  # Create variable for prime duration (in s)
  out$PD = 0
  out$PD[ out$Cnd > 4 & 
                          out$Cnd < 11 ] = .1
  out$PD[ out$Cnd > 10 ] = .8
  
  # Create meaningful label for correct and choice options
  out$CL = 'consonant'
  out$CL[ out$Co == 0 ] = 'vowel'
  out$Ch = 0 # Vowel
  sel = ( out$Ac == 1 & out$Co == 1 ) | 
    ( out$Ac == 0 & out$Co == 1 )
  out$Ch[sel] = 1 # Consonant
  out$RL = 'vowel'
  out$RL[ out$Ch == 1 ] = 'consonant'
  
  return( out )
}

quick_cnd_select = function( d, duration, type, flanker, correct ) {
  # Purpose: 
  # ... 
  # Arguments: 
  # ... 
  # Returns: 
  # ...
  
  # Select desired conditions
  v1 = exp_str$PD %in% duration
  v2 = exp_str$PT %in% type
  v3 = exp_str$F %in% flanker
  cnd = d$Cnd %in% exp_str$Cnd[ which( v1 + v2 + v3 == 3 ) ]
  
  # Select desired correct values
  ch = d$Co %in% correct
  
  # Return new dataset conditioned on variables
  return( d[ cnd & ch, ] )
}

# Lookup - 04
plot_jcdf_defaults = function( opt ) {
  # Purpose: 
  # ... 
  # Arguments: 
  # ... 
  # Returns: 
  # ...
  
  if ( is.null( opt$xl ) ) opt$xl = c(.2,1.2)
  if ( is.null( opt$yl ) ) opt$yl = c(0,1)
  if ( is.null( opt$lnSz ) ) opt$lnSz = 1.2
  if ( is.null( opt$txtA ) ) opt$txtA = 1.25
  if ( is.null( opt$txt ) ) opt$txt = 1.5
  if ( is.null( opt$ptSz ) ) opt$ptSz = 1.5
  if ( is.null( opt$uclr ) ) opt$uclr = rgb( .5, .5, .5, .3 )
  if ( is.null( opt$mrgn ) ) opt$mrgn = c( 4, 1, 3, 5 )
  
  return( opt )
}

# Lookup - 05
plot_jcdf = function( d, cnd, f = NULL, opt = list() ) {
  # Purpose: 
  # ... 
  # Arguments: 
  # ... 
  # Returns: 
  # ...
  
  # Generate default options
  opt = plot_jcdf_defaults( opt )
  
  # Extract plot characteristics
  lnSz = opt$lnSz
  ptSz = opt$ptSz
  txtA = opt$txtA
  txt = opt$txt
  xl = opt$xl
  yl = opt$yl
  uclr = opt$uclr
  mrgn = opt$mgrn
  
  # Generate blank plot
  par( opt$mrgn )
  blankPlot( xl, yl )
  
  # Add in x-axis information
  tck = seq( xl[1], xl[2], length = 5 )
  create_grid( xl, yl, .25, .25, col = 'grey80', lwd = lnSz )
  segments( xl[1], 0, xl[2], 0, lwd = lnSz )
  axis( 1, round( tck, 1 ), tick = F, line = -1.2, cex.axis = txtA )
  mtext( 'Time (s)', side = 1, line = 2, cex = txt )
  
  # Add in y-axis information
  segments( xl[2], yl[1], xl[2], yl[2], lwd = lnSz )
  axis( 4, c( -1, -.5, 0, .5, 1 ),
        c( 1, .5, 0, .5, 1 ), tick = F, line = -1.2, cex.axis = txtA )
  mtext( 'Joint CDF', side = 4, line = 2.5, cex = txt )
  
  for ( i in 1:length( cnd ) ) {
    
    # Extract data
    nd = quick_cnd_select( d, cnd[[i]]$PD, cnd[[i]]$PT, cnd[[i]]$FL, cnd[[i]]$Co )
    
    # Plotting options
    o = list( out = T, flip = cnd[[i]]$flip )
    
    # Consonant
    cv = tryCatch( cdf_curve( nd$RT, nd$Ch, grp = nd$S, lwd = lnSz, opt = o ),
                   error = function(e) NULL )
    if ( !is.null(cv) & !is.null(f) ) add_uncertainty( cv, f = f, border = NA, col = uclr )
    
    # Vowel
    cv = tryCatch( cdf_curve( nd$RT, nd$Ch, grp = nd$S, sel = 0, lwd = lnSz, lty= 2, opt = o ),
                   error = function(e) NULL )
    if ( !is.null(cv) & !is.null(f) ) add_uncertainty( cv, f = f, border = NA, col = uclr )
    
  }
}