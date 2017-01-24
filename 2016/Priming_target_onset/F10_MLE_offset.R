#------------------------------#
# Examination of offsets using #
# MLE                          #
# Kevin Potter                 #
# Updated 12/20/2016           #
#------------------------------#

# Initialize script
source('F4_starting_script.R')

# Indicate which code segments to run
runCode = c( F, F, T, T )

# Indicate if figures should be saved as PDF files
savePlot = T

# Load in data
setwd( 'Data' )
load( 'Priming_offset.RData' )
setwd( orig_dir )

# For easy manipulation
d = allDat
colnames( d ) = c('S','PD','PT','O','Co','RT',
                  'Ch','Ac','OT','DT','CDT','Cnd')

# Load in ML estimates and AIC values
if ( !runCode[1] ) {
  setwd( 'Data' )
  load( 'MLE_results.RData' )
  setwd( orig_dir )
}

# Define version of model to fit
type = 1

# Extract data and covariates
input = model_structure_create(type, d, Priors = NULL )

# Index
# Lookup - 01:  Define functions for current script
# Lookup - 02:  ML estimation per subject
# Lookup - 03:  Comparison of Akaike weights
# Lookup - 04:  Figures for model fit

###
### Define functions for current script
###
# Lookup - 01

# Load in package for numerical estimation of derivatives
library( numDeriv )

tran_par = function( par, min_rt, reverse = F ) {
  # Purpose:
  # Transforms or un-transforms parameter values.
  # Arguments:
  # par     - A vector of parameter values
  # min_rt  - The fastest response time value
  # reverse - A logical value, indicating whether a 
  #           transformation should be reversed or not
  # Returns:
  # A new vector of parameter values.
  
  if ( reverse ) {
    
    # Convert from positive-only
    par[ c( 1:2, 6:13 ) ] = log( par[ c( 1:2, 6:13 ) ] )
    # Convert from restricted range
    par[14] = logit( par[14]/min_rt )
    
  } else {
    
    # Restrict to be positive
    par[ c( 1:2, 6:13 ) ] = exp( par[ c( 1:2, 6:13 ) ] )
    # Restrict to be between a range
    par[14] = logistic( par[14] )*min_rt
    
  }
  
  return( par )
}

mle_fn = function( par, dat, priors = NULL, est = T ) {
  # Purpose:
  # Calculates the sum of the log-likelihoods given a vector of 
  # parameter values.
  # Arguments:
  # par    - A vector of parameter values
  # dat    - A list, the output of the function 'model_structure_create'
  # priors - An optional list of prior values for the parameters
  # est    - A logical value, if true, the sum of the log-likelihoods 
  #          is returned, else the vector of log-likelihoods is returned
  # Returns:
  # The sum of the log-likelihoods or the vector of log-likelihoods
  
  # Extract data
  No = dat$No[s]
  rt = dat$Y[s,1:No,1];
  ch = dat$Y[s,1:No,2];
  X = dat$X[s,,1:No]
  
  # Determine parameter estimates
  cf = tran_par( par, min(rt) )
  
  prm = param_est( X, cf, dat$fixed, dat$index, dat$parSel )
  
  # Extract offset times
  if ( yesOffset ) {
    # Adjustment for foil
    FO = d$OT[ d$S == s ]
    # Adjustment for target
    TO = -d$OT[ d$S == s ]
    
    Co = d$Co[ d$S == s ]
    
    ot1 = TO*Co + FO*(1-Co)
    ot0 = TO*(1-Co) + FO*Co
    
    prm[4,] = prm[4,] + ot1
    prm[8,] = prm[8,] + ot0
    
  }
  
  # Calculate likelihood
  ll = dwaldrace( rt, ch, prm[1,], prm[2,], prm[4,],
                  prm[5,], prm[6,], prm[8,], ln = 1 )
  sll = sum( ll )
  
  # Incorporate priors
  if ( length( priors ) > 0 ) {
    
    sll = sll + sum( dnorm( cf, priors[,1], priors[,2], log = T ) )
    
  }
  
  if ( is.na( sll ) ) sll = -Inf
  
  if ( est ) {
    return( sll )
  } else {
    return( ll )
  }
}

grad_fn = function( par, dat, priors = NULL ) {
  
  out = grad( mle_fn, par, dat = dat, priors = priors )
  
  return( out )
}

st_fn = function() {
  # Purpose:
  # Generates a set of starting values for the optimization routine.
  # Returns:
  # A set of parameter values.
  
  # If desired, generate starting values using a normal distribution 
  # centered as a set of pre-specified values
  if ( exists( 'startCenter' ) ) {
    
    par = rnorm( 14, startCenter[,1], startCenter[,2] )
    par[14] = par[14]*min_rt
    
  } else {
    
    par = runif( 14, 
                 c( rep( .5, 2 ), rep( -.2, 3 ), rep( .5, 8 ), 
                    .5*min_rt ),
                 c( rep( 1.5, 2 ), rep( .2, 3 ), rep( 4, 8 ), 
                    .9*min_rt ) )
  }
  
  par = tran_par( par, min_rt, reverse = T )
  
  return( par )
}

ui_1s = function( out ) {
  # Purpose:
  # Function to add uncertainty intervals to a RT plot
  # Arguments:
  # out - The output from the 'cdf_curve' function for a 
  #       single subject
  
  ui = rbind(
    qbinom( .16, out$v$n, out$pv$y )/out$v$n,
    qbinom( .84, out$v$n, out$pv$y )/out$v$n
  )
  polygon( c( out$pv$x, rev(out$pv$x) ),
           c( ui[1,], rev(ui[2,]) ), border = NA, 
           col = rgb( .5, .5, .5, .3 ) )
}

###
### ML estimation per subject
###
# Lookup - 02

if ( runCode[1] ) {
  
  if ( savePlot ) {
    setwd( 'Plots' )
    pdf( 'MLE_offsets_likelihoods.pdf',width=12)
    setwd( orig_dir )
  }
  
  # Obtain ML estimates, AICc values, and parameter values 
  # for each subject
  
  AICcNO = numeric( N )
  AICcYO = numeric( N )
  AICcYO2 = numeric( N )
  
  allPar = array( NA, dim = c( N, 14, 2 ) )
  rawPar = array( NA, dim = c( N, 14, 2 ) )
  
  # Create a progress bar using a base R function
  pb = txtProgressBar( min = 1, max = N, style = 3 )
  
  # Define a set of priors
  priors = cbind(
    
    c( 1.22, 1.14, .18, -.03, -.04,
       1.23, 2.03, 2.06, 1.44,
       1.49, .71, .78, 1.00,
       .27 ),
    
    c( .5, .5, .3, .3, .3,
       .8, .8, .8, .8,
       .8, .8, .8, .8,
       .2 )
  )
  
  for ( s in 1:N ) {
    
    startTime = Sys.time()
    
    # First estimate without the offset values
    yesOffset = F
    
    # Fastest RT
    min_rt = input$min_RT[ s ]
    
    # Carry out MLE
    
    # Nelder-Mead (Faster)
    # res = MLE( input, mle_fn, st_fn, nRep = 20,
    #            maxit = 10000 )
    
    # Conjugate gradient with numerical gradient function
    # (Very slow)
    res = MLE( input, mle_fn, st_fn, grad_fn = grad_fn,
               method = "CG", nRep = 4, maxit = 10000 )
    
    # Extract parameters
    rawPar[s,,1] = res$param
    allPar[s,,1] = tran_par( res$param, min_rt )
    
    # Calculate AICc
    AICcNO[s] = infoCrit( mle_fn( res$param, input ), 14, input$No[s] )
    # Calculate AICc after adjusting estimates by offsets
    yesOffset = T
    AICcYO[s] = infoCrit( mle_fn( res$param, input ), 14, input$No[s] )
    
    runTime = startTime - Sys.time()
    
    # Carry out MLE with offsets included
    
    # Nelder-Mead (Faster)
    # res2 = MLE( input, mle_fn, st_fn, nRep = 20, maxit = 10000 )
    
    # Conjugate gradient with numerical gradient function
    # (Very slow)
    res2 = MLE( input, mle_fn, st_fn, grad_fn = grad_fn,
               method = "CG", nRep = 4, maxit = 10000 )
    
    # Extract parameters
    rawPar[s,,2] = res2$param
    allPar[s,,2] = tran_par( res2$param, min_rt )
    
    # Calculate AICc
    AICcYO2[s] = infoCrit( mle_fn( res2$param, input ), 14, input$No[s] )
    
    if ( savePlot ) {
      
      rt = d$RT[ d$S == s ]
      ac = d$Ac[ d$S == s ]
      
      yesOffset = F
      ll = mle_fn( res$param, input, est = F )
      yesOffset = T
      ll2 = mle_fn( res$param, input, est = F )
      ll3 = mle_fn( res2$param, input, est = F )
      
      layout( cbind( 1, 2 ) )
      xl = lowerUpper( .2, rt )
      yl = lowerUpper( .1, exp( c(ll,ll2,ll3) ) )
      yl = c( -yl[2], yl[2] )
      
      blankPlot( xDim = xl, yDim = yl )
      abline( h = 0 )
      axis( 1, seq(xl[1],xl[2],.2) )
      mtext('Likelihood',side=2,line=2.5)
      mtext('RT (s)',side=1,line=2)
      title( 'MLE versus MLE + post-offset' )
      
      for ( i in 0:1 ) {
        
        sel = ac == i
        adj = 1; if ( i == 0 ) adj = -1
        
        points( rt[ sel ], adj*exp( ll[ sel ] ), pch = 19 )
        
        tst = exp( ll ) - exp( ll2 )
        
        segments( rt[ sel & tst < 0 ], adj*exp( ll[ sel & tst < 0 ] ),
                  rt[ sel & tst < 0 ], adj*exp( ll2[ sel & tst < 0 ] ),
                  col = 'red' )
        
        segments( rt[ sel & tst > 0 ], adj*exp( ll[ sel & tst > 0 ] ),
                  rt[ sel & tst > 0 ], adj*exp( ll2[ sel & tst > 0 ] ),
                  col = 'black' )
        
        points( rt[ sel ], adj*exp( ll2[ sel ] ), pch = 20, col = 'red' )
        
        if ( i == 1 ) {
          tmp = round( tran_par( res$param, min(rt) ), 2 )
          string = c( paste( tmp[1:5], collapse = '|' ),
                      paste( tmp[6:9], collapse = '|' ),
                      paste( tmp[10:13], collapse = '|' ),
                      as.character( tmp[14] ) )
          
          legend( 'topright', string, bty = 'n', text.col = 'black' )
        }
        
        if ( i == 0 ) {
          
          legend( 'bottomright', c('MLE','MLE+post-offset'),
                  fill = c('black','red'), bty = 'n' )
          
        }
        
      }
      
      blankPlot( xDim = xl, yDim = yl )
      abline( h = 0 )
      axis( 1, seq(xl[1],xl[2],.2) )
      mtext('Likelihood',side=2,line=2.5)
      mtext('RT (s)',side=1,line=2)
      title( 'MLE + post-offset versus MLE + pre-offset' )
      
      for ( i in 0:1 ) {
        
        sel = ac == i
        adj = 1; if ( i == 0 ) adj = -1
        
        points( rt[ sel ], adj*exp( ll2[ sel ] ), pch = 19, col = 'red' )
        
        tst = exp( ll ) - exp( ll2 )
        
        segments( rt[ sel & tst < 0 ], adj*exp( ll2[ sel & tst < 0 ] ),
                  rt[ sel & tst < 0 ], adj*exp( ll3[ sel & tst < 0 ] ),
                  col = 'green' )
        
        segments( rt[ sel & tst > 0 ], adj*exp( ll2[ sel & tst > 0 ] ),
                  rt[ sel & tst > 0 ], adj*exp( ll3[ sel & tst > 0 ] ),
                  col = 'red' )
        
        points( rt[ sel ], adj*exp( ll3[ sel ] ), pch = 20, col = 'green' )
        
        if ( i == 1 ) {
          tmp = round( tran_par( res2$param, min(rt) ), 2 )
          string = c( paste( tmp[1:5], collapse = '|' ),
                      paste( tmp[6:9], collapse = '|' ),
                      paste( tmp[10:13], collapse = '|' ),
                      as.character( tmp[14] ) )
          
          legend( 'topright', string, bty = 'n', text.col = 'black' )
        }
        
        if ( i == 0 ) {
          
          legend( 'bottomright', c('MLE+post-offset','MLE+pre-offset'),
                  fill = c('red','green'), bty = 'n' )
          
        }
        
      }
      
    }
    
    # Update the progress bar
    setTxtProgressBar(pb,s)
    
  }
  close(pb)
  if (savePlot) dev.off()
  
  MLE_res = list(
    AICcNO = AICcNO,
    AICcYO = AICcYO,
    AICcYO2 = AICcYO2,
    allPar = allPar )
  
  # Save MLE results
  setwd( 'Data' )
  save( MLE_res, file = 'MLE_results.RData' )
  setwd( orig_dir )
  
}

###
### Comparison of Akaike weights
###
# Lookup - 03

if ( runCode[2] ) {
  
  # Calculate Akaike weights
  w = cbind( MLE_res$AICcNO,
             MLE_res$AICcYO,
             MLE_res$AICcYO2 )
  w = w - apply( w, 1, min )
  w = exp( -.5*w )
  w = w/rowSums( w )
  colnames( w ) = c( 'NO','PO','AO' )
  
  if (savePlot) {
    setwd( 'Plots' )
    pdf( 'AICc_comparison_offset_MLE.pdf' )
  } else x11()
  
  lyt = matrix( 1, 10, 10 )
  lyt[1,] = 2
  layout( lyt )
  
  par( mar = c( 4, 5, 1, 1 ) )
  blankPlot( xDim = c(1,N), yDim = c(0,1) )
  
  axis( 2, seq( 0, 1, .25 ), cex.axis = 1.5,
        lwd = 2 )
  mtext( 'Akaike weights', side = 2, cex = 1.5, line = 2.5 )
  mtext( 'Subjects', side = 1, cex = 1.5, line = 2 )
  
  # Re-order subjects based on rankings for AIC values
  ord = rev( order( w[,1] ) )
  
  # Shade different regions of notable fit
  xa = c(.5, sum( w[,1] > .95 )+.5)
  polygon( c(xa,rev(xa)), c(-.2,-.2,1.2,1.2),
           col = 'grey80', border = NA )
  xa = c(N-sum( w[,3] > .95 )+.5,N+.5)
  polygon( c(xa,rev(xa)), c(-.2,-.2,1.2,1.2),
           col = 'grey90', border = NA )
  axis(3, c( sum( w[,1] > .95 )/2 + .5,
             N - sum( w[,3] > .95 )/2 + .5 ),
       paste( round( c( sum( w[,1] > .95 )/N, 
                        sum( w[,3] > .95 )/N ), 2 )*100, '%', sep = '' ), 
       tick = F, cex.axis = 1.5 )
  
  # Plot AKaike weights
  clr = c( 'blue', 'red', 'green' )
  for ( i in 1:3 ) points( 1:N, w[ord,i], pch = 19, col = clr[i],
                           cex = 1.5 )
  
  par( mar = c(0,4,0,0) )
  blankPlot()
  legend( 'left', c('No offset','Offset post-MLE','Offset pre-MLE'),
          fill = clr, bty = 'n', cex = 1.5, horiz = T )
  
  if (savePlot) { dev.off(); setwd( orig_dir ) }
}

###
### Figures for model fit (No offset)
###
# Lookup - 04

if ( runCode[3] ) {
  
  if ( savePlot ) {
    setwd( 'Plots' )
    pdf( 'MLE_fit_no_offset.pdf', width = 12 )
  }
  
  # Vector of residuals
  rsd = rep( NA, nrow(d) )
  
  for ( s in 1:N ) {
    
    # Extract parameter estimates
    curPar = MLE_res$allPar[s,,]
    prm = param_est( input$X_small,
                     curPar[,1], input$fixed, input$index,
                     input$parSel )
    
    ### Calculate residuals ###
    
    # Extract data
    rto = d$RT[ d$S == s ]
    cho = d$Ch[ d$S == s ]
    
    # Calculate observed and predicted joint CDFs
    y = numeric( length(rto) )
    yhat = numeric( length(rto) )
    for ( cnd in 1:8 ) {
      sel = d$CDT[ d$S == s ] == cnd
      
      rtC = rto[ sel ]
      chC = cho[ sel ]
      cury = numeric( length(rtC) )
      
      # Calculate predicted values
      yhat[ sel ] = pwaldrace( rtC, chC, 
                               prm[1,cnd], prm[2,cnd], prm[4,cnd],
                               prm[5,cnd], prm[6,cnd], prm[8,cnd] )
      
      for ( j in 0:1) cury[ chC == j ] = rank( rtC[ chC == j ] )/sum(sel)
      y[ sel ] = cury
      
    }
    
    # Save residuals
    rsd[ d$S == s ] = y - yhat
    
    ### Plot joint CDFs per subject ###
    
    if (!savePlot) x11( width = 12 )
    par( mar = c( 5, 4, 3, 1 ) )
    layout( rbind( c(1,3,5,7), c(2,4,6,8) ) )
    
    ttl = c( 'Target primed .05 s',
             'Target primed .4 s',
             'Foil primed .05 s',
             'Foil primed .4 s' )
    
    for ( cnd in 1:8 ) {
      
      blankRTplot( bty='l', tDim = c(.2, 2.8), 
                   cex.axis = 1.5, cex.lab = 1.5 )
      
      for (k in 1:4) { if ( cnd == c(1,3,5,7)[k] ) 
        title( ttl[k], cex = .9 ) }
      
      if ( cnd == 7 ) legend( 'topleft', 'Left correct',
                              bty = 'n' )
      if ( cnd == 8 ) legend( 'topleft', 'Right correct',
                              bty = 'n' )
      if ( cnd == 1 ) legend( 'topleft', paste( 'Subject', s ),
                              bty = 'n' )
      
      # Define test statistic
      T_x = function(x) quantile(x,seq(.1,.9,.2))
      
      rt = d$RT[ d$CDT == cnd & d$S == s ];
      ch = d$Ch[ d$CDT == cnd & d$S == s ];
      
      if ( mean( ch == 0 ) < 1 ) {
        out1 = cdf_curve( rt, ch, opt = list( out = T, draw = F ) )
        ui_1s( out1 )
        cdf_curve( rt, ch, lwd = 2 )
        add_points( out1, T_x = T_x, pch = 19, cex = 1.5 )
      }
      
      if ( mean( ch == 1 ) < 1 ) {
        
        out0 = cdf_curve( rt, ch, sel = 0, opt = list( out = T, draw = F ) )
        ui_1s( out0 )
        cdf_curve( rt, ch, sel = 0, lwd = 2, lty = 2 )
        add_points( out0, T_x = T_x, pch = 17, cex = 1.5 )
        
      }
      
      dist_plot( prm = prm[ -c(3,7), cnd ], 
                 rt = seq(.2,2.8,length=1000), 
                 lwd = 2, col = 'green' )
      dist_plot( prm = prm[ -c(3,7), cnd ], ch = 0, 
                 rt = seq(.2,2.8,length=1000), 
                 lwd = 2, lty = 2, col = 'red' )
      
    }
    
  }
  
  ### Plot residuals per condition ###
  if (!savePlot) x11( width = 12 )
  
  layout( rbind( c(1,3,5,7), c(2,4,6,8) ) )
  
  ag_rsd = aggregate( rsd, list( d$CDT, d$S ), quantile )
  colnames( ag_rsd ) = c( 'Cnd', 'S', 'Q' )
  m_rsd = aggregate( rsd, list( d$CDT ), mean )
  
  ttl = c( 'Target primed .05 s',
           'Target primed .4 s',
           'Foil primed .05 s',
           'Foil primed .4 s' )
  
  yl = lowerUpper( .2, as.vector( ag_rsd$Q ) )
  for ( cnd in 1:8 ) {
    
    blankPlot( c(1,N), yl  )
    axis( 2, seq(yl[1],yl[2],.2) )
    abline( h = 0, lwd = 2, col = 'grey' )
    abline( h = m_rsd$x[cnd], lwd = 2,lty = 2, col = 'blue' )
    
    for (k in 1:4) { if ( cnd == c(1,3,5,7)[k] ) 
      title( ttl[k], cex = .9 ) }
    if ( cnd == 7 ) legend( 'topleft', 'Left correct',
                            bty = 'n' )
    if ( cnd == 8 ) legend( 'topleft', 'Right correct',
                            bty = 'n' )
    
    sel = ag_rsd$Cnd == cnd
    
    segments( 1:N, ag_rsd$Q[ sel, 1 ],
              1:N, ag_rsd$Q[ sel, 5 ] )
    segments( 1:N, ag_rsd$Q[ sel, 2 ],
              1:N, ag_rsd$Q[ sel, 4 ], lwd = 3 )
    points( 1:N, ag_rsd$Q[ sel, 3 ],
            pch = 21, bg = 'white', cex = 1.5 )
    
  }
  mtext('Residuals',side=3,outer=T,line=-2.75,cex=1.5)
  
  if (savePlot) {
    dev.off()
    setwd( orig_dir )
  }
}

###
### Figures for model fit (Pre-MLE offsets)
###
# Lookup - 05

if ( runCode[4] ) {
  
  if ( savePlot ) {
    setwd( 'Plots' )
    pdf( 'MLE_fit_offset.pdf', width = 12 )
  }
  
  # Vector of residuals
  rsd = rep( NA, nrow(d) )
  
  for ( s in 1:N ) {
    
    # Extract parameter estimates
    curPar = MLE_res$allPar[s,,]
    
    ### Calculate residuals ###
    
    prmAll = param_est( input$X[s,,1:(input$No[s])],
                     curPar[,2], input$fixed, input$index,
                     input$parSel )
    # Adjustment for foil
    FO = d$OT[ d$S == s ]
    # Adjustment for target
    TO = -d$OT[ d$S == s ]
    
    Co = d$Co[ d$S == s ]
    
    ot1 = TO*Co + FO*(1-Co)
    ot0 = TO*(1-Co) + FO*Co
    
    prmAll[4,] = prmAll[4,] + ot1
    prmAll[8,] = prmAll[8,] + ot0
    
    prm = matrix( NA, 8, 8 )
    prm[c(3,7),] = 1
    
    for ( cnd in 1:8 ) {
      for ( i in c(1,2,4,5,6,8) ) {
        prm[ i, cnd ] = 
          mean( prmAll[ i, d$CDT[ d$S == s ] == cnd ] )
      }
    }
    
    # Extract data
    rto = d$RT[ d$S == s ]
    cho = d$Ch[ d$S == s ]
    
    # Calculate predicted values
    yhat = pwaldrace( rto, cho, 
                      prm[1,], prm[2,], prm[4,],
                      prm[5,], prm[6,], prm[8,] )
    
    # Calculate observed and predicted joint CDFs
    y = numeric( length(rto) )
    for ( cnd in 1:8 ) {
      sel = d$CDT[ d$S == s ] == cnd
      
      rtC = rto[ sel ]
      chC = cho[ sel ]
      cury = numeric( length(rtC) )
      
      for ( j in 0:1) cury[ chC == j ] = rank( rtC[ chC == j ] )/sum(sel)
      y[ sel ] = cury
      
    }
    
    # Save residuals
    rsd[ d$S == s ] = y - yhat
    
    ### Plot joint CDFs per subject ###
    
    if (!savePlot) x11( width = 12 )
    par( mar = c( 5, 4, 3, 1 ) )
    layout( rbind( c(1,3,5,7), c(2,4,6,8) ) )
    
    ttl = c( 'Target primed .05 s',
             'Target primed .4 s',
             'Foil primed .05 s',
             'Foil primed .4 s' )
    
    for ( cnd in 1:8 ) {
      
      blankRTplot( bty='l', tDim = c(.2, 2.8), 
                   cex.axis = 1.5, cex.lab = 1.5 )
      
      for (k in 1:4) { if ( cnd == c(1,3,5,7)[k] ) 
        title( ttl[k], cex = .9 ) }
      
      if ( cnd == 7 ) legend( 'topleft', 'Left correct',
                              bty = 'n' )
      if ( cnd == 8 ) legend( 'topleft', 'Right correct',
                              bty = 'n' )
      if ( cnd == 1 ) legend( 'topleft', paste( 'Subject', s ),
                              bty = 'n' )
      
      # Define test statistic
      T_x = function(x) quantile(x,seq(.1,.9,.2))
      
      rt = d$RT[ d$CDT == cnd & d$S == s ];
      ch = d$Ch[ d$CDT == cnd & d$S == s ];
      
      if ( mean( ch == 0 ) < 1 ) {
        out1 = cdf_curve( rt, ch, opt = list( out = T, draw = F ) )
        ui_1s( out1 )
        cdf_curve( rt, ch, lwd = 2 )
        add_points( out1, T_x = T_x, pch = 19, cex = 1.5 )
      }
      
      if ( mean( ch == 1 ) < 1 ) {
        
        out0 = cdf_curve( rt, ch, sel = 0, opt = list( out = T, draw = F ) )
        ui_1s( out0 )
        cdf_curve( rt, ch, sel = 0, lwd = 2, lty = 2 )
        add_points( out0, T_x = T_x, pch = 17, cex = 1.5 )
        
      }
      
      dist_plot( prm = prm[ -c(3,7), cnd ], 
                 rt = seq(.2,2.8,length=1000), 
                 lwd = 2, col = 'green' )
      dist_plot( prm = prm[ -c(3,7), cnd ], ch = 0, 
                 rt = seq(.2,2.8,length=1000), 
                 lwd = 2, lty = 2, col = 'red' )      
    }
    
  }
  
  ### Plot residuals per condition ###
  if (!savePlot) x11( width = 12 )
  
  layout( rbind( c(1,3,5,7), c(2,4,6,8) ) )
  
  ag_rsd = aggregate( rsd, list( d$CDT, d$S ), quantile )
  colnames( ag_rsd ) = c( 'Cnd', 'S', 'Q' )
  m_rsd = aggregate( rsd, list( d$CDT ), mean )
  
  ttl = c( 'Target primed .05 s',
           'Target primed .4 s',
           'Foil primed .05 s',
           'Foil primed .4 s' )
  
  yl = lowerUpper( .2, as.vector( ag_rsd$Q ) )
  for ( cnd in 1:8 ) {
    
    blankPlot( c(1,N), yl  )
    axis( 2, seq(yl[1],yl[2],.2) )
    abline( h = 0, lwd = 2, col = 'grey' )
    abline( h = m_rsd$x[cnd], lwd = 2,lty = 2, col = 'blue' )
    
    for (k in 1:4) { if ( cnd == c(1,3,5,7)[k] ) 
      title( ttl[k], cex = .9 ) }
    if ( cnd == 7 ) legend( 'topleft', 'Left correct',
                            bty = 'n' )
    if ( cnd == 8 ) legend( 'topleft', 'Right correct',
                            bty = 'n' )
    
    sel = ag_rsd$Cnd == cnd
    
    segments( 1:N, ag_rsd$Q[ sel, 1 ],
              1:N, ag_rsd$Q[ sel, 5 ] )
    segments( 1:N, ag_rsd$Q[ sel, 2 ],
              1:N, ag_rsd$Q[ sel, 4 ], lwd = 3 )
    points( 1:N, ag_rsd$Q[ sel, 3 ],
            pch = 21, bg = 'white', cex = 1.5 )
    
  }
  mtext('Residuals',side=3,outer=T,line=-2.75,cex=1.5)
  
  if (savePlot) {
    dev.off()
    setwd( orig_dir )
  }
}

setwd( orig_dir )