#---------------------------#
# Examination of posteriors #
# Kevin Potter              #
# Updated 12/15/2016        #
#---------------------------#

# Define model type
type = 4

# Load in estimation results for hierarchical model
source( 'F7_BE_MS.R' )

# Indicate which code segments to run
runCode = c(T, T, T )

# Indicate whether to save a pdf file
savePlot = T
if (savePlot) {
  setwd( 'Plots' )
  fName = paste('Model', type, 'fit_results.pdf', sep = '_' )
  pdf( fName, width = 12, height = 6 )
  rm( fName )
  setwd( orig_dir )
}

# Index
# Lookup - 01:  Convergence check
# Lookup - 02:  Retrodictive checks
# Lookup - 03:  Marginal posterior distributions for hyper-parameters

###
### Convergence check
###
# Lookup - 01

if ( runCode[1] ) {
  
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

###
### Retrodictive checks
###
# Lookup - 02

if ( runCode[2] ) {
  
  # Restrict to forced choice
  cD = d[ d$Ta == 0, ]
  
  # Extract the total number of posterior samples
  S = conv$totSampSize
  
  # Create a set of time-points
  tp = seq( 0, 2, length = 40 )
  
  # Create a list to store samples for the retrodictive checks
  RC = list()
  for ( cnd in 1:8 ) {
    RC[[cnd]] = list(
      ch1 = matrix( NA, length(tp), S ),
      ch0 = matrix( NA, length(tp), S )
    )
  }
  
  # Generate posterior retrodictive checks
  startTime = Sys.time()
  
  # Create a progress bar using a base R function
  pb = txtProgressBar( min = 1, max = S, style = 3 )
  
  # Create an array to store parameter estimates per condition and sample
  prm_all = array( NA, dim = c( 8, 8, S ) )
  
  # Simulate from the hierarchy and generate per condition estimates
  for (s in 1:S ) {
    prm_all[,,s] = gl_sim(s,post,input,type)
    
    # Update the progress bar
    setTxtProgressBar(pb,s)
    
  }
  close(pb)
  
  # Calculate the distribution function over samples
  
  # Create a progress bar using a base R function
  pb = txtProgressBar( min = 1, max = 8*2*length(tp), style = 3 )
  
  inc = 1
  for ( cnd in 1:8 ) {
    for ( ch in 1:0 ) {
      
      for (i in 1:length(tp)) {
        mp = pwaldrace( tp[i], ch,
                        prm_all[1,cnd,], 
                        prm_all[2,cnd,], 
                        prm_all[4,cnd,], 
                        prm_all[5,cnd,], 
                        prm_all[6,cnd,], 
                        prm_all[8,cnd,] )
        if ( ch == 1 ) RC[[cnd]]$ch1[i,] = mp else 
          RC[[cnd]]$ch0[i,] = mp
        inc = inc + 1
        # Update the progress bar
        setTxtProgressBar(pb,inc)
      }
      
    }
    
  }
  close(pb)
  runTime = Sys.time() - startTime
  print( runTime )
  
  if ( type == 1 ) {
    
    # Calculate residual latency
    theta = post$theta_alpha/(post$theta_alpha+post$theta_beta)
    
    # Determine MAP estimates for group-level
    cf = c(
      apply( post$kappa_mu, 2, findMode ), 
      apply( post$xi_mu, 2, findMode ),
      findMode( theta )*mean( input$min_RT ) )
    prm = param_est( input$X_small, cf, input$fixed, input$index,
                     input$parSel )
    
  }
  
  if ( type == 2 ) {
    
    # Calculate residual latency
    theta = post$theta_alpha/(post$theta_alpha+post$theta_beta)
    
    # Determine MAP estimates for group-level
    cf = c(
      findMode( post$kappa_mu ), 
      apply( post$xi_mu, 2, findMode ),
      findMode( theta )*mean( input$min_RT ) )
    prm = param_est( input$X_small, cf, input$fixed, input$index,
                     input$parSel )
    
  }
  
  if ( type == 3 ) {
    
    # Calculate residual latency
    theta = post$theta_alpha/(post$theta_alpha+post$theta_beta)
    
    # Determine MAP estimates for group-level
    cf = c(
      apply( post$kappa_mu, 2, findMode ), 
      apply( post$xi_mu, 2, findMode ),
      findMode( theta )*mean( input$min_RT ) )
    prm = param_est( input$X_small, cf, input$fixed, input$index,
                     input$parSel )
    
  }
  
  if ( type == 4 ) {
    
    # Calculate residual latency
    theta = post$theta_alpha/(post$theta_alpha+post$theta_beta)
    
    tmp = findMode( post$beta_mu[,1] ) + 
      findMode( post$beta_mu[,2] ) * colMeans( input$zIL )
    
    # Determine MAP estimates for group-level
    cf = c(
      apply( post$kappa_mu, 2, findMode ), 
      tmp,
      findMode( theta )*mean( input$min_RT ) )
    prm = param_est( input$X_small, cf, input$fixed, input$index,
                     input$parSel )
  }
  
  if (!savePlot) x11(width=12)
  layout( rbind( c(1,3,5,7), c(2,4,6,8) ) )
  
  ttl = c( 'Target primed .05 s',
           'Target primed .4 s',
           'Foil primed .05 s',
           'Foil primed .4 s' )
  
  for ( cnd in 1:8 ) {
    
    blankRTplot( bty='l', cex.axis = 1.5, cex.lab = 1.5 )
    
    for (k in 1:4) { if ( cnd == c(1,3,5,7)[k] ) title( ttl[k], cex = .9 ) }
    
    if ( cnd == 7 ) legend( 'topleft', 'Left correct',
                            bty = 'n' )
    if ( cnd == 8 ) legend( 'topleft', 'Right correct',
                            bty = 'n' )
    
    l = length( tp )
    ya = matrix( NA, 2, l )
    for (i in 1:l ) ya[,i] = quantile( RC[[cnd]]$ch1[i,], c(.16,.84) )
    polygon( c(tp,rev(tp)), c(ya[1,],rev(ya[2,])), col = rgb(0,1,0,.3),
             border = NA )
    for (i in 1:l ) ya[,i] = quantile( RC[[cnd]]$ch0[i,], c(.16,.84) )
    polygon( c(tp,rev(tp)), c(ya[1,],rev(ya[2,])), col = rgb(1,0,0,.3),
             border = NA )
    
    T_x = function(x) quantile(x,seq(.1,.9,.2))
    
    rt = cD$RT[ cD$Cnd == cnd ];
    ch = cD$Ch[ cD$Cnd == cnd ];
    grp = cD$S[ cD$Cnd == cnd ];
    out1 = cdf_curve( rt, ch, grp = grp, lwd = 2, 
                      opt = list( out = T ) )
    add_points( out1, T_x = T_x, pch = 19, cex = 1.5 )
    out0 = cdf_curve( rt, ch, sel = 0, grp = grp, lwd = 2, lty = 2, 
                      opt = list( out = T ) )
    add_points( out0, T_x = T_x, pch = 17, cex = 1.5 )
    
    dist_plot( prm[c(1,2,4,5,6,8),cnd], lwd = 2, col = 'green' )
    dist_plot( prm[c(1,2,4,5,6,8),cnd], ch = 0, lty = 2, lwd = 2, 
               col = 'red' )
    
  }
  
}

###
### Marginal posterior distributions for hyper-parameters
###
# Lookup - 03

if ( runCode[3] ) {
  
  # Define info for prior distributions
  if ( type == 1 ) {
    # Prior index
    priorIndex = list( 
      I = c( rep( 1, 5 ),
             rep( 2, 4 ),
             3 ),
      V = c( rep( 1, 5 ),
             rep( 1, 4 ),
             3 ) )
    lbls = c( 'Short', 'Long', 'Prime bias (short)', 'Prime bias (long)', 'Side bias', 
              'Target (.05 s)', 'Target (.4 s)',
              'Foil (.05 s)', 'Foil (.4 s)',
              ' ' )
  }
  
  if ( type == 2 ) {
    # Prior index
    priorIndex = list( 
      I = c( rep( 1, 1 ),
             rep( 2, 8 ),
             3 ),
      V = c( rep( 1, 1 ),
             rep( 1, 8 ),
             3 ) )
    
    lbls = c( ' ',
              'Target (STP)', 'Target (LTP)', 'Target (SFP)', 'Target (LFP)',
              'Foil (STP)', 'Foil (LTP)', 'Foil (SFP)', 'Foil (LFP)',
              ' ' )
    
  }
  
  if ( type == 3 ) {
    # Prior index
    priorIndex = list( 
      I = c( rep( 1, 5 ),
             rep( 2, 8 ),
             3 ),
      V = c( rep( 1, 5 ),
             rep( 1, 8 ),
             3 ) )
    
    lbls = c( 'Short', 'Long', 'Prime bias (short)', 'Prime bias (long)', 'Side bias', 
              'Target (STP)', 'Target (LTP)', 'Target (SFP)', 'Target (LFP)',
              'Foil (STP)', 'Foil (LTP)', 'Foil (SFP)', 'Foil (LFP)',
              ' ' )
    
  }
  
  # Define info for prior distributions
  if ( type == 4 ) {
    # Prior index
    priorIndex = list( 
      I = c( rep( 1, 5 ),
             rep( 2, 2 ),
             3 ),
      V = c( rep( 1, 5 ),
             rep( 1, 2 ),
             3 ) )
    lbls = c( 'Short', 'Long', 'Prime bias (short)', 'Prime bias (long)', 'Side bias', 
              'Intercept','Slope',
              ' ' )
  }
  
  
  ### Posterior distribution estimates ####
  
  lp = 1:3
  if ( type == 4 ) lp = c(1,4,3)
  inc = 0
  for (j in lp) {
    
    if ( j == 1 ) {
      pst = post$kappa_mu
      if ( type == 2 ) pst = cbind( as.numeric( pst ) )
      ttl = 'Threshold'
    }
    if ( j == 2 ) {
      pst = post$xi_mu
      ttl = 'Drift rates'
    }
    if ( j == 3 ) {
      pst = cbind( as.numeric( post$theta_alpha / 
                                 (post$theta_alpha+post$theta_beta) ) )
      ttl = 'Proportion of fastest RT'
    }
    if ( j == 4 ) {
      pst = post$beta_mu
      ttl = 'Regression coefficients'
    }
    
    if (!savePlot) x11(width=12)
    layout( 1 )
    yl = lowerUpper( .5, as.vector( pst ) )
    blankPlot( xDim=c( .5, .5+dim(pst)[2] ),
               yDim = yl )
    abline( h = 0, lwd = 2 )
    abline( v = .5, lwd = 2 )
    axis( 2, round( seq( yl[1], yl[2], length = 4 ), 2 ),
          tick = F, cex.axis = 1.5, line = -.5 )
    title(ttl,cex=1.5)
    
    # Determine heights of each marginal posterior density
    h = apply( pst, 2, function(x) max( density(x)$y ) )
    
    # Credible intervals
    for ( i in 1:ncol(pst) ) {
      
      inc = inc + 1
      
      # Determine prior distribution
      pSel = which( priorIndex$I == j )
      pCol = c(1,2); if ( j == 3 ) pCol = c(1,3)
      prs = prior_f( Priors[ pSel[i], pCol ], priorIndex$V[inc], 
                     v = seq(yl[1],yl[2],length=1000) )
      dn = .4*h[i]/max(h)
      p = .4*max( prs$y )/max(h)
      
      # Overlay prior distributions
      violinPlot( prs, i, scaleH = p, est = F, lty = 2, lwd = 2 )
      
      violinPlot( pst[,i], i, scaleH = dn, col = 'white', lwd = 2 )
      ci = quantile( pst[,i], c(.16,.84) )
      violinPlot( pst[,i], i, scaleH = dn,
                  interval = list( crit = ci ), 
                  col = 'grey', border = NA )
      violinPlot( pst[,i], i, scaleH = dn, lwd = 2 )
      
      axis( 1, i, lbls[inc], tick = F )
    }
    
    
  }
  
  # For model 3, assess whether the short/long threshold posteriors differ
  if ( type == 3 ) {
    
    pst = post$kappa_mu[,1] - post$kappa_mu[,2]
    if (!savePlot) x11(width=12);
    layout( cbind(1,2) )
    
    tmp = hist( pst, col = 'grey', border = 'white', breaks = 40, freq = F,
                xlab = 'Difference scores', main = 'Threshold: Short - Long' )
    
    # Calculate an uncertainty interval
    ui = quantile( pst, c(.16,.84) )
    h = c( max( which( tmp$mids <= ui[1] ) ),
           min( which( tmp$mids >= ui[2] ) ) )
    segments( ui, c(0,0), ui, tmp$density[h], lwd = 2 )
    blankPlot()
  }
  
}

if (savePlot) {
  setwd( orig_dir )
  setwd( 'Plots' )
  dev.off()
}

setwd( orig_dir )