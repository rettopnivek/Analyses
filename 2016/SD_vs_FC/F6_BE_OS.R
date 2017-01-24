#---------------------------#
# Sequential sampling model #
# using BE for one subject  #
# Kevin Potter              #
# 12/14/2016                #
#---------------------------#

# Initialize script
source('F4_starting_script.R')

# Load in data
setwd( 'Data' )
load( 'SD_v_FC.RData' )
setwd( orig_dir )

# For easy manipulation
d = allDat
colnames( d ) = c('S','Cnd','Ac','RT','E','PD',
                  'Ta','Co','PT','Ch','TL','FL','zTL','zFL')

# Indicate whether plots should be generated
plotYes = T
# Determine subject to fit
sbj = 1
# Determine type of model to fit
tp = 4

# Index
# Lookup - 01:  Define type of model to fit
#          01a: Model 1
#          01b: Model 2
#          01c: Model 3
#          01d: Model 4
# Lookup - 02:  Estimate model parameters for single subject
# Lookup - 03:  Plot model results

###
### Define type of model to fit
###
# Lookup - 01

# Define a function to check if initial values are valid
initCheck = function( curInit ) {
  
  nCh = length( curInit )
  for (nc in 1:nCh) {
    cf = as.vector( unlist( curInit[[nc]] ) )
    
    prm = param_est( input$X, cf, input$fixed, 
                     input$index, input$parSel )
    prm[4,] = prm[4,]*input$min_RT
    prm[8,] = prm[8,]*input$min_RT
    chk = dwaldrace( input$Y[,1], input$Y[,2],
                     prm[1,], prm[2,], prm[4,],
                     prm[5,], prm[6,], prm[8,], ln = 1 )
    print( sum( chk ) )
  }
  
}

### Model 1
# Lookup - 01a
if ( tp == 1 ) {
  
  # Define set of priors
  Priors = cbind(
    c( 1, 1, 0, 0, 0, # kappa
       rep(2.0,2), rep(1.4,2), # xi
       7 # theta
    ),
    c( .25, .25, .1, .1, .1, # kappa
       rep(0.5,2), rep(.5,2), # xi
       3 # theta
    )
  )
  
  # Prior index
  priorIndex = list( 
    I = c( rep( 1, 5 ),
    rep( 2, 4 ),
    3 ),
    V = c( rep( 1, 5 ),
       rep( 1, 4 ),
       3 ) )
  
  # Define function to initialize values
  init_f = function(nCh) {
    
    # Initialize output
    out = c()
    for (nc in 1:nCh) out = c( out, list() )
    
    # Generate initial values
    for (nc in 1:nCh) {
      
      out[[nc]] = list(
        kappa_S = runif( 1, .5, 2 ),
        kappa_L = runif( 1, .5, 2 ),
        kappa_SP = runif( 1, -.2, .2 ),
        kappa_LP = runif( 1, -.2, .2 ),
        kappa_B = runif( 1, -.2, .2 ),
        xi = runif( 4, .5, 4 ),
        theta = runif( 1, .5, .8 )
      )
      
    }
    
    return( out )
  }
  
}

### Model 2
# Lookup - 01b
if ( tp == 2 ) {
  
  # Define set of priors
  Priors = cbind(
    c( 1, # kappa
       rep(2.0,4), rep(1.4,4), # xi
       7 # theta
    ),
    c( .25, # kappa
       rep(0.5,4), rep(.5,4), # xi
       3 # theta
    )
  )
  
  # Prior index
  priorIndex = list( 
    I = c( rep( 1, 1 ),
           rep( 2, 8 ),
           3 ),
    V = c( rep( 1, 1 ),
           rep( 1, 8 ),
           3 ) )
  
  
  # Define function to initialize values
  init_f = function(nCh) {
    
    # Initialize output
    out = c()
    for (nc in 1:nCh) out = c( out, list() )
    
    # Generate initial values
    for (nc in 1:nCh) {
      
      out[[nc]] = list(
        kappa = runif( 1, .5, 1 ),
        xi = runif( 8, .5, 4 ),
        theta = runif( 1, .5, .8 )
      )
      
    }
    
    return( out )
  }
  
}

# Model 3
# Lookup - 01c
if ( tp == 3 ) {
  
  # Define set of priors
  Priors = cbind(
    c( 1, 1, 0, 0, 0, # kappa
       rep(2.0,4), rep(1.4,4), # xi
       7 # theta
    ),
    c( .25, .25, .1, .1, .1, # kappa
       rep(0.5,4), rep(.5,4), # xi
       3 # theta
    )
  )
  
  # Prior index
  priorIndex = list( 
    I = c( rep( 1, 5 ),
           rep( 2, 8 ),
           3 ),
    V = c( rep( 1, 5 ),
           rep( 1, 8 ),
           3 ) )
  
  # Define function to initialize values
  init_f = function(nCh) {
    
    # Initialize output
    out = c()
    for (nc in 1:nCh) out = c( out, list() )
    
    # Generate initial values
    for (nc in 1:nCh) {
      
      out[[nc]] = list(
        kappa_S = runif( 1, .5, 2 ),
        kappa_L = runif( 1, .5, 2 ),
        kappa_SP = runif( 1, -.2, .2 ),
        kappa_LP = runif( 1, -.2, .2 ),
        kappa_B = runif( 1, -.2, .2 ),
        xi = runif( 8, .5, 4 ),
        theta = runif( 1, .5, .8 )
      )
      
    }
    
    return( out )
  }
  
}

# Model 4
# Lookup - 01d
if ( tp == 4 ) {
  
  # Define set of priors
  Priors = cbind(
    c( 1, 1, 0, 0, 0, # kappa
       0, 0, 3, # Regression parameters
       7 # theta
    ),
    c( .25, .25, .1, .1, .1, # kappa
       .5, .5, 3, # Regression parameters
       3 # theta
    )
  )
  
  # Prior index
  priorIndex = list( 
    I = c( rep( 1, 5 ),
           c(1,1,2),
           3 ),
    V = c( rep( 1, 5 ),
           c(1,1,2),
           3 ) )
  
  # Define function to initialize values
  init_f = function(nCh) {
    
    # Initialize output
    out = c()
    for (nc in 1:nCh) out = c( out, list() )
    
    # Generate initial values
    for (nc in 1:nCh) {
      
      out[[nc]] = list(
        kappa_S = runif( 1, .5, 2 ),
        kappa_L = runif( 1, .5, 2 ),
        kappa_SP = runif( 1, -.2, .2 ),
        kappa_LP = runif( 1, -.2, .2 ),
        kappa_B = runif( 1, -.2, .2 ),
        xi = runif( 8, .5, 4 ),
        rho = runif( 1, -.5, .5 ),
        sigma_xi = runif( 1, .25, 2 ),
        theta = runif( 1, .5, .8 )
      )
      
    }
    
    return( out )
  }
  
}

###
### Estimate model parameters for single subject
###
# Lookup - 02

# Extract data and covariates
input = model_structure_create( tp, d, Priors, s = sbj )

# Define aspects of input to pass into Stan
listVal = 1:(length(input)-1)

warm = 250 # Warm-up
niter = 1250 # Number of samples to approximate posterior
chains = 8 # Number of chains to run

setwd('Stan_scripts')

curInit = init_f( chains )
startTime = Sys.time() # To assess run-time
fit = stan( file = input$mName, data = input[listVal], 
            warmup = warm, iter = warm+niter, 
            chains = chains,
            init = curInit,
            control = list( adapt_delta = .95 ) )

post = extract(fit)
# Report run time
runTime = Sys.time() - startTime
print( runTime )
rm( startTime )

###
### Plot model results
###
# Lookup - 03

if ( plotYes ) {
  
  if ( tp == 1 | tp == 3 | tp == 4 ) {
    # Collapse threshold parameters to single variable
    post$kappa = cbind( post$kappa_S, post$kappa_L, 
                        post$kappa_SP, post$kappa_LP,
                        post$kappa_B )
    post$theta = cbind( post$theta )
  }
  
  if ( tp == 2 ) {
    post$kappa = cbind( post$kappa )
    post$theta = cbind( post$theta )
  }
  
  ### Convergence check ###
  
  cnv = convergence_extract( fit, 'tau' )
  x11(width=12);
  layout(cbind(1,2))
  tmp = hist( cnv$Rhat, col='grey', border='white', 
              xlab = expression( hat(R) ),
              ylab = 'Frequency', bty = 'l', main = ' ',
              xlim = c(.95,1.15) );
  segments( tmp$mids, rep(0,length(tmp$mids)),
            tmp$mids, tmp$counts, col = 'grey', lwd = 3 )
  abline(v=1.1)
  hist( cnv$n_eff/cnv$totSampSize, col='grey', border='white', 
        xlab = 'Percentage of effective samples',
        ylab = 'Frequency', bty = 'l', main = ' ',
        xlim = c(0,1) )
  mtext( 'Convergence check', outer = T, line = -2, cex = 1.5 )
  
  ### Retrodictive checks ###
  
  retroCheck = T
  if ( retroCheck ) {
    
    S = 1000
    
    # Create a progress bar using a base R function
    pb = txtProgressBar( min = 1, max = S, style = 3 )
    
    # Create an array to store parameter estimates per condition and sample
    prm_all = array( NA, dim = c( 8, 8, S ) )
    
    ranP = sample( 1:cnv$totSampSize, S )
    
    # Draw samples from the posterior and generate per condition estimates
    for (s in 1:S ) {
      
      prm_all[,,s] = gl_sim(ranP[s],post,input,tp,H=F)
      
      # Update the progress bar
      setTxtProgressBar(pb,s)
      
    }
    close(pb)
    
    # Create a set of time-points
    ts = seq( 0, 4, length = 50 )
    
    # Create a list to store samples for the retrodictive checks
    RC = list()
    for ( cnd in 1:8 ) {
      RC[[cnd]] = list(
        ch1 = matrix( NA, length(ts), S ),
        ch0 = matrix( NA, length(ts), S )
      )
    }
    
    # Create a progress bar using a base R function
    pb = txtProgressBar( min = 1, max = 8*2*length(ts), style = 3 )
    
    inc = 1
    for ( cnd in 1:8 ) {
      for ( ch in 1:0 ) {
        
        for (i in 1:length(ts)) {
          mp = pwaldrace( ts[i], ch,
                          prm_all[1,cnd,], 
                          prm_all[2,cnd,], 
                          prm_all[4,cnd,], 
                          prm_all[5,cnd,], 
                          prm_all[6,cnd,], 
                          prm_all[8,cnd,],
                          1, 1, 1 )
          if ( ch == 1 ) RC[[cnd]]$ch1[i,] = mp else 
            RC[[cnd]]$ch0[i,] = mp
          inc = inc + 1
          # Update the progress bar
          setTxtProgressBar(pb,inc)
        }
        
      }
      
    }
    close(pb)
    
    cf = c( apply( post$kappa, 2, findMode ),
            apply( post$xi, 2, findMode ),
            apply( post$theta, 2, findMode )*input$min_RT )
    #if ( tp == 4 ) {
    #  cf[6:13] = findMode(post$beta[,1]) + findMode(post$beta[,2])*input$zIL
    #}
    map = param_est( input$X_small, cf, input$fixed, input$index, input$parSel )
    
    # Plot predicted against observed per condition
    x11(width=12)
    
    layout( rbind( c(1,3,5,7), c(2,4,6,8) ) )
    
    ttl = c( 'Target primed .05 s',
             'Target primed .4 s',
             'Foil primed .05 s',
             'Foil primed .4 s' )
    
    # Loop over conditions
    for ( cnd in 1:8 ) {
      
      par( mar = c( 5, 4, 1, 3 ) )
      
      xl = lowerUpper( .2, d$RT[ d$S == sbj ] )
      
      blankPlot(xDim=xl)
      
      if ( i == 1 | i == 5) {
        mtext('Cumulative probability',side=2,line=2,cex=1.3)
      }
      
      for (k in 1:4) { if ( cnd == c(1,3,5,7)[k] ) title( ttl[k], cex = .9 ) }
      
      if ( cnd == 7 ) legend( 'topleft', 'Left correct',
                              bty = 'n' )
      if ( cnd == 8 ) legend( 'topleft', 'Right correct',
                              bty = 'n' )
      
      mtext('RT (s)', side = 1, line = 2.5, cex = 1 )
      axis( 1, seq( ceiling(xl[1]),floor(xl[2]),.5), cex.axis = 1.25 )
      axis( 4, seq(0,1,.25), cex.axis = 1.25 )
      
      l = length( ts )
      ya = matrix( NA, 2, l )
      for (i in 1:l ) ya[,i] = quantile( RC[[cnd]]$ch1[i,], c(.16,.84) )
      polygon( c(ts,rev(ts)), c(ya[1,],rev(ya[2,])), col = rgb(0,1,0,.3),
               border = NA )
      for (i in 1:l ) ya[,i] = quantile( RC[[cnd]]$ch0[i,], c(.16,.84) )
      polygon( c(ts,rev(ts)), c(ya[1,],rev(ya[2,])), col = rgb(1,0,0,.3),
               border = NA )
      
      T_x = function(x) quantile(x,seq(.1,.9,.2))
      
      sel = d$Ta == 0 & d$S == sbj & d$Cnd == cnd
      rt = d$RT[sel]; ch = d$Ch[sel]
      # Observed data
      o1 = cdf_curve(rt,ch,opt=list(out=T),lwd=2);
      o0 = cdf_curve(rt,ch,sel=0,opt=list(out=T),lwd=2,lty=2);
      
      add_points( o1, T_x = T_x, pch = 19, cex = 1.5 )
      add_points( o0, T_x = T_x, pch = 17, cex = 1.5 )
      
      dist_plot( c( map[c(1,2,4,5,6,8,3,7),cnd], 1 ),
                 rt = seq(xl[1],xl[2],length=1000), 
                 lwd = 2, col = 'green' )
      dist_plot( c( map[c(1,2,4,5,6,8,3,7),cnd], 1 ),
                 rt = seq(xl[1],xl[2],length=1000), ch = 0, 
                 lwd = 2, col = 'red', lty = 2 )
      
      
    }
    
  }
  
  ### Posterior distribution estimates ####
  
  lp = 1:3
  for (j in lp) {
    
    if ( j == 1 ) {
      pst = post$kappa
      ttl = 'Threshold'
    }
    if ( j == 2 ) {
      pst = post$xi
      ttl = 'Drift rates'
    }
    if ( j == 3 ) {
      pst = cbind( as.numeric( post$theta ) )
      ttl = 'Proportion of fastest RT'
    }
    
    x11(width=12)
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
      
      # Determine prior distribution
      pSel = which( priorIndex$I == j )
      prs = prior_f( Priors[ pSel[i], ], priorIndex$V[i], 
                     v = seq(yl[1],yl[2],length=1000) )
      dn = .4*h[i]/max(h)
      p = .4*max( prs$y )/max(h)
      
      # Overlay prior distributions
      if ( tp < 4 | ( tp == 4 & j != 2 ) ) {
        violinPlot( prs, i, scaleH = p, est = F, lty = 2, lwd = 2 )
      }
      
      violinPlot( pst[,i], i, scaleH = dn, col = 'white', lwd = 2 )
      ci = quantile( pst[,i], c(.16,.84) )
      violinPlot( pst[,i], i, scaleH = dn,
                  interval = list( crit = ci ), 
                  col = 'grey', border = NA )
      violinPlot( pst[,i], i, scaleH = dn, lwd = 2 )
      
    }
    
    
  }
  
  x11(width=12)
  lyt = matrix( 1, 10, 18 )
  lyt[,10:18] = 3;
  lyt[1,] = rep( c(2,4), each = 9 )
  layout( lyt )
  
  # Extract maximum a posteriori estimates (posterior modes)
  if ( tp == 1 | tp == 3 ) {
    prm = c(
      findMode( post$kappa_S ),
      findMode( post$kappa_L ),
      findMode( post$kappa_SP ),
      findMode( post$kappa_LP ),
      findMode( post$kappa_B ),
      apply( post$xi, 2, findMode ),
      findMode( post$theta )*input$min_RT )
  }
  
  if ( tp == 2 ) {
    prm = c(
      findMode( post$kappa ),
      apply( post$xi, 2, findMode ),
      findMode( post$theta )*input$min_RT )
  }
  
  if ( tp == 4 ) {
    prm = c(
      findMode( post$kappa_S ),
      findMode( post$kappa_L ),
      findMode( post$kappa_SP ),
      findMode( post$kappa_LP ),
      findMode( post$kappa_B ),
      apply( post$xi, 2, findMode ),
      # findMode( post$beta[,1] ) + findMode( post$beta[,2] )*input$zIL, 
      findMode( post$theta )*input$min_RT )
  }
  
  # Compute estimates per condition
  tst = param_est( input$X_small, prm, input$fixed, input$index, input$parSel )
  
  # Plot threshold values per condition
  yl = lowerUpper( .2, as.vector( tst[c(1,5),] ) )
  par( mar = c( 5.5, 4, 0, 0 ) )
  blankPlot( c(1,8), yl )
  axis( 1, 1:8, rep( c('L','R'), 4 ), tick = F, line = 0, cex.axis = 1.5 )
  axis( 1, c(1.5,3.5,5.5,7.5), rep( c('Short','Long'), 2 ), 
        tick = F, line = 1.5, cex.axis = 1.5 )
  axis( 1, c(2.5,6.5), c('Target primed','Foil primed'), 
        tick = F, line = 3, cex.axis = 1.5 )
  axis( 2, seq( yl[1], yl[2], .2 ) )
  mtext( 'Threshold', side = 2, cex = 1.5, line = 2 )
  points( 1:8, tst[1,], pch = 21, bg = 'grey', cex = 1.5 )
  points( 1:8, tst[5,], pch = 24, bg = 'white', cex = 1.5 )
  par( mar = c(0,4,0,0) )
  blankPlot()
  legend( 'top', c( expression( kappa[R] ), expression( kappa[L] ) ),
          pt.bg = c('grey','white'), pch = c(21,24), cex = 1.5,
          horiz = T, bty = 'n' )
  
  # Plot threshold values per condition
  yl = lowerUpper( .2, as.vector( tst[c(2,6),] ) )
  par( mar = c( 5.5, 4, 0, 0 ) )
  blankPlot( c(1,8), yl )
  axis( 1, 1:8, rep( c('L','R'), 4 ), tick = F, line = 0, cex.axis = 1.5 )
  axis( 1, c(1.5,3.5,5.5,7.5), rep( c('Short','Long'), 2 ), 
        tick = F, line = 1.5, cex.axis = 1.5 )
  axis( 1, c(2.5,6.5), c('Target primed','Foil primed'), 
        tick = F, line = 3, cex.axis = 1.5 )
  axis( 2, seq( yl[1], yl[2], .2 ) )
  mtext( 'Drift rate', side = 2, cex = 1.5, line = 2 )
  points( 1:8, tst[2,], pch = 21, bg = 'grey', cex = 1.5 )
  points( 1:8, tst[6,], pch = 24, bg = 'white', cex = 1.5 )
  par( mar = c(0,4,0,0) )
  blankPlot()
  legend( 'top', c( expression( xi[R] ), expression( xi[L] ) ),
          pt.bg = c('grey','white'), pch = c(21,24), cex = 1.5,
          horiz = T, bty = 'n' )
  
}

setwd( orig_dir )