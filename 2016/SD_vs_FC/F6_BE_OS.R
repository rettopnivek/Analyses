#---------------------------#
# Sequential sampling model #
# using BE for one subject  #
# Kevin Potter              #
# 12/01/2016                #
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
sbj = 10
# Determine type of model to fit
tp = 1

# Index
# Lookup - 01:  Define type of model to fit
#          01a: Model 1
#          01b: Model 2
#          01c: Model 3
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

###
### Estimate model parameters for single subject
###
# Lookup - 02

# Extract data and covariates
input = model_structure_create( tp, d, Priors, s = sbj )

# Input for Stan
stanDat = list(
  N = input$N,
  V = input$V,
  K = input$K,
  U = input$U, 
  C = input$C, 
  X = input$X,
  fixed = input$fixed,
  index = input$index,
  parSel = input$parSel,
  Y = input$Y,
  min_RT = input$min_RT,
  Priors = input$Priors
)

warm = 250 # Warm-up
niter = 1250 # Number of samples to approximate posterior
chains = 8 # Number of chains to run

setwd('Stan_scripts')

curInit = init_f( chains )
startTime = Sys.time() # To assess run-time
fit = stan( file = input$mName, data = stanDat, 
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
  
  if ( tp == 1 | tp == 3 ) {
    # Collapse threshold parameters to single variable
    post$kappa = cbind( post$kappa_S, post$kappa_L, 
                        post$kappa_SP, post$kappa_LP,
                        post$kappa_B )
  }
  
  if ( tp == 2 ) {
    post$kappa = cbind( post$kappa )
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
    
    # Extract posterior samples for the set of 8 conditions
    RC = array( NA, dim = c( nrow(post$kappa), 8, 8 ) )
    for ( rc in 1:dim(RC)[1] ) {
      cf = c(
        post$kappa[rc,],
        post$xi[rc,],
        post$tau )
      RC[rc,,] = param_est(input$X_small,cf,input$fixed,input$index,input$parSel)
    }
    
    # Plot predicted against observed per condition
    x11(width=12)
    layout( rbind( 1:4, 5:8 ) )
    
    # Labels for conditions
    lbls = paste(
      rep(c('L','R'),4), # Left/right correct
      rep( rep(c('S','L'),each=2),2), # Short/long prime
      rep(c('T','F'),each=4),sep='') # Target/foil primed
    
    # Loop over conditions
    for ( i in 1:8 ) {
      
      par( mar = c( 5, 4, 1, 3 ) )
      blankPlot(xDim=c(0,2))
      if ( i == 1 | i == 5) {
        mtext('Cumulative probability',side=2,line=2,cex=1.3)
      }
      mtext('RT (s)', side = 1, line = 2.5, cex = 1.3 )
      axis( 1, seq(0,2,.5), cex.axis = 1.25 )
      axis( 4, seq(0,1,.25), cex.axis = 1.25 )
      
      cnd = d$Ta == 0 & d$S == sbj & d$Cnd == i
      rt = d$RT[cnd]; ch = d$Ch[cnd]
      # Observed data
      o1 = cdf_curve(rt,ch,opt=list(out=T),lwd=2);
      o0 = cdf_curve(rt,ch,sel=0,opt=list(out=T),lwd=2,lty=2);
      
      # Credible intervals
      wr_ui( RC, i, sort( o1$i[,1] ), 1, col=rgb(0,0,1,.2),border=NA )
      wr_ui( RC, i, sort( o0$i[,1] ), 0, col=rgb(0,0,1,.2),border=NA )
      
      # Maximum a posteriori estimates
      map = apply( RC[,,i], 2, findMode )
      if ( is.list( map ) ) { # Sometimes there are multiple modes
        tmp = map; map = numeric(8)
        for (j in 1:8) map[j] = tmp[[j]][1]
      }
      # Model estimates
      dist_plot( map[ c(1,2,4,5,6,8,3,7) ], rt = sort(o1$i[,1]), 
                 lwd = 2, col = 'blue' )
      dist_plot( map[ c(1,2,4,5,6,8,3,7) ], ch = 0, rt = sort(o0$i[,1]), 
                 lwd = 2, col = 'blue', lty = 2 )
      
      legend('topleft',lbls[i],bty='n',cex=1.3)
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
      violinPlot( prs, i, scaleH = p, est = F, lty = 2, lwd = 2 )
      
      violinPlot( pst[,i], i, scaleH = dn, col = 'white', lwd = 2 )
      ci = quantile( pst[,i], c(.16,.84) )
      violinPlot( pst[,i], i, scaleH = dn,
                  interval = list( crit = ci ), 
                  col = 'grey', border = NA )
      violinPlot( pst[,i], i, scaleH = dn, lwd = 2 )
      
    }
    
    
  }
  
}

setwd( orig_dir )