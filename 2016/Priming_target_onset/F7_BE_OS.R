#---------------------------#
# Sequential sampling model #
# using BE for one subject  #
# Kevin Potter              #
# 11/15/2016                #
#---------------------------#

# Initialize script
source('F5_starting_script.R')

# Load in data
setwd( 'Data' )
load( 'Priming_offset.RData' )
setwd( orig_dir )

# For easy manipulation
d = allDat
colnames( d ) = c('S','PD','PT','O','Co','RT',
                  'Ch','Ac','OT','DT','CDT','Cnd')

# Indicate whether plots should be generated
plotYes = T

###
### Estimate model parameters for single subject
###

# Define set of priors
Priors = cbind(
  c( rep(.85,2), rep(2.0,4), rep(1.4,4), rep(7,10) ),
  c( rep(.3,2), rep(.5,4), rep(.5,4), rep(3,10) )
)
# Define index for priors
priorIndex = c( rep(1,2), rep(2,8), rep(3,10) )
# Define distributions for priors for plotting purposes
prior_f = function( priors, type, v = seq(0,4,length=1000) ) {
  
  if ( type == 1 | type == 2 ) {
    out = list( x = v, y = dnorm( v, priors[1], priors[2] ) )
  }
  if ( type == 3 ) {
    if ( max(v) > 1 | min(v) < 0 ) v = seq(0,1,length=1000)
    out = list( x = v, y = dbeta( v, priors[1], priors[2] ) )
  }
  
  return( out )
}

# Extract data and covariates
input = model_structure_create( 1, d, Priors, s = 2 )

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

startTime = Sys.time() # To assess run-time
fit = stan(file = 'WR_OS.stan', data = stanDat, 
           warmup = warm, iter = warm+niter, 
           chains = chains)

post = extract(fit)
# Report run time
runTime = Sys.time() - startTime
print( runTime )
rm( startTime )

if ( plotYes ) {
  
  # Convergence check
  cnv = convergence_extract( fit, 'tau[2]' )
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
  
  for (j in 1:3) {
    
    if ( j == 1 ) {
      pst = post$kappa
      ttl = 'Threshold'
    }
    if ( j == 2 ) {
      pst = post$xi
      ttl = 'Drift rates'
    }
    if ( j == 3 ) {
      pst = post$theta
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
      pSel = which( priorIndex == j )
      prs = prior_f( Priors[ pSel[i], ], j, 
                     v = seq(yl[1],yl[2],length=1000) )
      d = .4*h[i]/max(h)
      p = .4*max( prs$y )/max(h)
      
      # Overlay prior distributions
      violinPlot( prs, i, scaleH = p, est = F, lty = 2, lwd = 2 )
      
      violinPlot( pst[,i], i, scaleH = d, col = 'white', lwd = 2 )
      ci = quantile( pst[,i], c(.16,.84) )
      violinPlot( pst[,i], i, scaleH = d,
                  interval = list( crit = ci ), 
                  col = 'grey', border = NA )
      violinPlot( pst[,i], i, scaleH = d, lwd = 2 )
      
    }
    
    
  }
  
}

setwd( orig_dir )

