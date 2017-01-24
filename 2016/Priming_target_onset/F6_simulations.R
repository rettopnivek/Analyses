#------------------------------------#
# Simulations and parameter recovery #
# Kevin Potter                       #
# Updated 12/01/2016                 #
#------------------------------------#

# Initialize script
source('F4_starting_script.R')

# Define segments of code to run
runCode = c( T, T, T, T, T )

# Indicate whether plots should be saved
savePlot = T
if (savePlot) {
  setwd( 'Plots' )
  pdf(file='Sim_results.pdf',width=12,height=6)
  setwd( orig_dir )
}

# Index
# Lookup - 01:  Test of likelihood and distribution functions
# Lookup - 02:  Parameter recovery for simulated data from 1 subject
# Lookup - 03:  Parameter recovery (over multiple repetitions with MLE)
# Lookup - 04:  Parameter recovery with guessing processes
# Lookup - 05:  Recovery with guessing over multiple repetitions (MLE)

###
### Test of likelihood and distribution functions
###
# Lookup - 01

if (runCode[1]) {
  
  # Define generating parameters
  gp = c( k1 = 1, xi1 = 3, tau1 = .2, 
          k0 = 1, xi0 = 1, tau0 = .15,
          rl = 1 )
  
  # Simulate values from two inverse gaussians
  n = 10000
  t1 = rinvgauss( n, gp[1], gp[2], 1 ) + gp[3]
  t0 = rinvgauss( n, gp[4], gp[5], 1 ) + gp[6]
  # Determine choice
  ch = as.numeric( t1 < t0 )
  # Determine response time
  rt = t1*ch + t0*(1-ch)
  
  if (!savePlot) x11(width=12)
  layout( cbind(1,2) )
  
  # Plot empirical PDF
  blankRTplot( ver = 'PDF', pDim = c(0,4), bty = 'l', 
               cex.axis = 1.5, cex.lab = 1.5 )
  out = pdf_curve( rt, ch, lwd = 2, opt = list( out = T ) )
  add_points( out, T_x = median, pch = 21, bg = 'white', cex  = 1.5 )
  out = pdf_curve( rt, ch, sel = 0, 
                   lwd = 2, opt = list( out = T ) )
  add_points( out, T_x = median, pch = 24, bg = 'white', cex  = 1.5 )
  
  # Plot generating PDF
  dist_plot( prm = gp, ver = 'PDF', 
             lty = 2, lwd = 2, col = 'blue' )
  dist_plot( prm = gp, ch = 0, ver = 'PDF', 
             lty = 2, lwd = 2, col = 'blue' )
  
  legend( 'topleft', c('Simulated','Predicted'),
          fill = c('black','blue'), cex = 1.5, 
          bty = 'n' )
  
  # Plot empirical CDF
  blankRTplot( ver = 'CDF', bty = 'l', cex.axis = 1.5, 
               cex.lab = 1.5 )
  out = cdf_curve( rt, ch, 
             lwd = 2, opt = list( out = T ) )
  add_points( out, T_x = median, pch = 21, bg = 'white', cex  = 1.5 )
  out = cdf_curve( rt, ch, sel = 0, 
             lwd = 2, opt = list( out = T ) )
  add_points( out, T_x = median, pch = 24, bg = 'white', cex  = 1.5 )
  
  dist_plot( prm = gp, lty = 2, lwd = 2, col = 'blue' )
  dist_plot( prm = gp, ch = 0, lty = 2, lwd = 2, col = 'blue' )
  
  mtext( 'Check of PDF/CDF functions', side = 3, outer = T, 
         cex = 1.5, line = -2 )
  
}

###
### Parameter recovery for simulated data from 1 subject
###
# Lookup - 02

if (runCode[2]) {
  
  ### Design ###
  
  # Thresholds           - Speed vs. accuracy parameters
  # Drift rates          - Correct versus incorrect parameters
  # Residual latencies   - Left versus right encoding
  # Coefficient of drift - Fixed to 1
  
  # Generating coefficients
  # k_S k_A (1:2)
  # x_C x_I (2:3)
  # tau_L tau_R (5:6)
  # sigma (Fixed to 1)
  gp = c( 0.8, 1.2, 2.5, 1.2, .27, .3 )
  
  # Define conditions for simulation
  n = 400
  SvA = rep( c(1,0), each = n/2 ) # Speed vs. accuracy conditions
  LvR = rep( c(1,0), each = n/2 ) # Position of correct answer
  LvR = sample( LvR ) # Randomly shuffle position
  int = rep( 1, n ) # Intercept
  
  # Assume a racer for the choice on the left and a racer for the 
  # choice on the right. This then requires matching the drift 
  # rate for the correct choice to the appropriate position.
  # Let 0 indicate picking left, and 1 indicate picking right.
  # We can then specify a design matrix that enables this switching.
  # Hence, the design matrix should have v = 12 variables.
  
  # Design matrix ( v by n )
  X = rbind( SvA, (1-SvA), # kappa (1)
             (1-LvR), LvR, # xi (1)
             int, # sigma (1)
             int, # tau (1)
             SvA, (1-SvA), # kappa (0)
             LvR, (1-LvR), # xi (0)
             int, # sigma (0)
             int  # tau (0)
  )
  
  # Next, we'll create an index for the parameter matrix 
  # rows and columns to use in the linear algebra function.
  # The columns of the parameter matrix are matched to the 
  # variables in the design matrix, while the rows are matched 
  # to the 8 parameters for the wald race model. However, rows 
  # 5 and 11 in the design matrix are for the coefficient of 
  # drift parameters, which are fixed. Fixed values are included 
  # last in the index.
  Clm = c( 1:2,  # Rows in design matrix for kappa(1)
           3:4,  # Rows in design matrix for xi(1)
           6,    # Rows in design matrix for tau(1)
           7:8,  # Rows in design matrix for kappa(0)
           9:10, # Rows in design matrix for xi(0)
           12    # Rows in design matrix for tau(0)
  )
  Rws = c( rep(1,2),  # kappa (1)
           rep(2,2),  # xi (1)
           4,         # tau (1)
           rep(5,2),  # kappa (0)
           rep(6,2),  # xi (0)
           8          # tau (0)
  )
  index = cbind( Rws, Clm )
  
  # Next, we specify the coefficients that correspond with 
  # each non-zero value in the parameter matrix. These indices 
  # are matched to the rows in the index variable.
  parSel = c( 1, 2, 3, 4, 5, 1, 2, 3, 4, 6 )
  
  # Coefficient of drift for both racers is fixed to 1
  fixed = c(1,1)
  index = rbind( index,
                 c(3,5), # sigma (1)
                 c(7,11) # sigma (0)
  )
  rm( Clm, Rws )
  
  # Finally, we include the dimensions (number of parameter P,
  # and the number of coefficients) at the end of the index
  index = rbind( index, c(8,length(parSel)) )
  
  # Generate parameter matrix for simulation purposes
  pm = param_est( X, gp, fixed, index, parSel )
  
  # Simulate data
  sim = rwaldrace( ncol( pm ), pm[1,], pm[2,], pm[4,],
                   pm[5,], pm[6,], pm[8,], 1, 1, rl = 1 )
  
  # Determine accuracy
  sim = cbind( sim, acc = sim[,2] == (1-LvR) )
  # Convert to data frame
  tmp = covCreate( cbind( SvA, LvR ) )
  sim = as.data.frame( cbind( sim, SvA = SvA, LvR = LvR, 
                              Cnd = tmp ) )
  
  if (!savePlot)  x11(width=12)
  layout( cbind(1,2) )
  
  blankRTplot(tDim = c(.4,.8),ver='PvT',bty='l',
              cex.lab = 1.5, cex.axis = 1.5 )
  abline( v = .5, col = 'grey', lwd = 2 )
  plt = list( pch = rep( c( 21, 22, 24, 25 ), 2 ),
              bg = rep( c('black','white'), each = 4 ) )
  pvt_points( sim$rt, sim$acc, sim$Cnd, plt = plt, cex = 1.5 )
  
  blankPlot()
  legend( 'topleft', c('Correct','Error'),
          pch = 21, pt.bg = c('white','black'), bty = 'n',
          cex = 1.5 )
  legend( 'left', c('Cond = L,S',
                       'Cond = L,A',
                       'Cond = R,S',
                       'Cond = R,A'),
          pch = c(21,22,24,25), pt.bg = 'black', bty = 'n',
          cex = 1.5 )
  mtext( 'MRT for simulated multi-condition data',
         side = 3, cex = 1.5, line = -.5 )
  
  # Plot simulation results
  if (!savePlot)  x11(width=12)
  
  layout( rbind( 1:2 ) )
  
  ttl = c('Speed | Left correct','Accuracy | Left correct',
          'Speed | Right correct','Accuracy | Right correct')
  
  for ( i in c(1,3) ) {
    
    par( mar=c(4,1,3,5) )
    blankPlot(xDim=c(.25,1.25),yDim=c(-1,1))
    axis( 4, seq(-1,1,.25), abs( seq(-1,1,.25) ),
          lwd = 2, cex.axis = 1.25 )
    axis( 1, seq(.25,1.25,.5), lwd = 2, cex.axis = 1.25 )
    
    sel = sim$Cnd == i
    out = cdf_curve( sim$rt[sel], sim$acc[sel],
                     opt = list( out = T ), lwd = 2 )
    add_segments( out, lwd = 2 )
    add_points( out, pch = 21, bg='black', cex = 1.5 )
    out = cdf_curve( sim$rt[sel], sim$acc[sel], sel = 0, 
                     opt = list( out = T ), lwd = 2 )
    add_segments( out, lwd = 2 )
    add_points( out, pch = 21, bg='black', cex = 1.5 )
    mtext(ttl[i],side=3,cex=1.5,line=-.5)
    
    sel = sim$Cnd == i+1
    out = cdf_curve( sim$rt[sel], sim$acc[sel],
                     opt = list( out = T, flip = T ), lwd = 2,
                     lty = 2)
    add_segments( out, lwd = 2, lty = 2 )
    add_points( out, pch = 21, bg='white', cex = 1.5 )
    out = cdf_curve( sim$rt[sel], sim$acc[sel], sel = 0, 
                     opt = list( out = T, flip = T ), lwd = 2,
                     lty = 2 )
    add_segments( out, lwd = 2, lty = 2 )
    add_points( out, pch = 21, bg='white', cex = 1.5 )
    mtext(ttl[i+1],side=1,cex=1.5,line=2)
    
    abline( h = 0, lwd = 2 )
  }
  
  ### Parameter estimation using Stan ###
  
  # Define set of priors
  Priors = cbind(
    c( rep(1,2), rep(2.0,2), rep(4, 2 ) ),
    c( rep(.5,2), rep(.5,2), rep(2, 2 ) )
  )
  
  # Input for Stan
  stanDat = list(
    N = ncol(X),
    V = nrow(X),
    K = length(parSel),
    U = length(fixed), 
    C = c(2,2,2), 
    X = X,
    fixed = fixed,
    index = index,
    parSel = parSel,
    Y = sim[,1:2],
    min_RT = rep( min(sim[,1]), 2 ),
    Priors = Priors
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
  
  # Plot of parameter recovery
  if (!savePlot) x11(width = 12)
  layout( cbind(1,2,3) )
  
  # Threshold recovery
  plot( c(.8,2.2), c(.4,1.8), type = 'n', xaxt = 'n',
        ylab = 'Threshold values', xlab = 'Condition',
        bty = 'l' )
  tmp = apply( post$kappa, 2, quantile, prob = c(.025,.5,.975) )
  segments( 1:2, tmp[1,], 1:2, tmp[3,] )
  points( 1:2, tmp[2,], pch = 15 )
  points( 1:2, gp[1:2], pch = 15, col='blue' )
  axis(1,1:2,c('Speed','Accuracy'),tick=F)
  
  legend('topleft',c('Generating','Estimated'),
         fill = c('blue','black'),bty='n')
  
  # Drift rate recovery
  plot( c(.8,2.2), c(.5,4), type = 'n', xaxt = 'n',
        ylab = 'Drift rate values', xlab = ' ',
        bty = 'l' )
  tmp = apply( post$xi, 2, quantile, prob = c(.025,.5,.975) )
  segments( 1:2, tmp[1,], 1:2, tmp[3,] )
  points( 1:2, tmp[2,], pch = 15 )
  points( 1:2, gp[3:4], pch = 15, col='blue' )
  axis(1,1:2,c('Correct','Incorrect'),tick=F)
  
  # Residual latency
  plot( c(.8,2.2), c(.5,1), type = 'n', xaxt = 'n',
        ylab = 'Proportion of shortest time', xlab = ' ',
        bty = 'l' )
  tmp = apply( post$theta, 2, quantile, prob = c(.025,.5,.975) )
  segments( 1:2, tmp[1,], 1:2, tmp[3,] )
  points( 1:2, tmp[2,], pch = 15 )
  points( 1:2, gp[5:6]/min(sim[,1]), pch = 15, col='blue' )
  axis(1,1:2,c('Left','Right'),tick=F)
  
  mtext( 'Parameter recovery', side = 3, outer = T, cex = 1.5,
         line = -2 )
  
  setwd( orig_dir )
}

###
### Parameter recovery (over multiple repetitions with MLE)
###
# Lookup - 03

if (runCode[3]) {
  
  ### Design ###
  
  # Thresholds           - Speed vs. accuracy parameters
  # Drift rates          - Correct versus incorrect parameters
  # Residual latencies   - Left versus right encoding
  # Coefficient of drift - Fixed to 1
  
  # Generating coefficients
  # k_S k_A (1:2)
  # x_C x_I (2:3)
  # tau_L tau_R (5:6)
  # sigma (Fixed to 1)
  gp = c( 0.8, 1.2, 2.5, 1.2, .27, .3 )
  
  # Define conditions for simulation
  n = 400
  SvA = rep( c(1,0), each = n/2 ) # Speed vs. accuracy conditions
  LvR = rep( c(1,0), each = n/2 ) # Position of correct answer
  LvR = sample( LvR ) # Randomly shuffle position
  int = rep( 1, n ) # Intercept
  
  # Assume a racer for the choice on the left and a racer for the 
  # choice on the right. This then requires matching the drift 
  # rate for the correct choice to the appropriate position.
  # Let 0 indicate picking left, and 1 indicate picking right.
  # We can then specify a design matrix that enables this switching.
  # Hence, the design matrix should have v = 12 variables.
  
  # Design matrix ( v by n )
  X = rbind( SvA, (1-SvA), (1-LvR), LvR, # kappa/xi (1)
             int, # sigma (1)
             int, # tau (1)
             SvA, (1-SvA), LvR, (1-LvR), # kappa/xi (0)
             int, # sigma (0)
             int  # tau (0)
  )
  
  # Next, we'll create an index for the parameter matrix 
  # rows and columns to use in the linear algebra function.
  # The columns of the parameter matrix are matched to the 
  # variables in the design matrix, while the rows are matched 
  # to the 8 parameters for the wald race model. However, rows 
  # 5 and 11 in the design matrix are for the coefficient of 
  # drift parameters, which are fixed. Fixed values are included 
  # last in the index.
  Clm = c( 1:2,  # Rows in design matrix for kappa(1)
           3:4,  # Rows in design matrix for xi(1)
           6,    # Rows in design matrix for tau(1)
           7:8,  # Rows in design matrix for kappa(0)
           9:10, # Rows in design matrix for xi(0)
           12    # Rows in design matrix for tau(0)
  )
  Rws = c( rep(1,2),  # kappa (1)
           rep(2, 2), # xi (1)
           4,         # tau (1)
           rep(5,2),  # kappa (0)
           rep(6, 2), # xi (0)
           8          # tau (0)
  )
  index = cbind( Rws, Clm )
  
  # Next, we specify the coefficients that correspond with 
  # each non-zero value in the parameter matrix. These indices 
  # are matched to the rows in the index variable.
  parSel = c( 1, 2, 3, 4, 5, 1, 2, 3, 4, 6 )
  
  # Coefficient of drift for both racers is fixed to 1
  fixed = c(1,1)
  index = rbind( index,
                 c(3,5), # sigma (1)
                 c(7,11) # sigma (0)
  )
  rm( Clm, Rws )
  
  # Finally, we include the dimensions (number of parameter P,
  # and the number of coefficients) at the end of the index
  index = rbind( index, c(8,length(parSel)) )
  
  # Generate parameter matrix for simulation purposes
  pm = param_est( X, gp, fixed, index, parSel )
  
  # Define function to transform/untransform parameters
  tran_par = function( prm, min_rt, reverse = F ) {
    
    if (reverse) {
      prm[c(1,2,3,4)] = log( prm[c(1,2,3,4)] )
      prm[ c(5,6) ] = logit( prm[ c(5,6) ]/min_rt )
    } else {
      prm[c(1,2,3,4)] = exp( prm[c(1,2,3,4)] )
      prm[ c(5,6) ] = logistic( prm[c(5,6)] )*min_rt
    }
    
    return( prm )
  }
  
  st_f = function() {
    
    v = runif( 6, c(.5,.5,.5,.5,.1,.1),
               c(2,2,4,4,.25,.25) )
    return( tran_par( v, min(sim[,1]), reverse = T ) )
  }
  
  # Define function for maximum likelihood estimation
  wr_mle = function( par, dat, priors ) {
    
    # Transform parameters
    pv = tran_par( par, min( dat[,1] ) )
    
    # Generate parameter matrix for simulation purposes
    pm = param_est( X, pv, fixed, index, parSel )
    
    # Trim out coefficients of drift
    pm = pm[ -c(3,7),]
    
    sll = sum(
      dwaldrace( dat[,1], dat[,2], pm[1,], pm[2,], pm[3,],
                 pm[4,], pm[5,], pm[6,], ln = 1 ) )
    if (is.na(sll)) sll = -Inf
    
    return( sll )
  }
  
  # Define number of simulation runs to conduct
  nRep = 20
  
  allPar = matrix( NA, nRep, length( gp ) )
  CI = array( NA, dim = c( nRep, 2, length(gp) ) )
  
  # Create a progress bar
  pb = txtProgressBar( min = 0, max = nRep, style = 3 )
  
  # Loop through repetitions
  for (nr in 1:nRep) {
    
    # Simulate data
    sim = rwaldrace( ncol( pm ), pm[1,], pm[2,], pm[4,],
                     pm[5,], pm[6,], pm[8,], 1, 1, rl = 1 )
    frt = min( sim[,1] )
    
    res = MLE( sim, wr_mle, st_f, maxit = 10000 )
    
    est = tran_par( res$param, frt )
    
    allPar[nr,] = est;
    CI[nr,1,] = tran_par( res$CI[1,], frt )
    CI[nr,2,] = tran_par( res$CI[2,], frt )
    
    # Update the progress bar
    setTxtProgressBar(pb,nr)
  }
  
  if (!savePlot) x11(width=12)
  layout( matrix(1:6,2,3,byrow=T) )
  
  lbls = c( expression(kappa[S]), expression(kappa[A]),
            expression(xi[C]), expression(xi[W]),
            expression(tau[L]), expression(tau[R]) )
  for (i in 1:6) {
    par( mar = c( 4, 5, 3, 1 ) )
    scl = .2; if (i > 4) scl = .05
    yl = lowerUpper(scl,as.vector( CI[,,i]) )
    plot( c(1,nRep), yl,
          type = 'n', xaxt = 'n', ylab = ' ', xlab = ' ',
          cex.axis = 1.5 )
    title( lbls[i], cex = 1.5 )
    abline( h = mean( allPar[,i] ), lwd = 2 )
    abline( h = gp[i], col = 'blue', lwd = 2, lty = 2 )
    segments( 1:nRep, CI[,1,i],
              1:nRep, CI[,2,i], lwd = 2 )
    points( 1:nRep, allPar[,i], pch = 19, cex = 1.5 )
    
    if (i==5) abline( h = gp[6], lty = 2, lwd = 2 )
    if (i==6) abline( h = gp[5], lty = 2, lwd = 2 )
    
    if (i==1) {
      par(xpd=TRUE)
      legend(1,yl[1] - diff(yl)*.1, c("True", "Estimate"), 
             lty = c(2,1), col = c('blue','black'),
             cex = 1.5, horiz = T, bty = 'n', lwd = 2 )
      par(xpd=FALSE)
    }
    
  }
  mtext( 'Parameter recovery', side = 1, outer = T,
         cex = 1.5, line = -2 )
  

}

###
### Parameter recovery with guessing processes
###
# Lookup - 04

if (runCode[4]) {
  
  ### Design ###
  
  # Thresholds           - Speed vs. accuracy parameters
  # Drift rates          - Correct versus incorrect parameters
  # Residual latencies   - Left versus right encoding
  # Coefficient of drift - Fixed to 1
  
  # Generating coefficients
  # k_S k_A (1:2)
  # x_C x_I (2:3)
  # tau_L tau_R (5:6)
  # sigma (Fixed to 1)
  gp = c( 0.8, 1.2, 2.5, 1.2, .27, .3 )
  
  # Define conditions for simulation
  n = 400
  SvA = rep( c(1,0), each = n/2 ) # Speed vs. accuracy conditions
  LvR = rep( c(1,0), each = n/2 ) # Position of correct answer
  LvR = sample( LvR ) # Randomly shuffle position
  int = rep( 1, n ) # Intercept
  
  # Assume a racer for the choice on the left and a racer for the 
  # choice on the right. This then requires matching the drift 
  # rate for the correct choice to the appropriate position.
  # Let 0 indicate picking left, and 1 indicate picking right.
  # We can then specify a design matrix that enables this switching.
  # Hence, the design matrix should have v = 12 variables.
  
  # Design matrix ( v by n )
  X = rbind( SvA, (1-SvA), (1-LvR), LvR, # kappa/xi (1)
             int, # sigma (1)
             int, # tau (1)
             SvA, (1-SvA), LvR, (1-LvR), # kappa/xi (0)
             int, # sigma (0)
             int  # tau (0)
  )
  
  # Next, we'll create an index for the parameter matrix 
  # rows and columns to use in the linear algebra function.
  # The columns of the parameter matrix are matched to the 
  # variables in the design matrix, while the rows are matched 
  # to the 8 parameters for the wald race model. However, rows 
  # 5 and 11 in the design matrix are for the coefficient of 
  # drift parameters, which are fixed. Fixed values are included 
  # last in the index.
  Clm = c( 1:2,  # Rows in design matrix for kappa(1)
           3:4,  # Rows in design matrix for xi(1)
           6,    # Rows in design matrix for tau(1)
           7:8,  # Rows in design matrix for kappa(0)
           9:10, # Rows in design matrix for xi(0)
           12    # Rows in design matrix for tau(0)
  )
  Rws = c( rep(1,2),  # kappa (1)
           rep(2, 2), # xi (1)
           4,         # tau (1)
           rep(5,2),  # kappa (0)
           rep(6, 2), # xi (0)
           8          # tau (0)
  )
  index = cbind( Rws, Clm )
  
  # Next, we specify the coefficients that correspond with 
  # each non-zero value in the parameter matrix. These indices 
  # are matched to the rows in the index variable.
  parSel = c( 1, 2, 3, 4, 5, 1, 2, 3, 4, 6 )
  
  # Coefficient of drift for both racers is fixed to 1
  fixed = c(1,1)
  index = rbind( index,
                 c(3,5), # sigma (1)
                 c(7,11) # sigma (0)
  )
  rm( Clm, Rws )
  
  # Finally, we include the dimensions (number of parameter P,
  # and the number of coefficients) at the end of the index
  index = rbind( index, c(8,length(parSel)) )
  
  # Generate parameter matrix for simulation purposes
  pm = param_est( X, gp, fixed, index, parSel )
  
  # Simulate data
  sim = rwaldrace( ncol( pm ), pm[1,], pm[2,], pm[4,],
                   pm[5,], pm[6,], pm[8,], 1, 1, rl = 1 )
  
  # Simulate a fast and slow guessing process
  guess1 = rexp( n, 6.5 ) + runif( n, -.05, .05 )
  summary( guess1 )
  sel = guess1 <= 0
  guess1[sel] = guess1[sel] + runif( sum(sel), .1, .4 )
  guess1 = cbind( guess1, rbinom( n, 1, .5 ) )
  
  guess2 = rlnorm( n, 1, .5 )
  sel = guess2 > 5
  guess2[sel] = 4.995
  guess2 = cbind( guess2, rbinom( n, 1, .5 ) )
  
  # Plot of guessing processes
  if (!savePlot) x11(width=12)
  
  layout( cbind(1,2) )
  
  # Fast responses
  blankRTplot(ver='PDF',pDim=c(0,4),
              bty = 'l', cex.axis = 1.5, cex.lab = 1.5 )
  
  pdf_curve( guess1[,1], guess1[,2], lwd = 2 )
  pdf_curve( guess1[,1], guess1[,2], sel = 0, lty = 2, lwd = 2 )
  
  legend( 'topright', c('Standard cut-offs'), cex = 1.5, 
          lty = 1, col= 'grey', lwd = 2, bty = 'n' )
  
  abline( v = c(.2,4), col = 'grey', lwd = 2 )
  
  # Slow responses
  blankRTplot(ver='PDF',pDim=c(0,4),tDim=c(0,5), bty = 'l',
              cex.axis = 1.5, cex.lab = 1.5 )
  
  pdf_curve( guess2[,1], guess2[,2], lwd = 2 )
  pdf_curve( guess2[,1], guess2[,2], sel = 0, lty = 2, lwd = 2 )
  
  abline( v = c(.2,4), col = 'grey', lwd = 2 )
  
  legend( 'top', c('Left','Right'), cex = 1.5, 
          lty = c(2,1), lwd = 2, bty = 'n' )
  
  mtext( 'Distribution of fast/slow guessing processes',
         side = 3, outer = T, cex = 1.5, line = -2 )
  
  # Set probability of attention failure
  alpha = runif(1,.9,.98) # Randomize
  Z = matrix( 0, n, 3 )
  Z[,1] = as.numeric( runif(n,0,1) < alpha )
  Z[ Z[,1] == 0,2] = rbinom( sum( Z[,1] == 0 ), 1, .5 )
  Z[ Z[,1] == 0, 3] = 1 - Z[ Z[,1] == 0,2]
  
  sim[,1] = sim[,1]*Z[,1] + guess1[,1]*Z[,2] + guess2[,1]*Z[,3]
  sim[,2] = sim[,2]*Z[,1] + guess1[,2]*Z[,2] + guess2[,2]*Z[,3]
  
  # Determine accuracy
  sim = cbind( sim, acc = sim[,2] == (1-LvR) )
  
  # Convert to data frame
  tmp = covCreate( cbind( SvA, LvR ) )
  sim = as.data.frame( cbind( sim, SvA = SvA, LvR = LvR, 
                              Cnd = tmp ) )
  
  # Pre-process data
  old_sim = sim
  trm = rtclean(sim[,1],nRep=20)
  sim = sim[ !trm$exclude_1, ]
  X = X[ , !trm$exclude_1 ]
  SvA = SvA[ !trm$exclude_1 ]
  LvR = LvR[ !trm$exclude_1 ]
  
  ### Parameter estimation using Stan ###
  
  # Define set of priors
  Priors = cbind(
    c( rep(1,2), rep(2.0,2), rep(4, 2 ) ),
    c( rep(.5,2), rep(.5,2), rep(2, 2 ) )
  )
  
  # Input for Stan
  stanDat = list(
    N = ncol(X),
    V = nrow(X),
    K = length(parSel),
    U = length(fixed), 
    C = c(2,2,2), 
    X = X,
    fixed = fixed,
    index = index,
    parSel = parSel,
    Y = sim[,1:2],
    min_RT = rep( min(sim[,1]), 2 ),
    Priors = Priors
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
  
  # Plot of parameter recovery
  if (!savePlot) x11(width=12)
  layout( cbind(1,2,3) )
  
  # Threshold recovery
  plot( c(.8,2.2), lowerUpper(.2,as.vector(post$kappa) ), 
        type = 'n', xaxt = 'n',
        ylab = 'Threshold values', xlab = 'Condition',
        bty = 'l', cex.lab = 1.5 )
  tmp = apply( post$kappa, 2, quantile, prob = c(.025,.5,.975) )
  segments( 1:2, tmp[1,], 1:2, tmp[3,] )
  points( 1:2, tmp[2,], pch = 15 )
  points( 1:2, gp[1:2], pch = 15, col='blue' )
  axis(1,1:2,c('Speed','Accuracy'),tick=F,cex.axis=1.5)
  
  legend('topleft',c('Generating','Estimated'),
         fill = c('blue','black'),bty='n',cex=1.5)
  
  # Drift rate recovery
  plot( c(.8,2.2), lowerUpper(.2,as.vector(post$xi) ), 
        type = 'n', xaxt = 'n',
        ylab = 'Drift rate values', xlab = ' ',
        bty = 'l', cex.lab = 1.5 )
  tmp = apply( post$xi, 2, quantile, prob = c(.025,.5,.975) )
  segments( 1:2, tmp[1,], 1:2, tmp[3,] )
  points( 1:2, tmp[2,], pch = 15 )
  points( 1:2, gp[3:4], pch = 15, col='blue' )
  axis(1,1:2,c('Correct','Incorrect'),tick=F,cex.axis=1.5)
  
  # Residual latency
  plot( c(.8,2.2), lowerUpper(.2,c(as.vector(post$tau),gp[5:6])), 
        type = 'n', xaxt = 'n',
        ylab = 'Residual latency', xlab = ' ',
        bty = 'l',cex.lab=1.5 )
  tmp = apply( post$tau, 2, quantile, prob = c(.025,.5,.975) )
  segments( 1:2, tmp[1,], 1:2, tmp[3,] )
  points( 1:2, tmp[2,], pch = 15 )
  points( 1:2, gp[5:6], pch = 15, col='blue' )
  axis(1,1:2,c('Left','Right'),tick=F,cex.axis=1.5)
  
  mtext( 'Parameter recovery', side = 3, cex = 1.5, 
         outer = T, line = -2 )
  
}

###
### Recovery with guessing over multiple repetitions (MLE)
###
# Lookup - 05

if ( runCode[5] ) {
  
  ### Design ###
  
  # Thresholds           - Speed vs. accuracy parameters
  # Drift rates          - Correct versus incorrect parameters
  # Residual latencies   - Left versus right encoding
  # Coefficient of drift - Fixed to 1
  
  # Generating coefficients
  # k_S k_A (1:2)
  # x_C x_I (2:3)
  # tau_L tau_R (5:6)
  # sigma (Fixed to 1)
  gp = c( 0.8, 1.2, 2.5, 1.2, .27, .3 )
  
  # Define conditions for simulation
  n = 400
  SvA = rep( c(1,0), each = n/2 ) # Speed vs. accuracy conditions
  LvR = rep( c(1,0), each = n/2 ) # Position of correct answer
  LvR = sample( LvR ) # Randomly shuffle position
  int = rep( 1, n ) # Intercept
  
  # Assume a racer for the choice on the left and a racer for the 
  # choice on the right. This then requires matching the drift 
  # rate for the correct choice to the appropriate position.
  # Let 0 indicate picking left, and 1 indicate picking right.
  # We can then specify a design matrix that enables this switching.
  # Hence, the design matrix should have v = 12 variables.
  
  # Design matrix ( v by n )
  X = rbind( SvA, (1-SvA), (1-LvR), LvR, # kappa/xi (1)
             int, # sigma (1)
             int, # tau (1)
             SvA, (1-SvA), LvR, (1-LvR), # kappa/xi (0)
             int, # sigma (0)
             int  # tau (0)
  )
  
  # Next, we'll create an index for the parameter matrix 
  # rows and columns to use in the linear algebra function.
  # The columns of the parameter matrix are matched to the 
  # variables in the design matrix, while the rows are matched 
  # to the 8 parameters for the wald race model. However, rows 
  # 5 and 11 in the design matrix are for the coefficient of 
  # drift parameters, which are fixed. Fixed values are included 
  # last in the index.
  Clm = c( 1:2,  # Rows in design matrix for kappa(1)
           3:4,  # Rows in design matrix for xi(1)
           6,    # Rows in design matrix for tau(1)
           7:8,  # Rows in design matrix for kappa(0)
           9:10, # Rows in design matrix for xi(0)
           12    # Rows in design matrix for tau(0)
  )
  Rws = c( rep(1,2),  # kappa (1)
           rep(2, 2), # xi (1)
           4,         # tau (1)
           rep(5,2),  # kappa (0)
           rep(6, 2), # xi (0)
           8          # tau (0)
  )
  index = cbind( Rws, Clm )
  
  # Next, we specify the coefficients that correspond with 
  # each non-zero value in the parameter matrix. These indices 
  # are matched to the rows in the index variable.
  parSel = c( 1, 2, 3, 4, 5, 1, 2, 3, 4, 6 )
  
  # Coefficient of drift for both racers is fixed to 1
  fixed = c(1,1)
  index = rbind( index,
                 c(3,5), # sigma (1)
                 c(7,11) # sigma (0)
  )
  rm( Clm, Rws )
  
  # Finally, we include the dimensions (number of parameter P,
  # and the number of coefficients) at the end of the index
  index = rbind( index, c(8,length(parSel)) )
  
  # Generate parameter matrix for simulation purposes
  pm = param_est( X, gp, fixed, index, parSel )
  
  # Maintain stable versions
  oX = X; oLvR = LvR; oSvA = SvA
  
  # Define number of repetitions
  nRep = 20
  
  # Define function to transform/untransform parameters
  tran_par = function( prm, min_rt, reverse = F ) {
    
    if (reverse) {
      prm[c(1,2,3,4)] = log( prm[c(1,2,3,4)] )
      prm[ c(5,6) ] = logit( prm[ c(5,6) ]/min_rt )
    } else {
      prm[c(1,2,3,4)] = exp( prm[c(1,2,3,4)] )
      prm[ c(5,6) ] = logistic( prm[c(5,6)] )*min_rt
    }
    
    return( prm )
  }
  
  st_f = function() {
    
    v = runif( 6, c(.5,.5,.5,.5,.1,.1),
               c(2,2,4,4,.25,.25) )
    return( tran_par( v, min(sim[,1]), reverse = T ) )
  }
  
  # Define function for maximum likelihood estimation
  wr_mle = function( par, dat, priors ) {
    
    # Transform parameters
    pv = tran_par( par, min( dat[,1] ) )
    
    # Generate parameter matrix for simulation purposes
    pm = param_est( X, pv, fixed, index, parSel )
    
    # Trim out coefficients of drift
    pm = pm[ -c(3,7),]
    
    sll = sum(
      dwaldrace( dat[,1], dat[,2], pm[1,], pm[2,], pm[3,],
                 pm[4,], pm[5,], pm[6,], ln = 1 ) )
    if (is.na(sll)) sll = -Inf
    
    return( sll )
  }
  
  # Define number of simulation runs to conduct
  nRep = 20
  
  allPar = matrix( NA, nRep, length( gp ) )
  CI = array( NA, dim = c( nRep, 2, length(gp) ) )
  
  # Create a progress bar
  pb = txtProgressBar( min = 0, max = nRep, style = 3 )
  
  # Loop through repetitions
  for (nr in 1:nRep) {
    
    # Simulate data
    sim = rwaldrace( ncol( pm ), pm[1,], pm[2,], pm[4,],
                     pm[5,], pm[6,], pm[8,], 1, 1, rl = 1 )
    
    # Simulate a fast and slow guessing process
    guess1 = rexp( n, 6.5 ) + runif( n, -.05, .05 )
    summary( guess1 )
    sel = guess1 <= 0
    guess1[sel] = guess1[sel] + runif( sum(sel), .1, .4 )
    guess1 = cbind( guess1, rbinom( n, 1, .5 ) )
    
    guess2 = rlnorm( n, 1, .5 )
    sel = guess2 > 5
    guess2[sel] = 4.995
    guess2 = cbind( guess2, rbinom( n, 1, .5 ) )
    
    # Set probability of attention failure
    alpha = runif(1,.9,.98) # Randomize
    Z = matrix( 0, n, 3 )
    Z[,1] = as.numeric( runif(n,0,1) < alpha )
    Z[ Z[,1] == 0,2] = rbinom( sum( Z[,1] == 0 ), 1, .5 )
    Z[ Z[,1] == 0, 3] = 1 - Z[ Z[,1] == 0,2]
    
    sim[,1] = sim[,1]*Z[,1] + guess1[,1]*Z[,2] + guess2[,1]*Z[,3]
    sim[,2] = sim[,2]*Z[,1] + guess1[,2]*Z[,2] + guess2[,2]*Z[,3]
    
    # Pre-process data
    old_sim = sim
    trm = rtclean(sim[,1],nRep=20)
    sim = sim[ !trm$exclude_1, ]
    X = oX[ , !trm$exclude_1 ]
    SvA = oSvA[ !trm$exclude_1 ]
    LvR = oLvR[ !trm$exclude_1 ]
    
    res = MLE( sim, wr_mle, st_f, maxit = 10000 )
    frt = min( sim[,1] )
    
    est = tran_par(res$param,frt)
    
    allPar[nr,] = est;
    CI[nr,1,] = tran_par( res$CI[1,], frt )
    CI[nr,2,] = tran_par( res$CI[2,], frt )
    
    # Update the progress bar
    setTxtProgressBar(pb,nr)
  }
  
  if (!savePlot) x11(width=12)
  layout( matrix(1:6,2,3,byrow=T) )
  
  lbls = c( expression(kappa[S]), expression(kappa[A]),
            expression(xi[C]), expression(xi[W]),
            expression(tau[L]), expression(tau[R]) )
  for (i in 1:6) {
    par( mar = c( 4, 5, 3, 1 ) )
    scl = .2; if (i > 4) scl = .05
    yl = lowerUpper(scl,as.vector( CI[,,i]) )
    plot( c(1,nRep), yl,
          type = 'n', xaxt = 'n', ylab = ' ', xlab = ' ',
          cex.axis = 1.5 )
    title( lbls[i], cex = 1.5 )
    abline( h = mean( allPar[,i] ), lwd = 2 )
    abline( h = gp[i], col = 'blue', lwd = 2, lty = 2 )
    segments( 1:nRep, CI[,1,i],
              1:nRep, CI[,2,i], lwd = 2 )
    points( 1:nRep, allPar[,i], pch = 19, cex = 1.5 )
    
    if (i==5) abline( h = gp[6], lty = 2, lwd = 2 )
    if (i==6) abline( h = gp[5], lty = 2, lwd = 2 )
    
    if (i==1) {
      par(xpd=TRUE)
      legend(1,yl[1] - diff(yl)*.1, c("True", "Estimate"), 
             lty = c(2,1), col = c('blue','black'),
             cex = 1.5, horiz = T, bty = 'n', lwd = 2 )
      par(xpd=FALSE)
    }
    
  }
  mtext( 'Parameter recovery', side = 1, outer = T,
         cex = 1.5, line = -2 )
  
}

if (savePlot) {
  setwd( orig_dir )
  setwd( 'Plots' )
  dev.off()
}

setwd( orig_dir )