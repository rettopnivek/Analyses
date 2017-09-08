#--------------------#
# Manuscript figures #
# Kevin Potter       #
# Updated 08/21/2017 #
#--------------------#

# Initialize script
source('F3_starting_script.R')

# Load in package for nROUSE model
library( nROUSE )

# Indicate which code segments to run
runCode = c( T, T, T, T, T, T, T, T, T, F )

# Change directory to where figures will be stored
setwd( 'Manuscript_figures' )

# Indicate whether to save figures as PDFs
savePlot = T

# Index
# Lookup - 01:  Figure 1
# Lookup - 02:  Figure 3
# Lookup - 03:  Figure 4
# Lookup - 04:  Figure 6
# Lookup - 05:  Figure 7
# Lookup - 06:  Figure 8
# Lookup - 07:  Figure 9
# Lookup - 08:  Figure 13
# Lookup - 09:  Figure 14
# Lookup - 10:  Figures B.1, B.2, and B.3

###
### Figure 1
###
# Lookup - 01

if ( runCode[1] ) {
  
  # Data from Huber et al. (2008)
  ExDat = cbind( 
    Accuracy = c( .921, .875, .776, .594, .550, 
                  .616, .393, .488, .809, .832 ),
    Duration = rep( c(.017,.05,.15,.4,2), 2 ), 
    Prime = rep( c(1,0), each = 5 ) )
  ExDat = as.data.frame( ExDat )
  
  # Plotting characteristics
  lnSz = 2
  txtSz = 1.5
  ptSz = 1.5
  
  if ( savePlot ) 
    pdf('Example_data_Huber_2008.pdf',width=6,height=6)
  else x11()
  
  # Create a blank plot
  xl = c( -4.6, .7 )
  yl = c( 0, 1 )
  blankPlot( xl, yl )
  
  # Axes
  customAxes( xl, yl, lnSz = lnSz )
  axis( 1, 
        log( c( .01, .1, 1 ) ), 
        c('10','100','1000'), 
        cex.axis = txtSz,
        tick = F, 
        line = -1 )
  axis( 2, 
        c(0,.2,.4,.6,.8,1), 
        tick = F, line = -1, 
        cex.axis = txtSz )
  segments( xl[1], .5, xl[2], .5, 
            col = 'grey', lwd = lnSz )
  mtext('Prime duration (ms - log scale)', 
        side = 1, cex = txtSz, line = 2.25 )
  mtext('Proportion correct', 
        side = 2, cex = txtSz, line = 2.75 )
  legend( 'bottomleft', c('Target primed','Foil primed'),
          pch = c(21,22), pt.bg = c('black','white'), bty = 'n',
          cex = txtSz, inset = c(.1,.1) )
  
  pts = c( 21, 22 ); clr = c( 'black', 'white' )
  for ( i in 1:0 ) {
    sel = ExDat$Prime == i
    lines( log( ExDat$Duration[sel] ), 
           ExDat$Accuracy[sel], lwd = lnSz )
    points( log( ExDat$Duration[sel] ), 
            ExDat$Accuracy[sel], pch = pts[2-i], cex = ptSz,
            bg = clr[2-i] )
  }
  if ( savePlot ) dev.off()
  
}

###
### Figure 3
###
# Lookup - 02

if ( runCode[2] ) {
  
  # Default parameters for nROUSE
  param = c( Fe = .25, N = .0302, L = .15, 
             D = .324, R = .022, I = .9844, 
             Th = .15, Ta = 1.0, SV = .0294, 
             SO = .0609, SS = .015 )
  
  # Initialize variables for simulation
  presentations = c( 400, # Prime duration
                     50, # Target duration
                     450, # Mask duration
                     500 # On-screen choices duration
                     )
  # Prime types
  TargetPrimed = c(2,0)
  FoilPrimed = c(0,2)
  
  # Initialize variables for plotting
  shft = .2 * (3:0)
  shft2 = c( 400, 50, 400, 50 )
  input = TargetPrimed
  loop_setup = function(cnd) {
    if ( cnd == 1 ) {
      input <<- TargetPrimed
      presentations[1] <<- 50
    }
    if ( cnd == 2 ) {
      input <<- TargetPrimed
      presentations[1] <<- 400
    }
    if ( cnd == 3 ) {
      input <<- FoilPrimed
      presentations[1] <<- 50
    }
    if ( cnd == 4 ) {
      input <<- FoilPrimed
      presentations[1] <<- 400
    }
  }
  
  # Create blank plot
  if ( savePlot ) 
    pdf('nROUSE_activation.pdf',width=6.8,height=6.8)
  else x11()
  xl = c( 0, sum( presentations ) + 50 )
  yl = c( 0, .2*4 )
  par( mar = c( 4, 4, 3, 1.5 ) )
  blankPlot( xl, yl )
  
  # Fill in target flash duration
  for ( cnd in 1:4 ) {
    
    loop_setup(cnd)
    
    # Indicate prime durations
    x = rep( c( presentations[1] + shft2[cnd], 
                cumsum( presentations )[2] + shft2[cnd] ),
             each = 2 )
    y = c( 0, .2, .2, 0 ) + shft[cnd]
    polygon( x, y, border = NA, col = 'grey80' )
    
  }
  
  # Add in final plotting elements
  for ( cnd in 1:4 ) {
    
    loop_setup(cnd)
    
    total = sum( presentations )
    
    # Simulate activation
    sim = simulate_nROUSE( presentations, input, param )
    sel = (sum( presentations[1:3] ) + 1 ):sum( presentations )
    act = sim$Activation[ sel, ]
    
    lines( 1:total + shft2[cnd], 
           sim$Activation[,'Target'] + shft[cnd], lwd = 2 )
    lines( 1:total + shft2[cnd], 
           sim$Activation[,'Foil'] + shft[cnd], 
           lwd = 2, lty = 2 )
    
    # Axis on left
    segments( shft2[cnd], 0 + shft[cnd],
              shft2[cnd], .2 + shft[cnd],
              lwd = 2 )
    # Top and bottom axes
    segments( shft2[cnd], shft[cnd], 
              xl[2], shft[cnd], lwd = 2 )
    segments( shft2[cnd], .2 + shft[cnd], 
              xl[2], .2 + shft[cnd], lwd = 2 )
    # Choice segment
    segments( cumsum( presentations )[3] + shft2[cnd],
              0 + shft[cnd],
              cumsum( presentations )[3] + shft2[cnd],
              .2 + shft[cnd], lwd = 2 )
    
    # Peak activation
    sel = 1:total > cumsum( presentations )[3] + shft2[cnd] - 1
    sel2 = ( nrow( sim$Activation ) - 
               presentations[4] ):nrow( sim$Activation )
    x = 950 + sim$Latencies[1:2]
    y = sim$Activation[ sel2, ][ sim$Latencies, ]
    y = c( y[1,1], y[2,2] ) + shft[cnd]
    segments( x, rep( 0 + shft[cnd], 2 ),
              x, y, lwd = 2, lty = 1:2 )
              
    points( x, y, pch = c(21,22), bg = c('black','white'),
            cex = 1.2 )
    
    axis( 4, shft[cnd] + .1, 
          paste( round( sim$Latencies[3] * 100 ), '%',
          sep = '' ), cex.axis = 1.4, tick = F, line = -1.5 )
    
  }
  customAxes( xl, yl, pos = 4 )
  # Axis labels
  axis( 3, 475, 'Target flash', cex.axis = 1.4, 
        tick = F, line = -1.5 )
  axis( 1, c(50,500,950), c('Prime','Mask','Choice'), cex.axis = 1.4, 
        tick = F, line = -1.5, hadj = 0 )
  axis( 2, yl[2]/2, 'Lexical-semantic activation', 
        cex.axis = 1.4, tick = F, line = 1.75 )
  axis( 2, c( yl[2]/4, 3*yl[2]/4 ),
        c( 'Foil primed', 'Target primed' ),
        tick = F, cex.axis = 1.4, line = -1 )
  
  axis( 2, shft[c(2,4)] + .1, c( '400', '400' ),
        tick = F, line = -2.8, cex.axis = 1.4 )
  axis( 2, shft[c(1,3)] + .1, c( '50', '50' ),
        tick = F, line = -9.2, cex.axis = 1.4 )
  
  mtext( 'Time (ms)', side = 1, line = 1, cex = 1.4 )
  
  par( xpd = T )
  legend( 'bottom', c( 'Target representation',
                       'Foil representation' ),
          lty = 1:2, lwd = 2, bty = 'n', horiz = T,
          cex = 1.4, inset = -.15 )
  par( xpd = F )
  
  if (savePlot) dev.off()
}

###
### Figure 4
###
# Lookup - 03

if ( runCode[3] ) {
  
  if ( savePlot ) pdf('Seq_samp_ex.pdf',width=6,height=6) else 
    x11(width=7,height=7);
  
  # Generate two diffusion processes
  
  sim.diff.proc = function(N=80,mu=c(.5,.45),k=c(25,25)) {
    
    evA = cumsum(rnorm(N,mu[1],1.0))
    evB = cumsum(rnorm(N,mu[2],1.0))
    
    evA[evA<=0]=0 # Restrict values to be positive
    evB[evB<=0]=0
    
    kA = k[1] # Threshold values
    kB = k[2]
    
    if (max(evA) > kA) {
      sel = min( which(evA>=kA ) )
      if (sel < N) evA[(sel+1):N] = NA
    }
    if (max(evB) > kB) {
      sel = min( which(evB>=kB ) )
      if (sel < N) evB[(sel+1):N] = NA
    }
    
    # Finishing times
    FT = c( which( na.omit(evA)==max(na.omit(evA)) ),
            which( na.omit(evB)==max(na.omit(evB)) ) )
    return( list( evA = evA, evB = evB, FT = FT,
                  k = k, N = N) )
  }
  
  # Simulate from the process
  "
  tmp = sim.diff.proc(mu=c(.7,.5))
  while (tmp$FT[1] > tmp$FT[2]) {
  tmp = sim.diff.proc(mu=c(.7,.5))
  }
  v1 = tmp
  "
  
  # Save the best examples of evidence accumulation
  # save(v1,file='Example_diffusion_race_process.RData')
  setwd( orig_dir ); setwd( 'Data' )
  load( 'Example_diffusion_race_process.RData' )
  setwd( orig_dir ); setwd( 'Manuscript_figures' )
  
  par(mar=c(0,0,0,0),oma=c(4,4,3,4))
  
  txtsz = 1.3
  cax = 1.5
  
  # Framework for plot
  part = 1
  if (part==1) {
    xl = c(0,100)
    plot( xl, c(0, 40), type='n',xaxt='n',
          yaxt='n', bty='n',ylab = ' ')
    axis(1,seq(0,100,20),seq(0,1,.2),cex.axis=cax,
         lwd=2)
    axis(2,25,expression( kappa[1]==kappa[0] ),
         line=-1,tick=F,cex.axis=cax)
    segments(0,-2,0,40,lwd=2)
    mtext('Evidence for a choice',side=2,line=2,cex=1.5)
    mtext('RT (s)',side=1,cex=1.5,line=3)
    
    legend(0,40,
           c( expression( paste(
             '( ',kappa[1],', ',kappa[0],' ) = Thresholds',
             sep='') ),
             expression( paste(
               '( ',xi[1],', ',xi[0],' ) = Rates of evidence accumulation',
               sep='') ),
             expression( paste(
               tau,' = Non-decision time',sep='') ) ), 
           bty='n',cex=1.2)
  }
  
  # Plot the evidence accumulation and label
  part = 2
  if (part==2) {
    
    lines(21:(v1$FT[1]+20),
          c(v1$evA[1:(v1$FT[1]-1)],v1$k[1]),lwd=2)
    lines(21:(v1$FT[2]+20),
          c(v1$evB[1:(v1$FT[2]-1)],v1$k[2]),
          lty=2,lwd=2)
    abline(h=25,lwd=2)
    
    segments( c(1,1,9.5,19),c(0,2,2.5,0),
              c(1,19,9.5,19),c(2,2,2,2), lwd=2)
    text(9.5,2.5,expression(tau),pos=3,cex=1.4)
    
    points(v1$FT[1]+20,25,pch=19,cex=1.5)
    text(v1$FT[1]+20,25,'Observed RT',pos=3,cex=1.2)
    text(v1$FT[1]+20,28,'Choose choice 1',pos=3,cex=1.2)
    
  }
  
  part = 3
  if (part==3) {
    # Arrows denoting drift rate values
    text(37,17.5,expression(xi[1]),pos=2,cex=1.3)
    arrows(34,15,40,20,length=.1,lwd=2)
    
    text(57,11,expression(xi[0]),pos=4,cex=1.3)
    arrows(50,10,60,13.5,length=.1, lwd=2)
  }
  
  # Clean up workspace
  rm( xl, v1 )
  
  if ( savePlot ) dev.off()
  
}

###
### Figure 6
###
# Lookup - 04

if (runCode[4]) {
  
  ### Accuracy and choice ###
  
  ds = aggregate( d$Ac, list( d$PDL, d$PTL, d$TaL), mean )
  colnames( ds ) = c('Duration','Prime_type','Task','P')
  
  lnSz = 2
  ptSz = 1.5
  txtSz = 1.5
  lgnd = .26
  
  if ( savePlot ) pdf('Accuracy_and_choice.pdf',width=12,height=6) else x11(width=12)
  layout( cbind(1,2,3) )
  
  par(mar=c(3,5,1,1))
  xl = c(.8,2.2); yl = c(0,1); blankPlot( xl, yl )
  axis(1,1:2,c('50','400'),cex.axis = txtSz*1.25, 
       lwd = lnSz, tick = F, line = -1.5 )
  axis(2,seq(0,1,.25),cex.axis = txtSz*1.25, lwd = lnSz )
  abline(h=0,lwd=lnSz)
  mtext('Accuracy',side=2,line=3,cex=txtSz)
  
  pts = c(21,24,21,24)
  clr = rep(c('black','white'),2)
  ltp = c(1,1,2,2)
  for (i in 1:4) {
    sel = 1:2 + (i-1)*2
    lines( 1:2, ds$P[sel], lty = ltp[i], lwd = lnSz )
    points( 1:2, ds$P[sel], pch = pts[i], bg = clr[i], cex = ptSz )
  }
  
  legend(.8,lgnd,c( 'Forced-choice','Same-different',
                    'Foil primed','Target primed' ),
         pch = c(NA,NA,21,24), 
         pt.bg = c(NA,NA,'black','white'), 
         lty = c(1,2,NA,NA), 
         lwd = c(2,2,1,1), 
         bty = 'n',
         cex = txtSz*1.25 )
  
  
  ### Choice ###
  
  ds = aggregate( d$Ac, list( d$PDL, d$PTL, d$Co, d$TaL), mean )
  colnames( ds ) = c('Duration','Prime_type','Correct','Task','P')
  
  ### Forced choice ###
  
  par(mar=c(3,5,1,1))
  xl = c(.8,2.2); yl = c(0,1); blankPlot( xl, yl )
  axis(1,1:2,c('50','400'),cex.axis = txtSz*1.25, 
       lwd = lnSz, tick = F, line = -1.5 )
  axis(2,seq(0,1,.25),cex.axis = txtSz*1.25, lwd = lnSz )
  abline(h=0,lwd=lnSz)
  
  pts = c(21,24,21,24)
  clr = rep(c('black','white'),2)
  ltp = c(1,1,2,2)
  for (i in 1:4) {
    sel = 1:2 + (i-1)*2
    lines( 1:2, ds$P[sel], lty = ltp[i], lwd = lnSz )
    points( 1:2, ds$P[sel], pch = pts[i], bg = clr[i], cex = ptSz )
  }
  
  legend(.8,lgnd,c( 'Left','Right' ),
         lty = c(1,2), 
         lwd = c(2,2), 
         bty = 'n',
         cex = txtSz*1.25 )
  
  ### Same-different ###
  
  par(mar=c(3,5,1,1))
  xl = c(.8,2.2); yl = c(0,1); blankPlot( xl, yl )
  axis( 1,1:2,c('50','400'),cex.axis = txtSz*1.25, 
        lwd = lnSz, tick = F, line = -1.5 )
  axis(2,seq(0,1,.25),cex.axis = txtSz*1.25, lwd = lnSz )
  abline(h=0,lwd=lnSz)
  mtext('Accuracy',side=2,line=3,cex=txtSz)
  
  pts = c(21,24,21,24)
  clr = rep(c('black','white'),2)
  ltp = c(1,1,2,2)
  for (i in 1:4) {
    sel = 1:2 + (4+i-1)*2
    lines( 1:2, ds$P[sel], lty = ltp[i], lwd = lnSz )
    points( 1:2, ds$P[sel], pch = pts[i], bg = clr[i], cex = ptSz )
  }
  
  legend(.8,lgnd,c( 'Same','Different' ),
         lty = c(1,2), 
         lwd = c(2,2), 
         bty = 'n',
         cex = txtSz*1.25 )
  
  mtext('Prime duration (ms)',side=1,line=-1,cex=txtSz,
        outer = T )
  if (savePlot) dev.off()
  
}

###
### Figure 7
###
# Lookup - 05

if ( runCode[5] ) {
  
  # Median RTs
  tmp = aggregate( d$RT, list( d$PDL, d$PTL, d$Co, 
                               d$Ac,
                               d$TaL, d$S ), median )
  
  # Aggregate over subjects
  ds = aggregate( tmp$x, list( tmp[,1], tmp[,2], tmp[,3], 
                               tmp[,4], tmp[,5] ), mean )
  colnames( ds ) = c('Duration','Prime_type','Correct',
                     'Accuracy','Task','M')
  
  if ( savePlot ) pdf('Median_RT.pdf',width=12,height=6) else 
    x11(width=12)
  layout( cbind(1,2) )
  
  lnSz = 2
  ptSz = 1.5
  txtSz = 1.5
  lgnd = .26
  
  ### Forced choice ###
  
  par(mar=c(3,5,2,1))
  xl = c(.5,4.5); yl = c(.4,1.0); blankPlot( xl, yl )
  axis(1,1:4,rep(c('50','400'),2),cex.axis = txtSz*1, 
       lwd = lnSz, tick = F, line = -1.25 )
  axis(3,c(1.5,3.5),c('Errors','Correct'),tick=F,
       cex.axis=txtSz*1, line = -1 )
  axis(2,seq(yl[1],yl[2],.2),cex.axis = txtSz*1, lwd = lnSz )
  abline(h=yl[1],lwd=lnSz)
  segments(2.5,yl[1],2.5,yl[2],lwd=lnSz)
  mtext('Median RT (s)',side=2,line=3,cex=txtSz)
  
  pts = c(21,24,21,24)
  clr = rep(c('black','white'),2)
  ltp = c(1,1,2,2)
  for (i in 1:4) {
    sel = 1:2 + (i-1)*2
    lines( 1:2, ds$M[sel], lty = ltp[i], lwd = lnSz )
    points( 1:2, ds$M[sel], pch = pts[i], bg = clr[i], cex = ptSz )
    
    sel = 1:2 + (4+i-1)*2
    lines( 3:4, ds$M[sel], lty = ltp[i], lwd = lnSz )
    points( 3:4, ds$M[sel], pch = pts[i], bg = clr[i], cex = ptSz )
  }
  
  legend('topleft',c( 'Left','Right',
                      'Foil primed','Target primed' ),
         pch = c(NA,NA,21,24), 
         pt.bg = c(NA,NA,'black','white'), 
         lty = c(1,2,NA,NA), 
         lwd = c(2,2,1,1), 
         bty = 'n',
         cex = txtSz*.9 )
  
  ### Same-different ###
  
  par(mar=c(3,5,2,1))
  xl = c(.5,4.5); c(.4,1.0); blankPlot( xl, yl )
  axis(1,1:4,rep(c('50','400'),2),cex.axis = txtSz*1, 
       lwd = lnSz, tick = F, line = -1.25 )
  axis(3,c(1.5,3.5),c('Errors','Correct'),tick=F,
       cex.axis=txtSz*1, line = -1 )
  axis(2,seq(yl[1],yl[2],.2),cex.axis = txtSz*1, lwd = lnSz )
  abline(h=yl[1],lwd=lnSz)
  segments(2.5,yl[1],2.5,yl[2],lwd=lnSz)
  
  pts = c(21,24,21,24)
  clr = rep(c('black','white'),2)
  ltp = c(1,1,2,2)
  for (i in 1:4) {
    sel = 1:2 + (8+i-1)*2
    lines( 1:2, ds$M[sel], lty = ltp[i], lwd = lnSz )
    points( 1:2, ds$M[sel], pch = pts[i], bg = clr[i], cex = ptSz )
    
    sel = 1:2 + (12+i-1)*2
    lines( 3:4, ds$M[sel], lty = ltp[i], lwd = lnSz )
    points( 3:4, ds$M[sel], pch = pts[i], bg = clr[i], cex = ptSz )
  }
  
  legend('topleft',c( 'Same','Different' ),
         lty = c(1,2), 
         lwd = c(2,2), 
         bty = 'n',
         cex = txtSz*.9 )
  
  mtext('Prime duration (ms)',side=1,line=-1,cex=txtSz,
        outer = T )
  if ( savePlot ) dev.off()
  
}

###
### Figure 8
###
# Lookup - 06

if ( runCode[6] ) {
  
  # Extract relevant plotting variables
  
  # Target durations
  setwd( orig_dir )
  setwd( 'Data' )
  load( 'nROUSE_est_from_Dave.RData' )
  setwd( orig_dir )
  setwd( 'Manuscript_figures' )
  tarDur = nROUSE_prev_est$target_duration
  
  # Accuracy over prime type and duration
  cD = allDat[ allDat$TaskLab == 'Forced-choice', ]
  cD = aggregate( cD$Accuracy, list( cD$PrimeDurLab,
                                     cD$PrimeTypeLab,
                                     cD$Subject ), mean )
  cD = cbind(
    SFP = cD$x[ cD[,1] == .05 & cD[,2] == 'Foil-primed' ],
    LFP = cD$x[ cD[,1] == .4 & cD[,2] == 'Foil-primed' ],
    STP = cD$x[ cD[,1] == .05 & cD[,2] == 'Target-primed' ],
    LTP = cD$x[ cD[,1] == .4 & cD[,2] == 'Target-primed' ] )
  cD = as.data.frame( cD )
  cD$S = 1:N
  
  # Parameter estimates
  NM = nROUSE_res$NoiseMult
  I = nROUSE_res$Inhibit
  TA = nROUSE_res$TemporalAtten
  
  ### Create scatterplots ###
  
  txtSz = 1.5
  
  if ( !savePlot ) x11()
  if ( savePlot ) {
    pdf(file='nROUSE_data_correlations.pdf',width=6,height=6)
  }
  layout( rbind(  c( 3, 3, 4, 4 ),
                  c( 1, 5, 5, 2 ) ) )
  blankPlot()
  blankPlot()
  
  # Function to extract boundaries for plot
  match_xy = function(x,y) {
    
    sx = scale(x)
    sy = scale(y)
    
    lmnts = rbind( c( min(sx), max(sx) ), 
                   c( min(sy), max(sy) ) )
    lmnts = c( min(lmnts[,1]), max(lmnts[,2]) )
    lmnts = lowerUpper( .5, lmnts )
    
    x_lmnts = lmnts*sd(x) + mean(x)
    y_lmnts = lmnts*sd(y) + mean(y)
    
    return( list( x_lmnts, y_lmnts ) )
  }
  
  x = TA; y = tarDur
  tmp = match_xy( x, y )
  xl = tmp[[1]]; yl = tmp[[2]]
  
  par( mar = c( 5, 5, 1, 1 ) )
  plot( xl, yl,
        type = 'n', xaxt = 'n', yaxt = 'n', bty = 'n',
        xlab = ' ', ylab = ' ' )
  axis( 1, round( seq( xl[1], xl[2], length = 4 ), 1 ), 
        cex.axis = txtSz )
  mtext('Temporal attention',side=1,cex=txtSz*.9,line = 3 )
  axis( 2, round( seq( yl[1], yl[2], length = 4 ) ), 
        cex.axis = txtSz )
  mtext('Target duration (ms)',side=2,cex=txtSz*.9,line = 3)
  
  points( x, y, pch = 19, cex = txtSz )
  legend( 'topright', paste('R =', round( cor(x,y), 2 ) ),
          bty = 'n', cex = txtSz*1.2 )
  legend( 'bottomleft', '(A)',
          bty = 'n', cex = txtSz*1.2 )
  
  
  x = NM
  y = cD$STP - cD$SFP
  tmp = match_xy( x, y )
  xl = tmp[[1]]; yl = tmp[[2]]
  
  par( mar = c( 5, 5, 1, 1 ) )
  plot( xl, yl,
        type = 'n', xaxt = 'n', yaxt = 'n', bty = 'n',
        xlab = ' ', ylab = ' ' )
  axis( 1, c( .024, .028, .032 ), 
        cex.axis = txtSz )
  mtext('Noise multiplier',side=1,cex=txtSz*.9,line = 3 )
  axis( 2, round( seq( yl[1], yl[2], length = 4 ), 2 ), 
        cex.axis = txtSz )
  mtext('Short target - foil',side=2,cex=txtSz*.9,line = 3)
  
  points( x, y, pch = 19, cex = txtSz )
  legend( 'topright', paste('R =', round( cor(x,y), 2 ) ),
          bty = 'n', cex = txtSz*1.2 )
  legend( 'bottomleft', '(B)',
          bty = 'n', cex = txtSz*1.2 )
  
  x = I
  y = cD$LTP - cD$LFP
  tmp = match_xy( x, y )
  xl = tmp[[1]]; yl = tmp[[2]]
  
  par( mar = c( 5, 5, 1, 1 ) )
  plot( xl, yl,
        type = 'n', xaxt = 'n', yaxt = 'n', bty = 'n',
        xlab = ' ', ylab = ' ' )
  axis( 1, round( seq( xl[1], xl[2], length = 4 ), 2 ), 
        cex.axis = txtSz )
  mtext('Inhibition',side=1,cex=txtSz*.9,line = 3 )
  axis( 2, round( seq( yl[1], yl[2], length = 4 ), 2 ), 
        cex.axis = txtSz )
  mtext('Long target - foil',side=2,cex=txtSz*.9,line = 3)
  
  points( x, y, pch = 19, cex = txtSz )
  legend( 'topleft', paste('R =', round( cor(x,y), 2 ) ),
          bty = 'n', cex = txtSz*1.2 )
  legend( 'bottomright', '(C)',
          bty = 'n', cex = txtSz*1.2 )
  if (savePlot) dev.off()
  
}

###
### Figure 9
###
# Lookup - 07

if ( runCode[7] ) {
  
  id = d[ 1:(N*20), ]
  for ( s in 1:N ) {
    
    sta = 1 + 20 * ( s - 1 )
    sto = 20 + 20 * ( s - 1 )
    
    tmp = d[ d$S == s, ]
    id[ sta:sto, ] = tmp[ 1:20, ]
    
  }
  rm( sta, sto, tmp, s )
  
  ds = aggregate( id$Ac, list( id$PDL, id$PTL, id$TaL), mean )
  colnames( ds ) = c('Duration','Prime_type','Task','P')
  
  lnSz = 2
  ptSz = 1.5
  txtSz = 1.5
  lgnd = .51
  
  if ( savePlot ) pdf('Initial_accuracy.pdf',width=6,height=6) else x11()
  
  par(mar=c(3,5,1,1))
  xl = c(.8,2.2); yl = c(.25,1); blankPlot( xl, yl )
  axis(1,1:2,c('50','400'),cex.axis = txtSz*1.25, 
       lwd = lnSz, tick = F, line = -1.5 )
  axis(2,seq(0,1,.25),cex.axis = txtSz*1.25, lwd = lnSz )
  abline(h=yl[1],lwd=lnSz)
  mtext('Accuracy',side=2,line=3,cex=txtSz)
  
  pts = c(21,24,21,24)
  clr = rep(c('black','white'),2)
  ltp = c(1,1,2,2)
  for (i in 1:4) {
    sel = 1:2 + (i-1)*2
    lines( 1:2, ds$P[sel], lty = ltp[i], lwd = lnSz )
    points( 1:2, ds$P[sel], pch = pts[i], bg = clr[i], cex = ptSz )
  }
  
  legend(.8,lgnd,c( 'Foil primed','Target primed' ),
         pch = c(21,24), 
         pt.bg = c('black','white'),
         bty = 'n',
         cex = txtSz*1.25 )
  
  if ( savePlot ) dev.off()
  
}


###
### Figure 13
###
# Lookup - 08

if ( runCode[8] ) {
  
  quick_plot = function( lbl, pos ) {
    # Purpose: 
    # ... 
    # Arguments: 
    # ... 
    # Returns: 
    # ...
    
    vrb = paste( lbl, '_xi', sep = '' )
    
    dn = density( cor_val[,vrb] )
    trm = dn$x > 1 | dn$x < -1;
    dn$x = dn$x[ !trm ];
    dn$x = c( dn$x[1], dn$x, dn$x[length(dn$x)] )
    dn$y = dn$y[ !trm ];
    dn$y = c( 0, dn$y, 0 )
    pts = densityPoints( cor_val[,vrb] )
    scl = max( dn$y )
    polygon( .4 * -dn$y/scl + pos, dn$x, col = 'grey50', lwd = 2 )
    points( .4 * -pts$y/scl + pos, pts$x, pch = 19 )
    
    #vrb = paste( lbl, '_k', sep = '' )
    #dn = density( cor_val[,vrb] )
    #trm = dn$x > 1 | dn$x < -1;
    #dn$x = dn$x[ !trm ];
    #dn$x = c( dn$x[1], dn$x, dn$x[length(dn$x)] )
    #dn$y = dn$y[ !trm ];
    #dn$y = c( 0, dn$y, 0 )
    #pts = densityPoints( cor_val[,vrb] )
    #scl = max( dn$y )
    #polygon( .4 * dn$y/scl + pos, dn$x, col = 'grey80', lwd = 2 )
    #points( .4 * pts$y/scl + pos, pts$x, pch = 19 )
    
  }
  
  setwd( orig_dir )
  setwd( 'RT_and_choice_MLE' )
  setwd( 'Estimation_results' )
  
  # Load in saturated model
  load( 'DRM_SD_M6_Sat_results.RData' )
  setwd( orig_dir )
  setwd( 'Manuscript_figures' )
  
  # Extract drift rates
  xi = all_prm[,grep('xi',colnames(raw_prm))]
  # Extract thresholds
  kappa = all_prm[,grep('kappa',colnames(raw_prm))]
  # Extract latencies
  Lat = nROUSE_res[,1:8+3]
  
  # Create data-frame with data to plot
  dtbp = data.frame( 
    # Subject index
    S = rep( 1:N, each = 8 ), 
    # Correct response
    R = rep( rep( c('Target','Foil'), each = 4 ), N ), 
    # Prime type
    PT = rep( rep( rep( c('Foil','Target'), each = 2 ), 2 ), N ), 
    # Prime duration
    PD = rep( rep( c('50','400'), 4 ), N ), 
    # Label
    L = rep( paste( 'xi-', 1:8, sep = '' ), N ) )
  
  # Initialize variables for nROUSE latencies
  dtbp$lat = 0
  # Initialize variables for drift rates
  dtbp$xi = 0
  # Initialize variables for thresholds
  dtbp$kappa = 0
  
  # Loop over subjects
  for ( n in 1:N ) {
    
    # Extract drift rates
    dtbp$xi[ dtbp$S == n ] = xi[n,1:8]
    # Extract thresholds
    dtbp$kappa[ dtbp$S == n ] = kappa[n,1:8]
    
    sel = d$S == n & d$TaL == 'Same-different'
    tmp = by( d$TL[sel], list( d$PDL[sel], d$PTL[sel], d$CoL[sel] ), 
              unique )
    tmp = as.numeric( unlist( tmp ) )
    dtbp$lat[ dtbp$S == n & dtbp$R == 'Target' ] = tmp[1:4]
    tmp = by( d$FL[sel], list( d$PDL[sel], d$PTL[sel], d$CoL[sel] ), 
              unique )
    tmp = as.numeric( unlist( tmp ) )
    dtbp$lat[ dtbp$S == n & dtbp$R == 'Foil' ] = tmp[1:4]
  }
  
  # Data frame for correlations
  cor_val = data.frame( S = 1:N, 
                        R_xi = rep( NA, N ),
                        Tau_xi = rep( NA, N ),
                        Rho_xi = rep( NA, N ),
                        R_k = rep( NA, N ),
                        Tau_k = rep( NA, N ),
                        Rho_k = rep( NA, N ) )
  
  ### Individual correlations with threshold ###
  # Drift rates
  for ( s in 1:N ) {
    sel = dtbp$S == s
    xa = scale( 1/dtbp$lat[sel] )
    ya = scale( dtbp$xi[sel] )
    cor_val$R_xi[s] = cor( xa, ya )
    cor_val$Tau_xi[s] = cor( xa, ya, method = 'kendall' )
    cor_val$Rho_xi[s] = cor( xa, ya, method = 'spearman' )
  }
  
  # Threshold
  for ( s in 1:N ) {
    sel = dtbp$S == s
    xa = scale( 1/dtbp$lat[sel] )
    ya = scale( dtbp$kappa[sel] )
    cor_val$R_k[s] = cor( xa, ya )
    cor_val$Tau_k[s] = cor( xa, ya, method = 'kendall' )
    cor_val$Rho_k[s] = cor( xa, ya, method = 'spearman' )
  }
  
  # Reset margins
  par( mar = c( 4, 5, 3, .5 ) )
  
  if ( savePlot ) 
    pdf('Correlation_distributions.pdf',width=6,height=6)
  else x11()
  layout( cbind( 1 ) )
  
  xl = c( .5, 3.5 ); yl = c( 0, 1 )
  blankPlot( xl, yl )
  
  segments( rep( xl[1], 3 ), c(.25,.5,.75), 
            rep( xl[2], 3 ), c(.25,.5,.75), col = 'grey70', lwd = 2 )
  
  quick_plot( 'R', 1 )
  quick_plot( 'Tau', 2 )
  quick_plot( 'Rho', 3 )
  
  customAxes( xl, yl )
  axis( 2, seq(0,1,.25), tick = F, line = -1.75, cex = 1 )
  axis( 1, 1:3, c( "Pearson's R", "Kendall's Tau",
                   "Spearman's Rho" ), tick = F, line = -1 )
  mtext( 'Density', side = 1, line = 1.5 )
  mtext( 'Correlation with nROUSE latencies', side = 2, line = 1.5 )
  
  #par( xpd = NA )
  #legend( 'top', c('Drift rates', 'Thresholds'),
  #        fill = c('grey50','grey80'), bty = 'n',
  #        horiz = T, inset = c(0,-.1) )
  #par( xpd = T )
  
  if ( savePlot ) dev.off()
}

###
### Figure 14
###
# Lookup - 09

if ( runCode[9] ) {
  
  # Change directory to where model results are saved
  setwd( orig_dir )
  setwd( 'RT_and_choice_MLE' )
  setwd( 'Estimation_results' )
  
  # Load in reference data
  load( 'DRM_SD_M7_PM_results.RData' )
  setwd( orig_dir )
  setwd( 'Manuscript_figures' )
  
  if ( savePlot ) 
    pdf('Parameter_estimates.pdf',width=12,height=6) else 
      x11( width = 12 )
  layout( cbind( 1, 2 ) )
  
  # Drift rates
  xl = c( .2, 8.8 )
  yl = c( -.04, 4 )
  blankPlot( xl, yl )
  
  segments( 2.5, yl[1], 2.5, yl[2], col = 'grey80', lwd = 2 )
  segments( 4.5, yl[1], 4.5, yl[2], col = 'grey80', lwd = 2 )
  segments( 6.5, yl[1], 6.5, yl[2], col = 'grey80', lwd = 2 )
  
  xi = all_prm[ , grep( 'xi', colnames( all_prm ) ) ]
  
  f = function(x) {
    q = quantile( x, prob = c(0,.25,.5,.75,1) )
    return( q )
  }
  # Different
  tmp = apply( xi[,c(5:8,1:4)+8], 2, f )
  segments( 1:8 + .1, tmp[2,], 1:8 + .1, tmp[4,],
            lwd = 2, col = 'black' )
  for ( j in 1:4 ) {
    lines( 1:2 + 2*(j-1) + .1, tmp[3,1:2 + 2*(j-1)], 
           type = 'b',lwd = 2, pch = 21, bg = 'black', 
           cex = 1.2, lty = 1 )
  }
  
  # Same
  tmp = apply( xi[,1:8], 2, f )
  segments( 1:8 - .1, tmp[2,], 1:8 - .1, tmp[4,],
            lwd = 2, col = 'black' )
  for ( j in 1:4 ) {
    lines( 1:2 + 2*(j-1) - .1, tmp[3,1:2 + 2*(j-1)], 
           type = 'b',lwd = 2, pch = 22, bg = 'white', 
           cex = 1.2, lty = 2 )
  }
  
  # Axes
  axis( 1, 1:8, rep( c(.05, .4 ), 4 ), 
        tick = F, cex = 1, line = -1.5 )
  axis( 1, seq(1,7,2)+.5, rep( c('Foil', 'Target' ), 2 ), 
        tick = F, cex = 1, line = -.5 )
  axis( 1, seq(2,6,4)+.5, c('Same','Different'), 
        tick = F, cex = 1, line = .5 )
  axis( 2, seq(0,4,1), tick = F, line = -1.5, cex = 1.2 )
  mtext( 'Median drift rate', side = 2, line = 2, cex = 1.5 )
  
  par(xpd=T)
  # segments( 4.5, -.03, 4.5, -.75, lwd = 2, col = 'grey80' )
  # segments( 2.5, -.03, 2.5, -.4, lwd = 2, col = 'grey80' )
  # segments( 6.5, -.03, 6.5, -.4, lwd = 2, col = 'grey80' )
  legend( 'top', c( 'Different', 'Same' ),
          pt.bg = c( 'black', 'white' ),
          pch = c( 21, 22 ), horiz = T, bty = 'n',
          inset = c(0,-.1), cex = 1.25 )
  par(xpd=F)
  customAxes( xl, yl )
  
  # Thresholds
  xl = c( .2, 4.8 )
  yl = c( .5, 2 )
  blankPlot( xl, yl )
  
  segments( 2.5, yl[1], 2.5, yl[2], col = 'grey80', lwd = 2 )
  segments( 4.5, yl[1], 4.5, yl[2], col = 'grey80', lwd = 2 )
  segments( 6.5, yl[1], 6.5, yl[2], col = 'grey80', lwd = 2 )
  
  kappa = all_prm[ , grep( 'kappa', colnames( all_prm ) ) ]
  
  f = function(x) {
    q = quantile( x, prob = c(0,.25,.5,.75,1) )
    return( q )
  }
  
  # Same
  tmp = apply( kappa[,1:4], 2, f )
  segments( 1:4 - .1, tmp[2,], 1:4 - .1, tmp[4,],
            lwd = 2, col = 'black' )
  for ( j in 1:2 ) {
    lines( 1:2 + 2*(j-1) - .1, tmp[3,1:2 + 2*(j-1)], 
           type = 'b',lwd = 2, pch = 22, bg = 'white', 
           cex = 1.2, lty = 2 )
  }
  
  # Different
  tmp = apply( kappa[,c(7:8,5:6)], 2, f )
  segments( 1:4 + .1, tmp[2,], 1:4 + .1, tmp[4,],
            lwd = 2, col = 'black' )
  for ( j in 1:2 ) {
    lines( 1:2 + 2*(j-1) + .1, tmp[3,1:2 + 2*(j-1)], 
           type = 'b',lwd = 2, pch = 21, bg = 'black', 
           cex = 1.2, lty = 1 )
  }
  
  # Axes
  axis( 1, 1:8, rep( c(.05, .4 ), 4 ), 
        tick = F, cex = 1, line = -1.5 )
  axis( 1, c(1.5,3.5), c('Unprimed', 'Primed' ), 
        tick = F, cex = 1, line = -.5 )
  axis( 2, seq(0,2,.5), tick = F, line = -1.5, cex = 1.2 )
  mtext( 'Median threshold', side = 2, line = 2, cex = 1.5 )
  
  par(xpd=T)
  # segments( 4.5, 0, 4.5, -.32, lwd = 2, col = 'grey80' )
  # segments( 2.5, 0, 2.5, -.17, lwd = 2, col = 'grey80' )
  # segments( 6.5, 0, 6.5, -.17, lwd = 2, col = 'grey80' )
  par(xpd=F)
  customAxes( xl, yl )
  if ( savePlot ) dev.off()
  
}

###
### Figures B.1, B.2, and B.3
###
# Lookup - 10

if (runCode[10]) {
  
  # Conditions
  # 2AFC
  # 0 - Target primed for .05s; Left;
  # 1 - Target primed for .05s; Right
  # 2 - Target primed for .4s; Left
  # 3 - Target primed for .4s; Right
  # 4 - Foil primed for .05s; Left;
  # 5 - Foil primed for .05s; Right
  # 6 - Foil primed for .4s; Left
  # 7 - Foil primed for .4s; Right
  
  # S/D
  # 8  - Target primed for .05s; Same
  # 9  - Target primed for .05s; Diff
  # 10 - Target primed for .4s; Same
  # 11 - Target primed for .4s; Diff
  # 12 - Foil primed for .05s; Same
  # 13 - Foil primed for .05s; Diff
  # 14 - Foil primed for .4s; Same
  # 15 - Foil primed for .4s; Diff
  
  ### Boxplot of accuracy ###
  
  # Extract accuracy over conditions
  tmp = aggregate( d$Ac, list(d$Cnd,d$S), mean )
  colnames(tmp) = c('Condition','Subject','P')
  prp = aggregate( tmp$P, list(tmp[,1]), quantile, prob = c(.05,.25,.5,.75,.95) )
  colnames(prp) = c('Condition','P')
  # Determine outliers
  outliers = c()
  for (i in 1:nrow(prp)) {
    x = tmp$P[ tmp$Condition == prp$Condition[i] ]
    x = x[ x < prp$P[i,1] | x > prp$P[i,5] ]
    outliers = c( outliers, list(x) )
  }
  rm(tmp)
  
  # Set width of lines and text size
  lW = 2.5
  txtS = 1.5
  
  if (savePlot) pdf('Desc_accuracy.pdf',width=6,height=6) else x11()
  par(mar=c(3,5,1,.5))
  # Create blank plot
  xl = c(.5,16.5); yl = c(0,1); blankPlot( xl, yl )
  # Define color for ???
  clr = rep( c('grey40','grey70'), 8)
  # Indicate same-different conditions
  polygon( c(4.5,4.5,4.5+8,4.5+8),
           c(0,1,1,0), col = 'grey90', border = NA )
  # Add boxplots and points for outliers
  inc = 1
  for (i in c(1:8,8:1+8)) {
    bplt( x = as.numeric(prp$P[i,]), pos = inc,
          f = function(x) return(x),
          lwd = lW, col = clr[i] )
    points( rep(inc,length(outliers[[i]])), outliers[[i]], pch = 19,
            cex = .5 )
    inc = inc + 1
  }
  segments(8.5,0,8.5,1,lwd=lW)
  axis(2,seq(0,1,.25),lwd=lW,cex.axis=txtS)
  axis(1,c(4.5,4.5+8),c('2AFC','S/D'),tick=F,lwd=lW,cex.axis=txtS)
  abline(h=0,lwd=lW)
  axis(1,c(2.5,6.5,10.5,14.5),
       c('Target','Foil','Foil','Target'),
       cex.axis=txtS,tick=F,line=-1.6)
  legend(-.5,.2,c('Left','Right'),fill=c('grey40','grey70'),
         bty='n',cex=txtS*.9)
  legend(17.5,.2,c('Same','Different'),fill=c('grey40','grey70'),
         bty='n',cex=txtS*.9,xjust=1)
  mtext('Proportion correct',2,cex=txtS,line=3)
  
  axis( 3,seq(1.5,15.5,2),
        c( rep( c('50','400'), 2 ),
           rep( c('400','50'), 2 ) ), tick = F, cex.axis = txtS*.9,
        line = -2 )
  
  if (savePlot) dev.off()
  
  # Extract median RT
  tmp = aggregate( d$RT, list(d$Cnd,d$Ac,d$S), median )
  colnames(tmp) = c('Condition','Accuracy','Subject','RT')
  rt = aggregate( tmp$RT, list(tmp[,1],tmp[,2]), quantile, prob=c(.05,.25,.5,.75,.95))
  colnames(rt) = c('Condition','Accuracy','RT')
  # Determine outliers
  outliers = c()
  for (i in 1:nrow(rt)) {
    x = tmp$RT[ tmp$Condition == rt$Condition[i] & tmp$Accuracy==rt$Accuracy[i] ]
    x = x[ x < rt$RT[i,1] | x > rt$RT[i,5] ]
    outliers = c( outliers, list(x) )
  }
  rm(tmp)
  
  
  ### Boxplot of median RTs for 2AFC ###
  
  lW = 3
  txtS = 1.5
  
  
  if (savePlot) pdf('Desc_RT_2AFC.pdf',width=6,height=6) else x11()
  par(mar=c(3,5,1,.5))
  xl = c(.5,16.5); yl = c(0,1.2); blankPlot( xl, yl )
  
  polygon( c(4.5,4.5,4.5+8,4.5+8),
           c(0,1.2,1.2,0), col='grey90', border=NA )
  
  # Accurate
  ps = c( seq(1,7,2), seq(10,16,2) )
  sel = c(1,2,3,4,8,7,6,5)
  for (i in 1:8) {
    ya = as.numeric(rt$RT[sel[i]+16,])
    bplt( x = ya, pos = ps[i],
          f = function(x) return(x),
          lwd = lW, col = 'grey40' )
    points( rep(ps[i],length(outliers[[sel[i]+16]])), 
            outliers[[sel[i]+16]], pch = 19,
            cex = .5 )
    inc = inc + 1
  }
  
  # Errors
  ps = c( seq(2,8,2), seq(9,15,2) )
  sel = c(1,2,3,4,8,7,6,5)
  for (i in 1:8) {
    ya = as.numeric(rt$RT[sel[i],])
    bplt( x = ya, pos = ps[i],
          f = function(x) return(x),
          lwd = lW, col = 'grey70' )
    points( rep(ps[i],length(outliers[[sel[i]]])), 
            outliers[[sel[i]]], pch = 19,
            cex = .5 )
    inc = inc + 1
  }
  
  segments(8.5,0,8.5,1.2,lwd=lW)
  
  axis(2,seq(0,1.2,.3),lwd=lW,cex.axis=txtS)
  axis(1,c(4.5,4.5+8),c('Target primed','Foil primed'),tick=F,lwd=lW,cex.axis=txtS)
  abline(h=0,lwd=lW)
  axis(1,c(2.5,6.5,10.5,14.5),
       c('50 ms','400 ms','400 ms','50 ms'),
       cex.axis=txtS,tick=F,line=-1.6)
  legend(-.5,.2,c('Correct','Error'),fill=c('grey40','grey70'),
         bty='n',cex=txtS*.9)
  mtext('Median response time (s)',2,cex=txtS,line=3)
  
  axis( 3,seq(1.5,15.5,2),
        c( rep( c('Lt','Rt'), 2 ),
           rep( c('Rt','Lt'), 2 ) ), tick = F, cex.axis = txtS,
        line = -2 )
  
  if ( savePlot ) dev.off()
  
  ### Same/different ###
  
  if ( savePlot ) pdf('Desc_RT_SD.pdf',width=6,height=6) else x11()
  
  par(mar=c(3,5,1,.5))
  xl = c(.5,16.5); yl = c(0,1.2); blankPlot( xl, yl )
  
  polygon( c(4.5,4.5,4.5+8,4.5+8),
           c(0,1.2,1.2,0), col='grey90', border=NA )
  
  # Accurate
  ps = c( seq(1,7,2), seq(10,16,2) )
  sel = c(1,2,3,4,8,7,6,5)
  for (i in 1:8) {
    ya = as.numeric(rt$RT[sel[i]+16+8,])
    bplt( x = ya, pos = ps[i],
          f = function(x) return(x),
          lwd = lW, col = 'grey40' )
    points( rep(ps[i],length(outliers[[sel[i]+16+8]])), 
            outliers[[sel[i]+16+8]], pch = 19,
            cex = .5 )
    inc = inc + 1
  }
  
  # Errors
  ps = c( seq(2,8,2), seq(9,15,2) )
  sel = c(1,2,3,4,8,7,6,5)
  for (i in 1:8) {
    ya = as.numeric(rt$RT[sel[i]+8,])
    bplt( x = ya, pos = ps[i],
          f = function(x) return(x),
          lwd = lW, col = 'grey70' )
    points( rep(ps[i],length(outliers[[sel[i]+8]])),
            outliers[[sel[i]+8]], pch = 19,
            cex = .5 )
    inc = inc + 1
  }
  
  segments(8.5,0,8.5,1.2,lwd=lW)
  
  axis(2,seq(0,1.2,.3),lwd=lW,cex.axis=txtS)
  axis(1,c(4.5,4.5+8),c('Target primed','Foil primed'),tick=F,lwd=lW,cex.axis=txtS)
  abline(h=0,lwd=lW)
  axis(1,c(2.5,6.5,10.5,14.5),
       c('50 ms','400 ms','400 ms','50 ms'),
       cex.axis=txtS,tick=F,line=-1.6)
  legend(-.5,.2,c('Correct','Error'),fill=c('grey40','grey70'),
         bty='n',cex=txtS*.9)
  mtext('Median response time (s)',2,cex=txtS,line=3)
  
  axis( 3,seq(1.5,15.5,2),
        c( rep( c('S','D'), 2 ),
           rep( c('D','S'), 2 ) ), tick = F, cex.axis = txtS,
        line = -2 )
  
  if (savePlot) dev.off()
}

setwd( orig_dir )