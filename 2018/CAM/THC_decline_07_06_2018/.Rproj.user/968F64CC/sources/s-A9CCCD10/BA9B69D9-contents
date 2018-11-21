
setwd( fig_dir )
pdf( 'Scatter_plot_of_age_variables.pdf' )

# x11()
lyt = rbind(
  c( 1, 7, 8 ),
  c( 4, 2, 9 ),
  c( 5, 6, 3 ) )
layout( lyt )

ttlSz = 1.5

# Variables

# par( mar = c( 3, 4, 3, 1 )

blankPlot()
legend( 'center', 
        'Overall age',
        bty = 'n',
        cex = ttlSz )
blankPlot()
legend( 'center', 
        'Age of initiation',
        bty = 'n',
        cex = ttlSz )
blankPlot()
legend( 'center', 
        'Years of use',
        bty = 'n',
        cex = ttlSz )

# Scatter plots

qsp = function( x, y ) {
  zx = my_standardize( x )
  zy = my_standardize( y )
  
  ll = max( abs( zx ), abs( zy ) )
  ll = c( -ll, ll )
  
  blankPlot( ll, ll )
  points( zx, zy, pch = 19 )
  customAxes( ll, ll )
  
  lmr = lm( zy ~ 0 + zx )
  sm = summary( lmr )
  p_val = sm$coefficients[1,4]
  if ( p_val < .001 ) {
    p_val_sr = ', p < .001'
  } else {
    p_val_sr = paste(', p =', round( p_val, 3 ) )
  }
  abline( lmr )
  
  title( paste( 'R = ', round( coef( lmr ), 2 ), p_val_sr, sep = '' ) )
  
}

qsp( dtbf$Age, dtbf$Age_first_used_MJ )
qsp( dtbf$Age, dtbf$Years_of_MJ_use )
qsp( dtbf$Age_first_used_MJ, dtbf$Years_of_MJ_use )

# Histograms
hist( dtbf$Age, col = 'grey', border = 'white',
      xlab = ' ', main = ' ' )
hist( dtbf$Age_first_used_MJ, col = 'grey', border = 'white',
      xlab = ' ', main = ' ' )
hist( dtbf$Years_of_MJ_use, col = 'grey', border = 'white',
      xlab = ' ', main = ' ' )

dev.off()
setwd( R_dir )

