# SEM using summed scores
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directory if you 
# have any questions or comments
# Last updated 2018-04-26

# Table of contents
# 1) Initial setup
# 2) SEM following dimension reduction via PCA
#   2.1) PCA
#   2.2) Lavaan syntax
#   2.3) Model estimation with Lavaan

###
### 1) Initial setup
###

# Indicate which code segments 
# should be run
run_code = c(
  T,
  F
)

# Save current working directory
orig_dir = getwd()

# Navigate to project directory
setwd( '..' ); proj_dir = getwd();

# Indicate whether to save figures
savePlot = T
if ( savePlot ) {
  setwd( 'Plots' )
  pdf( 'SEM_summed_scores.pdf' )
  setwd( proj_dir )
}

# Load in useful packages

# Collection of useful functions
# devtools::install_github("rettopnivek/utilityf")
library( utilityf )

# Packages for easy manipulating of data frames
# install.packages( 'dplyr' )
library(dplyr)

# Package for SEM
# install.packages( 'lavaan' )
library( lavaan )

# Package for PCA analysis
# install.packages( "FactoMineR" )
library( FactoMineR )

# Package for plotting SEM path diagrams
# install.packages( 'semPlot' )
library( semPlot )

# Load in data
setwd( 'Data' )
load( 'Peer_influence_alcohol.RData' )
setwd( proj_dir )

# Standardize variables
cn = colnames( SC )
sel = c( grep( 'MISS', cn ),
         grep( 'UPPS', cn ),
         grep( 'DMQ', cn ),
         grep( 'BCEOA', cn ), 
         grep( 'Total_days', cn ),
         grep( 'Total_drinks', cn )
)

# Standardize all summary scores
my_standardize = function( x ) {
  out = ( x - mean( x, na.rm = T ) ) / sd( x, na.rm = T )
  return( out )
}
tmp = apply( SC[,sel], 2, my_standardize )
SC[,sel] = tmp

###
### 2) SEM following dimension reduction via PCA
###

if ( run_code[1] ) {
  
  # 2.1) PCA
  
  # Extract data to use with factor analysis
  sel = sel = c( grep( 'DMQ', cn ),
                 grep( 'BCEOA', cn ), 
                 grep( 'Total_days', cn ),
                 grep( 'Total_drinks', cn )
  )
  M = as.matrix( SC[,sel] )
  no_na = apply( M, 1, function(x) any( is.na( x ) ) )
  M = M[!no_na,]
  No = nrow( M )
  
  # Carry out principal components analysis 
  # to reduce variables to smaller set of 
  # components
  # Components are linear combinations of 
  # original variables, weighted by 
  # their contributions for explaining variance 
  # of a particular orthogonal dimension
  pca = PCA( M, ncp = 4, graph = FALSE )
  
  # Bootstrapped uncertainty intervals 
  # for loadings
  PCA_loadings_bootstrap = function( nr ) {
    
    # Resampling
    ord = sample( 1:nrow( M ), replace = T )
    
    NM = M[ord,]
    pca = PCA( NM, ncp = 4, graph = FALSE )
    # Extract loadings
    out = pca$var$cor
    
    return( out )
  }
  nRep = 10000
  tst = lapply( 1:nRep, PCA_loadings_bootstrap )
  BSL = array( NA, dim = c( 10, 4, nRep ) )
  for ( nr in 1:nRep ) {
    BSL[,,nr] = tst[[nr]]
  }
  
  loading_p_val = matrix( NA, 10, 4 )
  rownames( loading_p_val ) = rownames( pca$var$cor )
  for ( i in 1:10 ) {
    for ( j in 1:4 ) {
      loading_p_val[i,j] = sum( BSL[i,j,] > 0 )/nRep
      if ( loading_p_val[i,j] > .5 ) {
        loading_p_val[i,j] = 1 - loading_p_val[i,j]
      }
      
    }
  }
  
  # Adjust for multiple comparisons
  loading_p_val_adj = matrix(
    p.adjust( as.vector( loading_p_val ), "BH" ),
    10, 4, byrow = F )
  rownames( loading_p_val_adj ) = rownames( pca$var$cor )
  colnames( loading_p_val_adj ) = c( 'C1', 'C2', 'C3', 'C4' )
  
  # Matrix with eigenvalues and percent variance explained
  # pca$eig
  # First 4 components 
  
  # Determine which components to largest 
  # shared variance between different variables
  # apply( pca$var$cor, 1, which.max )
  # Note there is no guarantee that the 
  # components are interpretable...
  
  # Component 1 = DMQ
  #                 -Social motives, Coping
  #               B-CEOA
  #                 - Tension reduction, 
  #                   Liquid courage, Sexuality
  # Component 2 = Total days spent drinking,
  #               Total number of drinks
  #               DMQ
  #                 - Enhancement
  # Component 3 = B-CEOA
  #                 - Self perception
  # Component 4 = DMQ
  #                 - Conformity
  
  ### Color map of loadings
  
  CM = list(
    LM = pca$var$cor,
    cmap = list(
      z = seq( -1, 1, length = 12 ),
      clr = c(
    rgb( seq( 0, .8, length = 5 ), 1, 1 ),
    rgb( 1, 1, 1 ),
    rev( rgb( 1, seq( 0, .8, length = 5 ), 1 ) ) ) ),
    x = cbind(
      rep( 0:3, each = 10 ),
      rep( 0:3, each = 10 ),
      rep( 1:4, each = 10 ),
      rep( 1:4, each = 10 ) ),
    y = cbind(
      rep( 0:9, 4 ),
      rep( 1:10, 4 ),
      rep( 1:10, 4 ),
      rep( 0:9, 4 ) )
    )
  CM$clr = data.frame( clr = rep( rgb( 0,0,0 ), 40 ),
                       x = NA, y = NA,
                       l = NA, 
                       z = NA, 
                       stringsAsFactors = F )
  inc = 1
  for ( i in 1:4 ) {
    for ( j in 10:1 ) {
      sel = sum( CM$LM[j,i] >= CM$cmap$z[ -12 ] )
      CM$clr$clr[inc] = 
        CM$cmap$clr[sel]
      CM$clr$x[inc] = i; CM$clr$y[inc] = j
      CM$clr$l[inc] = round( CM$LM[j,i] * 100 )
      CM$clr$z[inc] = CM$cmap$z[sel]
      inc = inc + 1
    }
  }
  
  if ( !savePlot ) x11();
  par( mar = c( 5, 9, 3, 3 ) )
  blankPlot( c(0,4), c(0,10) )
  for ( i in 1:40 ) {
    polygon( CM$x[i,],
             CM$y[i,],
             border = NA, col = CM$clr$clr[i] )
  }
  
  ya = 1:10 - .5
  lbl = c( 'Total drinks', 'Total days', 'Tension red.',
     'Sexuality', 'Self-percept.', 'Liquid courage',
     'Conformity', 'Enhancement', 'Coping', 'Social motives' )
  axis( 2, ya, lbl, tick = F, las = 1,
        line = -.5 )
  axis( 2, c( 4, 8 ), c( 'B-CEOA', 'DMQ' ),
        tick = F, line = 6.5 )
  axis( 2, c( 2.25, 6.25 ), c( '_______', '_______' ),
        tick = F, line = .5, las = 1 )
  axis( 1, 1:4 - .5,
        c( 'Dim. 1', 'Dim. 2', 'Dim. 3', 'Dim. 4' ),
        tick = F, line = 0 )
  
  legend( 4, 7,
          as.character( seq( -1, 1, .2 ) ), 
          fill = CM$cmap$clr,
          xpd = TRUE, bty = 'n' )
  
  title( 'PCA component loadings' )

  # 2.2) Lavaan syntax
  
  # Define data to be fitted
  dtbf = SC[ !no_na, ]
  pc = pca$ind$coord[,1:4]
  dtbf$PCA1 = pc[,1]
  dtbf$PCA2 = pc[,2]
  dtbf$PCA3 = pc[,3]
  dtbf$PCA4 = pc[,4]
  
  full_path_syntax = "
  # Measurement models
  
  # MISS
  Suggest =~ MISS_CS + MISS_P + MISS_PS + 
             MISS_PR + MISS_PC
  # UPPS-P
  Rash_emotions =~ UPPS_NU + UPPS_PU
  Deficit_diligence =~ UPPS_Pr + UPPS_Pe
  
  # Regressions
  Suggest ~ Rash_emotions + Deficit_diligence + UPPS_SS
  PCA1 ~ Suggest + Rash_emotions + Deficit_diligence + UPPS_SS
  PCA2 ~ Suggest + Rash_emotions + Deficit_diligence + UPPS_SS
  PCA3 ~ Suggest + Rash_emotions + Deficit_diligence + UPPS_SS
  PCA4 ~ Suggest + Rash_emotions + Deficit_diligence + UPPS_SS
  
  "
  
  # 2.3) Model estimation with Lavaan
  
  full_path = sem( model = full_path_syntax, 
                   data = dtbf )
  
  # x11()
  # semPaths( full_path, what = "std", 
  #           residuals = FALSE, edge.label.cex = 1.25 )
  
}

if ( savePlot ) dev.off()
setwd( orig_dir )