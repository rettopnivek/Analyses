#--------------------#
# Useful functions   #
# Kevin Potter       #
# Updated 12/01/2016 #
#--------------------#

# Index
# Lookup - 01:  param_est
# Lookup - 02:  p_f
# Lookup - 03:  model_structure_create
# Lookup - 04:  convergence_extract
# Lookup - 05:  find_dec
# Lookup - 06:  wr_ui
# Lookup - 07:  gl_sim
# Lookup - 08:  prior_f

# Lookup - 01
param_est = function( X, coef, fixed, index, 
                      parSel ) {
  # Purpose:
  # A function to calculate the weighted sum of a design matrix  
  # and a set of coefficients to produce a matrix of parameters 
  # by observations. Fixed coefficient values can also be specified.
  # Arguments:
  # X      - A design matrix
  # coef   - A vector of C coefficients
  # fixed  - A vector of F values for the fixed coefficients
  # index  - A matrix with two columns, the first giving the 
  #          row positions for the free coefficients (followed by the 
  #          row positions for any fixed values), the second giving the 
  #          column positions for the free coefficients (followed by 
  #          the positions for the fixed values). The first column of 
  #          the final row of the matrix indicates the total number of 
  #          desired parameters so that the parameter matrix can be 
  #          created.
  # parSel - A vector giving the indices mapping the values in the 
  #          coefficient vector to the rows of the 'index' matrix, 
  #          thereby allowing different conditions to be constrained 
  #          to have the same free coefficient.
  # Details:
  # Given N observations and P desired parameters, the goal is to produce 
  # a P x N matrix given a set of V covariates, C coefficients, and F 
  # fixed values. To do so, a P x V parameter matrix M is specified, and 
  # the 'index' matrix along with the 'parSel' vector are used to fill 
  # the positions of the P x V matrix. Fixed values from the 'fixed' 
  # vector are additionally included in the matrix M. To produce the 
  # desired P x N output matrix, the P x V matrix M is multiplied by 
  # the V x N design matrix X.
  # Returns:
  # A P x N matrix giving the set of parameter values for each of the 
  # N observations.
  
  N <- ncol(X)
  V <- nrow(X)
  C <- length( parSel )
  L <- dim(index)[1]
  P <- index[ L, 1 ]
  finish <- L - 1 - C
  
  pm <- matrix( 0.0, P, V )
  
  for (i in 1:C ) {
    pm[ index[i,1], index[i,2] ] <- coef[ parSel[i] ]
  }
  
  if (finish > 0) {
    for (i in 1:finish) {
      pm[ index[i+C,1], 
          index[i+C,2] ] <- fixed[i]
    }
  }
  
  return( pm %*% X )
}

# Lookup - 02
p_f = function(x) {
  # Purpose;
  # A function to calculate choice probabilities.
  # Arguments:
  # x - A vector of binary choices (0 or 1)
  # Returns:
  # A vector of two values, the proportions for 
  # picking 0 and 1 respectively.
  
  out = c(0,0); names(out) = c("0","1")
  p = table(x);
  out[ names(p) ] = as.vector( p )
  return( out/sum(out) )
}

# Lookup - 03
model_structure_create = function( type, cD, Priors, s = NULL ) {
  # Purpose:
  # A function to create the necessary input to fit a 
  # Bayesian response time model in Stan (either a 
  # single subject or hierarchical).
  # Arguments:
  # type   - The desired type of model structure and inputs to 
  #          generate, where...
  #          Forthcoming
  # cD     - A data frame with the covariates and dependent variables
  # Priors - A matrix with the desired parameter values for the priors
  # s      - An optional parameter indicating the subset of data 
  #          for a particular subject to fit (if NULL, the 
  #          hierarchical model is fit)
  # Returns:
  # A list of the inputs.
  
  # There are a total of 5 model types of interest:
  # Model 1: The effects of interest are determined by 
  #          thresholds primarily
  # Model 2: The effects of interest are determined by 
  #          drift rates primarily
  # Model 3: The effects of interest are determined by 
  #          by a hybrid of drift rates and thresholds
  # Model 4: The effects of interest are determined by 
  #          a hybrid of thresholds and drift rates 
  #          constrained by nROUSE inverse latencies
  # Model 5: A over-parameterized model
  
  # Restrict to forced-choice only
  cD = cD[ cD$Ta == 0, ]
  
  # Number of observations by subject
  No = aggregate( rep(1,nrow(cD)),list(cD$S), sum )$x
  No_max = max( No ) # Largest number of observations
  Ns = length( No ) # Number of subjects
  allS = cD$S
  sbj = unique( allS ) # Subject identifiers
  
  # Fastest RTs by subject
  min_RT = aggregate( cD$RT, list( allS ), min )$x
  
  # Abbreviated subscripts key
  # S - short
  # L - long
  # T - target
  # F - foil
  # P - Primed choice
  # B - Bias
  # I - Intercept
  
  # Model 1
  if ( type == 1 ) {
    
    # Coefficients
    # kappa     -> S, L, PS, PL, B (1:5)
    # xi target -> S, L            (6:7)
    # xi foil   -> S, L            (8:9)
    # tau       -> I               (10)
    Cf = 10 # 10 coefficients in total
    
    # Create index for linear algebra function
    Clm = c( 1:5,   # Rows in design matrix for kappa(1)
             6:9,   # Rows in design matrix for xi(1)
             11,    # Rows in design matrix for tau(1)
             12:16, # Rows in design matrix for kappa(0)
             17:20, # Rows in design matrix for xi(0)
             22     # Rows in design matrix for tau(0)
    )
    Rws = c( rep( 1, 5 ),   # kappa (1)
             rep( 2, 4 ),   # xi (1)
             rep( 4, 1 ),   # tau (1)
             rep( 5, 5 ),   # kappa (0)
             rep( 6, 4 ),   # xi (0)
             rep( 8, 1 )    # tau (0)
    )
    # Create index for parameter selection
    parSel = c( 1:5, 6:9, 10, 1:5, 6:9, 10 )
    
    index = cbind( Rws, Clm )
    # Fixed values
    fixed = c(1,1)
    index = rbind( index,
                   c(3,10), # sigma (1)
                   c(7,21) # sigma (0)
    )
    rm( Clm, Rws )
    # Dimensions
    index = rbind( index, c(8,length(parSel)) )
    
    if ( length(s) == 0 ) {
      
      # Define set of arrays for design matrix and 
      # for matrix of response pairs
      X = array( 0, dim = c( Ns, length( parSel ) + 
                               length( fixed ), No_max ) )
      Y = array( 0, dim = c( Ns, No_max, 2 ) )
      
      # Create a progress bar using a base R function
      pb = txtProgressBar( min = 1, max = Ns, style = 3 )
      
      for ( s in 1:Ns ) {
        
        # Extract data
        cS = cD[ allS == sbj[s], ]
        
        # Extract response pairs
        Y[s, 1:No[s], ] = cbind( cS$RT, cS$Ch )
        
        # Define covariates
        Int = rep( 1, nrow( cS ) ) # Intercept
        Co = cS$Co # Correct (Left/Right)
        PT = cS$PT # Prime type (Foil/Target)
        PD = cS$PD # Duration (Short/Long)
        PL = as.numeric( PT != Co ) # Left choice primed
        PR = as.numeric( PT == Co ) # Right choice primed
        # Priming duration x type
        DvT = Int;
        DvT[ cS$PD == 1 & cS$PT == 1 ] = 2;
        DvT[ cS$PD == 0 & cS$PT == 0 ] = 3;
        DvT[ cS$PD == 1 & cS$PT == 0 ] = 4;
        X_x = designCoding( DvT, Levels = 1:4, 
                            type = 'Intercept' )
        
        # Desired output
        # 8 x N matrix
        # Necessary input
        # 8 x row(X) -> parameter matrix
        # row(X) x N -> Design matrix
        
        # Design matrix
        X[ s, , 1:No[s] ] = 
          t( cbind( 
            (1-PD), PD, PR*(1-PD), PR*PD, Int,
            PD*Co, (1-PD)*Co, PD*(1-Co), (1-PD)*(1-Co),
            Int,
            Int,
            (1-PD), PD, PL*(1-PD), PL*PD, -Int,
            PD*(1-Co), (1-PD)*(1-Co), PD*Co, (1-PD)*Co,
            Int,
            Int
          ) )
        
        if ( No[s] == No_max ) {
          
          # Create small design matrix for later plotting etc...
          nCnd = length( sort( unique( cS$Cnd ) ) )
          X_small = matrix( NA, dim( X )[2], nCnd )
          for ( nc in 1:nCnd ) {
            X_small[,nc] = X[ s, , cS$Cnd == nc ][,1]
          }
          curCnd = sort( unique( cS$Cnd ) )
          
        }
        
        # Update the progress bar
        setTxtProgressBar(pb,s)
      }
      close(pb)
      
      # Redefine fastest RT object
      min_RT = as.array( min_RT )
      
      # Return results
      return( list(
        Ns = Ns, 
        No = No, 
        No_max = No_max, 
        V = dim(X)[2], 
        K = length( parSel ),
        U = length( fixed ), 
        C = Cf,
        X = X, 
        fixed = fixed, 
        index = index, 
        parSel = parSel, 
        Y = Y, 
        min_RT = min_RT, 
        Priors = Priors, 
        X_small = X_small,
        curCnd = curCnd,
        mName = 'WR_MS_M1.stan') )
      
    } else {
      
      # Extract data
      cS = cD[ allS == sbj[s], ]
      
      # Extract response pairs
      Y = cbind( cS$RT, cS$Ch )
      
      # Define covariates
      Int = rep( 1, nrow( cS ) ) # Intercept
      Co = cS$Co # Correct (Left/Right)
      PT = cS$PT # Prime type (Foil/Target)
      PD = cS$PD # Duration (Short/Long)
      PL = as.numeric( PT != Co ) # Left choice primed
      PR = as.numeric( PT == Co ) # Right choice primed
      # Priming duration x type
      DvT = Int;
      DvT[ cS$PD == 1 & cS$PT == 1 ] = 2;
      DvT[ cS$PD == 0 & cS$PT == 0 ] = 3;
      DvT[ cS$PD == 1 & cS$PT == 0 ] = 4;
      X_x = designCoding( DvT, Levels = 1:4, 
                          type = 'Intercept' )
      
      # Desired output
      # 8 x N matrix
      # Necessary input
      # 8 x row(X) -> parameter matrix
      # row(X) x N -> Design matrix
      
      # Design matrix
      X = 
        t( cbind( 
          (1-PD), PD, PR*(1-PD), PR*PD, Int,
          PD*Co, (1-PD)*Co, PD*(1-Co), (1-PD)*(1-Co),
          Int,
          Int,
          (1-PD), PD, PL*(1-PD), PL*PD, -Int,
          PD*(1-Co), (1-PD)*(1-Co), PD*Co, (1-PD)*Co,
          Int,
          Int
        ) )
      
      # Redefine minimum RT
      min_RT = min_RT[ sbj[s] ]
      
      # Create small design matrix for later plotting etc...
      nCnd = length( sort( unique( cS$Cnd ) ) )
      cndVal = sort( unique( cS$Cnd ) )
      X_small = matrix( NA, dim( X )[1], nCnd )
      for ( nc in 1:nCnd ) {
        X_small[,nc] = X[ , cS$Cnd == cndVal[nc] ][,1]
      }
      curCnd = sort( unique( cS$Cnd ) )
      
      # Return results
      return( list(
        N = ncol(X), 
        V = dim(X)[1], 
        K = length( parSel ),
        U = length( fixed ), 
        C = Cf,
        X = X, 
        fixed = fixed, 
        index = index, 
        parSel = parSel, 
        Y = Y, 
        min_RT = min_RT, 
        Priors = Priors, 
        X_small = X_small,
        curCnd = curCnd,
        mName = 'WR_OS_M1.stan') )
      
    }
    
  }
  
  # Model 2
  if ( type == 2 ) {
    
    # Coefficients
    # kappa     -> I                  (1)
    # xi target -> STP, LTP, SFP, LFP (2:5)
    # xi foil   -> STP, LTP, SFP, LFP (7:9)
    # tau       -> I                  (10)
    Cf = 10 # 10 coefficients in total
    
    # Create index for linear algebra function
    Clm = c( 1,     # Rows in design matrix for kappa(1)
             2:9,   # Rows in design matrix for xi(1)
             11,    # Rows in design matrix for tau(1)
             12,    # Rows in design matrix for kappa(0)
             13:20, # Rows in design matrix for xi(0)
             22     # Rows in design matrix for tau(0)
    )
    Rws = c( rep( 1, 1 ),   # kappa (1)
             rep( 2, 8 ),   # xi (1)
             rep( 4, 1 ),   # tau (1)
             rep( 5, 1 ),   # kappa (0)
             rep( 6, 8 ),   # xi (0)
             rep( 8, 1 )    # tau (0)
    )
    # Create index for parameter selection
    parSel = c( 1, 2:9, 10, 1, 2:9, 10 )
    
    index = cbind( Rws, Clm )
    # Fixed values
    fixed = c(1,1)
    index = rbind( index,
                   c(3,10), # sigma (1)
                   c(7,21) # sigma (0)
    )
    rm( Clm, Rws )
    # Dimensions
    index = rbind( index, c(8,length(parSel)) )
    
    if ( length(s) == 0 ) {
      
      # Define set of arrays for design matrix and 
      # for matrix of response pairs
      X = array( 0, dim = c( Ns, length( parSel ) + 
                               length( fixed ), No_max ) )
      Y = array( 0, dim = c( Ns, No_max, 2 ) )
      
      # Create a progress bar using a base R function
      pb = txtProgressBar( min = 1, max = Ns, style = 3 )
      
      for ( s in 1:Ns ) {
        
        # Extract data
        cS = cD[ allS == sbj[s], ]
        
        # Extract response pairs
        Y[s, 1:No[s], ] = cbind( cS$RT, cS$Ch )
        
        # Define covariates
        Int = rep( 1, nrow( cS ) ) # Intercept
        Co = cS$Co # Correct (Left/Right)
        PT = cS$PT # Prime type (Foil/Target)
        PD = cS$PD # Duration (Short/Long)
        PL = as.numeric( PT != Co ) # Left choice primed
        PR = as.numeric( PT == Co ) # Right choice primed
        # Priming duration x type
        DvT = Int;
        DvT[ cS$PD == 1 & cS$PT == 1 ] = 2;
        DvT[ cS$PD == 0 & cS$PT == 0 ] = 3;
        DvT[ cS$PD == 1 & cS$PT == 0 ] = 4;
        X_x = designCoding( DvT, Levels = 1:4, 
                            type = 'Intercept' )
        
        # Desired output
        # 8 x N matrix
        # Necessary input
        # 8 x row(X) -> parameter matrix
        # row(X) x N -> Design matrix
        
        # Design matrix
        X[ s, , 1:No[s] ] = 
          t( cbind( 
            Int,
            X_x*Co, X_x*(1-Co),
            Int,
            Int,
            Int,
            X_x*(1-Co), X_x*Co,
            Int,
            Int
          ) )
        
        if ( No[s] == No_max ) {
          
          # Create small design matrix for later plotting etc...
          nCnd = length( sort( unique( cS$Cnd ) ) )
          X_small = matrix( NA, dim( X )[2], nCnd )
          for ( nc in 1:nCnd ) {
            X_small[,nc] = X[ s, , cS$Cnd == nc ][,1]
          }
          curCnd = sort( unique( cS$Cnd ) )
          
        }
        
        # Update the progress bar
        setTxtProgressBar(pb,s)
      }
      close(pb)
      
      # Redefine fastest RT object
      min_RT = as.array( min_RT )
      
      # Return results
      return( list(
        Ns = Ns, 
        No = No, 
        No_max = No_max, 
        V = dim(X)[2], 
        K = length( parSel ),
        U = length( fixed ), 
        C = Cf,
        X = X, 
        fixed = fixed, 
        index = index, 
        parSel = parSel, 
        Y = Y, 
        min_RT = min_RT, 
        Priors = Priors, 
        X_small = X_small,
        curCnd = curCnd,
        mName = 'WR_MS_M2.stan') )
      
    } else {
      
      # Extract data
      cS = cD[ allS == sbj[s], ]
      
      # Extract response pairs
      Y = cbind( cS$RT, cS$Ch )
      
      # Define covariates
      Int = rep( 1, nrow( cS ) ) # Intercept
      Co = cS$Co # Correct (Left/Right)
      PT = cS$PT # Prime type (Foil/Target)
      PD = cS$PD # Duration (Short/Long)
      PL = as.numeric( PT != Co ) # Left choice primed
      PR = as.numeric( PT == Co ) # Right choice primed
      # Priming duration x type
      DvT = Int;
      DvT[ cS$PD == 1 & cS$PT == 1 ] = 2;
      DvT[ cS$PD == 0 & cS$PT == 0 ] = 3;
      DvT[ cS$PD == 1 & cS$PT == 0 ] = 4;
      X_x = designCoding( DvT, Levels = 1:4, 
                          type = 'Intercept' )
      
      # Desired output
      # 8 x N matrix
      # Necessary input
      # 8 x row(X) -> parameter matrix
      # row(X) x N -> Design matrix
      
      # Design matrix
      X = 
        t( cbind( 
          Int,
          X_x*Co, X_x*(1-Co),
          Int,
          Int,
          Int,
          X_x*(1-Co), X_x*Co,
          Int,
          Int
        ) )
      
      # Redefine minimum RT
      min_RT = min_RT[ sbj[s] ]
      
      # Create small design matrix for later plotting etc...
      nCnd = length( sort( unique( cS$Cnd ) ) )
      cndVal = sort( unique( cS$Cnd ) )
      X_small = matrix( NA, dim( X )[1], nCnd )
      for ( nc in 1:nCnd ) {
        X_small[,nc] = X[ , cS$Cnd == cndVal[nc] ][,1]
      }
      curCnd = sort( unique( cS$Cnd ) )
      
      # Return results
      return( list(
        N = ncol(X), 
        V = dim(X)[1], 
        K = length( parSel ),
        U = length( fixed ), 
        C = Cf,
        X = X, 
        fixed = fixed, 
        index = index, 
        parSel = parSel, 
        Y = Y, 
        min_RT = min_RT, 
        Priors = Priors, 
        X_small = X_small,
        curCnd = curCnd,
        mName = 'WR_OS_M2.stan') )
      
    }
    
  }
  
  # Model 3
  if ( type == 3 ) {
    
    # Coefficients
    # kappa     -> S, L, PS, PL, B    (1:5)
    # xi target -> STP, LTP, SFP, LFP (6:9)
    # xi foil   -> STP, LTP, SFP, LFP (10:13)
    # tau       -> I                  (14)
    Cf = 14 # 10 coefficients in total
    
    # Create index for linear algebra function
    Clm = c( 1:5,   # Rows in design matrix for kappa(1)
             6:13,  # Rows in design matrix for xi(1)
             15,    # Rows in design matrix for tau(1)
             16:20, # Rows in design matrix for kappa(0)
             21:28, # Rows in design matrix for xi(0)
             30     # Rows in design matrix for tau(0)
    )
    Rws = c( rep( 1, 5 ),   # kappa (1)
             rep( 2, 8 ),   # xi (1)
             rep( 4, 1 ),   # tau (1)
             rep( 5, 5 ),   # kappa (0)
             rep( 6, 8 ),   # xi (0)
             rep( 8, 1 )    # tau (0)
    )
    # Create index for parameter selection
    parSel = c( 1:5, 6:13, 14, 1:5, 6:13, 14 )
    
    index = cbind( Rws, Clm )
    # Fixed values
    fixed = c(1,1)
    index = rbind( index,
                   c(3,14), # sigma (1)
                   c(7,29) # sigma (0)
    )
    rm( Clm, Rws )
    # Dimensions
    index = rbind( index, c(8,length(parSel)) )
    
    if ( length(s) == 0 ) {
      
      # Define set of arrays for design matrix and 
      # for matrix of response pairs
      X = array( 0, dim = c( Ns, length( parSel ) + 
                               length( fixed ), No_max ) )
      Y = array( 0, dim = c( Ns, No_max, 2 ) )
      
      # Create a progress bar using a base R function
      pb = txtProgressBar( min = 1, max = Ns, style = 3 )
      
      for ( s in 1:Ns ) {
        
        # Extract data
        cS = cD[ allS == sbj[s], ]
        
        # Extract response pairs
        Y[s, 1:No[s], ] = cbind( cS$RT, cS$Ch )
        
        # Define covariates
        Int = rep( 1, nrow( cS ) ) # Intercept
        Co = cS$Co # Correct (Left/Right)
        PT = cS$PT # Prime type (Foil/Target)
        PD = cS$PD # Duration (Short/Long)
        PL = as.numeric( PT != Co ) # Left choice primed
        PR = as.numeric( PT == Co ) # Right choice primed
        # Priming duration x type
        DvT = Int;
        DvT[ cS$PD == 1 & cS$PT == 1 ] = 2;
        DvT[ cS$PD == 0 & cS$PT == 0 ] = 3;
        DvT[ cS$PD == 1 & cS$PT == 0 ] = 4;
        X_x = designCoding( DvT, Levels = 1:4, 
                            type = 'Intercept' )
        
        # Desired output
        # 8 x N matrix
        # Necessary input
        # 8 x row(X) -> parameter matrix
        # row(X) x N -> Design matrix
        
        # Design matrix
        X[ s, , 1:No[s] ] = 
          t( cbind( 
            (1-PD), PD, PR*(1-PD), PR*PD, Int,
            X_x*Co, X_x*(1-Co),
            Int,
            Int,
            (1-PD), PD, PR*(1-PD), PR*PD, Int,
            X_x*(1-Co), X_x*Co,
            Int,
            Int
          ) )
        
        if ( No[s] == No_max ) {
          
          # Create small design matrix for later plotting etc...
          nCnd = length( sort( unique( cS$Cnd ) ) )
          X_small = matrix( NA, dim( X )[2], nCnd )
          for ( nc in 1:nCnd ) {
            X_small[,nc] = X[ s, , cS$Cnd == nc ][,1]
          }
          curCnd = sort( unique( cS$Cnd ) )
          
        }
        
        # Update the progress bar
        setTxtProgressBar(pb,s)
      }
      close(pb)
      
      # Redefine fastest RT object
      min_RT = as.array( min_RT )
      
      # Return results
      return( list(
        Ns = Ns, 
        No = No, 
        No_max = No_max, 
        V = dim(X)[2], 
        K = length( parSel ),
        U = length( fixed ), 
        C = Cf,
        X = X, 
        fixed = fixed, 
        index = index, 
        parSel = parSel, 
        Y = Y, 
        min_RT = min_RT, 
        Priors = Priors, 
        X_small = X_small,
        curCnd = curCnd,
        mName = 'WR_MS_M3.stan') )
      
    } else {
      
      # Extract data
      cS = cD[ allS == sbj[s], ]
      
      # Extract response pairs
      Y = cbind( cS$RT, cS$Ch )
      
      # Define covariates
      Int = rep( 1, nrow( cS ) ) # Intercept
      Co = cS$Co # Correct (Left/Right)
      PT = cS$PT # Prime type (Foil/Target)
      PD = cS$PD # Duration (Short/Long)
      PL = as.numeric( PT != Co ) # Left choice primed
      PR = as.numeric( PT == Co ) # Right choice primed
      # Priming duration x type
      DvT = Int;
      DvT[ cS$PD == 1 & cS$PT == 1 ] = 2;
      DvT[ cS$PD == 0 & cS$PT == 0 ] = 3;
      DvT[ cS$PD == 1 & cS$PT == 0 ] = 4;
      X_x = designCoding( DvT, Levels = 1:4, 
                          type = 'Intercept' )
      
      # Desired output
      # 8 x N matrix
      # Necessary input
      # 8 x row(X) -> parameter matrix
      # row(X) x N -> Design matrix
      
      # Design matrix
      X = 
        t( cbind( 
          (1-PD), PD, PR*(1-PD), PR*PD, Int,
          X_x*Co, X_x*(1-Co),
          Int,
          Int,
          (1-PD), PD, PR*(1-PD), PR*PD, Int,
          X_x*(1-Co), X_x*Co,
          Int,
          Int
        ) )
      
      # Redefine minimum RT
      min_RT = min_RT[ sbj[s] ]
      
      # Create small design matrix for later plotting etc...
      nCnd = length( sort( unique( cS$Cnd ) ) )
      cndVal = sort( unique( cS$Cnd ) )
      X_small = matrix( NA, dim( X )[1], nCnd )
      for ( nc in 1:nCnd ) {
        X_small[,nc] = X[ , cS$Cnd == cndVal[nc] ][,1]
      }
      curCnd = sort( unique( cS$Cnd ) )
      
      # Return results
      return( list(
        N = ncol(X), 
        V = dim(X)[1], 
        K = length( parSel ),
        U = length( fixed ), 
        C = Cf,
        X = X, 
        fixed = fixed, 
        index = index, 
        parSel = parSel, 
        Y = Y, 
        min_RT = min_RT, 
        Priors = Priors, 
        X_small = X_small,
        curCnd = curCnd,
        mName = 'WR_OS_M3.stan') )
      
    }
    
  }
  
  # Model 4
  # Forthcoming
  
  # Model 5
  # Forthcoming
  
}

# Lookup - 04
convergence_extract = function( fit, par_name = NULL ) {
  # Purpose:
  # Extract convergence diagnostics from a Stan fit object.
  # Arguments:
  # fit      - A Stan fit object
  # par_name - An optional string giving the final parameter label 
  #            of the subset of the output to include
  # Notes:
  # Extracting the convergence statistics can be slow, especially 
  # when a large number of parameters were stored.
  # Returns:
  # A list with the Gelman-Rubin convergence statistic, the 
  # effective number of samples, and the total number of samples.
  
  Rhat = summary(fit)$summary[,"Rhat"]
  n_eff = summary(fit)$summary[,"n_eff"]
  totSampSize = length(extract(fit, pars = "lp__")[[1]])
  # We're only interested in a subset of parameters
  if ( length( par_name ) == 0 ) 
    par_name = names( Rhat )[ 
      which( names( Rhat ) == "log_lik[1]" ) - 1 ]
  sel = 1:which( names(Rhat) == par_name )
  Rhat = Rhat[sel]; n_eff = n_eff[sel];
  
  return( list( Rhat = Rhat, n_eff = n_eff, totSampSize = totSampSize ) )
}

# Lookup - 05
find_dec = function( x, spacing = 10 ) {
  # Purpose:
  # Determines the rounded leading digit and the 
  # number of trailing zeros for a number or the 
  # number of decimal places.
  # Arguments:
  # x       - A vector of values
  # spacing - The value whose exponent should be increased
  # Returns:
  # A vector giving the leading digit, the number of 
  # trailing/leading zeros, the same but in scientific 
  # notation, and 1 if it's trailing zeros, -1 if it's 
  # decimal places.
  
  mx = max( x )
  rnd = mx
  
  if ( round( mx ) > 1 ) {
    
    inc = 0;
    
    while ( rnd > 1 ) {
      inc = inc + 1;
      rnd = round( mx/( spacing^inc ) )
    }
    
    v = round( mx/spacing^(inc-1) )
    f = spacing^(inc-1)
    out = c( v,f,inc-1,1)
  }
  
  if ( round( mx ) == 1 ) {
    
    out = c( 1, 1, 1, 0 )
    
  }
  
  if ( round( mx ) == 0 ) {
    
    inc = 0;
    
    while ( rnd < 1 ) {
      inc = inc + 1;
      rnd = round( mx*( spacing^inc ) )
    }
    
    v = round( mx*spacing^(inc) )
    f = spacing^(inc)
    out = c( v,f,inc,-1)
    
  }
  
  return( out )
}

# Lookup - 06
wr_ui = function( RC, clm, t, ch = 1, prb = c(.025,.975), ... ) {
  # Purpose:
  # Given an array of posterior samples for the wald race parameters 
  # over a set of conditions, calculates the lower and upper boundaries 
  # for the associated credible intervals and adds them to an existing plot.
  # Arguments:
  # RC  - A S x 8 x C array, where S is the number of posterior samples 
  #       and C is the number of conditions
  # clm - The slice (condition) of the array to plot the credible intervals
  # t   - The vector of sorted response times to plot over
  # ch  - The choice to condition on
  # prb - The lower and upper boundaries of the uncertainty interval
  # ... - Additional plotting parameters for the 'polygon' function
  
  # Determine uncertainty intervals
  ui = apply( RC[,,clm], 2, quantile, prob = prb )
  
  lb = pwaldrace( t, ch, ui[1,1], ui[1,2], ui[1,4],
                  ui[1,5], ui[1,6], ui[1,8], ui[1,3], ui[1,7] )
  ub = pwaldrace( t, ch, ui[2,1], ui[2,2], ui[2,4],
                  ui[2,5], ui[2,6], ui[2,8], ui[2,3], ui[2,7] )
  polygon( c(t,rev(t)), c(lb,rev(ub)), ... )
}

# Lookup - 07
gl_sim = function( s, post, input, type ) {
  # Purpose:
  # Given a list of hierarchical posterior samples and an index number,
  # simulates subject-level parameters and arranges them by condition.
  # Arguments:
  # s     - The particular sample index
  # post  - The list of posterior samples
  # input - The original input for Stan
  # type  - The model type (1-5)
  # Returns:
  # An 8 x C x N array, where C is the number of conditions and N is the 
  # number of subjects, giving the associated 8 wald race parameters per 
  # each.
  
  # Determine the average fastest RT
  avgfrt = mean( input$min_RT )
  
  # Model 1
  if ( type == 1 ) {
    
    mu = c( post$kappa_mu[s,], 0.0, 
            post$xi_mu[s,], post$theta_alpha[s] )
    sig = c( post$kappa_sig[s,], 0.0, 
             post$xi_sig[s,], post$theta_beta[s] )
    
  }
  
  # Model 1
  if ( type == 2 ) {
    
    mu = c( post$kappa_mu[s], 
            post$xi_mu[s,], post$theta_alpha[s] )
    sig = c( post$kappa_sig[s], 
             post$xi_sig[s,], post$theta_beta[s] )
    
  }
  
  # Model 3
  if ( type == 3 ) {
    
    tmp = apply( post$kappa_mu, 2, findMode )
    mu = c( tmp[1:2], 0.0, tmp[3], 0.0, 
            post$xi_mu[s,], post$theta_alpha[s] )
    sig = c( post$kappa_sig[s,], 
             post$xi_sig[s,], post$theta_beta[s] )
    
  }
  
  
  # Simulate subject-level parameters
  badVal = T
  while ( badVal ) {
    
    nP = input$C
    cf = c( rnorm( nP-1, mu[-nP], sig[-nP] ), rbeta( 1, mu[nP], sig[nP])*avgfrt )
    
    # Calculate corresponding parameter estimates per condition
    prm = param_est( input$X_small, cf, input$fixed, input$index,
                     input$parSel )
    if ( sum( prm < 0 ) == 0 ) badVal = F
    
  }
  
  # Return results
  return( prm )
}

# Lookup - 08
prior_f = function( priors, ver, v = seq(0,4,length=1000) ) {
  # Purpose:
  # Calculates density values for a particular set of distributions
  # to be passed into the 'violinPlot' function.
  # Arguments:
  # priors - A vector of parameter values
  # ver    - The type of prior distribution to use
  # v      - The sequence of values over which calculate the likelihood
  # Returns:
  # The x and y-axis values to use for plotting.
  
  if ( ver == 1 ) {
    out = list( x = v, y = dnorm( v, priors[1], priors[2] ) )
  }
  if ( ver == 2 ) {
    out = list( x = v, y = dgamma( v, priors[1], priors[2] ) )
  }
  if ( ver == 3 ) {
    if ( max(v) > 1 | min(v) < 0 ) v = seq(0,1,length=1000)
    out = list( x = v, y = dbeta( v, priors[1], priors[2] ) )
  }
  
  return( out )
}
