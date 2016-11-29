#--------------------#
# Useful functions   #
# Kevin Potter       #
# Updated 11/26/2016 #
#--------------------#

# Index
# Lookup - 01:  param_est
# Lookup - 02:  p_f
# Lookup - 03:  model_structure_create
# Lookup - 04:  convergence_extract
# Lookup - 05:  find_dec

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
model_structure_create = function( type, cD, Priors, s = NULL,
                                   cod = F ) {
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
  # cod    - A logical value, indicating if the coefficient of drift 
  #          should be freely estimated.
  # Returns:
  # A list of the inputs.
  
  # Restrict to forced-choice only
  if ( type == 1 | type == 2 ) {
    cD = cD[ cD$Ta == 0, ]
  }
  
  # Number of observations by subject
  No = aggregate( rep(1,nrow(cD)),list(cD$S), sum )$x
  No_max = max( No ) # Largest number of observations
  Ns = length( No ) # Number of subjects
  allS = cD$S
  sbj = unique( allS ) # Subject identifiers
  
  if ( cod ) {
    
    # Select type of model to fit
    if ( type == 1) {
      
      # Coefficients
      # kappa       -> SPC, SUC, LPC, LUC
      # xi target   -> STP,LTP,SFP,LFP
      # xi foil     -> STP,LTP,SFP,LFP
      # tau         -> tau (Single residual latency)
      # sigma       -> F (Estimate for foil racer, target fixed to 1)
      # 14 parameters per subject
      
      # Coefficients
      # k     ( 1:4 )
      # xi_T  ( 5:8 )
      # xi_F  ( 9:12 )
      # sigma ( 13 )
      # tau   ( 14 )
      
      # Define the number of thresholds, drift rates, 
      # residual latencies, and coefficients of drift
      Cf = c( 4, 8, 1, 1 )
      
      # Create index for linear algebra function
      Clm = c( 1:4,   # Rows in design matrix for kappa(1)
               5:12,  # Rows in design matrix for xi(1)
               14,    # Rows in design matrix for sigma(1)
               15,    # Rows in design matrix for tau(1)
               16:19, # Rows in design matrix for kappa(0)
               20:27, # Rows in design matrix for xi(0)
               29,    # Rows in design matrix for sigma(0)
               30     # Rows in design matrix for tau(0)
      )
      
      Rws = c( rep( 1, 4 ),   # kappa (1)
               rep( 2, 8 ),   # xi (1)
               rep( 3, 1 ),   # sigma (1)
               rep( 4, 1 ),   # tau (1)
               rep( 5, 4 ),   # kappa (0)
               rep( 6, 8 ),   # xi (0)
               rep( 7, 1 ),   # sigma (0)
               rep( 8, 1 )    # tau (0)
      )
      # Create index for parameter selection
      parSel = c( 1:4, 5:12, 13, 14, 1:4, 5:12, 13, 14 )
      
      index = cbind( Rws, Clm )
      # Fixed values
      fixed = c(1,1)
      index = rbind( index,
                     c(3,13), # sigma (1)
                     c(7,28) # sigma (0)
      )
      rm( Clm, Rws )
      # Dimensions
      index = rbind( index, c(8,length(parSel)) )
      
      # Fastest RTs by subject
      min_RT = aggregate( cD$RT, list( allS ), min )$x
      
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
          
          # Intercept
          Int = rep( 1, nrow( cS ) )
          
          # Position of correct answer
          Co = cS$Co
          # Define covariate indicating when a choice on the 
          # left was primed
          PL = numeric( nrow( cS ) );
          PL[ (cS$PT == 1 & Co == 0) | 
                (cS$PT == 0 & Co == 1) ] = 1
          # Define covariate indicating when a choice on the 
          # right was primed
          PR = numeric( nrow( cS ) );
          PR[ (cS$PT == 1 & Co == 1) | 
                (cS$PT == 0 & Co == 0) ] = 1
          # Priming duration
          PD = cS$PD
          # Priming duration x type
          tmp = covCreate( cbind( cS$PD, cS$PT ) )
          X_x = designCoding( cS$Cnd, Levels = c(2,1,4,3), 
                              type = 'Intercept' )
          
          # Desired output
          # 8 x N matrix
          # Necessary input
          # 8 x row(X) -> parameter matrix
          # row(X) x N -> Design matrix
          
          # Design matrix
          X[ s, , 1:No[s] ] = 
            t( cbind( 
              PR*PD, PR*(1-PD), PL*PD, PL*(1-PD), # kappa (1)
              X_x*Co, X_x*(1-Co), # xi (1)
              Co, 1-Co, # sigma (1)
              Int, # tau (1)
              PL*PD, PL*(1-PD), PR*PD, PR*(1-PD), # kappa (0)
              X_x*(1-Co), X_x*Co, # xi (0)
              Co, 1-Co, # sigma (1)
              Int  # tau (0)
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
        min_RT = matrix( rep( min_RT, each = Cf[4] ), Ns, Cf[4],
                         byrow = T )
        
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
          curCnd = curCnd ) )
        
      } else {
        
        # Extract data
        cS = cD[ allS == sbj[s], ]
        
        # Extract response pairs
        Y = cbind( cS$RT, cS$Ch )
        
        # Define covariates
        
        # Intercept
        Int = rep( 1, nrow( cS ) )
        # Position of correct answer
        Co = cS$Co
        # Define covariate indicating when a choice on the 
        # left was primed
        PL = numeric( nrow( cS ) );
        PL[ (cS$PT == 1 & Co == 0) | 
              (cS$PT == 0 & Co == 1) ] = 1
        # Define covariate indicating when a choice on the 
        # right was primed
        PR = numeric( nrow( cS ) );
        PR[ (cS$PT == 1 & Co == 1) | 
              (cS$PT == 0 & Co == 0) ] = 1
        # Priming duration
        PD = cS$PD
        # Priming duration x type
        tmp = covCreate( cbind( cS$PD, cS$PT ) )
        X_x = designCoding( tmp, Levels = 1:4, 
                            type = 'Intercept' )
        
        # Desired output
        # 8 x N matrix
        # Necessary input
        # 8 x row(X) -> parameter matrix
        # row(X) x N -> Design matrix
        
        # Design matrix
        X = 
          t( cbind( 
            PR*PD, PR*(1-PD), PL*PD, PL*(1-PD), # kappa (1)
            X_x*Co, X_x*(1-Co), # xi (1)
            Co, 1-Co, # sigma (1)
            Int, # tau (1)
            PL*PD, PL*(1-PD), PR*PD, PR*(1-PD), # kappa (0)
            X_x*(1-Co), X_x*Co, # xi (0)
            Co, 1-Co, # sigma (1)
            Int  # tau (0)
          ) )
        
        # Redefine minimum RT
        min_RT = array( rep( min_RT[ sbj[s] ], Cf[4] ), dim = c(Cf[4]) )
        
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
          curCnd = curCnd ) )
        
      }
      
    }
    
    # Select type of model to fit
    if ( type == 2) {
      
      # Coefficients
      # kappa       -> SL,SR,LL,LR
      # xi target   -> GM, zTL
      # xi foil     -> GM, zFL
      # tau         -> tau (Single residual latency)
      # sigma       -> F (Estimate for foil racer, target fixed to 1)
      # 10 parameters per subject
      
      # Coefficients
      # k     ( 1:4 )
      # xi_T  ( 5:6 )
      # xi_F  ( 7:8 )
      # sigma ( 9 )
      # tau   ( 10 )
      
      # Define the number of thresholds, drift rates, 
      # residual latencies, and coefficients of drift
      Cf = c( 4, 4, 1, 1 )
      
      # Create index for linear algebra function
      Clm = c( 1:4,   # Rows in design matrix for kappa(1)
               5:8,   # Rows in design matrix for xi(1)
               10,    # Rows in design matrix for sigma(1)
               11,    # Rows in design matrix for tau(1)
               12:15, # Rows in design matrix for kappa(0)
               16:19, # Rows in design matrix for xi(0)
               21,    # Rows in design matrix for sigma(0)
               22     # Rows in design matrix for tau(0)
      )
      
      Rws = c( rep( 1, 4 ),   # kappa (1)
               rep( 2, 4 ),   # xi (1)
               rep( 3, 1 ),   # sigma (1)
               rep( 4, 1 ),   # tau (1)
               rep( 5, 4 ),   # kappa (0)
               rep( 6, 4 ),   # xi (0)
               rep( 7, 1 ),   # sigma (0)
               rep( 8, 1 )    # tau (0)
      )
      # Create index for parameter selection
      parSel = c( 1:4, 5:8, 9, 10, 1:4, 5:8, 9, 10 )
      
      index = cbind( Rws, Clm )
      # Fixed values
      fixed = c(1,1)
      index = rbind( index,
                     c(3,9), # sigma (1)
                     c(7,20) # sigma (0)
      )
      rm( Clm, Rws )
      # Dimensions
      index = rbind( index, c(8,length(parSel)) )
      
      # Fastest RTs by subject
      min_RT = aggregate( cD$RT, list( allS ), min )$x
      
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
          
          # Intercept
          Int = rep( 1, nrow( cS ) )
          # Position of correct answer
          Co = cS$Co
          # Define covariate indicating when a choice on the 
          # left was primed
          PL = numeric( nrow( cS ) );
          PL[ (cS$PT == 1 & Co == 0) | 
                (cS$PT == 0 & Co == 1) ] = 1
          # Define covariate indicating when a choice on the 
          # right was primed
          PR = numeric( nrow( cS ) );
          PR[ (cS$PT == 1 & Co == 1) | 
                (cS$PT == 0 & Co == 0) ] = 1
          # Priming duration
          PD = cS$PD
          # Priming duration x type
          X_x = cbind( Int, cS$zTL, Int, cS$zFL )
          
          # Desired output
          # 8 x N matrix
          # Necessary input
          # 8 x row(X) -> parameter matrix
          # row(X) x N -> Design matrix
          
          # Design matrix
          X[ s, , 1:No[s] ] = 
            t( cbind( 
              PR*PD, PR*(1-PD), PL*PD, PL*(1-PD), # kappa (1)
              X_x[,1:2]*Co, X_x[,3:4]*(1-Co), # xi (1)
              Co, 1-Co, # sigma (1)
              Int, # tau (1)
              PL*PD, PL*(1-PD), PR*PD, PR*(1-PD), # kappa (0)
              X_x[,1:2]*(1-Co), X_x[,3:4]*Co, # xi (0)
              Co, 1-Co, # sigma (0)
              Int  # tau (0)
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
        min_RT = matrix( rep( min_RT, each = Cf[4] ), Ns, Cf[4],
                         byrow = T )
        
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
          curCnd = curCnd ) )
        
      } else {
        
        # Extract data
        cS = cD[ allS == sbj[s], ]
        
        # Extract response pairs
        Y = cbind( cS$RT, cS$Ch )
        
        # Define covariates
        
        # Intercept
        Int = rep( 1, nrow( cS ) )
        # Position of correct answer
        Co = cS$Co
        # Define covariate indicating when a choice on the 
        # left was primed
        PL = numeric( nrow( cS ) );
        PL[ (cS$PT == 1 & Co == 0) | 
              (cS$PT == 0 & Co == 1) ] = 1
        # Define covariate indicating when a choice on the 
        # right was primed
        PR = numeric( nrow( cS ) );
        PR[ (cS$PT == 1 & Co == 1) | 
              (cS$PT == 0 & Co == 0) ] = 1
        # Priming duration
        PD = cS$PD
        # Priming duration x type
        X_x = cbind( Int, cS$zTL, Int, cS$zFL )
        
        # Desired output
        # 8 x N matrix
        # Necessary input
        # 8 x row(X) -> parameter matrix
        # row(X) x N -> Design matrix
        
        # Design matrix
        X = 
          t( cbind( 
            PR*PD, PR*(1-PD), PL*PD, PL*(1-PD), # kappa (1)
            X_x[,1:2]*Co, X_x[,3:4]*(1-Co), # xi (1)
            Co, 1-Co, # sigma (1)
            Int, # tau (1)
            PL*PD, PL*(1-PD), PR*PD, PR*(1-PD), # kappa (0)
            X_x[,1:2]*(1-Co), X_x[,3:4]*Co, # xi (0)
            Co, 1-Co, # sigma (0)
            Int  # tau (0)
          ) )
        
        # Redefine minimum RT
        min_RT = array( rep( min_RT[ sbj[s] ], Cf[4] ), dim = c(Cf[4]) )
        
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
          curCnd = curCnd ) )
        
      }
      
    }
    
  } else {
    
    # Select type of model to fit
    if ( type == 1) {
      
      # Coefficients
      # kappa       -> SPC, SUC, LPC, LUC
      # xi target   -> STP,LTP,SFP,LFP
      # xi foil     -> STP,LTP,SFP,LFP
      # tau         -> tau (Single residual latency)
      # 13 parameters per subject
      
      # Coefficients
      # k     ( 1:4 )
      # xi_T  ( 5:8 )
      # xi_F  ( 9:12 )
      # tau   ( 13 )
      
      # Define the number of thresholds, drift rates, 
      # and residual latencies
      Cf = c( 4, 8, 1 )
      
      # Create index for linear algebra function
      Clm = c( 1:4,   # Rows in design matrix for kappa(1)
               5:12,  # Rows in design matrix for xi(1)
               14,    # Rows in design matrix for tau(1)
               15:18, # Rows in design matrix for kappa(0)
               19:26, # Rows in design matrix for xi(0)
               28     # Rows in design matrix for tau(0)
      )
      
      Rws = c( rep( 1, 4 ),   # kappa (1)
               rep( 2, 8 ),   # xi (1)
               rep( 4, 1 ),   # tau (1)
               rep( 5, 4 ),   # kappa (0)
               rep( 6, 8 ),   # xi (0)
               rep( 8, 1 )    # tau (0)
      )
      # Create index for parameter selection
      parSel = c( 1:4, 5:12, 13, 1:4, 5:12, 13 )
      
      index = cbind( Rws, Clm )
      # Fixed values
      fixed = c(1,1)
      index = rbind( index,
                     c(3,13), # sigma (1)
                     c(7,27) # sigma (0)
      )
      rm( Clm, Rws )
      # Dimensions
      index = rbind( index, c(8,length(parSel)) )
      
      # Fastest RTs by subject
      min_RT = aggregate( cD$RT, list( allS ), min )$x
      
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
          
          # Intercept
          Int = rep( 1, nrow( cS ) )
          
          # Position of correct answer
          Co = cS$Co
          
          # Define covariate indicating when a choice on the 
          # left was primed
          PL = numeric( nrow( cS ) );
          PL[ (cS$PT == 1 & Co == 0) | 
                (cS$PT == 0 & Co == 1) ] = 1
          # Define covariate indicating when a choice on the 
          # right was primed
          PR = numeric( nrow( cS ) );
          PR[ (cS$PT == 1 & Co == 1) | 
                (cS$PT == 0 & Co == 0) ] = 1
          
          # Priming duration
          PD = cS$PD
          
          # Covariate for prime duration x type
          dxt = rep( 1, nrow(cS) )
          dxt[ cS$PT == 1 & cS$PD == 1 ] = 2;
          dxt[ cS$PT == 0 & cS$PD == 0 ] = 3;
          dxt[ cS$PT == 0 & cS$PD == 1 ] = 4;
          
          # Priming duration x type
          X_x = designCoding( dxt, Levels = 1:4, 
                              type = 'Intercept' )
          
          # Desired output
          # 8 x N matrix
          # Necessary input
          # 8 x row(X) -> parameter matrix
          # row(X) x N -> Design matrix
          
          # Design matrix
          X[ s, , 1:No[s] ] = 
            t( cbind( 
              PR*PD, PR*(1-PD), PL*PD, PL*(1-PD), # kappa (1)
              X_x*Co, X_x*(1-Co), # xi (1)
              Int, # sigma (1)
              Int, # tau (1)
              PL*PD, PL*(1-PD), PR*PD, PR*(1-PD), # kappa (0)
              X_x*(1-Co), X_x*Co, # xi (0)
              Int, # sigma (0)
              Int  # tau (0)
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
        min_RT = matrix( rep( min_RT, each = Cf[3] ), Ns, Cf[3],
                         byrow = T )
        
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
          curCnd = curCnd ) )
        
      } else {
        
        # Extract data
        cS = cD[ allS == sbj[s], ]
        
        # Extract response pairs
        Y = cbind( cS$RT, cS$Ch )
        
        # Define covariates
        
        # Intercept
        Int = rep( 1, nrow( cS ) )
        
        # Position of correct answer
        Co = cS$Co
        
        # Define covariate indicating when a choice on the 
        # left was primed
        PL = numeric( nrow( cS ) );
        PL[ (cS$PT == 1 & Co == 0) | 
              (cS$PT == 0 & Co == 1) ] = 1
        # Define covariate indicating when a choice on the 
        # right was primed
        PR = numeric( nrow( cS ) );
        PR[ (cS$PT == 1 & Co == 1) | 
              (cS$PT == 0 & Co == 0) ] = 1
        
        # Priming duration
        PD = cS$PD
        
        # Covariate for prime duration x type
        dxt = rep( 1, nrow(cS) )
        dxt[ cS$PT == 1 & cS$PD == 1 ] = 2;
        dxt[ cS$PT == 0 & cS$PD == 0 ] = 3;
        dxt[ cS$PT == 0 & cS$PD == 1 ] = 4;
        
        # Priming duration x type
        X_x = designCoding( dxt, Levels = 1:4, 
                            type = 'Intercept' )
        
        # Desired output
        # 8 x N matrix
        # Necessary input
        # 8 x row(X) -> parameter matrix
        # row(X) x N -> Design matrix
        
        # Design matrix
        X = 
          t( cbind( 
            PR*PD, PR*(1-PD), PL*PD, PL*(1-PD), # kappa (1)
            X_x*Co, X_x*(1-Co), # xi (1)
            Int, # sigma (1)
            Int, # tau (1)
            PL*PD, PL*(1-PD), PR*PD, PR*(1-PD), # kappa (0)
            X_x*(1-Co), X_x*Co, # xi (0)
            Int, # sigma (0)
            Int  # tau (0)
          ) )
        
        # Redefine minimum RT
        min_RT = array( rep( min_RT[ sbj[s] ], Cf[3] ), dim = c(Cf[3]) )
        
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
          curCnd = curCnd ) )
        
      }
      
    }
    
    
  }
  
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
