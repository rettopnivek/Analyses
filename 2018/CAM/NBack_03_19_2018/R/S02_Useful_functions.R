# Useful functions
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you 
# have any questions or comments
# Last updated 2018-04-30

# Table of contents
# 1) tableContents
# 2) cbbPalette
# 3) sdt_prob
# 4) quick_desc
# 5) tableContents
# 6) sdt_calc_binary
# 7) create_design_mat

# 1)
tableContents = function( string_vector, txtSz = 1.5, 
                          shift = NULL, ... ) {
  # Purpose: 
  # Creates a table of contents at the start of a 
  # PDF file with figures.
  # Arguments: 
  # string_vector - A vector of each line for the table 
  #                 table of contents
  # txtSz         - Text size
  # shift         - Optional shift parameter to page numbers
  # ...           - Additional parameters for the legend
  
  blankPlot()
  
  nLines = length( string_vector )
  
  if ( is.null( shift ) ) shift = rep( 0, nLines )
  
  shift = shift + 1;
  
  all_text = paste( 'Page ', 1:nLines + shift, '.......  ', 
                    string_vector,
                    sep = '' )
  
  title( 'Table of contents', cex = txtSz * 1.2 )
  legend( 0, 1, all_text, bty = 'n', cex = txtSz, ... )
  
}

# 2)
# Colorblind palette
cbbPalette = c( "#000000", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )

# 3)
sdt_prob = function( dp, crt, co ) {
  # Purpose:
  # Computes the probability of indicating the 
  # presence of a target based on whether 
  # a target is present or absent.
  # Arguments:
  # dp  - Discriminability parameter (higher values 
  #       indicate a greater likelihood of hits 
  #       and correct rejections)
  # crt - Criterion parameter (positive values 
  #       indicate a bias against saying the 
  #       target is present)
  # co  - Whether the target is present or not
  # Returns:
  # A vector of probabilities.
  
  # Probability of a hit
  prob_H = 1 - pnorm( crt, dp/2 )
  # Probability of a false alarm
  prob_FA = 1 - pnorm( crt, -dp/2 )
  
  # Probability of indicating the presence of a target
  theta = co * prob_H + (1-co) * prob_FA
  
  return( theta )
}

# 4)
quick_desc = function( x,
                       prob = c(.025,.16,.84,.975) ) {
  # Purpose:
  # A convenience function to quickly 
  # compute a set of descriptive 
  # statistics.
  # Arguments:
  # x    - A vector of values
  # prob - A vector giving the cumulative probabilities 
  #        for the quantile function
  # Returns:
  # A vector with the mode, mean, median, 
  # standard deviation, and the quantiles for 
  # the cumulative probabilities.
  
  out = c(
    findMode(x),
    mean(x),
    median(x), 
    sd(x), 
    quantile(x, prob) )
  names( out ) = c( 'Mode', 'Mean', 'Median', 'SD', 
                  paste( 'Q', prob*100, '%', 
                         sep = '' ) )
  
  return(out)
}

# 5)
tableContents = function( string_vector, txtSz = 1.5, 
                          shift = NULL, ... ) {
  # Purpose: 
  # Creates a table of contents at the start of a 
  # PDF file with figures.
  # Arguments: 
  # string_vector - A vector of each line for the table 
  #                 table of contents
  # txtSz         - Text size
  # shift         - Optional shift parameter to page numbers
  # ...           - Additional parameters for the legend
  
  blankPlot()
  
  nLines = length( string_vector )
  
  if ( is.null( shift ) ) shift = rep( 0, nLines )
  
  shift = shift + 1;
  
  all_text = paste( 'Page ', 1:nLines + shift, '.......  ', 
                    string_vector,
                    sep = '' )
  
  title( 'Table of contents', cex = txtSz * 1.2 )
  legend( 0, 1, all_text, bty = 'n', cex = txtSz, ... )
  
}

# 6)
sdt_calc_binary = function( fH, Np, fFA, Nn, correct = 0,
                            dp = TRUE ) {
  # Purpose:
  # Calculates estimates for d' and criterion values of  
  # equal-variance SDT for binary data using the algebraic method.
  # Arguments:
  # fH       - Frequency of hits
  # Np       - Number of positive trials
  # fFA      - Frequency of false alarms
  # Nn       - Number of negative trials
  # correct  - The type of correction to use, where...
  #            0 = none
  #            1 = The log-linear approach, add .5 to the hits and
  #                false alarm frequencies, then add 1 to the 
  #                total number of trials ( Hautus, 1995 )
  #            2 = The conditional approach, where only proportions 
  #                equal to 0 or 1 are adjusted by .5/N or (N-.5)/N
  #                respectively, where N is the associated number of 
  #                total trials for the given proportion 
  #                ( Macmillan & Kaplan, 1985 )
  # dp       - If TRUE, returns the d' estimate; otherwise, 
  #            returns the criterion estimate
  # Returns:
  # A vector giving the estimates for d' and the criterion.
  
  H = fH/Np; FA = fFA/Nn;
  
  # Apply corrections if indicated
  if ( correct == 1 ) {
    H = (fH+.5)/(Np+1)
    FA = (fFA+.5)/(Nn+1)
  }
  if ( correct == 2 ) {
    if ( H == 0 ) H = .5/Np
    if ( FA == 0 ) FA = .5/Nn
    if ( H == 1 ) H = (Ns-.5)/Np
    if ( FA == 1 ) FA = (Nn-.5)/Nn
  }
  
  # Obtain estimate of criterion
  crt_est = -.5*( qnorm( H ) + qnorm( FA ) )
  
  # Obtain estimate of d'
  dp_est = 2*( qnorm( H ) + crt_est )
  
  if ( dp ) {
    return( dp_est )
  } else {
    return( crt_est )
  }
  
}

# 7)
create_design_mat = function( df, 
                              X = NULL,
                              index = NULL ) {
  # Purpose:
  # Creates indicator variables for experimental 
  # design and adds additional predictors, or 
  # checks whether a design matrix has been 
  # correctly specified.
  # Arguments:
  # df    - The data frame to use when creating 
  #         or checking the design matrices
  # X     - The design matrix to check (Optional)
  # index - The columns of the design matrix 
  #         to check (Optional)
  # Returns:
  # A data frame with the design matrix or a 
  # check of the design matrix over the 
  # variables of interest.
  
  # Create design matrix
  if ( is.null( X ) ) {
    
    # Empty vector
    empty_vec = rep( 0, nrow( df ) )
    
    # Initialize data frame with variables
    X = data.frame(
      Intercept = rep( 1, nrow( df ) ),
      NBack_0 = empty_vec,
      NBack_2 = empty_vec,
      Placebo = empty_vec,
      Drug = empty_vec,
      T1 = empty_vec,
      T2 = empty_vec,
      T3 = empty_vec,
      T1_0 = empty_vec,
      T2_0 = empty_vec,
      T3_0 = empty_vec,
      T1_2 = empty_vec,
      T2_2 = empty_vec,
      T3_2 = empty_vec,
      T1_0_D = empty_vec,
      T2_0_D = empty_vec,
      T3_0_D = empty_vec,
      T1_2_D = empty_vec,
      T2_2_D = empty_vec,
      T3_2_D = empty_vec,
      T1_0_P = empty_vec,
      T2_0_P = empty_vec,
      T3_0_P = empty_vec,
      T1_2_P = empty_vec,
      T2_2_P = empty_vec,
      T3_2_P = empty_vec,
      O = empty_vec
    )
    
    # Task
    sel = df$Task == 'NBack_0'
    X$NBack_0[ sel ] = 1
    X$NBack_2[ !sel ] = 1
    
    # Condition
    sel = df$Condition == 'Placebo'
    X$Placebo[ sel ] = 1
    X$Drug[ !sel ] = 1
    
    # Timepoints
    tm = paste( 'T', 1:3, sep = '' )
    
    for ( i in 1:3 ) {
      
      if ( i == 1 )
        val = paste( tm[i], 'Pre_drug', sep = '_' ) else 
          val = paste( tm[i], 'Post_drug', sep = '_' ) 
        
        sel = df$Timepoints == val
        X[sel,tm[i]] = 1
        
    }
    
    # Timepoints by task and condition
    for ( i in 1:3 ) {
      
      for ( j in 1:2 ) {
        if ( j == 1 ) {
          stub = '0'
          adj = X$NBack_0
        } else {
          stub = '2'
          adj = X$NBack_2
        }
        
        vn = paste( tm[i], stub, sep = '_' )
        X[,vn] = as.numeric( X[,tm[i]] * adj )
        
        vn = paste( vn, 'D', sep = '_' )
        X[,vn] = as.numeric( X[,tm[i]] * adj * X$Drug )
        
        vn = paste( tm[i], stub, sep = '_' )
        vn = paste( vn, 'P', sep = '_' )
        X[,vn] = as.numeric( X[,tm[i]] * adj * X$Placebo )
        
      }
    }
    
    # Visit order
    X$O[ df$Visit == 2 ] = 1
    
    # Include rescaled self-report on high
    X$SR = as.numeric( df$Self_report_on_high )/ 100
    
    tmp = X$SR
    tmp[ tmp == 0 ] = .005
    tmp[ tmp == 1 ] = .995
    X$lSR = logit( tmp )
    
    # Variable names for ROI
    roi = c( 'R_DLPFC',
             'L_DLPFC',
             'MPFC',
             'R_VLPFC',
             'L_VLPFC' )
    
    # Initialize values
    X$R_DLPFC = NA
    X$L_DLPFC = NA
    X$MPFC = NA
    X$R_VLPFC = NA
    X$L_VLPFC = NA
    
    # Mean and standard deviation over all estimates
    BOLD_ds = c(
      m = mean( as.vector( unlist( df[,roi] ) ), na.rm = T ),
      sd = sd( as.vector( unlist( df[,roi] ) ), na.rm = T )
    )
    
    # Standardize BOLD estimates
    BOLD = df[,roi]
    BOLD = BOLD %>% 
      mutate(
        R_DLPFC = ( R_DLPFC - BOLD_ds[1] )/BOLD_ds[2],
        L_DLPFC = ( L_DLPFC - BOLD_ds[1] )/BOLD_ds[2],
        MPFC = ( MPFC - BOLD_ds[1] )/BOLD_ds[2],
        R_VLPFC = ( R_VLPFC - BOLD_ds[1] )/BOLD_ds[2],
        L_VLPFC = ( L_VLPFC - BOLD_ds[1] )/BOLD_ds[2]
      )
    
    # Find subjects with complete pre-post pairs 
    # for the 2-back task
    pre = 'T1_Pre_drug'
    po = c( 'T2_Post_drug', 'T3_Post_drug' )
    tsk = df$Task == 'NBack_2'
    
    # Create difference scores for BOLD estimates
    for ( i in 1:2 ) {
      
      sel0 = df$Timepoints == pre & 
        tsk
      sel1 = df$Timepoints == po[i] & 
        tsk
      s0 = df$Subject[ sel0 ]
      s1 = df$Subject[ sel1 ]
      
      if ( length( s0 ) > length( s1 ) ) {
        subjects = sort( unique( s0[ s0 %in% s1 ] ) )
      } else {
        subjects = sort( unique( s1[ s1 %in% s0 ] ) )
      }
      
      # Loop through subjects
      for ( s in 1:length( subjects ) ) {
        
        sel_pre = df$Task == 'NBack_2' & 
          df$Subject == subjects[s] & 
          df$Timepoints == pre
        sel_po = df$Task == 'NBack_2' & 
          df$Subject == subjects[s] & 
          df$Timepoints == po[i]
        
        if ( sum( sel_pre ) == sum( sel_po ) ) {
          X[sel_po,roi] = 
            BOLD[sel_po,] - BOLD[sel_pre,]
        } else {
          sel_co = df$Condition[sel_pre] %in% df$Condition[sel_po]
          X[sel_po,roi] = 
            BOLD[sel_po,] - BOLD[sel_pre,][sel_co,]
        }
        
      }
      
      X[ df$Task == 'NBack_2' & 
           df$Timepoints == 'T1_Pre_drug', roi ] = 0
      
      
    }
    
    return( X )
    
  } else {
    
    # List of conditions
    cnd = list( df$Timepoints, df$Condition, df$Task )
    
    tmp = aggregate( X, cnd, 
                     function(x) mean(x,na.rm = T ) )
    
    check = data.frame(
      Task = tmp[,3],
      Condition = tmp[,2],
      Timepoints = tmp[,1]
    )
    check$X = tmp[,-(1:3)]
    colnames( check$X ) = colnames(X)
    
    if ( !is.null(index) ) {
      check$X = check$X[,index]
    }
    
    return( check )
  }
  
}

