# 
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you 
# have any questions or comments
# Last updated 2019-01-14

# Table of counts
# 1) Initial setup

###
### 1) Initial setup
###

library( shiny )

###
### 2) 
###

# 2.1) 
blm = function( y, X, mu_0, Lambda_0, a_0, b_0 ) {
  # Purpose:
  # Fits a Bayesian linear regression model to a 
  # set of observations given a matrix of predictors 
  # and priors for the regression coefficients and 
  # residual standard deviation.
  # Arguments:
  # y        - A column vector of observed values to fit
  # X        - A design matrix with the predictors
  # mu_0     - A column vector with the means for the 
  #            conditional normal prior on the coefficients
  # Lambda_0 - A precision matrix (inverse of the covariance 
  #            matrix for the coefficients)
  # a_0      - The shape parameter for the inverse gamma 
  #            prior on the residual standard deviation
  #            (equivalent to the shape parameter for 
  #            the gamma prior on the precision)
  # b_0      - The scale parameter for the inverse gamma 
  #            prior on the residual standard deviation 
  #            (equivalent to  1 over the rate or inverse scale 
  #            parameter for the gamma prior on precision)
  # Returns:
  # 
  
  # Ensure dependent variable is a column vector
  y = as.matrix( y )
  # Ensure predictors are in matrix form
  if ( is.null( dim( X ) ) ) {
    X = rbind( X )
  }
  
  # Sample size
  n = nrow( y )
  
  # Check that inputs are correctly specified
  issues = c(
    paste( 'Number of prior means for coefficients', 
           'must match number of predictors' ),
    paste( 'Prior precision matrix must be',
           'positive-definite square matrix' ),
    paste( 'Number of predictors must',
           'match dimensions of prior precision matrix' ),
    paste( 'Shape and scale parameters must be positive' )
  )
  
  checks = c(
    nrow( mu_0 ) == ncol( X ),
    nrow( Lambda_0 ) == ncol( Lambda_0 ) & 
      det( Lambda_0 ) > 0,
    nrow( Lambda_0 ) == ncol( X ),
    a_0 > 0 & b_0 > 0
  )
  
  # If not, return error message
  if ( !all( checks ) ) {
    
    val = sum( !checks )
    if ( val > 1 ) {
      string = paste(
        issues[ !check ][ -val ],
        collapse = ', ' )
      string = paste( string, ', and ', 
                      issues[ !check ][val], sep = '' )
    } else {
      string = issues[ !check ]
    }
    
    stop( string, call. = F )
  }
  
  # Update prior parameters based on current data
  
  # Posterior precision matrix for coefficients
  tXX = t( X ) %*% X
  Lambda_n = tXX + Lambda_0
  # New covariance matrix
  Sigma_n = solve( Lambda_n )
  
  # Posterior means for coefficients
  mu_n = solve( Lambda_n ) %*% ( Lambda_0 %*% mu_0 + t(X) %*% y )
  
  # Posterior shape parameter
  a_n = a_0 + n/2
  # Posterior scale parameter
  b_n = ( t(y) %*% y ) + ( t(mu_0) %*% Lambda_0 %*% mu_0 )
  b_n = b_n - ( t( mu_n ) %*% Lambda_n %*% mu_n )
  b_n = b_0 + .5*b_n
  
  # Create output
  out = list(
    mu_n = mu_n,
    Lambda_n = Lambda_n,
    Sigma_n = Sigma_n,
    a_n = a_n,
    b_n = b_n
  )
  
  return( out )
}

# 2.2) 
comma_sep_to_num = function( string ) {
  # Purpose:
  # ...
  # Arguments:
  # string - 
  # Returns:
  # ...
  
  out = as.numeric( strsplit( string, split = "," )[[1]] )
  
  return( out )
}

# 2.3) 
create_CI_int = function( m, s, alpha ) {
  # Purpose:
  # ...
  # Arguments:
  # m     - 
  # s     - 
  # alpha - 
  # Returns:
  # ...
  
  interval = numeric(2)
  interval[1] = (1-alpha)/2
  interval[2] = interval[1] + alpha
  
  string = paste(
    round( qnorm( interval[1], m, s ), 2 ),
    " - ",
    round( qnorm( interval[2], m, s ), 2 ),
    sep = "" )
  
  return( string )
}

# 2.4) 
create_yhat = function( mu, sigma, X ) {
  # Purpose:
  # ...
  # Arguments:
  # mu    - 
  # sigma - 
  # X     - 
  # Returns:
  # ...
  
  log_yhat = X %*% mu
  yhat = exp( log_yhat + (sigma^2)/2 )
  
  return( list( yhat = yhat[,1], log_yhat = log_yhat[,1] ) )
}

# 2.5) 
compute_ci_or_pi = function( X, param, 
                             alpha = .95 ) {
  
  mu = param[[ grep( 'mu', names( param ) ) ]]
  Lambda = param[[ grep( 'Lambda', names( param ) ) ]]
  a = param[[ max( grep( 'a', names( param ) ) ) ]]
  b = param[[ max( grep( 'b', names( param ) ) ) ]]
  
  sigma = b/( -1 + a )
  sigma = as.vector( sigma )
  
  # Width for credible interval
  interval = numeric(2)
  interval[1] = (1-alpha)/2
  interval[2] = interval[1] + alpha
  
  # Convert to covariance matrix
  Sigma = solve( Lambda )
  
  # Credible interval
  
  
  # mu = X + bY + cZ
  
  # mu = E( X + bY + cZ )
  mu_ci = X %*% mu
  
  # sigma^2 = Var( X + bY + cZ )
  #         = Var(X) + b^2 Var(Y) + c^2 Var(Z) + 
  #           2b Cov(X,Y) + 2c Cov(X,Z) + 2bc Cov( Y, Z )
  sigma_ci = sqrt(
    Sigma[1,1] + 
    pow( X[,2], 2 )*Sigma[2,2] + 
    pow( X[,3], 2 )*Sigma[3,3] + 
    2*X[,2]*Sigma[2,1] + 
    2*X[,3]*Sigma[3,1] + 
    2*X[,2]*X[,3]*Sigma[3,2]
  )
  
  # Compute boundaries
  out = data.frame(
    lLB = qnorm( interval[1], mu_ci, sigma_ci ),
    lUB = qnorm( interval[2], mu_ci, sigma_ci )
  )
  out$LB = exp( out$lLB + pow(sigma,2)/2 )
  out$UB = exp( out$lUB + pow(sigma,2)/2 )
  
  return( out )
}

# 2.6)
pow = function( x, a ) {
  # Purpose:
  # ...
  # Arguments:
  # ...
  # Returns:
  # ...
  
  return( x^a )
}

# 2.7) 
text_for_details = list(
  '$$\\alpha^2$$'
)

###
### 3) User interface
###

# Define UI for app
ui <- fluidPage(
  
  # App title
  titlePanel("CN-THCCOOH Elimination"),
  
  # Vertical layout
  verticalLayout(
    
    # 3.1) Output: Tabset w/ plot, 
    #                        summary, 
    #                        explanation, 
    #                        details, 
    #                        internal settings 
    tabsetPanel( type = "tabs",
                 tabPanel("Explanation",
                          includeMarkdown( "Explanation_for_app.md" ) ), 
                 tabPanel("Plot", plotOutput("plot")),
                 tabPanel("Summary", tableOutput("summary")),
                 tabPanel("Internal settings",
                          fluidPage(
                            textInput( inputId = "prior_mean", 
                                       label = paste( "Priors for",
                                                      "posterior",
                                                      "means",
                                                      "(separated",
                                                      "by commas)" ),
                                       value = "4.62, .61, .3" ),
                            textInput( inputId = "prior_var", 
                                       label = paste( "Priors for",
                                                      "posterior",
                                                      "variances",
                                                      "(separated",
                                                      "by commas)" ),
                                       value = "0.758, .011, .131" ),
                            textInput( inputId = "prior_cor", 
                                       label = paste( "Priors for",
                                                      "posterior",
                                                      "correlations",
                                                      "(separated",
                                                      "by commas)" ),
                                       value = "-0.013, .364, -.008" ),
                            textInput( inputId = "prior_res", 
                                       label = paste( "Priors for",
                                                      "residual SD",
                                                      "(shape and",
                                                      "scale)" ),
                                       value = "146.8, 23.4" )
                            ) ) ),
    
    # 3.2) Panel for inputs
    wellPanel(
      
      
      # 3.2.1) String with CN-THCCOOH values
      textInput( inputId = "CN_THCCOOH", 
                 label = "CN-THCCOOH values (separated by commas)",
                 value = "101, 5" ),
      
      # 3.2.2) String with days at which specimens were collected
      textInput( inputId = "Days_collected", 
                 label = paste( "Days at which specimens were collected",
                                "(separated by commas)" ),
                 value = "0, 14" ),
      
      # 3.2.3) Numeric input for number of days spent using
      numericInput( inputId = "Level_of_use",
                    label = paste( "Number of days spent using",
                                   "during last 30" ),
                    value = 20,
                    min = 0,
                    max = 30 )
      
    )
    
  )
  
)

###
### 4) Server code
###

server <- function(input, output) {
  
  extract_priors = reactive({
    
    # Extract inputs for priors
    tmp = matrix( comma_sep_to_num( input$prior_mean ), 3, 1 )
    priors = list(
      mu_0 = tmp
    )
    V = matrix( 0, 3, 3 )
    diag( V ) = comma_sep_to_num( input$prior_var )
    O = matrix( 1, 3, 3 )
    tmp = comma_sep_to_num( input$prior_cor )
    O[2,1] = tmp[1]
    O[3,1] = tmp[2]
    O[3,2] = tmp[3]
    priors$Sigma_0 = V %*% O %*% V
    priors$Lambda_0 = solve( priors$Sigma_0 )
    priors$a_0 = comma_sep_to_num( input$prior_res )[1]
    priors$b_0 = comma_sep_to_num( input$prior_res )[1]
    
    priors
  })
  
  if ( FALSE ) {
    priors = list(
      mu_0 = matrix( c( 4.62, 0.61, 0.3 ), 3, 1 ),
      Lambda_0 = rbind(
        c(  1.52,  0.13, -1.33 ),
        c(  0.13, 93.03,  0.09 ),
        c( -1.33,  0.09,  8.80 ) ),
      a_0 = 146.8,
      b_0 = 23.4
    )
  }
  
  mLoU = 20.2
  sLoU = 8.21381
  
  # 4.1)
  model_est = reactive({
    
    # Extract priors
    priors = extract_priors()
    
    # Extract inputs
    raw_y = comma_sep_to_num( input$CN_THCCOOH )
    raw_x = comma_sep_to_num( input$Days_collected )
    cv = input$Level_of_use
    
    # Transform values
    y = log( raw_y )
    x = -raw_x
    no_na = !is.na( y )
    y = y[ no_na ]
    x = x[ no_na ]
    zcv = ( cv - mLoU )/sLoU
    
    # Initialize output
    out = list(
      raw_y = raw_y,
      raw_x = raw_x,
      y = y,
      x = x,
      cv = cv,
      zcv = zcv,
      posterior = NULL
    )
    
    if ( length( y ) > 0 ) {
      
      # Create design matrix
      X = matrix( 1, length( y ), 3 )
      X[,2] = zcv
      X[,3] = x
      # Convert dependent variable to column vector
      y = matrix( y, length( y ), 1 )
      
      # Fit Bayesian linear regression model
      post = blm( y, X, 
                  mu_0 = priors$mu_0,
                  Lambda_0 = priors$Lambda_0,
                  a_0 = priors$a_0,
                  b_0 = priors$b_0
      )
      
      out$posterior = post
    }
    
    out
  })
  
  # 4.2) 
  output$plot <- renderPlot({
    
    # Extract priors
    priors = extract_priors()
    
    # Extract model results
    cur_res = model_est()
    
    # Create two plotting panels
    layout( cbind( 1, 2 ) )
    
    # Specify x-axis and 
    # y-axis limits
    yl = c( 0, 400 )
    lyl = c( 0, 6 )
    xl = c( 0, 31 )
    
    # Check if there is data to plot
    check = 
      length( cur_res$raw_x ) > 0 & 
      length( cur_res$raw_y ) > 0 & 
      length( cur_res$raw_x ) == length( cur_res$raw_y )
    
    # Adjust based on current data
    if ( check ) {
      
      # Shift y-axis based on observed CN-THCCOOH
      if ( max( cur_res$raw_y ) > yl[2] ) {
        yl = c( 0, round( max( cur_res$raw_y )*1.1 ) )
        lyl = c( 0, round( log( max( cur_res$raw_y )*1.1 ), 2 ) )
      }
      # Shift x-axis based on final collection day
      if ( max( cur_res$raw_x ) > xl[2] ) {
        xl = c( 0, round( max( cur_res$raw_x )*1.1 ) )
      }
    }
    
    # Create a initial plot
    plot( xl,
          yl,
          bty = 'n',
          type = 'n',
          xlab = 'Specimen collection - days',
          ylab = 'CN-THCCOOH - ng/mg',
          xlim = xl,
          ylim = yl )
    
    # Add legend
    legend( 'topright',
            c( 'Observations',
               'Predictions',
               'Prior' ),
            pch = c(
              19,
              NA,
              NA
            ),
            lty = c(
              NA,
              1,
              2
            ),
            bty = 'n'
    )
            
    
    # Create range of points to plot 
    # estimates over
    rxa = seq( 0, xl[2], length = 100 )
    xa = -rxa
    
    # Create design matrix
    X = matrix( 1, length( rxa ), 3 )
    # Standardized level of use
    X[,2] = cur_res$zcv;
    # Negative of days for elimination rate
    X[,3] = xa
    
    # If able, add information from model fit
    if ( !is.null( cur_res$posterior ) ) {
      
      # Extract posterior
      post = cur_res$posterior
      
      # Posterior mean for elimination
      yhat = create_yhat( post$mu_n,
                          as.vector( post$b_n/( -1 + post$a_n ) ),
                          X )
      # Credible interval around path
      ui_n = compute_ci_or_pi( X, post )
      
      # Add to figure
      polygon( c( rxa, rev( rxa ) ),
               c( ui_n$LB, rev( ui_n$UB ) ),
               col = 'grey70', border = NA )
      lines( rxa, yhat$yhat )
      
    }
    
    # If able, add observed data
    if ( check ) {
      points(
        cur_res$raw_x,
        cur_res$raw_y,
        pch = 19 )
    }
    
    # Add prior information
    yhat_0 = create_yhat( priors$mu_0, 
                          priors$b_0/( -1 + priors$a_0 ), 
                          X )
    ui = compute_ci_or_pi( X, priors )
    
    # Prior credible interval around path
    lines( rxa, ui$LB, lty = 3 )
    lines( rxa, ui$UB, lty = 3 )
    lines( rxa, yhat_0$yhat, lty = 2 )
    
    plot( xl,
          lyl,
          bty = 'n',
          type = 'n', 
          xlab = 'Specimen collection - days',
          ylab = 'CN-THCCOOH - log( ng/mg )',
          xlim = xl,
          ylim = lyl )
    
    # Add legend
    legend( 'topright',
            '95% prediction interval',
            fill = 'grey',
            bty = 'n'
    )
    
    # If able, add information from model fit
    if ( !is.null( cur_res$posterior ) ) {
      
      # Add to figure
      polygon( c( rxa, rev( rxa ) ),
               c( ui_n$lLB, rev( ui_n$lUB ) ),
               col = 'grey70', border = NA )
      lines( rxa, yhat$log_yhat )
      
    }
    
    # If able, add observed data
    if ( check ) {
      points(
        cur_res$raw_x,
        cur_res$y,
        pch = 19 )
    }
    
    # Add prior information
    lines( rxa, ui$lLB, lty = 3 )
    lines( rxa, ui$lUB, lty = 3 )
    lines( rxa, yhat_0$log_yhat, lty = 2 )
    
  }, height = 400, width = 800 )
  
  # 4.3) Generate a summary of the data
  output$summary <- renderTable({
    
    # Extract model results
    cur_res = model_est()
    
    dt = data.frame(
      Variable = c( 
        "Level of log(CN-THCCOOH) at ingestion",
        "Influence of level of use",
        "Elimination rate",
        "Half-life",
        "Window detection" ),
      Status = "",
      Mean = "",
      CI_50 = "",
      CI_80 = "",
      CI_90 = "",
      CI_95 = "",
      CI_99 = "",
      stringsAsFactors = F
    )
    alpha = c( .5, .8, .9, .95, .99 )
    
    dt$Status[1] = "Could not fit model"
    
    if ( !is.null( cur_res$posterior ) ) {
      
      dt$Status[1] = ""
      
      # Extract parameters summarizing posterior distribution
      mu_n = cur_res$posterior$mu_n
      Lambda_n = cur_res$posterior$Lambda_n
      Sigma_n = solve( Lambda_n )
      a_n = cur_res$posterior$a_n
      b_n = cur_res$posterior$b_n
      
      # Regression coefficients
      for ( i in 1:3 ) {
        dt$Mean[i] = as.character( round( mu_n[i,1], 2 ) )
        for ( j in 1:5 ) {
          dt[[j+3]][i] = create_CI_int(
            mu_n[i,1], sqrt( Sigma_n[i,i] ), alpha[j] )
        }
      }
      
      # Sample draws from posterior
      draws = MASS::mvrnorm( 8e+05, 
                             mu_n,
                             Sigma_n )
      
      # Half-life
      
      # Transform to half-life
      hl = log(2)/draws[,3]
      
      # Compute statistics
      
      # Compute median (mean undefined)
      dt$Mean[4] = round( median( hl ), 1 )
      for ( j in 1:5 ) {
        interval = c( (1-alpha[j])/2, alpha[j] + (1-alpha[j])/2 )
        tmp = quantile( hl, prob = interval )
        tmp = paste( round( tmp[1], 1 ), '-', round( tmp[2], 1 ) )
        dt[[j+3]][4] = tmp
      }
      
      # Window of detection
      
      wod = log(5) - draws[,1]/-draws[,3]
      
      # Compute median (mean undefined)
      dt$Mean[5] = round( median( wod ), 1 )
      for ( j in 1:5 ) {
        interval = c( (1-alpha[j])/2, alpha[j] + (1-alpha[j])/2 )
        tmp = quantile( wod, prob = interval )
        tmp = paste( round( tmp[1], 1 ), '-', round( tmp[2], 1 ) )
        dt[[j+3]][5] = tmp
      }
      
      
    }
    
    colnames( dt )[4:(4+4)] = paste( alpha*100, '% CI', sep = '' )
    dt
  })
  
}

# Return a Shiny app object
shinyApp(ui = ui, server = server)

