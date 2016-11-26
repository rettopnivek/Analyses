# Define a function to generate starving values
init_f = function(nChains) {
  
  allOut = c()
  
  for (nc in 1:nChains) {
    
    out = list( )
    
    # Initialize subject-level variables
    out$kappa = array( runif( Ns*Cf,
                              .5, 2 ), dim = 
                         c( input$Ns, input$C[1] ) )
    out$xi = array( runif( input$Ns*input$C[2],
                           .5, 2.5 ), dim = 
                      c( input$Ns, input$C[2] ) )
    out$theta = array( runif( input$Ns*input$C[3],
                              .4, .9 ), dim = 
                         c( input$Ns, input$C[3] ) )
    
    # Initialize group-level variables
    out$kappa_mu = array( runif( input$C[1],
                                 .5, 2 ), dim = 
                            c( input$C[1] ) )
    out$kappa_sig = array( runif( input$C[1],
                                  .25, 1 ), dim = 
                             c( input$C[1] ) )
    out$xi_mu = array( runif( input$C[2],
                              .5, 2.5 ), dim = 
                         c( input$C[2] ) )
    out$xi_sig = array( runif( input$C[2],
                               .25, 1 ), dim = 
                          c( input$C[2] ) )
    out$theta_alpha = array( runif( input$C[2],
                                    6, 9 ), dim = 
                               c( input$C[2] ) )
    out$theta_beta = array( runif( input$C[2],
                                   1, 4 ), dim = 
                              c( input$C[2] ) )
    
    allOut = c( allOut, list( out ) )
  }
  
  return( allOut )
}
