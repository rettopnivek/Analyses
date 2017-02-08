#--------------------#
# Useful functions   #
# Kevin Potter       #
# Updated 02/05/2017 #
#--------------------#

# Lookup - 01
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
