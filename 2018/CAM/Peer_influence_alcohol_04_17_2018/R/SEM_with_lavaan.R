# 
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directory if you 
# have any questions or comments
# Last updated 2018-04-23

# Table of contents
# 1) Initial setup
# 

###
### 1) Initial setup
###

# Indicate whether models should be estimated
run_model = c(
  T,
  F,
  F,
  F
)

# Save current working directory
orig_dir = getwd()

# Navigate to project directory
setwd( '..' ); proj_dir = getwd();

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

# Load in data
setwd( 'Data' )
load( 'Peer_influence_alcohol.RData' )
setwd( proj_dir )

# Create standardized scores for drinking and days
DR$Scores$zTD = DR$Scores$Total_drinks
DR$Scores$zTD = ( DR$Scores$zTD - mean( DR$Scores$zTD ) ) / 
  sd( DR$Scores$zTD )
# Quick access for variable names
DR$Variables = colnames( DR$Scores )

###
### 2) Syntax for Lavaan
###

# 2.1) MISS scale
# sel = DR$MISS$subscale_indices$consumer_suggestibility
# sel = DR$MISS$subscale_indices$persuadibility
# sel = DR$MISS$subscale_indices$physiological_suggestibility
# sel = DR$MISS$subscale_indices$physiological_reactivity
# sel = DR$MISS$subscale_indices$peer_conformity
# sel = DR$MISS$subscale_indices$mental_control
# sel = DR$MISS$subscale_indices$unpersuadibility
# paste( colnames( DR$MISS$scores )[sel], collapse = ' + ' )

Suggestibility = "
  
  # MISS
  Consumer_suggest =~ MISS_Q2 + MISS_Q10 + MISS_Q14 + 
                      MISS_Q20 + MISS_Q24 + MISS_Q32 + 
                      MISS_Q45 + MISS_Q51 + MISS_Q57 + 
                      MISS_Q63 + MISS_Q70
          Persuade =~ MISS_Q35 + MISS_Q5 + MISS_Q89 + 
                      MISS_Q69 + MISS_Q47 + MISS_Q62 + 
                      MISS_Q75 + MISS_Q44 + MISS_Q88 + 
                      MISS_Q13 + MISS_Q1 + MISS_Q22 + 
                      MISS_Q82 + MISS_Q76
      Phys_suggest =~ MISS_Q64 + MISS_Q94 + MISS_Q71 + 
                      MISS_Q25 + MISS_Q52 + MISS_Q15 + 
                      MISS_Q58 + MISS_Q66 + MISS_Q33 + 
                      MISS_Q68 + MISS_Q11 + MISS_Q77
        Phys_react =~ MISS_Q12 + MISS_Q40 + MISS_Q85 + 
                      MISS_Q91 + MISS_Q31 + MISS_Q43 + 
                      MISS_Q73 + MISS_Q27 + MISS_Q50 + 
                      MISS_Q17 + MISS_Q21 + MISS_Q3 + 
                      MISS_Q60
      Peer_conform =~ MISS_Q4 + MISS_Q16 + MISS_Q29 + 
                      MISS_Q34 + MISS_Q39 + MISS_Q46 + 
                      MISS_Q53 + MISS_Q59 + MISS_Q65 + 
                      MISS_Q72 + MISS_Q78 + MISS_Q84 + 
                      MISS_Q95 + MISS_Q90
  
  # Higher order factor
  Suggestibility =~ Consumer_suggest + Persuade + 
                    Phys_suggest + Phys_react + 
                    Peer_conform
  
"

# Companion scales for suggestibility (Not included)
b("
    Mental_control =~ MISS_Q4 + MISS_Q16 + MISS_Q29 + 
                      MISS_Q34 + MISS_Q39 + MISS_Q46 + 
                      MISS_Q53 + MISS_Q59 + MISS_Q65 + 
                      MISS_Q72 + MISS_Q78 + MISS_Q84 + 
                      MISS_Q95 + MISS_Q90
          Stubborn =~ MISS_Q7 + MISS_Q9 + MISS_Q19 + 
                      MISS_Q37 + MISS_Q38 + MISS_Q42 + 
                      MISS_Q49 + MISS_Q87 + MISS_Q81 + 
                      MISS_Q86 + MISS_Q93 + MISS_Q26 + 
                      MISS_Q56 + MISS_Q30 + MISS_Q54 + 
                      MISS_Q61
")

# 2.2) UPPS-P
# sel = DR$UPPS$subscale_indices$negative_urgency
# sel = DR$UPPS$subscale_indices$premediation
# sel = DR$UPPS$subscale_indices$perseverance
# sel = DR$UPPS$subscale_indices$sensation_seeking
# sel = DR$UPPS$subscale_indices$positive_urgency
# paste( colnames( DR$UPPS$scores )[sel], collapse = ' + ' )

Impulsivity = "

  # UPPS
  Negative_urgency =~ UPPS_Q2 + UPPS_Q7 + UPPS_Q12 + UPPS_Q17 + 
                      UPPS_Q22 + UPPS_Q29 + UPPS_Q34 + UPPS_Q39 + 
                      UPPS_Q44 + UPPS_Q50 + UPPS_Q53 + UPPS_Q58
      Premediation =~ UPPS_Q1 + UPPS_Q6 + UPPS_Q11 + UPPS_Q16 + 
                      UPPS_Q21 + UPPS_Q28 + UPPS_Q33 + UPPS_Q38 + 
                      UPPS_Q43 + UPPS_Q48 + UPPS_Q55
     Perserverance =~ UPPS_Q4 + UPPS_Q9 + UPPS_Q14 + UPPS_Q19 + 
                      UPPS_Q24 + UPPS_Q27 + UPPS_Q32 + UPPS_Q37 + 
                      UPPS_Q42 + UPPS_Q47
    Sensation_seek =~ UPPS_Q3 + UPPS_Q8 + UPPS_Q13 + UPPS_Q18 + 
                      UPPS_Q23 + UPPS_Q26 + UPPS_Q31 + UPPS_Q36 + 
                      UPPS_Q41 + UPPS_Q46 + UPPS_Q51 + UPPS_Q56
  Positive_urgency =~ UPPS_Q5 + UPPS_Q10 + UPPS_Q15 + UPPS_Q20 + 
                      UPPS_Q25 + UPPS_Q30 + UPPS_Q35 + UPPS_Q40 + 
                      UPPS_Q45 + UPPS_Q49 + UPPS_Q52 + UPPS_Q54 + 
                      UPPS_Q57 + UPPS_Q59
  
  # Higher order factor structure
        Rash_action =~ Negative_urgency + Premediation
  Deficit_deligence =~ Sensation_seek + Positive_urgency
  
"


# 2.3) DMQ
# sel = DR$DMQ$subscale_indices$social_motives
# sel = DR$DMQ$subscale_indices$coping
# sel = DR$DMQ$subscale_indices$enhancement
# sel = DR$DMQ$subscale_indices$conformity
# paste( colnames( DR$DMQ$scores )[sel], collapse = ' + ' )

Drinking_motivations = "
  
  # DMQ
  Social_motives =~ DMQ_Q3 + DMQ_Q5 + DMQ_Q11 + DMQ_Q14 + 
                    DMQ_Q16
          Coping =~ DMQ_Q3 + DMQ_Q5 + DMQ_Q11 + DMQ_Q14 + 
                    DMQ_Q16
     Enhancement =~ DMQ_Q7 + DMQ_Q9 + DMQ_Q10 + DMQ_Q13 + 
                    DMQ_Q18
      Conformity =~ DMQ_Q2 + DMQ_Q8 + DMQ_Q12 + DMQ_Q19 + 
                    DMQ_Q20
  
"


###
### 3) Model 1
###

model_1_syntax = paste(
  Suggestibility,
  Drinking_motivations,
  "
  # Regressions
  zTD ~ Suggestibility
  ", sep = "\n" )

if ( run_model[1] ) {
  
  m1 = sem( model = model_1_syntax, 
            data = DR$Scores,
            ordered = c(
              grep( 'MISS', DR$Variables ),
              grep( 'DMQ', DR$Variables ) ) )
  
}

setwd( orig_dir )