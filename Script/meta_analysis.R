####################################################################################
# TITLE #
####################################################################################

##################
# R script to generate the analysis
# Analysis performed with R (v. XXX) and Rstudio (v. XXX)
# Created by Stefano Mammola
# Initiated 16/12/2025
##################

####################################################################################
# Preparation ----------------------------------------------------------------------
####################################################################################

# Loading R packages ------------------------------------------------------

if(!require("pacman")) {install.packages("pacman")}
pacman::p_load("dplyr","ggplot2","metafor","readxl")

# Settings ----------------------------------------------------------------

# Plotting parameters
theme_set(theme_bw())
theme_update(
  legend.background = element_blank(), #No background (legend)
  plot.background = element_blank(), #No background
  panel.grid = element_blank(), #No gridlines
  axis.text  = element_text(size = 10, colour = "grey10"), #Size and color of text
  axis.title = element_text(size = 12, colour = "grey10") #Size and color of text
)

# Load data ---------------------------------------------------------------
db <- read_excel("Data/Database_meta_analysis_bias.xlsx", na = "NA") |>
  dplyr::mutate_if(is.character, as.factor) |> 
  dplyr::arrange(Pearson_r) 
str(db)

####################################################################################
# Data exploration -----------------------------------------------------------------
####################################################################################

#Number of papers
nlevels(db$ID) #54

#Number of estimates
nrow(db) #881

#Number of journals
nlevels(db$Journal) #78

#Summary of estimates in each combination of response and predictor
(table_tot <- table(db$Dependent_zoomed,db$Indipendent_macro))

# Derive Fischer's Z and its variance -------------------------------------
db <- escalc(measure = "ZCOR", ri = Pearson_r, ni = N, data = db)



