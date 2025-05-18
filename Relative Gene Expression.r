## Relative Expression Calculations for Expeirmental Samples ##

rm(list=ls())
dev.off()

# Load packages  
library(RQdeltaCT)
library(dplyr)
library(tidyr)
library(ggsignif)

# Raw Ct Data for plates 2 and 3 
setwd("C:/Users/Ana Mogrovejo/Desktop/Thesis/qPCR/Trial 1")

plate23 <- read_Ct_long(
  path = "plates23.csv",
  sep = ";",
  dec = ",",
  skip = 1, 
  column.Sample = 3, 
  column.Gene = 2, 
  column.Ct = 5, 
  column.Group = 4, 
)

plate23$Ct < as.numeric(plate23$Ct)  # Convert Ct values to numerical data 
plate23$Gene <- trimws(plate23$Gene) # Trim spaces in gene names 


# Average Ct value for samples  
collapsed23 <- make_Ct_ready(
  plate23,
  imput.by.mean.within.groups = FALSE,
  save.to.txt = TRUE,
  name.txt = "collapsed"
)

# dCt + 2^-(Ct)
delta.Ct23 <- delta_Ct(
  collapsed23,
  ref = "HPRT",
  normalise = TRUE,
  transform = TRUE,
  save.to.txt = FALSE
)

delta.Ct23$MUC2 <- NULL # Unsuccesful primers 

# Wide to long format  
delta.Ct.long23 <- delta.Ct23 %>%
  pivot_longer(
    cols = -c(Sample, Group),  # Protect Sample + Group
    names_to = "Gene",
    values_to = "value"
  )


delta.Ct.long23$value <- format(delta.Ct.long23$value, scientific = FALSE)

