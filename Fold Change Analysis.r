## Fold Change Calculations for Positive Control Samples ##
rm(list=ls())

# Load RQdeltaCT 
library(RQdeltaCT)
library(dplyr)
library(tidyr)
library(ggsignif)

# Make df for DE + HG genes (FoxA2, VIM, CDX2, SOX) 
setwd("C:/Users/Ana Mogrovejo/Desktop/Thesis/qPCR/Trial 1")

plate <- read_Ct_long(
  path = "plate.1.csv",
  sep = ";",
  dec = ",",
  skip = 1, 
  column.Sample = 3, 
  column.Gene = 2, 
  column.Ct = 5, 
  column.Group = 4, 
)

# Format data correctly 
plate$Ct < as.numeric(plate$Ct) # Convert Ct values to numerical data 
plate$Ct[is.na(plate$Ct)] <- 48 # Convert NA to 48 to indicate low exp. 

# Average Ct values (combine technical replicates) 
collapsed <- make_Ct_ready(
  plate,
  imput.by.mean.within.groups = FALSE,
)

# Calculate Delta Ct 
delta.Ct <- delta_Ct(
  collapsed,
  ref = "HPRT",
  normalise = TRUE,
  transform = FALSE,
  save.to.txt = FALSE
)

# Definitive endoderm fold change
ddCt.DE <- RQ_ddCt(delta.Ct,
                   group.study = "Definitive Endoderm",
                   group.ref = "Control (H9 cells)",
                   do.tests = TRUE,
                   pairwise = FALSE,
                   alternative = "two.sided",
                   p.adjust.method = "BH")

ddCt.DE$log10FCh <- format(ddCt.DE$log10FCh, scientific = FALSE)

# Hindgut fold change
ddCt.HG <- RQ_ddCt(
  delta.Ct,
  group.study = "Hindgut",
  group.ref = "Control (H9 cells)",
  do.tests = TRUE,
  pairwise = FALSE,
  alternative = "two.sided",
  p.adjust.method = "BH",
)

ddCt.HG$log10FCh <- format(ddCt.HG$log10FCh, scientific = FALSE)  

## Intestinal organoid fold change ##

# Make df for positive control IO + EB for experimental genes  
IO.plate <- read_Ct_long(
  path = "IO.plate.csv",
  sep = ",",
  dec = ".",
  skip = 1, 
  column.Sample = 3, 
  column.Gene = 2, 
  column.Ct = 5, 
  column.Group = 4, 
)
IO <- rbind(plate23, IO.plate) # Combine IO qPCR and EB files  

# Format data 
IO$Ct < as.numeric(IO$Ct) # Convert Ct values to numerical data 
IO$Ct[is.na(IO$Ct)] <- 48 # Convert NA to 48 to indicate low exp.
IO$Gene <- trimws(IO$Gene)

# Average Ct values (combine technical replicates) 
IO.collapsed <- make_Ct_ready(
  IO,
  imput.by.mean.within.groups = FALSE,
)

IO.collapsed$MUC2 <- NULL

# Calculate Delta Ct 
IO.delta.Ct <- delta_Ct(
  IO.collapsed,
  ref = "HPRT",
  normalise = TRUE,
  transform = FALSE,
  save.to.txt = FALSE
)

# Intestinal organoid fold change
ddCt.IO <- RQ_ddCt(IO.delta.Ct,
                   group.study = "IO",
                   group.ref = "Control (H9 cells)",
                   do.tests = TRUE,
                   pairwise = FALSE,
                   alternative = "two.sided",
                   p.adjust.method = "BH")

ddCt.IO$log10FCh <- format(ddCt.IO$log10FCh, scientific = FALSE)

