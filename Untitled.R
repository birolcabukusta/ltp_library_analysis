library(readxl)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggalluvial)


'~/Desktop/GenomeComposition.xlsx' %>% 
   read_xlsx() %>% 
   ggplot()
