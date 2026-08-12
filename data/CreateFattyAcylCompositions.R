
# LOAD FILES  ####
library(reshape2)
library(stringr)
library(readr)
library(dplyr)


#setwd('~/Desktop/AnalysisLTPlibrary/')
ltp.string <- readxl::read_xlsx('~/Desktop/AnalysisLTPlibrary/_lbp_library_full2.xlsx',1)
path = '~/Desktop/AnalysisLTPlibrary/data/'

fa1 <- read_csv(paste0(path, 'LipidSpeciesConcentration_NB_2021_021.csv'))
fa2 <- read_csv(paste0(path, 'LipidSpeciesConcentration_NB_2021_034.csv'))
fa3 <- read_csv(paste0(path, 'LIpidSpeciesConcentration_NB_2021_042.csv'))


fa1 <- melt(fa1[-1], id.vars = c('SampleID'),variable.name = 'Type', value.name = 'nmolg')
fa2 <- melt(fa2[-1], id.vars = c('SampleID'),variable.name = 'Type', value.name = 'nmolg')
fa3 <- melt(fa3[-1], id.vars = c('SampleID'),variable.name = 'Type', value.name = 'nmolg')

quad1 = rbind(fa1, fa2, fa3)

quad1 = quad1%>% 
   mutate(Extraction = substring(SampleID, 1,1),
          Roster = substring(SampleID, 2,3) %>% as.numeric) %>%
   filter(Extraction != 'M')

## MERGE WITH LTP INFO ####
quad2 = merge(quad1, ltp.string[1:4], by= 'Roster', all.x = T)
#write.csv(unique(quad2$Type), 'txt.txt')

head(quad2)
## fix strings ####
quad3 = quad2%>%
   #filter(str_detect(Type, 'd18')) %>% 
   mutate(Type2 = str_replace(Type, 'Cer.d18.0/','DCER')) %>% 
   mutate(Type2 = str_replace(Type2, "d18:1/","")) %>% 
   mutate(Type2 = str_replace(Type2, "TG","TAG")) %>% 
   mutate(Type2 = str_replace_all(Type2, 'DG','DAG')) %>% 
   mutate(Type2 = str_replace_all(Type2, 'Cer','CER')) %>% 
   mutate(Class = gsub("[^a-zA-Z]", "", Type2)) %>% 
   mutate(Class = str_replace_all(Class, 'TAGFA','TAG')) %>% 
   mutate(Class = str_replace_all(Class, '^FA','FFA')) %>% 
   mutate(Tails = gsub("[^0-9]", "", Type2)) 

   
#View(quad3)
quad3$Type2 %>% unique()
quad3$Class %>% unique()




##fix tail length annotations ####
library(plyr)
quad4 = rbind.fill(
   quad3 %>% 
      filter(str_detect(Class, 'CE')) %>% 
      mutate(Tail_1 = Tails) %>%  
      mutate(Tail_1_length = substr(Tail_1, 1,2),
             Tail_1_bend = substr(Tail_1, 3,3)) %>% 
      mutate(Total_tail_length = Tail_1_length,
             Total_tail_bend = Tail_1_bend) ,
   
   quad3 %>% 
      filter(Class %in% c("FFA", 'PA','LPC','LPE','SM')) %>% 
      mutate(Tail_1 = Tails) %>%  
      mutate(Tail_1_length = substr(Tail_1, 1,2),
             Tail_1_bend = substr(Tail_1, 3,3)) %>% 
      mutate(Total_tail_length = Tail_1_length,
             Total_tail_bend = Tail_1_bend) ,
   
   
   quad3 %>% 
      filter(!str_detect(Class, 'TAG')) %>% 
      filter(!str_detect(Class, 'CE')) %>% 
      filter(!Class %in% c("FFA", 'PA','LPC','LPE','SM')) %>% 
      mutate(Tail_1 = substr(Tails, 1,3),
             Tail_2 = substr(Tails, 4,6)) %>%  
      mutate(Tail_1_length = substr(Tail_1, 1,2) %>% as.numeric(),
             Tail_1_bend = substr(Tail_1, 3,3) %>% as.numeric(),
             Tail_2_length = substr(Tail_2, 1,2) %>% as.numeric(),
             Tail_2_bend = substr(Tail_2, 3,3) %>% as.numeric() ) %>% 
      mutate(Total_tail_length = Tail_1_length + Tail_2_length,
             Total_tail_bend = Tail_1_bend + Tail_2_bend),
   
   quad3 %>% 
      filter(str_detect(Class, 'TAG')) %>% 
      filter(nchar(Tails)==6) %>% 
      
      mutate(Tail_1 = substr(Tails, 4,6)) %>%  
      mutate(Tail_1_length = substr(Tail_1, 1,2),
             Tail_1_bend = substr(Tail_1, 3,3)) %>% 
      mutate(Total_tail_length = substr(Tails, 1,2),
             Total_tail_bend = substr(Tails, 3,3)) ,
   
   quad3 %>% 
      filter(str_detect(Class, 'TAG')) %>% 
      filter(nchar(Tails)==7) %>% 
      
      mutate(Tail_1 = substr(Tails, 5,7)) %>%  
      mutate(Tail_1_length = substr(Tail_1, 1,2),
             Tail_1_bend = substr(Tail_1, 3,3)) %>% 
      mutate(Total_tail_length = substr(Tails, 1,2),
             Total_tail_bend = substr(Tails, 3,4)) 
   
   
)





#View(quad4)

head(quad4)
tail(quad4)
colnames(quad4)

quad5 = quad4 %>% 
   select(
      LBP_Family, 
      Target.Gene.Symbol, 
      Roster, 
      SampleID, 
      Extraction, 
      Type2, 
      Class,
      Tails,
      nmolg, 
      Tail_1,
      Tail_1_length,
      Tail_1_bend,
      Tail_2,
      Tail_2_length,
      Tail_2_bend,
      Total_tail_length,
      Total_tail_bend
   )

View(quad5)
#write.csv(quad5, '~/Desktop/GitProjects/ltp_library_analysis/data/LTP_acyl_concentrations_raw1.csv')







## fix tail length annotations ####
# indentTails = function(x) {
#    if(nchar(x)==7)
#       return(paste0(substr(x,1,2), ":",
#                     substr(x, 3,4), "-",
#                     substr(x, 5,6), ":",
#                     substr(x, 7,7)))
#    else if(nchar(x) == 6)
#       return(paste0(substr(x,1,2), ":",
#                     substr(x, 3,3), "-",
#                     substr(x, 4,5),":",
#                     substr(x, 6,6)))
#    else if(nchar(x) == 3)
#       return(paste0(substr(x,1,2), ":",substr(x, 3,3)))
# }
# 
# 
# quad2$Tails2 = sapply(quad2[,'Tail'], indentTails)

head(quad3)
colnames(quad3)
quad4 = quad3 %>% 
   select(
      LBP_Family, 
      Target.Gene.Symbol, 
      Roster, 
      SampleID, 
      Extraction, 
      Type2, 
      Class,
      Tails,
      nmolg, 
      Tail_1,
      Tail_1_length,
      Tail_1_bend,
      Tail_2,
      Tail_2_length,
      Tail_2_bend,
      Total_tail_length,
      Total_tail_bend
   )

View(quad4)
#write.csv(quad4, '__git_GitHub/ltp_library_2024/data/LTP_acyl_concentrations_raw1.csv')




