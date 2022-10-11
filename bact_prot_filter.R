#Load packages and proteinprophet data
library(tidyverse)

setwd("")
raw_data <- list.files(pattern="*.tsv") %>%
  map_df(~read_tsv(.))

#filter data 
filter_data <- raw_data %>%
  filter(!str_detect(protein, "HUMAN")) %>%   #remove human entries
  filter(`protein probability`>= 0.9) %>%  #retain entries with prob>0.9
  filter(!`protein probability`> 1) %>%   
  filter(!str_count(`protein`, "SEQ")>1)   #remove bact entries with multiple protein mappings
filter_data$`num peptides` <- str_count(filter_data$peptides,"\\+")+1
filter_data <- filter_data %>% relocate (`num peptides`, .after = peptides)


#Aggregate 1 concatenates and collapses the libra data
aggregate_data_1 <- aggregate(x=filter_data,by=list(cross_check=filter_data$protein), mean, na.rm=TRUE) %>%
select(., c("cross_check", 
            "C3N-03620_H","C3L-00995_H","C3N-03008_H",
            "C3N-01944_H","C3N-03028_H","C3L-03378_H",
            "C3N-03933_H","C3L-00999_H","C3N-03045_H",
            "C3N-01758_H","C3N-03619_H","C3L-01237_H",
            "C3N-03027_H","C3N-03015_H","C3N-01947_H",
            "C3N-04279_H","C3N-02714_H","C3N-03013_H",
            "C3N-01754_H","C3N-03928_H","C3L-00994_H",
            "C3N-00829_H","C3N-00498_H",
            "C3N-03620_T","C3L-00995_T","C3N-03008_T",
            "C3N-01944_T","C3N-03028_T","C3L-03378_T",
            "C3N-03933_T","C3L-00999_T","C3N-03045_T",
            "C3N-01758_T","C3N-03619_T","C3L-01237_T",
            "C3N-03027_T","C3N-03015_T","C3N-01947_T",
            "C3N-04279_T","C3N-02714_T","C3N-03013_T",
            "C3N-01754_T","C3N-03928_T","C3L-00994_T",
            "C3N-00829_T","C3N-00498_T", "C3N-03042_T*",
            "C3N-03664_TO", "C3N-00825_TO", "C3N-02694_TO", 
            "C3N-01948_TO", "C3N-03458_TO", "C3N-00846_TO", 
            "C3N-03226_TO", "C3N-02695_TO", "C3N-03456_TO",
            "C3L-00977_TO", "C3N-00871_TO", "C3N-03433_TO", 
            "C3N-03889_TO", "C3N-03783_TO", "C3N-03009_TO", 
            "C3N-03457_TO", "C3N-02730_TO", "C3N-03785_TO",
            "C3N-03487_TO", "C3L-00987_TO", "C3L-04025_TO",
            "C3N-02925_TO", "C3N-01752_TO", "C3L-04849_TO",
            "C3N-03782_TO"))

#Aggregate 2 maintains all the other columns and produces means for values such as probability
aggregate_data_2 <- filter_data %>%  group_by(protein) %>% 
  summarise("protein" = paste(`protein`, collapse = "~"),"mean protein probability"= mean(`protein probability`), 
            "protein description" = paste(`protein description`, collapse = "~"),
            "protein length"=paste(`protein length`, collapse = "+"), "mean percent coverage"=mean(`percent coverage`),
            "mean tot indep spectra"=mean(`tot indep spectra`), "mean percent share of spectrum ids"=mean(`percent share of spectrum ids`),
            "peptides"=paste(`peptides`, collapse = "|"), "mean peptides" = mean(`num peptides`), 
            "max_peptides"= max(`num peptides`)) %>%
  mutate_at("protein", str_replace,"\\~.*","") %>%
  mutate_at("protein description", str_replace,"\\~.*","") %>%
  mutate("entries collapsed"= 1+(str_count(`protein length`,"\\+"))) %>%
  mutate_at("protein length", str_replace,"\\+.*","")
  
  
#Merge the two dataframes
aggregate_data <- bind_cols(aggregate_data_2,aggregate_data_1) %>% relocate(cross_check,.after=protein) %>%            
             relocate(any_of(c("C3N-03620_H","C3L-00995_H","C3N-03008_H",
                               "C3N-01944_H","C3N-03028_H","C3L-03378_H",
                               "C3N-03933_H","C3L-00999_H","C3N-03045_H",
                               "C3N-01758_H","C3N-03619_H","C3L-01237_H",
                               "C3N-03027_H","C3N-03015_H","C3N-01947_H",
                               "C3N-04279_H","C3N-02714_H","C3N-03013_H",
                               "C3N-01754_H","C3N-03928_H","C3L-00994_H",
                               "C3N-00829_H","C3N-00498_H",
                               "C3N-03620_T","C3L-00995_T","C3N-03008_T",
                               "C3N-01944_T","C3N-03028_T","C3L-03378_T",
                               "C3N-03933_T","C3L-00999_T","C3N-03045_T",
                               "C3N-01758_T","C3N-03619_T","C3L-01237_T",
                               "C3N-03027_T","C3N-03015_T","C3N-01947_T",
                               "C3N-04279_T","C3N-02714_T","C3N-03013_T",
                               "C3N-01754_T","C3N-03928_T","C3L-00994_T",
                               "C3N-00829_T","C3N-00498_T", "C3N-03042_T*",
                               "C3N-03664_TO", "C3N-00825_TO", "C3N-02694_TO", 
                               "C3N-01948_TO", "C3N-03458_TO", "C3N-00846_TO", 
                               "C3N-03226_TO", "C3N-02695_TO", "C3N-03456_TO",
                               "C3L-00977_TO", "C3N-00871_TO", "C3N-03433_TO", 
                               "C3N-03889_TO", "C3N-03783_TO", "C3N-03009_TO", 
                               "C3N-03457_TO", "C3N-02730_TO", "C3N-03785_TO",
                               "C3N-03487_TO", "C3L-00987_TO", "C3L-04025_TO",
                               "C3N-02925_TO", "C3N-01752_TO", "C3L-04849_TO",
                               "C3N-03782_TO")),.after=`entries collapsed`)
#Make obs count and species columns
mutate_aggregate_data <- aggregate_data 
mutate_aggregate_data[mutate_aggregate_data==Inf] <- NA
mutate_aggregate_data <- mutate_aggregate_data%>%
mutate(species=str_extract(`protein description`, "(?<=HMT-... ).+(?=\\])")) %>%
  mutate(HMT_ID=str_extract(`protein description`, "(?<=HMT-)[[:digit:]]+(?=)")) %>%
  relocate(any_of(c("species", "HMT_ID")), .after = `protein description`) 
mutate_aggregate_data <- mutate_aggregate_data %>%
  mutate(healthy_obs=rowSums(!is.na(mutate_aggregate_data[,c(15:37)]))) %>%
  mutate(tumour_obs=rowSums(!is.na(mutate_aggregate_data[,c(38:86)]))) 

write.csv(mutate_aggregate_data, row.names = FALSE, file = "output.csv")

