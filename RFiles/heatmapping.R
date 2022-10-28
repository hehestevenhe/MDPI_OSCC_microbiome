#Load packages and data
library(tidyverse)
library(gplots)

raw_data <- read_csv("output.csv")

pepcut_data <- raw_data %>% filter(max_peptides>=2)

## Create reduced dataframe
reduced_data <- select(pepcut_data,-c(2:16))

## Create dataframes with obs cutoffs
twentyfour_obs <- reduced_data %>% filter(tumour_obs>=24) %>% select(1:73) %>% column_to_rownames(var ="accession")

#Colour pallettes
log_colour <- colorRampPalette(c("blue","grey", "red"))(n=300)
sample.types <- c(rep("Darkblue",23),rep("Darkred",72))

#Make individual patient matrix
twentyfour_obs_patient <- data.matrix(twentyfour_obs)
twentyfour_obs_patient_log <- log10(twentyfour_obs_patient)
twentyfour_obs_patient_log[twentyfour_obs_patient_log==-Inf] <- 0
#Make heatmap
png(filename = "log10intensity_individual_peptide_heatmap.png", width = 3200, height = 3200, res = 180)
heatmap.2(twentyfour_obs_patient_log, col=log_colour, na.color = "black", colCol = sample.types ,main = "Log10(protein TMT intensity)\n in patient tissue \n twentyfour_obs", 
          trace = "none", margins = c(7,12) ,cexRow = 0.8, cexCol = 0.8)
dev.off()
