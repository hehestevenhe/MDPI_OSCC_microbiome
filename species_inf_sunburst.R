#Load packages
library(tidyverse)
library(plotly)
library(data.table)

### Raw data upload and processing
raw_data <- read.csv(file="output.csv")
colnames(raw_data)[colnames(raw_data)=="species"] <- "tax_name"
trim_data <- raw_data %>% select(-c(24:95))

### Database upload, processing and merge
db <- read_tsv("rankedlineage.dmp") %>% 
  select(seq(1,by=2,len=10))
colnames(db) <- c("tax_id", "tax_name", "species", "genus", "family", "order", "class",
                  "phylum", "kingdom", "superkingdom")
merge_data <- merge(trim_data, db, by = "tax_id")

colnames(merge_data)[colnames(merge_data)=="tax_name.y"] <- "tax_name"
merge_data <- merge_data %>% select(-c("tax_name.x"))
merge_table <- as.data.table(merge_data)
merge_table[,count := .N, by = .(tax_name)]
merge_Libra_mean <- merge_table %>% group_by(tax_name) %>% 
  summarise(mean_log2_FC = mean(log2foldchange))
merge_table <- merge(merge_table, merge_Libra_mean, by = "tax_name")

### Define function
as.sunburstDF <- function(DF, value_column = NULL, add_root = FALSE){
  require(data.table)
  
  colNamesDF <- names(DF)
  
  if(is.data.table(DF)){
    DT <- copy(DF)
  } else {
    DT <- data.table(DF, stringsAsFactors = FALSE)
  }
  
  if(add_root){
    DT[, root := "Total"]  
  }
  
  colNamesDT <- names(DT)
  hierarchy_columns <- setdiff(colNamesDT, value_column)
  DT[, (hierarchy_columns) := lapply(.SD, as.factor), .SDcols = hierarchy_columns]
  
  if(is.null(value_column) && add_root){
    setcolorder(DT, c("root", colNamesDF))
  } else if(!is.null(value_column) && !add_root) {
    setnames(DT, value_column, "values", skip_absent=TRUE)
    setcolorder(DT, c(setdiff(colNamesDF, value_column), "values"))
  } else if(!is.null(value_column) && add_root) {
    setnames(DT, value_column, "values", skip_absent=TRUE)
    setcolorder(DT, c("root", setdiff(colNamesDF, value_column), "values"))
  }
  
  hierarchyList <- list()
  
  for(i in seq_along(hierarchy_columns)){
    current_columns <- colNamesDT[1:i]
    if(is.null(value_column)){
      currentDT <- unique(DT[, ..current_columns][, values := .N, by = current_columns], by = current_columns)
    } else {
      currentDT <- DT[, lapply(.SD, sum, na.rm = TRUE), by=current_columns, .SDcols = "values"]
    }
    setnames(currentDT, length(current_columns), "labels")
    hierarchyList[[i]] <- currentDT
  }
  
  hierarchyDT <- rbindlist(hierarchyList, use.names = TRUE, fill = TRUE)
  
  parent_columns <- setdiff(names(hierarchyDT), c("labels", "values", value_column))
  hierarchyDT[, parents := apply(.SD, 1, function(x){fifelse(all(is.na(x)), yes = NA_character_, no = paste(x[!is.na(x)], sep = ":", collapse = " - "))}), .SDcols = parent_columns]
  hierarchyDT[, ids := apply(.SD, 1, function(x){paste(x[!is.na(x)], collapse = " - ")}), .SDcols = c("parents", "labels")]
  hierarchyDT[, c(parent_columns) := NULL]
  return(hierarchyDT)
}

#### OSCC
OSCC <- merge_table %>% filter(mean_log2_FC>-1)
OSCC_taxa <- OSCC %>% select("phylum","class","order","family","genus","species","tax_name","count")
OSCC_taxa[(is.na(OSCC_taxa))] <- "Unclassified"
# OSCC >=2 proteins for inference
OSCC_multi <- filter(OSCC_taxa,count>=2)
OSCC_multi <- unique(OSCC_multi)
OSCC_multi_sunburst <- as.sunburstDF(OSCC_multi)
plot_ly(data = OSCC_multi_sunburst, ids = ~ids, labels= ~labels, parents = ~parents, 
        values= ~values, type='sunburst', branchvalues = 'total', maxdepth=3,
        textinfo='label+percent root', marker = list(colors = OSCC_multi_sunburst$col))


#### NAT
NAT <- merge_table %>% filter(mean_log2_FC<1)
NAT_taxa <- NAT %>% select("phylum","class","order","family","genus","species","tax_name","count")
NAT_taxa[(is.na(NAT_taxa))] <- "Unclassified"
# NAT >=2 proteins for inference
NAT_multi <- filter(NAT_taxa,count>=2)
NAT_multi <- unique(NAT_multi)
NAT_multi_sunburst <- as.sunburstDF(NAT_multi)
plot_ly(data = NAT_multi_sunburst, ids = ~ids, labels= ~labels, parents = ~parents, 
        values= ~values, type='sunburst', branchvalues = 'total', maxdepth=3,
        textinfo='label+percent root', marker = list(colors = NAT_multi_sunburst$col))

#Venn
library(eulerr)

venn <- filter(merge_table, merge_table$count>=2)
venn <- venn[!duplicated(venn$tax_name),]

sum(venn$mean_log2_FC>=1)
sum(venn$mean_log2_FC<=-1)
sum(venn$mean_log2_FC>=-1 & venn$mean_log2_FC<=1)

fungal_venn <- euler(c("NAT"=4, "OSCC"=70, "NAT&OSCC"=154)) #fungal numbers used here
plot(fungal_venn,
     fills = c("#B3E5FC", "#F7Dc64"),
     labels = list(cex=2),
     quantities = list(cex=2.5))

write.csv(venn, row.names = FALSE, file = "venn_list.csv")
