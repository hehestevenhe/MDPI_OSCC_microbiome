# Load packages and data
library(tidyverse)
library(ggrepel)

raw_data <- read_csv(file="output.csv")

pepcut_data <- raw_data %>% filter(max_peptides>=2)

# Modification of data into an amenable format
col_data <-  select(pepcut_data, -c(1,3:16,90:93))
mod_data <- col_data %>% relocate(healthy_obs,.before = "C3N-03620_H")
long_data <- gather(mod_data,patient,libra_intensity,"C3N-03620_H":"C3N-03782_TO")
long_data <- filter(long_data, healthy_obs>=2)
test_data <- long_data %>% 
  mutate(condition=patient) %>%
  mutate(condition=str_replace(condition, ".*_H","NAT")) %>%
  mutate(condition=str_replace(condition, ".*_T.*","OSCC")) %>%
  relocate(condition,.before = libra_intensity)
test_data[test_data==Inf] <- NA

# Preparation of new table with summary statistics
summaryall <- na.exclude(test_data) %>% 
  group_by(protein, condition) %>%
  summarise(mean=mean(libra_intensity, na.rm=TRUE),
            sd=sd(libra_intensity, na.rm=TRUE),
            len=n())

# Use of loop to perform t.tests for every protein here !!!! THIS FOLD CHANGE IS NOT CORRECT 
protlist <- unique(test_data$protein)

allresult <- data.frame(protein = NULL,    #Create table that will record all results
                        foldchange = NULL,
                        tstat = NULL,
                        pvalue = NULL)
for (i in 1:length(protlist)){
  sub <- test_data %>% filter(protein %in% protlist[i]) %>%
    filter(!is.na(libra_intensity))
  result <- t.test(data = sub, libra_intensity~condition)
  tmp <- data.frame(protein = unique(sub$protein),
                    foldchange = result$estimate[2]/result$estimate[1], #[1] is NAT, [2] is OSCC
                    tstat = result$statistic,
                    pvalue = result$p.value)
  allresult <- rbind(allresult,tmp)
  
  print(i)
}
# Adjust p values
allresult$adj.pvalue <- p.adjust(allresult$pvalue, method = "BH")
allresult <- mutate(allresult,significance=ifelse(allresult$adj.pvalue<0.05,"P<.05", "P>.05"))
allresult <-  mutate(allresult, log2foldchange=log2(allresult$foldchange))
allresult[allresult==-Inf] <- NaN
allresult[allresult==Inf] <- NaN
allresult_graph <- allresult %>% filter(!is.na(log2foldchange))
allresult_graph$accession <- str_extract(allresult_graph$protein, "\\|.*\\|") %>% str_remove_all("\\|")
# Volcano plot
allresult_graph %>%
  ggplot(aes(x=log2foldchange, y=-log10(adj.pvalue))) +
  geom_point(aes(color = significance), size=1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dotted") +
  geom_vline(xintercept = -1, linetype = "dotted") +
  labs(x= expression(paste("Fold change"~~"(",log[2]~frac(tumour,healthy),")")),
       y= expression(-log[10]~adjusted~pvalue),
       title= "Volcano plot of NAT protein fold change compared to OSCC") +
  theme_classic() +
  scale_x_continuous(limits = c(-5,5), breaks = seq(-5, 5, by = 1)) +
  scale_color_manual(values = c("P<.05" = "deeppink3", "P>.05" = "springgreen3")) +
  geom_label_repel(data=filter(allresult_graph, adj.pvalue<=1.904884e-05), mapping = aes(label=accession),size=3,max.overlaps = 20) +
  theme(legend.background = element_rect(fill = 'white' , colour = 'black')) +
  theme(panel.background = element_rect(fill = 'lemonchiffon', colour =)) +
  theme(legend.position = c(.9,.8))