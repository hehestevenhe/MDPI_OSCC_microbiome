# MDPI_OSCC_microbiome
R scripts used for processing and secondary analysis of OSCC MS data.

*_prot_filter.R scripts used the TPP ProteinProphet .tsv outputs as input to perform file concatenation and applying filtering criteria for the identified
proteins.

species_inf_sunburst.R: R script for species inference based on log-fold change cutoffs of identified proteins, and visualisation as sunburst plots.

BH_ttesting_volcano_plot.R: R script for performing Benjamini-Hochberg t-testing to identify sig. differentially abundant proteins between NAT and OSCC samples
The results of this were visualised as volcano plots

heatmapping.R: R script for the production of heatmaps based on the protein log-fold change data
