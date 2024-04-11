#Indicator Species Analysis

#Load packages and import data:
library(phyloseq)
library(indicspecies)
library(tidyverse)

# function for extracting phyloseq objects into table
dephyloseq <- function(phylo_obj){ 
	## get the metadata
	meta = as.data.frame(as.matrix(phylo_obj@sam_data)) 
	## how many metadata columns you have
	metacols = ncol(meta)+1
	## get out the otu table
	## if your metadta is empty after running this, you need to use 
	## #otu = as.data.frame(t(as.matrix(phylo_obj@otu_table)))
	otu = as.data.frame(t(as.matrix(phylo_obj@otu_table)))
	## merge the metadata and otu table by the rownames (sample ids from 
	## the Illumina sequence)
	mo = merge(meta, otu, by=0) 
	## get out the taxonomy file
	tax = as.data.frame(phylo_obj@tax_table)
	## get the ASV ID out. This the matches the placeholder ASV ID in the OTU 
	## table
	tax = tax %>% rownames_to_column(var="asv_id")
	## pivot longer to be able to match the ASVs in the OTU table to the taxonomy 
	## table
	mo = mo %>% pivot_longer(cols = -c(1:metacols), 
													 names_to = "asv_id", 
													 values_to="asv_abundance") 
	## Join the metadata and otu table with the taoxnomy table
	mot = full_join(mo, tax)
	## Specify the output for the dephyloseq funciton
	output = mot 
}

starch <- readRDS("../21_phyloseq/starch_filtered_phyloseq.RDS")

# separate genus and family
starch_genus <- tax_glom(starch, taxrank = "Genus")
starch_family <- tax_glom(starch, taxrank = "Family")


# genus -------------------------------------------------------------------

#Transform phyloseq object to relative abundance
starch_ra_genus <- transform_sample_counts(starch_genus, fun=function(x) x/sum(x))

#Run ISA using group variable
starch_isa_output_genus <- multipatt(t(otu_table(starch_ra_genus)), cluster = sample_data(starch_ra_genus)$group)

summary(starch_isa_output_genus, indvalcomp=TRUE)

#Extract the taxonomy table
starch_taxtable <- tax_table(starch_ra_genus) %>% as.data.frame() %>% rownames_to_column(var="ASV")

#Merge taxonomy table with phyloseq object and filter by significant p-value
starch_res_genus <- starch_isa_output_genus$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(starch_taxtable) %>%
  filter(p.value<0.05)

#View results
genus_filtered <- starch_res_genus %>%
	filter(stat >= 0.8) 

genusdf <- dephyloseq(starch_genus)

## keep only relevant data
groupdf <- inner_join(genusdf, genus_filtered)

## calcuate the relative abundance of each ASV in each sample
groupdf$relative_abundance <- as.numeric(groupdf$asv_abundance)/as.numeric(groupdf$sample_sums)
## make a presence/abscence variable for each asv in each sample
groupdf$presabs <- if_else(groupdf$relative_abundance == 0, "asb.", "pres.")
## use this to make the 0s white later
groupdf$relative_abundance <- ifelse(groupdf$relative_abundance==0,NA,groupdf$relative_abundance)



## make the bubble plot
groupdf %>% 
	ggplot(aes(x=as.character(Row.names), y=Genus, size=relative_abundance))+
	geom_point()+
	theme_bw() +
	theme(panel.grid = element_blank(),
				axis.text.y = element_text(size = 8, colour = "black", face="bold"), 
				axis.text.x = element_blank(),
				axis.title = element_text(size=15),
				strip.background =element_rect(fill="grey95"),
				strip.text = element_text(color="black", size=12), 
				legend.text=element_text(size=12),
				strip.text.y = element_text(angle = 0, size=8))+ 
	facet_grid(Family~group, space="free", scales="free")+ 
	labs(x="Samples", y = "Genus", size="relative abundance")

ggsave("../40_core-microbiome/indval.png", units = "px", width = 3270, height = 1418)

# family ------------------------------------------------------------------

#Transform phyloseq object to relative abundance
starch_ra_family <- transform_sample_counts(starch_family, fun=function(x) x/sum(x))

#Run ISA using group variable
starch_isa_output_family <- multipatt(t(otu_table(starch_ra_family)), cluster = sample_data(starch_ra_family)$group)

summary(starch_isa_output_family)

#Extract the taxonomy table
starch_taxtable <- tax_table(starch_ra_family) %>% as.data.frame() %>% rownames_to_column(var="ASV")

#Merge taxonomy table with phyloseq object and filter by significant p-value
starch_res_family <- starch_isa_output_family$sign %>%
	rownames_to_column(var="ASV") %>%
	left_join(starch_taxtable) %>%
	filter(p.value<0.05)

#View results
view(starch_res_family)
