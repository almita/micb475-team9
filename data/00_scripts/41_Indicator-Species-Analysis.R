#Indicator Species Analysis

#Load packages and import data:
library(phyloseq)
library(indicspecies)

starch_filtered <- readRDS("starch_filtered_phyloseq.RDS")

#Transform phyloseq object to relative abundance
starch_relative_abundance <- transform_sample_counts(starch_filtered, fun=function(x) x/sum(x))

#Run ISA using group variable
starch_isa_output <- multipatt(t(otu_table(starch_relative_abundance)), cluster = sample_data(starch_relative_abundance)$group)

summary(starch_isa_output)

#Extract the taxonomy table
starch_taxtable <- tax_table(starch_relative_abundance) %>% as.data.frame() %>% rownames_to_column(var="ASV")

#Merge taxonomy table with phyloseq object and filter by significant p-value
starch_res <- starch_isa_output$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(starch_taxtable) %>%
  filter(p.value<0.05)

#View results
view(starch_res)
