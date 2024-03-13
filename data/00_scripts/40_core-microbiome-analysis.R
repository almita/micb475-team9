#Core Microbiome Analysis

#Load packages and import data:
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

starch_filtered <- readRDS("starch_filtered_phyloseq.RDS")

#Transform phyloseq object to relative abundance
starch_relative_abundance <- transform_sample_counts(starch_filtered, fun=function(x) x/sum(x))

#Create groups
mouse_control_phyloseq <- subset_samples(starch_relative_abundance, group=="mouse-C")
mouse_RS_phyloseq <- subset_samples(starch_relative_abundance, group=="mouse-RS")
human_control_phyloseq <- subset_samples(starch_relative_abundance, group=="human-C")
human_RS_phyloseq <- subset_samples(starch_relative_abundance, group=="human-RS")

#Core microbiome analysis
mouse_control_core <- core_members(mouse_control_phyloseq, detection=0, prevalence=0.8)
mouse_RS_core <- core_members(mouse_RS_phyloseq, detection=0, prevalence=0.8)
human_control_core <- core_members(human_control_phyloseq,detection=0, prevalence=0.8)
human_RS_core <- core_members(human_RS_phyloseq,detection=0, prevalence=0.8)

#Venn diagrams
##mouse-C vs mouse-RS
mouse_groups_comparison <- ggVennDiagram(x =list(mouse_control_core, mouse_RS_core), set_labels(mouse_groups_comparison, label_list = c("Mouse Control", "Mouse RS"))
)
mouse_groups_comparison <- set_labels(ggVennDiagram(x = list(mouse_control_core, mouse_RS_core)), label_list = c("Mouse Control", "Mouse RS"))                                                     
mouse_groups_comparison
                                                     
##human-C vs human-RS
human_groups_comparison <- ggVennDiagram(x =list(human_control_core, human_RS_core))
human_groups_comparison
##mouse-C vs human-C
control_cohort_comparison <- ggVennDiagram(x = list(mouse_control_core, human_control_core))
control_cohort_comparison
##mouse-RS vs human-RS
RS_cohort_comparison <- ggVennDiagram(x= list(mouse_RS_core, human_RS_core))
RS_cohort_comparison
##4 way comparison
complete_digaram <- ggVennDiagram(x = list(mouse_control_core, human_control_core, mouse_RS_core, human_RS_core))
complete_digaram
