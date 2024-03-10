library(tidyverse)
library(phyloseq)
library(DESeq2)



#### Load data ####
starch <- read_rds("starch_filtered_phyloseq.RDS")



stach_plus1 <- transform_sample_counts(starch, function(x) x+1)
starch_deseq <- phyloseq_to_deseq2(stach_plus1, ~`group`)
starch_DESEQ <- DESeq(starch_deseq)

#comparison between human control and human resistant starch
res_human <- results(starch_DESEQ, tidy=TRUE,
               contrast = c("group", "human-C", "human-RS"))


#comparison between human control and mouse control
res_humanC_mouseC <- results(starch_DESEQ, tidy=TRUE,
                contrast = c("group", "human-C", "mouse-C"))


#comparison between human control and mouse resistant starch
res_humanC_mouseRS <- results(starch_DESEQ, tidy=TRUE,
                           contrast = c("group", "human-C", "mouse-RS"))


#comparison between mouse control and mouse resistant starch
res_mouse <- results(starch_DESEQ, tidy=TRUE,
                           contrast = c("group", "mouse-C", "mouse-RS"))


#comparison between mouse control and human resistant starch
res_mouseC_humanRS <- results(starch_DESEQ, tidy=TRUE,
                     contrast = c("group", "mouse-C", "human-RS"))



#comparison between human resistant starch and mouse resistant starch
res_humanRS_mouseRS <- results(starch_DESEQ, tidy=TRUE,
                               contrast = c("group", "human-RS", "mouse-RS"))

## Volcano plot: effect size VS significance

vol_plot_human <- res_human %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

ggsave("vol_plot_human.png",vol_plot_human)


vol_plot_humanC_mouseC <- res_humanC_mouseC %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

ggsave("vol_plot_humanC_mouseC.png",vol_plot_humanC_mouseC)


vol_plot_humanC_mouseR <- res_humanC_mouseRS %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

ggsave("vol_plot_humanC_mouseR.png",vol_plot_humanC_mouseR)


vol_plot_mouse <- res_mouse %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

ggsave("vol_plot_mouse.png",vol_plot_mouse)


vol_plot_mouseC_humanRS <- res_mouseC_humanRS %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

ggsave("vol_plot_mouseC_humanRS.png",vol_plot_mouseC_humanRS)


vol_plot_humanRS_mouseRS <- res_humanRS_mouseRS %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

ggsave("vol_plot_humanRS_mouseRS.png",vol_plot_humanRS_mouseRS)



# Create bar plot
# To get table of results
sigASVs_human <- res_human %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

sigASVs_humanC_mouseC <- res_humanC_mouseC %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

sigASVs_humanC_mouseRS <- res_humanC_mouseRS %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

sigASVs_mouse <- res_mouse %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

sigASVs_mouseC_humanRS <- res_mouseC_humanRS %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

sigASVs_humanRS_mouseRS <- res_humanRS_mouseRS %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

# Get a vector of ASV names
sigASVs_human_vec <- sigASVs_human %>%
  pull(ASV)

sigASVs_humanC_mouseC_vec <- sigASVs_humanC_mouseC %>%
  pull(ASV)

sigASVs_humanC_mouseRS_vec <- sigASVs_humanC_mouseRS %>%
  pull(ASV)

sigASVs_mouse_vec <- sigASVs_mouse %>%
  pull(ASV)

sigASVs_mouseC_humanRS_vec <- sigASVs_mouseC_humanRS %>%
  pull(ASV)

sigASVs_humanRS_mouseRS_vec <- sigASVs_humanRS_mouseRS %>%
  pull(ASV)

# Prune phyloseq file
human_filt <- prune_taxa(sigASVs_human_vec,starch)
humanC_mouseC_filt <- prune_taxa(sigASVs_humanC_mouseC_vec,starch)
humanC_mouseRS_filt <- prune_taxa(sigASVs_humanC_mouseRS_vec,starch)
mouse_filt <- prune_taxa(sigASVs_mouse_vec,starch)
mouseC_humanRS_filt <- prune_taxa(sigASVs_mouseC_humanRS_vec,starch)
humanRS_mouseRS_filt <- prune_taxa(sigASVs_humanRS_mouseRS_vec,starch)


# Add taxonomy onto DESeq results table
merged_results_human <- tax_table(human_filt) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_human) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))


merged_results_humanC_mouseC <- tax_table(humanC_mouseC_filt) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_humanC_mouseC) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))


merged_results_humanC_mouseRS <- tax_table(humanC_mouseRS_filt) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_humanC_mouseRS) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))


merged_results_mouse <- tax_table(mouse_filt) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_mouse) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))


merged_results_mouseC_humanRS <- tax_table(mouseC_humanRS_filt) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_mouseC_humanRS) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))


merged_results_humanRS_mouseRS <- tax_table(humanRS_mouseRS_filt) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_humanRS_mouseRS) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))


# Make DESeq plot
bar_plot_human <- ggplot(merged_results_human) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

ggsave("bar_plot_human.png", bar_plot_human)


bar_plot_humanC_mouseC <- ggplot(merged_results_humanC_mouseC) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

ggsave("bar_plot_humanC_mouseC.png", bar_plot_humanC_mouseC, width = 20)


bar_plot_humanC_mouseRS <- ggplot(merged_results_humanC_mouseRS) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

ggsave("bar_plot_humanC_mouseRS.png", bar_plot_humanC_mouseRS, width = 20)


bar_plot_mouse <- ggplot(merged_results_mouse) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

ggsave("bar_plot_mouse.png", bar_plot_mouse)


bar_plot_mouseC_humanRS <- ggplot(merged_results_mouseC_humanRS) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

ggsave("bar_plot_mouseC_humanRS.png", bar_plot_mouseC_humanRS, width = 20)

bar_plot_humanRS_mouseRS <- ggplot(merged_results_humanRS_mouseRS) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

ggsave("bar_plot_humanRS_mouseRS.png", bar_plot_humanRS_mouseRS, width = 20)