## PICRUST ##

# setup -------------------------------------------------------------------
# Load all necessary libraries
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(DESeq2)
library(ggh4x)

# import data
abundance_file <- "../60_picrust2/pathway_abundance.tsv"
ko_file <- "../60_picrust2/kometa.tsv"

abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE, skip = 1) 
ko_data <- read_delim(ko_file, delim = "\t", col_names = TRUE, trim_ws = TRUE, skip = 1)
metadata <- read_delim("../10_metadata-and-manifest/metadata.tsv") 

# change sample-id name
metadata <- metadata %>% as_tibble() %>% select(sample_name = `sample-id`, trt, gender, cohort, group)

# change ko to kegg terms
kegg <- ko2kegg_abundance(data = ko_data)

# create function to count the number of non-zero abundances for each feature
nonzero_count <- function(x){return(sum(x>0))}

# keep only features that have more than 3 non-zero abundances
kegg_highfreq <- kegg %>%
	mutate(non_zero_count = apply(kegg, 1, nonzero_count)) %>%
	filter(non_zero_count > 3) %>%
	select(-non_zero_count)

# by group analysis ------------------------------------------------------------------
# differential abundance analysis using LinDa method with Bonferroni p-value adjustment
daa <- pathway_daa(abundance = kegg_highfreq, 
										 metadata = metadata, 
										 group = "group", 
										 daa_method = "LinDA",
										 p.adjust = "bonferroni",
										 reference = "human-C"
										 )

# filter to keep only significantly different pathways
signif_daa <- daa %>% filter(p_adjust < 0.001)

rows = nrow(signif_daa)
half1 = ceiling(rows / 2)
half2 = rows - half1

# divide signif_daa into two to run pathway_annotation function
path_annot1 <- signif_daa %>% 
	slice_head(n = half1) %>% 
	pathway_annotation(daa_results_df = .,
																 pathway = "KO",
																 ko_to_kegg = TRUE)

path_annot2 <- signif_daa %>%
	slice_tail(n = half2) %>%
	pathway_annotation(daa_results_df = .,
										 pathway = "KO",
										 ko_to_kegg = TRUE)

# for some reason this column outputs as a list, so we'll turn it back to a character for merging
path_annot2$pathway_name <- as.character(path_annot2$pathway_name)

# merge annotated tables and remove paths without a name
path_annot <- full_join(path_annot1, path_annot2) %>% filter(!is.na(pathway_name))

# keep only significant ko in abundance table then change rownames to pathway names
annotated_abundance <- kegg_highfreq[rownames(kegg_highfreq) %in% path_annot$feature,] %>%
	rownames_to_column("feature") %>%
	right_join(path_annot %>% select(c("feature","pathway_name")),
						 by = "feature") %>%
	distinct() %>%
	column_to_rownames("pathway_name") %>%
	select(-feature)

# generate heatmap
pathway_heatmap(abundance = annotated_abundance, 
									metadata = metadata, 
									group = "group")

ggsave("../60_picrust2/group_heatmap.png", units = "px", width = 4000, height = 3200)

# generate pca plot
pathway_pca(abundance = annotated_abundance, 
						metadata = metadata, 
						group = "group")

ggsave("../60_picrust2/group_pca.png", units = "px", width = 2400, height = 2000)


# human comparisons -------------------------------------------------------
human_metadata <- metadata %>% filter(cohort == "human")
human_kegg <- kegg_highfreq[,colnames(kegg_highfreq) %in% human_metadata$sample_name]

# differential abundance analysis using LinDa method with Bonferroni p-value adjustment
human_daa <- pathway_daa(abundance = human_kegg, 
									 metadata = human_metadata, 
									 group = "trt", 
									 daa_method = "LinDA",
									 p.adjust = "bonferroni"
)

# filter to keep only significantly different pathways
human_signif_daa <- human_daa %>% filter(p_adjust < 0.05)

## THERE IS NO SIGNIFICANT DIFFERENTIAL ABUNDANCE KO

# mouse comparisons -------------------------------------------------------

mouse_metadata <- metadata %>% filter(cohort == "mouse")
mouse_kegg <- kegg_highfreq[,colnames(kegg_highfreq) %in% mouse_metadata$sample_name]

# differential abundance analysis using LinDa method with Bonferroni p-value adjustment
mouse_daa <- pathway_daa(abundance = mouse_kegg, 
												 metadata = mouse_metadata, 
												 group = "trt", 
												 daa_method = "LinDA",
												 p.adjust = "bonferroni"
)

# filter to keep only significantly different pathways
mouse_signif_daa <- mouse_daa %>% filter(p_adjust < 0.05)

# annotate pathways
mouse_path_annot <- mouse_signif_daa %>% 
	pathway_annotation(daa_results_df = .,
										 pathway = "KO",
										 ko_to_kegg = TRUE)

# keep only significant ko in abundance table then change rownames to pathway names
mouse_annotated_abundance <- mouse_kegg[rownames(mouse_kegg) %in% mouse_path_annot$feature,] %>%
	rownames_to_column("feature") %>%
	right_join(mouse_path_annot %>% select(c("feature","pathway_name")),
						 by = "feature") %>%
	distinct() %>%
	column_to_rownames("pathway_name") %>%
	select(-feature)

pathway_heatmap(abundance = mouse_annotated_abundance, 
								metadata = mouse_metadata, 
								group = "trt")

ggsave("../60_picrust2/mouse_heatmap.png", units = "px", width = 4000, height = 3200)

pathway_pca(abundance = mouse_annotated_abundance, 
						metadata = mouse_metadata, 
						group = "trt")

ggsave("../60_picrust2/mouse_pca.png", units = "px", width = 2400, height = 2000)


# human mouse RS ----------------------------------------------------------
RS_metadata <- metadata %>% filter(trt == "RS")
RS_kegg <- kegg_highfreq[,colnames(kegg_highfreq) %in% RS_metadata$sample_name]

# differential abundance analysis using LinDa method with Bonferroni p-value adjustment
RS_daa <- pathway_daa(abundance = RS_kegg, 
												 metadata = RS_metadata, 
												 group = "cohort", 
												 daa_method = "LinDA",
												 p.adjust = "bonferroni"
					)

# filter to keep only significantly different pathways
RS_signif_daa <- RS_daa %>% filter(p_adjust < 0.01)


rows = nrow(RS_signif_daa)
half1 = ceiling(rows / 2)
half2 = rows - half1

# divide signif_daa into two to run pathway_annotation function
RS_path_annot1 <- RS_signif_daa %>% 
	slice_head(n = half1) %>% 
	pathway_annotation(daa_results_df = .,
										 pathway = "KO",
										 ko_to_kegg = TRUE)

RS_path_annot2 <- RS_signif_daa %>%
	slice_tail(n = half2) %>%
	pathway_annotation(daa_results_df = .,
										 pathway = "KO",
										 ko_to_kegg = TRUE)

# for some reason this column outputs as a list, so we'll turn it back to a character for merging
RS_path_annot2$pathway_name <- as.character(RS_path_annot2$pathway_name)

# merge annotated tables and remove paths without a name
RS_path_annot <- full_join(RS_path_annot1, RS_path_annot2) %>% filter(!is.na(pathway_name))

# keep only significant ko in abundance table then change rownames to pathway names
RS_annotated_abundance <- RS_kegg[rownames(RS_kegg) %in% RS_path_annot$feature,] %>%
	rownames_to_column("feature") %>%
	right_join(RS_path_annot %>% select(c("feature","pathway_name")),
						 by = "feature") %>%
	distinct() %>%
	column_to_rownames("pathway_name") %>%
	select(-feature)

pathway_heatmap(abundance = RS_annotated_abundance, 
								metadata = RS_metadata, 
								group = "cohort")

ggsave("../60_picrust2/RS_heatmap.png", units = "px", width = 4000, height = 4000)


pathway_pca(abundance = RS_annotated_abundance, 
						metadata = RS_metadata, 
						group = "cohort")

ggsave("../60_picrust2/RS_pca.png", units = "px", width = 2400, height = 2000)


# human mouse C -----------------------------------------------------------

C_metadata <- metadata %>% filter(trt == "C")
C_kegg <- kegg_highfreq[,colnames(kegg_highfreq) %in% C_metadata$sample_name]

# differential abundance analysis using LinDa method with Bonferroni p-value adjustment
C_daa <- pathway_daa(abundance = C_kegg, 
											metadata = C_metadata, 
											group = "cohort", 
											daa_method = "LinDA",
											p.adjust = "bonferroni"
)

# filter to keep only significantly different pathways
C_signif_daa <- C_daa %>% filter(p_adjust < 0.001)

# annotate pathways
C_path_annot <- C_signif_daa %>% 
	pathway_annotation(daa_results_df = .,
										 pathway = "KO",
										 ko_to_kegg = TRUE) %>% 
	filter(!is.na(pathway_name))

# keep only significant ko in abundance table then change rownames to pathway names
C_annotated_abundance <- C_kegg[rownames(C_kegg) %in% C_path_annot$feature,] %>%
	rownames_to_column("feature") %>%
	right_join(C_path_annot %>% select(c("feature","pathway_name")),
						 by = "feature") %>%
	distinct() %>%
	column_to_rownames("pathway_name") %>%
	select(-feature)

pathway_heatmap(abundance = C_annotated_abundance, 
								metadata = C_metadata, 
								group = "cohort")

ggsave("../60_picrust2/C_heatmap.png", units = "px", width = 3200, height = 3000)


pathway_pca(abundance = C_annotated_abundance, 
						metadata = C_metadata, 
						group = "cohort")

ggsave("../60_picrust2/C_pca.png", units = "px", width = 2400, height = 2000)

# set group pairs
#groups <- list(c("human-C", "human-RS"), c("mouse-C", "mouse-RS"), c("human-C", "mouse-C"), c("human-RS", "mouse-RS"))

#human <- metadata %>% filter(group %in% group1)


