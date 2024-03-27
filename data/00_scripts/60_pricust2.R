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
metadata <- read_delim("../10_metadata-and-manifest/metadata.tsv") 

# change sample-id name
metadata <- metadata %>% as_tibble() %>% select(sample_name = `sample-id`, trt, gender, cohort, group)

# remove pathways with less than 3 non-zero abundances (otherwise it can cause problems later on)
abundance_data_filtered <-  abundance_data %>% filter(rowSums(. != 0) > 3)

# Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(abundance_data_filtered) <- NULL

# verify samples in metadata match samples in abundance_data
samples <- rownames(t(abundance_data_filtered[,-1])) #Getting a list of the sample names in the newly filtered abundance data
metadata <- metadata[metadata$sample_name %in% samples,] #making sure the filtered metadata only includes these samples

# by group analysis ------------------------------------------------------------------
# differential abundance analysis using LinDa method with Bonferroni p-value adjustment
daa <- pathway_daa(abundance = abundance_data_filtered %>% column_to_rownames("#OTU ID"), 
										 metadata = metadata, 
										 group = "group", 
										 daa_method = "LinDA",
										 p.adjust = "bonferroni",
										 reference = "human-C"
										 )

# filter to keep only significantly different pathways 
notsignif_daa <- daa %>% filter(p_adjust > 0.95)
signif_daa <- daa %>% filter(p_adjust <= 0.0001)

rows = nrow(signif_daa)
half1 = ceiling(rows / 2)
half2 = rows - half1

# divide signif_daa into two to run pathway_annotation function
path_annot1 <- signif_daa %>% 
	slice_head(n = half1) %>% 
	pathway_annotation(daa_results_df = .,
																 pathway = "MetaCyc",
																 ko_to_kegg = FALSE)

path_annot2 <- signif_daa %>%
	slice_tail(n = half2) %>%
	pathway_annotation(daa_results_df = .,
										 pathway = "MetaCyc",
										 ko_to_kegg = FALSE)

# merge annotated tables and remove paths without a name
path_annot <- full_join(path_annot1, path_annot2) #%>% filter(!is.na(pathway_name))

# change the pathway column to description in the abundance table 
annotated_abundance <- abundance_data_filtered %>% filter(`#OTU ID` %in% path_annot$feature) %>%
	dplyr::rename(feature = `#OTU ID`) %>%
	left_join(path_annot %>% select(c("feature","description")),
						by = "feature") %>%
	distinct() %>%
	column_to_rownames("description") %>%
	select(-feature)

# generate heatmap
pathway_heatmap(abundance = annotated_abundance, 
									metadata = metadata, 
									group = "group")

ggsave("../60_picrust2/group_heatmap.png", units = "px", width = 4500, height = 4500)

# generate not significant heatmap
notsignif_annot <- pathway_annotation(daa_results_df = notsignif_daa,
									 pathway = "MetaCyc",
									 ko_to_kegg = FALSE)

notsignif_annotated_abundance <- abundance_data_filtered %>% filter(`#OTU ID` %in% notsignif_annot$feature) %>%
	dplyr::rename(feature = `#OTU ID`) %>%
	left_join(notsignif_annot %>% select(c("feature","description")),
						by = "feature") %>%
	distinct() %>%
	column_to_rownames("description") %>%
	select(-feature)

pathway_heatmap(abundance = notsignif_annotated_abundance, 
								metadata = metadata, 
								group = "group")

ggsave("../60_picrust2/group_heatmap_notsignif.png", units = "px", width = 4500, height = 4500)

# generate pca plot
pathway_pca(abundance = annotated_abundance, 
						metadata = metadata, 
						group = "group")

ggsave("../60_picrust2/group_pca.png", units = "px", width = 2400, height = 2000)


# human comparisons -------------------------------------------------------
human_metadata <- metadata %>% filter(cohort == "human")
human_abund <- abundance_data_filtered[,colnames(abundance_data_filtered) %in% c("#OTU ID",human_metadata$sample_name)]

# differential abundance analysis using LinDa method with Bonferroni p-value adjustment
human_daa <- pathway_daa(abundance = human_abund %>% column_to_rownames("#OTU ID"), 
									 metadata = human_metadata, 
									 group = "trt", 
									 daa_method = "LinDA",
									 p.adjust = "bonferroni"
)

# filter to keep only significantly different pathways 
human_notsignif_daa <- human_daa %>% filter(p_adjust > 0.95)
#human_signif_daa <- human_daa %>% filter(p_adjust < 0.1)
# THERE ARE NO PATHWAYS WITH SIGNIFICANT DIFFERENTIAL ABUNDANCE

human_path_annot <- pathway_annotation(daa_results_df = human_daa,
																			 pathway = "MetaCyc",
																			 ko_to_kegg = FALSE)

# change the pathway column to description for the results 
human_annotated_abundance <- human_abund %>% filter(`#OTU ID` %in% human_path_annot$feature) %>%
	dplyr::rename(feature = `#OTU ID`) %>%
	left_join(human_path_annot %>% select(c("feature","description")),
						by = "feature") %>%
	distinct() %>%
	column_to_rownames("description") %>%
	select(-feature)

# generate heatmap
pathway_heatmap(abundance = human_annotated_abundance %>% filter(rowSums(. != 0) > 0), 
								metadata = human_metadata, 
								group = "trt")

ggsave("../60_picrust2/human_heatmap_all.png", units = "px", width = 4500, height = 4500)

# generate not significant heatmap
human_notsignif_annot <- pathway_annotation(daa_results_df = human_notsignif_daa %>% filter(rowSums(. != 0) > 0),
																			pathway = "MetaCyc",
																			ko_to_kegg = FALSE)

human_notsignif_annotated_abundance <- human_abund %>% filter(`#OTU ID` %in% human_notsignif_annot$feature) %>%
	dplyr::rename(feature = `#OTU ID`) %>%
	left_join(human_notsignif_annot %>% select(c("feature","description")),
						by = "feature") %>%
	distinct() %>%
	column_to_rownames("description") %>%
	select(-feature)

pathway_heatmap(abundance = human_notsignif_annotated_abundance, 
								metadata = human_metadata, 
								group = "trt")

ggsave("../60_picrust2/human_heatmap_notsignif.png", units = "px", width = 4500, height = 4500)

# generate pca plot
pathway_pca(abundance = human_annotated_abundance %>% filter(rowSums(. != 0) > 0), 
						metadata = human_metadata, 
						group = "trt")

ggsave("../60_picrust2/human_pca.png", units = "px", width = 2400, height = 2000)


# mouse comparison --------------------------------------------------------
mouse_metadata <- metadata %>% filter(cohort == "mouse")
mouse_abund <- abundance_data_filtered[ ,colnames(abundance_data_filtered) %in% c("#OTU ID", mouse_metadata$sample_name)]

# differential abundance analysis using LinDa method with Bonferroni p-value adjustment
mouse_daa <- pathway_daa(abundance = mouse_abund %>% column_to_rownames("#OTU ID"), 
												 metadata = metadata, 
												 group = "trt", 
												 daa_method = "LinDA",
												 p.adjust = "bonferroni"
)

# filter to keep only significantly different pathways 
mouse_notsignif_daa <- mouse_daa %>% filter(p_adjust == 1) 
mouse_signif_daa <- mouse_daa %>% filter(p_adjust <= 0.05)

# annotate all mouse pathways
mouse_path_annot <- pathway_annotation(daa_results_df = mouse_daa,
																 pathway = "MetaCyc",
																 ko_to_kegg = FALSE)

mouse_notsignif_annot <- mouse_path_annot %>% filter(feature %in% mouse_notsignif_daa$feature)

# change the pathway column to description for the results 
mouse_annotated_abundance <- mouse_abund %>% filter(`#OTU ID` %in% mouse_path_annot$feature) %>%
	dplyr::rename(feature = `#OTU ID`) %>%
	left_join(mouse_path_annot %>% select(c("feature","description")),
						by = "feature") %>%
	distinct() %>%
	column_to_rownames("description") %>%
	select(-feature)

# generate heatmap
pathway_heatmap(abundance = mouse_annotated_abundance %>% filter(rowSums(. != 0) > 0), 
								metadata = mouse_metadata, 
								group = "trt")

ggsave("../60_picrust2/mouse_heatmap_all.png", units = "px", width = 4500, height = 4500)

# generate most not significant heatmap
pathway_heatmap(abundance = mouse_annotated_abundance[rownames(mouse_annotated_abundance) %in% mouse_notsignif_annot$description, ] %>% filter(rowSums(. != 0) > 0),
								metadata = mouse_metadata, 
								group = "trt")

ggsave("../60_picrust2/mouse_heatmap_notsignifp1.png", units = "px", width = 4500, height = 4500)


# generate pca plot
pathway_pca(abundance = mouse_annotated_abundance %>% filter(rowSums(. != 0) > 0), 
						metadata = mouse_metadata, 
						group = "trt")

ggsave("../60_picrust2/mouse_pca.png", units = "px", width = 2400, height = 2000)

# human mouse RS ----------------------------------------------------------
RS_metadata <- metadata %>% filter(trt == "RS")
RS_abund <- abundance_data_filtered[ ,colnames(abundance_data_filtered) %in% c("#OTU ID", RS_metadata$sample_name)]

# differential abundance analysis using LinDa method with Bonferroni p-value adjustment
RS_daa <- pathway_daa(abundance = RS_abund %>% column_to_rownames("#OTU ID"), 
												 metadata = metadata, 
												 group = "cohort", 
												 daa_method = "LinDA",
												 p.adjust = "bonferroni"
)

# filter to keep only significantly different pathways 
RS_notsignif_daa <- RS_daa %>% filter(p_adjust < 0.0001)
RS_signif_daa <- RS_daa %>% filter(p_adjust < 0.0001)

# divide signif_daa into two to run pathway_annotation function
rows = nrow(RS_signif_daa)
half1 = ceiling(rows / 2)
half2 = rows - half1

RS_path_annot1 <- RS_signif_daa %>% 
	slice_head(n = half1) %>% 
	pathway_annotation(daa_results_df = .,
										 pathway = "MetaCyc",
										 ko_to_kegg = FALSE)

RS_path_annot2 <- RS_signif_daa %>%
	slice_tail(n = half2) %>%
	pathway_annotation(daa_results_df = .,
										 pathway = "MetaCyc",
										 ko_to_kegg = FALSE)

# merge annotated tables and remove paths without a name
RS_path_annot <- full_join(RS_path_annot1, RS_path_annot2) #%>% filter(!is.na(pathway_name))

# change the pathway column to description in the abundance table 
RS_annotated_abundance <- RS_abund %>% filter(`#OTU ID` %in% RS_path_annot$feature) %>%
	dplyr::rename(feature = `#OTU ID`) %>%
	left_join(RS_path_annot %>% select(c("feature","description")),
						by = "feature") %>%
	distinct() %>%
	column_to_rownames("description") %>%
	select(-feature)

# generate heatmap
pathway_heatmap(abundance = RS_annotated_abundance, 
								metadata = RS_metadata, 
								group = "cohort")

ggsave("../60_picrust2/RS_heatmap.png", units = "px", width = 4500, height = 4500)

# generate not significant heatmap
RS_notsignif_annot <- pathway_annotation(daa_results_df = RS_notsignif_daa,
																			pathway = "MetaCyc",
																			ko_to_kegg = FALSE)

RS_notsignif_annotated_abundance <- RS_abund %>% filter(`#OTU ID` %in% RS_notsignif_annot$feature) %>%
	dplyr::rename(feature = `#OTU ID`) %>%
	left_join(RS_notsignif_annot %>% select(c("feature","description")),
						by = "feature") %>%
	distinct() %>%
	column_to_rownames("description") %>%
	select(-feature)

pathway_heatmap(abundance = RS_notsignif_annotated_abundance, 
								metadata = RS_metadata, 
								group = "cohort")

ggsave("../60_picrust2/RS_heatmap_notsignif.png", units = "px", width = 4500, height = 4500)

# generate pca plot
pathway_pca(abundance = RS_annotated_abundance %>% filter(rowSums(. != 0) > 0), 
						metadata = RS_metadata, 
						group = "cohort")

ggsave("../60_picrust2/RS_pca.png", units = "px", width = 2400, height = 2000)


# human mouse C -----------------------------------------------------------
C_metadata <- metadata %>% filter(trt == "C")
C_abund <- abundance_data_filtered[ ,colnames(abundance_data_filtered) %in% c("#OTU ID", C_metadata$sample_name)]

# differential abundance analysis using LinDa method with Bonferroni p-value adjustment
C_daa <- pathway_daa(abundance = C_abund %>% column_to_rownames("#OTU ID"), 
											metadata = C_metadata, 
											group = "cohort", 
											daa_method = "LinDA",
											p.adjust = "bonferroni"
)

# filter to keep only significantly different pathways 
C_notsignif_daa <- C_daa %>% filter(p_adjust == 1)
C_signif_daa <- C_daa %>% filter(p_adjust < 0.0001)

# divide signif_daa into two to run pathway_annotation function
rows = nrow(C_signif_daa)
half1 = ceiling(rows / 2)
half2 = rows - half1

C_path_annot1 <- C_signif_daa %>% 
	slice_head(n = half1) %>% 
	pathway_annotation(daa_results_df = .,
										 pathway = "MetaCyc",
										 ko_to_kegg = FALSE)

C_path_annot2 <- C_signif_daa %>%
	slice_tail(n = half2) %>%
	pathway_annotation(daa_results_df = .,
										 pathway = "MetaCyc",
										 ko_to_kegg = FALSE)

# merge annotated tables and remove paths without a name
C_path_annot <- full_join(C_path_annot1, C_path_annot2)

# change the pathway column to description in the abundance table 
C_annotated_abundance <- C_abund %>% filter(`#OTU ID` %in% C_path_annot$feature) %>%
	dplyr::rename(feature = `#OTU ID`) %>%
	left_join(C_path_annot %>% select(c("feature","description")),
						by = "feature") %>%
	distinct() %>%
	column_to_rownames("description") %>%
	select(-feature)

# generate heatmap
pathway_heatmap(abundance = C_annotated_abundance, 
								metadata = C_metadata, 
								group = "cohort")

ggsave("../60_picrust2/C_heatmap.png", units = "px", width = 4500, height = 4500)

# generate not significant heatmap
C_notsignif_annot <- pathway_annotation(daa_results_df = C_notsignif_daa,
																				 pathway = "MetaCyc",
																				 ko_to_kegg = FALSE)

C_notsignif_annotated_abundance <- C_abund %>% filter(`#OTU ID` %in% C_notsignif_annot$feature) %>%
	dplyr::rename(feature = `#OTU ID`) %>%
	left_join(C_notsignif_annot %>% select(c("feature","description")),
						by = "feature") %>%
	distinct() %>%
	column_to_rownames("description") %>%
	select(-feature)

pathway_heatmap(abundance = C_notsignif_annotated_abundance, 
								metadata = C_metadata, 
								group = "cohort")

ggsave("../60_picrust2/C_heatmap_notsignif.png", units = "px", width = 4500, height = 4500)

# generate pca plot
pathway_pca(abundance = C_annotated_abundance %>% filter(rowSums(. != 0) > 0), 
						metadata = C_metadata, 
						group = "cohort")

ggsave("../60_picrust2/C_pca.png", units = "px", width = 2400, height = 2000)


