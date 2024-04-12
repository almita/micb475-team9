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

# by cohort 1 ------------------------------------------------------------------
# differential abundance analysis using LinDa method with Bonferroni p-value adjustment
daa <- pathway_daa(abundance = abundance_data_filtered %>% column_to_rownames("#OTU ID"), 
										 metadata = metadata, 
										 group = "cohort", 
										 daa_method = "LinDA",
										 p.adjust = "bonferroni"
										 )

# filter to keep only significantly different pathways 
notsignif_daa <- daa %>% filter(p_adjust >= 0.05)
signif_daa <- daa %>% filter(p_adjust < 0.05)

# divide signif_daa into two to run pathway_annotation function
all_path_annot <- pathway_annotation(daa_results_df = daa,
																		 pathway = "MetaCyc",
																		 ko_to_kegg = FALSE)

signif_path_annot <- pathway_annotation(daa_results_df = signif_daa,
																 pathway = "MetaCyc",
																 ko_to_kegg = FALSE)

notsignif_path_annot <- pathway_annotation(daa_results_df = notsignif_daa,
										 pathway = "MetaCyc",
										 ko_to_kegg = FALSE)

# create annotated abundances tables
annotated_abundance <- abundance_data_filtered %>% 
	dplyr::rename(feature = `#OTU ID`) %>%
	inner_join(all_path_annot %>% select(c("feature","description"))) %>%
	distinct() %>%
	column_to_rownames("description") %>%
	select(-feature)

annotated_abundance_signif <- abundance_data_filtered %>% 
	dplyr::rename(feature = `#OTU ID`) %>%
	inner_join(signif_path_annot %>% select(c("feature","description"))) %>%
	distinct() %>%
	column_to_rownames("description") %>%
	select(-feature)

annotated_abundance_notsignif <- abundance_data_filtered %>% 
	dplyr::rename(feature = `#OTU ID`) %>%
	inner_join(notsignif_path_annot %>% select(c("feature","description"))) %>%
	distinct() %>%
	column_to_rownames("description") %>%
	select(-feature)

# generate PCA plot
pathway_pca(abundance = annotated_abundance, 
						metadata = metadata, 
						group = "group")

ggsave("../60_picrust2/group_pca_all.pdf", units = "px", width = 2400, height = 2000)

pathway_pca(abundance = annotated_abundance_signif, 
						metadata = metadata, 
						group = "group")

ggsave("../60_picrust2/group_pca_signif.pdf", units = "px", width = 2400, height = 2000)

pathway_pca(abundance = annotated_abundance_notsignif, 
						metadata = metadata, 
						group = "group")

ggsave("../60_picrust2/group_pca_notsignif.pdf", units = "px", width = 2400, height = 2000)


# pathway barplots --------------------------------------------------------

pathways_relab <- annotated_abundance %>%
	mutate(across(everything(), ~ .x / sum(.x))) %>%
	rownames_to_column(var = "feature") %>%
	pivot_longer(cols = all_of(samples),
							 names_to = "sample_name",
							 values_to = "abundance") %>%
	left_join(metadata %>% select(sample_name, cohort, trt, group))


pathways_ab <- annotated_abundance %>%
	rownames_to_column(var = "feature") %>%
	pivot_longer(cols = all_of(samples),
							 names_to = "sample_name",
							 values_to = "abundance") %>%
	left_join(metadata) %>%
	select(-c("gender"))

# relab
pathways_top <- pathways_relab %>%
	group_by(group) %>%
	slice_max(abundance, n = 10) 

pathways <- pathways_relab %>%
	filter(feature %in% pathways_top$feature)

# abundance
pathways_top <- pathways_ab %>%
	group_by(group) %>%
	slice_max(abundance, n = 10) 

pathways <- pathways_ab %>%
	filter(feature %in% pathways_top$feature)

pal <- c("#FBE183FF", "#F4C40FFF", "#FE9B00FF", "#D8443CFF", "#9B3441FF", "#DE597CFF", "#E87B89FF", "#E6A2A6FF", "#AA7AA1FF", "#9F5691FF", "#633372FF", "#1F6E9CFF", "#2B9B81FF", "#92C051FF", "#c8e077", "#FBE183FF", "#F4C40FFF")

pathways %>% ggplot(aes(x = forcats::fct_reorder(sample_name,abundance), y = abundance, fill = forcats::fct_reorder(feature, abundance, .desc = T))) +
	geom_col(position = "fill") +
	ylab("relative abundance") +
	scale_fill_manual(values = pal) +
	scale_y_continuous(expand = expansion(add = c(0,0))) +
	facet_grid(.~ group, space = "free", scales = "free") +
	theme_bw(base_size = 14) +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		strip.background =element_rect(fill="white"),
		panel.spacing = unit(2, "pt"),
		#strip.text = element_text(size = 11),
		axis.title.x = element_blank(),
		axis.ticks.x = element_blank(),
		axis.text.x = element_blank()
	)

#ggsave("../60_picrust2/pathway_barplot_legend.pdf", units = "px", width = 2954*2, height = 894*2, dpi = 320)
ggsave("../60_picrust2/pathway_barplot.pdf", units = "px", width = 2502*2, height = 1200*2, dpi = 320)


# by cohort ---------------------------------------------------------------
daa <- pathway_daa(abundance = abundance_data_filtered %>% column_to_rownames("#OTU ID"), 
									 metadata = metadata, 
									 group = "cohort", 
									 daa_method = "LinDA",
									 p.adjust = "bonferroni"
)

# filter to keep only significantly different pathways 
signif_daa <- daa %>% filter(p_adjust <= 0.05)
notsignif_daa <- daa %>% filter(p_adjust > 0.05)

# divide signif_daa into two to run pathway_annotation function
path_annot <- signif_daa %>% 
	pathway_annotation(daa_results_df = .,
										 pathway = "MetaCyc",
										 ko_to_kegg = FALSE)

notsignif_path_annot <- notsignif_daa %>%
	pathway_annotation(daa_results_df = .,
										 pathway = "MetaCyc",
										 ko_to_kegg = FALSE)

# change the pathway column to description in the abundance table 
signif_annotated_abundance <- abundance_data_filtered %>% filter(`#OTU ID` %in% path_annot$feature) %>%
	dplyr::rename(feature = `#OTU ID`) %>%
	left_join(path_annot %>% select(c("feature","description")),
						by = "feature") %>%
	distinct() %>%
	column_to_rownames("description") %>%
	select(-feature)

# generate heatmap
pathway_heatmap(abundance = signif_annotated_abundance, 
								metadata = metadata, 
								group = "cohort")

#ggsave("../60_picrust2/group_heatmap.png", units = "px", width = 4500, height = 4500)

# generate not significant heatmap
notsignif_annotated_abundance <- abundance_data %>% filter(`#OTU ID` %in% notsignif_path_annot$feature) %>%
	dplyr::rename(feature = `#OTU ID`) %>%
	left_join(notsignif_path_annot %>% select(c("feature","description")),
						by = "feature") %>%
	distinct() %>%
	column_to_rownames("description") %>%
	select(-feature)

pathway_heatmap(abundance = notsignif_annotated_abundance, 
								metadata = metadata, 
								group = "cohort")

#ggsave("../60_picrust2/group_heatmap_notsignif.png", units = "px", width = 4500, height = 4500)

# generate pca plot
pathway_pca(abundance = signif_annotated_abundance, 
						metadata = metadata, 
						group = "group")

#ggsave("../60_picrust2/cohort_pca.png", units = "px", width = 2400, height = 2000)

pathway_pca(abundance = notsignif_annotated_abundance, 
						metadata = metadata, 
						group = "group")

# annotate all and plot
all_annot <- pathway_annotation(daa_results_df = daa,
																pathway = "MetaCyc",
																ko_to_kegg = FALSE)

all_annot_abundance <- abundance_data_filtered %>% filter(`#OTU ID` %in% all_annot$feature) %>%
	dplyr::rename(feature = `#OTU ID`) %>%
	left_join(all_annot %>% select(c("feature","description")),
						by = "feature") %>%
	distinct() %>%
	column_to_rownames("description") %>%
	select(-feature)

pathway_pca(abundance = all_annot_abundance, 
						metadata = metadata, 
						group = "group")


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
RS_notsignif_daa <- RS_daa %>% filter(p_adjust > 0.05)
RS_signif_daa <- RS_daa %>% filter(p_adjust < 0.05)

# divide signif_daa into two to run pathway_annotation function

RS_signif_path_annot <- RS_signif_daa %>% 
	pathway_annotation(daa_results_df = .,
										 pathway = "MetaCyc",
										 ko_to_kegg = FALSE)

RS_notsignif_path_annot <- RS_notsignif_daa %>%
	pathway_annotation(daa_results_df = .,
										 pathway = "MetaCyc",
										 ko_to_kegg = FALSE)


# change the pathway column to description in the abundance table 
RS_annotated_abundance <- RS_abund %>% filter(`#OTU ID` %in% RS_signif_path_annot$feature) %>%
	dplyr::rename(feature = `#OTU ID`) %>%
	left_join(RS_signif_path_annot %>% select(c("feature","description")),
						by = "feature") %>%
	distinct() %>%
	column_to_rownames("description") %>%
	select(-feature)

# significant PCA plot
pathway_pca(abundance = RS_annotated_abundance %>% filter(rowSums(. != 0) > 0), 
						metadata = RS_metadata, 
						group = "cohort")

ggsave("../60_picrust2/RS_pca.png", units = "px", width = 2400, height = 2000)

# generate not significant PCA plot
RS_notsignif_annotated_abundance <- RS_abund %>% filter(`#OTU ID` %in% RS_notsignif_path_annot$feature) %>%
	dplyr::rename(feature = `#OTU ID`) %>%
	left_join(RS_notsignif_path_annot %>% select(c("feature","description")),
						by = "feature") %>%
	distinct() %>%
	column_to_rownames("description") %>%
	select(-feature)

# generate pca plot
pathway_pca(abundance = RS_notsignif_annotated_abundance %>% filter(rowSums(. != 0) > 0), 
						metadata = RS_metadata, 
						group = "cohort")

#ggsave("../60_picrust2/RS_pca.png", units = "px", width = 2400, height = 2000)


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
C_notsignif_daa <- C_daa %>% filter(p_adjust > 0.05)
C_signif_daa <- C_daa %>% filter(p_adjust < 0.05)

# run pathway_annotation function
C_signif_path_annot <- C_signif_daa %>% 
	pathway_annotation(daa_results_df = .,
										 pathway = "MetaCyc",
										 ko_to_kegg = FALSE)

C_notsignif_path_annot <- C_notsignif_daa %>%
	pathway_annotation(daa_results_df = .,
										 pathway = "MetaCyc",
										 ko_to_kegg = FALSE)

C_all_path_annot <- C_daa %>%
	pathway_annotation(daa_results_df = .,
										 pathway = "MetaCyc",
										 ko_to_kegg = FALSE)

# change the pathway column to description in the abundance table 
C_annotated_abundance <- C_abund %>% filter(`#OTU ID` %in% C_signif_path_annot$feature) %>%
	dplyr::rename(feature = `#OTU ID`) %>%
	left_join(C_signif_path_annot %>% select(c("feature","description")),
						by = "feature") %>%
	distinct() %>%
	column_to_rownames("description") %>%
	select(-feature)

C_notsignif_annotated_abundance <- C_abund %>% filter(`#OTU ID` %in% C_notsignif_path_annot$feature) %>%
	dplyr::rename(feature = `#OTU ID`) %>%
	left_join(C_notsignif_path_annot %>% select(c("feature","description")),
						by = "feature") %>%
	distinct() %>%
	column_to_rownames("description") %>%
	select(-feature)

C_all_annot_abundance <- C_abund %>% filter(`#OTU ID` %in% C_all_path_annot$feature) %>%
	dplyr::rename(feature = `#OTU ID`) %>%
	left_join(C_all_path_annot %>% select(c("feature","description")),
						by = "feature") %>%
	distinct() %>%
	column_to_rownames("description") %>%
	select(-feature)

# generate PCAs
pathway_pca(abundance = C_annotated_abundance %>% filter(rowSums(. != 0) > 0), 
						metadata = C_metadata, 
						group = "cohort")

ggsave("../60_picrust2/C_pca.png", units = "px", width = 2400, height = 2000)

pathway_pca(abundance = C_notsignif_annotated_abundance %>% filter(rowSums(. != 0) > 0), 
						metadata = C_metadata, 
						group = "cohort")
