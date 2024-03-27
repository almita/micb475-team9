## DIFFERENTIAL ABUNDANCE ##

# load dependencies
library(tidyverse)
library(phyloseq)
library(DESeq2)

# read data
starch <- read_rds("starch_filtered_phyloseq.RDS")

# adjust for zeroes
starch_1 <- transform_sample_counts(starch, function(x) x+1)

# create deseq objects

starch_dsq_grp <- phyloseq_to_deseq2(starch_1, ~`group`)
starch_DESEQ_grp <- DESeq(starch_dsq_grp)

starch_dsq_gnd <- phyloseq_to_deseq2(starch_1, ~`gender`)
starch_DESEQ_gnd <- DESeq(starch_dsq_gnd)

starch_dsq <- phyloseq_to_deseq2(starch, ~group)
starch_sizefactors <- estimateSizeFactors(starch_dsq, type="poscounts")
starch_DESEQ <- DESeq(starch_sizefactors,
											test="Wald", 
											fitType="parametric")

# deseq results

res_human <- results(starch_DESEQ, tidy=TRUE,
										 contrast = c("group", "human-RS", "human-C"))

alpha = 0.05 
res = results(starch_DESEQ, alpha=alpha)
res = res[order(res$padj, na.last=NA), ]


tmp %>%
	ggplot(aes(x=padj)) +
	geom_density()
