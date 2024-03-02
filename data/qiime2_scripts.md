# project 2: qiime2 steps

## importing

### human:

```bash
qiime tools import \
--type "SampleData[SequencesWithQuality]" \
--input-format SingleEndFastqManifestPhred33V2 \
--input-path /home/qiime2/data/metadata/human_manifest.tsv \
--output-path demux_human.qza
```

Visualizing demux file:

```bash
qiime demux summarize \
--i-data demux_human.qza \
--o-visualization demux_human.qzv
```

### mouse:

```bash
qiime tools import \
--type "SampleData[SequencesWithQuality]" \
--input-format SingleEndFastqManifestPhred33V2 \
--input-path /home/qiime2/data/metadata/mouse_manifest_corrected.tsv \
--output-path demux_mouse.qza
```

Visualizing demux file:

```bash
qiime demux summarize \
--i-data demux_mouse.qza \
--o-visualization demux_mouse.qzv
```

## denoising

### human:

```bash
qiime dada2 denoise-single \
--i-demultiplexed-seqs /home/qiime2/data/output/demux_human.qza \
--p-trim-left 0 \
--p-trunc-len 260 \
--o-representative-sequences rep-seqs_human.qza \
--o-table table_human.qza \
--o-denoising-stats stats_human.qza
```

#### denoising visualizations

ASVs

```bash
qiime feature-table tabulate-seqs \
  --i-data rep-seqs_human.qza \
  --o-visualization rep-seqs_human.qzv
```

DADA2 stats

```bash
qiime metadata tabulate \
  --m-input-file stats_human.qza  \
  --o-visualization stats_human.qzv
```

Feature table

```bash
qiime feature-table summarize \
  --i-table table_human.qza \
  --o-visualization table_human.qzv \
  --m-sample-metadata-file /home/qiime2/data/metadata/metadata.tsv
```

### mouse:

```bash
qiime dada2 denoise-single \
--i-demultiplexed-seqs demux_mouse.qza \
--p-trim-left 0 \
--p-trunc-len 260 \
--o-representative-sequences rep-seqs_mouse.qza \
--o-table table_mouse.qza \
--o-denoising-stats stats_mouse.qza
```

#### denoising visualizations

ASVs

```bash
qiime feature-table tabulate-seqs \
  --i-data rep-seqs_mouse.qza \
  --o-visualization rep-seqs_mouse.qzv
```

DADA2 stats

```bash
qiime metadata tabulate \
  --m-input-file stats_mouse.qza  \
  --o-visualization stats_mouse.qzv
```

Feature table

```bash
qiime feature-table summarize \
  --i-table table_mouse.qza \
  --o-visualization table_mouse.qzv \
  --m-sample-metadata-file /home/qiime2/data/metadata/metadata.tsv
```

## merging datasets

Feature tables

```bash
qiime feature-table merge \
 --i-tables table_human.qza \
 --i-tables table_mouse.qza \
 --o-merged-table merged_table.qza
```

Representative sequences

```bash
qiime feature-table merge-seqs \
 --i-data rep-seqs_human.qza \
 --i-data rep-seqs_mouse.qza \
 --o-merged-data merged_rep-seqs.qza
```

Visualization
```bash
qiime feature-table tabulate-seqs \
  --i-data merged_rep-seqs.qza \
  --o-visualization merged_rep-seqs.qzv
```

```bash
qiime feature-table summarize \
  --i-table merged_table.qza \
  --o-visualization merged_table.qzv \
  --m-sample-metadata-file /home/qiime2/data/metadata/metadata.tsv
```
## taxonomic classification

Taxonomic analysis
(Use the complete SILVA classifier, not the region specific ones.)
merged
```bash
qiime feature-classifier classify-sklearn \
  --i-classifier /mnt/datasets/classifiers/silva-138-99-nb-classifier.qza \
  --i-reads merged_rep-seqs.qza \
  --o-classification taxonomy.qza
```

human
```bash
qiime feature-classifier classify-sklearn \
  --i-classifier /mnt/datasets/classifiers/silva-138-99-nb-classifier.qza \
  --i-reads rep-seqs_human.qza \
  --o-classification taxonomy_human.qza
```

Taxonomy barplots
Merged
```bash
qiime taxa barplot \
  --i-table merged_table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file /home/qiime2/data/metadata/metadata.tsv \
  --o-visualization taxa-bar-plots.qzv
```

human
```bash
qiime taxa barplot \
  --i-table table_human.qza \
  --i-taxonomy taxonomy_human.qza \
  --m-metadata-file /home/qiime2/data/metadata/metadata.tsv \
  --o-visualization taxa-bar-plots.qzv
```

Visualization
```bash
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
```
## filter

taxonomy based filtering
```bash
qiime taxa filter-table \
  --i-table merged_table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table merged_table-no-mitochondria-no-chloroplast.qza
```
Visualization
```bash
qiime feature-table summarize \
  --i-table merged_table-no-mitochondria-no-chloroplast.qza \
  --o-visualization merged_table-no-mitochondria-no-chloroplast.qzv \
  --m-sample-metadata-file /home/qiime2/data/metadata/metadata.tsv
```

## phylogenetic tree

```bash
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences merged_rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
```

## rarefaction

alpha-rarefaction curves

```bash
qiime diversity alpha-rarefaction \
  --i-table merged_table-no-mitochondria-no-chloroplast.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 120000 \
  --m-metadata-file /home/qiime2/data/metadata/metadata.tsv \
  --o-visualization alpha-rarefaction.qzv
```

```bash
qiime diversity alpha-rarefaction \
  --i-table merged_table-no-mitochondria-no-chloroplast.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 80000 \
  --m-metadata-file /home/qiime2/data/metadata/metadata.tsv \
  --o-visualization alpha-rarefaction_2.qzv
```

## diversity metrics
`--p-sampling-depth` is the rarefaction parameter

```bash
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table merged_table-no-mitochondria-no-chloroplast.qza \
  --p-sampling-depth 8193 \
  --m-metadata-file /home/qiime2/data/metadata/metadata.tsv \
  --output-dir core-metrics-results
```

## alpha diversity analysis


## beta diversity analysis


## taxonomic classification

Use the complete SILVA classifier, not the region specific ones.

```bash
qiime feature-classifier classify-sklearn \
  --i-classifier silva-v4-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza
```
