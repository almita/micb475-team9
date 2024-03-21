# PICRUST QIIME SECTION

## Filter the table.qza to remove all the features with 5 or lower counts. 
## (Quality control) This will make the PICRUSt2 analysis much faster.

qiime feature-table filter-features \
  --i-table merged_table-no-mitochondria-no-chloroplast-no-unassigned.qza \
  --p-min-frequency 5 \
  --o-filtered-table feature-frequency-filtered-table.qza

## Running the picrust2 qiime2 plugin (takes 15-20 mins) do on detached screen
  
qiime picrust2 full-pipeline \
  --i-table feature-frequency-filtered-table.qza \
  --i-seq merged_rep-seqs.qza \
  --output-dir q2-picrust2_output \
  --p-placement-tool sepp \
  --p-hsp-method pic \
  --p-max-nsti 2 \
  --verboseqiime tools export \
  --input-path ec_metagenome.qza \
  --output-path ecmeta_exported

## Create export files for R analysis downstream

qiime tools export \
  --input-path ec_metagenome.qza \
  --output-path ecmeta_exported

qiime tools export \
  --input-path ko_metagenome.qza \
  --output-path kometa_exported

qiime tools export \
  --input-path pathway_abundance.qza \
  --output-path pathabun_exported

## Convert feature table to tsv
biom convert \
  -i ecmeta_exported/feature-table.biom \
  -o ecmeta_exported/ecmeta.tsv \
  --to-tsv

biom convert \
  -i kometa_exported/feature-table.biom \
  -o kometa_exported/kometa.tsv \
  --to-tsv

biom convert \
  -i pathabun_exported/feature-table.biom \
  -o pathabun_exported/pathway_abundance.tsv \
  --to-tsv

## Transfer to Local Computer
scp -r root@10.19.139.174:~/data/20_qiime/q2-picrust2_output .
