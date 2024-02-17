# micb475-team9

## meetings

### Meeting 1 - 01/FEB/2024
#### notes
- Discuss datasets that other groups are looking into
- Want to explore dataset that is not heavily studied in MICB 475
- Considering the Starch Resistant Dataset → received research paper from Evelyn from last term
- Considering Wetland dataset → no research groups have looked into
- For next meeting: select dataset, compose potential research questions

### Meeting 2 - 08/FEB/2024
#### agenda
- marine microbiome datasets
- diabetic and children gut microbiome datasets
- decide datasets and research questions

#### notes
- validate mouse microbiome model vs human microbiome:
	- compare first within cohort then between cohorts
	- a good model would cluster very similarly
	- if the datasets are different (paired vs single end reads) we need to process them separately until we get to the table.qzv
		- qiime has a merging function, Dr. Evelyn will send us a tutorial for it
	- if we find that the model is bad we need to reach out to the researchers just to let them know that we're going to be "criticizing" their paper
- record everything we do to clean up the metadata
- get a head start on the metadata
- research on strains of the bacteria we found and the relation with human health


### Meeting 3 - 15/FEB/2024
#### agenda
- metadata cleanup files
- proposal aims
- divide the work

#### notes
- take out mnt files from server for both datasets reconsile with metadata (only same datafiles stay)
- use new manifest file when running the qiime in parallel
- merge meta datafiles
- proposal aims:
- 	aim 1 - data wrangling
-	aim 2 - qiime2 processing, done in parallel on two datasets - check what kind of data is used in mouse and human dataset, create a detached screen when running in parallel
- 	aim 3 - alpha beta diversity metrics to compare cohorts (qiime and R) - PCA plot to see if points are clustered together
- 	aim 4 - core microbiome to compare
- 	aim 5 - differential abundance
- sample sizes:
- 	human dataset - 31
- 	mouse dataset - 78 (15 of them are non-resistant starch, 63 are resistant starch)
- 	add new column for both called cohort
- 	when merging, create new column on metadata that combines organism and starch resistance (human-starch, human-no starch, mouse-starch, mouse-no starch)
- send email of quality plot or demux(qzv files) before trimming data, also send alpha-rarefaction data, combine datasets after denoising steps
- have alpha and beta diversity metrics ready after reading week
- for agendas, include summary of progress along with what we want to discuss
