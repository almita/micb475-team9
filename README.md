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

### Meeting with group members only - 10/FEB/2024
#### notes
- Compare metadata from both sets to find comparing variables
- Only found treatment variable to compare
- Group together RS treatments of mouse to 1 treatment group
- Changed column names in both datasets 

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
-	aim 2 - qiime2 processing, done in parallel on two datasets - check what kind of data is used in mouse and human dataset, create a detached screen when running in parallel, use the whole SILVA database for taxonomic classification
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


### Meeting 4 - 22/FEB/2024
#### Agenda
- ask about paired-end reads for mouse dataset - how do we distinguish forward & reverse reads?
- got the quality plot for human dataset, unsure if demultiplexed or not? - need a barcode if we want to demultiplex ourselves
- review team proposal

#### notes
- Initial analysis was not as expected
- Write up proposal as is and redo meta data analysis afterwards
- take a look at different resistant starch in the mouse dataset
- See greater seperation if PcoA plot did not have the control
- Do analysis with positive control to see if there is seperation
	- include Inulin sample in mouse data set
- Filter out more sampples to control for more variables in human like weight to get more groups
- Potentially add or explore new data set after proposal
- use mice data for PcoA for proposal
- make sure all bioinformatic tools are cited Qiime2 and dada2
- figure legends title for overall takeaway for the figure and detail explanation for each figure


### Meeting 5 - 29/Feb/2024
#### Agenda
- identify next steps for analyses
- update proposal if necessary

### Notes
- PcoA plot:
  	- there are some groups that are formed, but metadata is not comprehensive enough to show significant data
  	- for both mice and human, there is no significance between control and starch groups
- indicator taxa analyses between human and mice data
  	- use merged dataset
- all analyses (everything that will be done in R) will be continued using the merged dataset, except for the PcoA plot which will require qiime2 processing separately
- mouse data does not specify which diet types include type 2 resistant starch, just mentions resistant starch as a whole
  	- we can do our own research to see which diet types refer to a certain resistant starch type and perform comparisons from there
  	- difference between type 2 resistant starch and other types is based on differences in processing
  	- can potentially use "resistant starch" as an umbrella term from now to be inclusive of all resistant starch types that may have used in the mouse 	and human data sets
 
- work on aims 3 and 4 for next week

### Meeting 5 - 29/Feb/2024
#### Progress/Agenda
- processed human and mouse files individually through Qiime for core metrics
- generated diversity plots and statistical analysis 

#### Agenda
- look over progress (qiime core metrics and diversity plots and statisical analysis)

