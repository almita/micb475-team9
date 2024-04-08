# micb475-team9

## Meetings

### Meeting 1 - 01/FEB/2024
#### Notes
- Discuss datasets that other groups are looking into
- Want to explore dataset that is not heavily studied in MICB 475
- Considering the Starch Resistant Dataset → received research paper from Evelyn from last term
- Considering Wetland dataset → no research groups have looked into
- For next meeting: select dataset, compose potential research questions

### Meeting 2 - 08/FEB/2024
#### Agenda
- marine microbiome datasets
- diabetic and children gut microbiome datasets
- decide datasets and research questions

#### Notes
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
#### Notes
- Compare metadata from both sets to find comparing variables
- Only found treatment variable to compare
- Group together RS treatments of mouse to 1 treatment group
- Changed column names in both datasets 

### Meeting 3 - 15/FEB/2024
#### Agenda
- metadata cleanup files
- proposal aims
- divide the work

#### Notes
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

#### Notes
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

### Meeting 6 - 07/Mar/2024
#### Agenda
- look over progress (qiime core metrics and diversity plots and statisical analysis)

### Notes
- Plan to make corrections for the proposal (due March 12th)
- Meeting goals: discuss the figures and review scripts
- Looked at Kruskal-Wallis --> 2 groups to look at human vs mice
- Looking at our bar graph: evident difference between human and mice
- For alpha diversity: not completed yet
- For Beta diversiy: have completed
- PiCRUST and Indicator Test to be done
- Functional analysis to be done
- Differential abundance to be done
- After filtering, we eneded up having 5 and 5 for the 2 human groups, but did not change the mice data
- Can look at sex in human data and in mouse data
- Can make 2 more metadata category: RS type & Sex


### Meeting 7 - 14/Mar/2024
#### Agenda
- discuss timepoint 4 of the human paper and see if it should be included to boost number of samples
- go over core microbiome, indicator taxa, differential abundance analyses

### Notes
- PERMANOVA - use beta-diversity and look for statistical significance
- can incease sample size of human data by considering timepoint 4
- differential abundance:
- 	Bacteriodes genus had multiple bars on the bar plots (1,2,3,4)
- 	volcano plot - annotate more by labelling the significant plots with genus
- 	shows significance, which is good
- redo metadata, qiime processing, and other analyses with the additional human samples for next week

### Meeting 8 - 21/Mar/2024
#### Notes
- Chao 1 index - compare with original paper (trimming parameters)
- alpha diversity - significance in human cohort in shannon_group between RS2 and control
- beta diversity - no significance in human treatments, significance seen in mice cohort
- PERMANOVA - significance seen in mice but not human - human cohort had control group being fed some resistant starch, just not as much, but for mice, the control group was not fed anything
- core microbiome - no overlap found - donors of human gut microbiome for the mouse cohort study may have been from different from the human cohort - could go a bit higher than genus (maybe family) to see 
  if there is increased overlap between the two groups
- differential abundance - keep barplot for within mouse and within human groups and don't need a cutoff for either within groups - get rid of barplots for the between cohort analyses - also move up a classification level to family
- perform functional analysis and clean up above analyses for next week


### Meeting 9 - 28/Mar/2024
#### Notes
- abundance is driving the difference in RS in alpha and beta diversity
- core microbiome still at genus level and not changed to family
- lots of significant pathways in PiCRUST analysis - groups are very different
Story:
- humanized mouse model is not good for studying starch resistance in this specific mice and human dataset
- Figure 1: Alpha diversity all 3 panels - overall mention diversity were the same -
- Digure 2: Beta Diversity - only Bray-Curtis
- Figure 3: DESeq volcano plots only - human only, mouse only (shows within species effect), human control-mouse control, human RS (4 panels total) 
- Figure 4: Core Microbiome - 4 way Venn Diagram - 2 panels (Genus and Family)
- Figure 5 - PiCRUST Heat map - don't present (maybe just pcoA plot), but put in manuscript - show complete unique taxonomic groups and different metabolic pathways
- emphasize in paper that the two papers are very different studies - humans used as donors in mice study may have been very different than human cohort used in human study

### Meeting 10 - 04/04/2024
#### Notes
- Presenting slides:
- Title: has to be a conclusion: Mouse model was not an adequate comparison/model for human dataset
- show animations on text heavy sides
- show actual dataset rather than findings of the studies - mentioned humanized mice and how humanized mice are done
- experimental aims - only talk about 3-6
- alpha diversity - title of results slide should be overall finding - make fonts for graphs bigger - can get rid of text and put on presenter's notes - only keep Shannon's and Faith's
- differential abundance - title: resistant starch does affect abundace, but more profound in mouse - label graphs with the comparison rather than having it in text on the side
- Core microbiome title: Few Shared taxa between the mouse and human cohorts - remove text and blow up Venn diagrams - keep Genus diagram only, not family
- Conclusion: go through each analysis and mention what the main finding was
- Future Directions: good
- Narrative - we initially thought they would be similiar, we did these tests to take a look, we ended up seeing drastic differences
