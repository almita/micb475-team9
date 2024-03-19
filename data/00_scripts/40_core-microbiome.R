## CORE MICROBIOME

# dependencies
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggvenn)
library(sf)

# data
starch <- read_rds("../21_phyloseq/starch_filtered_phyloseq.RDS")

# color palettes
pal <- c("#DD4124FF", "#faa478", "#0269a1FF", "#6ac1de")
pal3 <- c("#FDAE61FF", "#ABDDA4FF", "#FDAE61FF", "#ABDDA4FF")

# convert to relative abundance
starch_ra <- transform_sample_counts(starch, fun=function(x) x/sum(x))

# get groups
mouse_control <- subset_samples(starch_ra, group=="mouse-C")
mouse_RS <- subset_samples(starch_ra, group=="mouse-RS")
human_control <- subset_samples(starch_ra, group=="human-C")
human_RS <- subset_samples(starch_ra, group=="human-RS")

# get genders by group
mouse_f <- subset_samples(starch_ra, cohort=="mouse" | gender=="F")
mouse_m <- subset_samples(starch_ra, cohort=="mouse" | gender=="M")
human_f <- subset_samples(starch_ra, cohort=="human" | gender=="F")
human_m <- subset_samples(starch_ra, cohort=="human" | gender=="M")

# get genders
f <- subset_samples(starch_ra, gender=="F")
m <- subset_samples(starch_ra, gender=="M")

# core microbiome
## parameters
det_v <- c(0, 1/1000, 1/100, 1/100)
prev_v <- c(1/100, 0.7, 0.6, 0.5)
names <- c("_0-01", "_1000-70", "_100_60", "_100_50")

# image parameters
w <- 2000
h <- 2000
path <- "../40_core-microbiome/venn/"
ext <- ".png"

# venn diagram loop
for (i in 1:length(det_v)){
	det <- det_v[i]
	prev <- prev_v[i]
	name <- names[i]

	mouse_control_core <- core_members(mouse_control, detection=det, prevalence=prev)
	mouse_RS_core <- core_members(mouse_RS, detection=det, prevalence=prev)
	human_control_core <- core_members(human_control,detection=det, prevalence=prev)
	human_RS_core <- core_members(human_RS,detection=det, prevalence=prev)
	
	f_core <- core_members(f, detection=det, prevalence=prev)
	m_core <- core_members(m, detection=det, prevalence=prev)
	
	mouse_f_core <- core_members(mouse_f, detection=det, prevalence=prev)
	mouse_m_core <- core_members(mouse_m, detection=det, prevalence=prev)
	human_f_core <- core_members(human_f,detection=det, prevalence=prev)
	human_m_core <- core_members(human_m,detection=det, prevalence=prev)
	
	core_group <- list(human_control_core, human_RS_core, mouse_control_core, mouse_RS_core)
	group_names <- c("human-C", "human-RS", "mouse-C", "mouse-RS")
	names(core_group) <- group_names
	
	core_gender <- list(f_core, m_core)
	gender_names <- c("F", "M")
	names(core_gender) <- gender_names
	
	core_group_gender <- list(human_f_core, human_m_core, mouse_f_core, mouse_m_core)
	group_gender_names <- c("human-F", "human-M", "mouse-F", "mouse-M")
	names(core_group_gender) <- group_gender_names
	
	core <- append(core_group, core_gender)
	core_names <- append(group_names, gender_names)
	names(core) <- core_names
	
	ggvenn(core_group, 
				 fill_color = pal, 
				 stroke_size = 0.5, 
				 set_name_size = 4,
				 show_percentage = FALSE)
	
	ggsave(paste0(path,"group",name,ext), 
				 units = "px", width = w, height = h, 
				 create.dir = TRUE)
	
	ggvenn(core_gender, 
				 fill_color = pal3, 
				 stroke_size = 0.5, 
				 set_name_size = 4,
				 show_percentage = FALSE)
	
	ggsave(paste0(path,"gender",name,ext), units = "px", width = w, height = h)
	
	ggvenn(core_group_gender, 
				 fill_color = pal3, 
				 stroke_size = 0.5, 
				 set_name_size = 4,
				 show_percentage = FALSE)
	
	ggsave(paste0(path,"group_gender",name,ext), 
				 units = "px", width = w, height = h)
	
}


