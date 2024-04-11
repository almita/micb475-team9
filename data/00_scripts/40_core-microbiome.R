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

# separate genus and family
starch_genus <- tax_glom(starch, taxrank = "Genus")
starch_family <- tax_glom(starch, taxrank = "Family")



# genus -------------------------------------------------------------------
# convert to relative abundance
starch_ra <- transform_sample_counts(starch_genus, fun=function(x) x/sum(x))

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
det <- 0.001
prev <- 0.7
name <- "_det0.001_prev0.7_"

group_names <- c("human-C", "human-RS", "mouse-C", "mouse-RS")
gender_names <- c("F", "M")
group_gender_names <- c("human-F", "human-M", "mouse-F", "mouse-M")


# image parameters
w <- 2000
h <- 2000
path <- "../40_core-microbiome/venn/"
ext <- ".png"

# venn diagrams 

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
names(core_group) <- group_names

core_gender <- list(f_core, m_core)
names(core_gender) <- gender_names

core_group_gender <- list(human_f_core, human_m_core, mouse_f_core, mouse_m_core)
names(core_group_gender) <- group_gender_names

ggvenn(core_group, 
			 fill_color = pal, 
			 stroke_size = 0.5, 
			 set_name_size = 4,
			 show_percentage = TRUE)

ggsave(paste0(path,"group",name,"genus",ext), 
			 units = "px", width = w, height = h, 
			 create.dir = TRUE)

ggvenn(core_gender, 
			 fill_color = pal3, 
			 stroke_size = 0.5, 
			 set_name_size = 4,
			 show_percentage = TRUE)

ggsave(paste0(path,"gender",name,"genus",ext), units = "px", width = w, height = h)

ggvenn(core_group_gender, 
			 fill_color = pal3, 
			 stroke_size = 0.5, 
			 set_name_size = 4,
			 show_percentage = TRUE)

ggsave(paste0(path,"group_gender",name,"genus",ext), 
			 units = "px", width = w, height = h)

core_genus <- list_to_data_frame(core_group)
tax_table_genus <- tax_table(starch_genus) %>% 
	as.data.frame() %>% 
	rownames_to_column(var = "key") %>% 
	select(c(key, Genus))

core_genus <- inner_join(tax_table_genus, core_genus)

write.csv(core_genus, "../40_core-microbiome/shared_genus.csv")

print("Shared genera between all groups")

core_genus %>% filter(`human-C`== TRUE &
									 `human-RS`== TRUE &
									 `mouse-C`== TRUE &
									 `mouse-RS` == TRUE) %>%
	pull(Genus) %>% 
	print()

print("Shared genera between all except mouse-RS")

core_genus %>% filter(`human-C`== TRUE &
												`human-RS`== TRUE &
												`mouse-C`== TRUE &
												`mouse-RS` == FALSE) %>%
	pull(Genus) %>% 
	print()

print("Human specific genera")

core_genus %>% filter(`human-C`== TRUE &
												`human-RS`== TRUE &
												`mouse-C`== FALSE &
												`mouse-RS` == FALSE) %>%
	pull(Genus) %>% 
	print()

print("Human C genera")

core_genus %>% filter(`human-C`== TRUE &
											 	`human-RS`== FALSE &
											 	`mouse-C`== FALSE &
											 	`mouse-RS` == FALSE) %>%
	pull(Genus) %>% 
	print()

print("Human RS genera")

core_genus %>% filter(`human-C`== F &
											 	`human-RS`== T &
											 	`mouse-C`== FALSE &
											 	`mouse-RS` == FALSE) %>%
	pull(Genus) %>% 
	print()

print("Mouse specific genera")

core_genus %>% filter(`human-C`== F &
											 	`human-RS`== F &
											 	`mouse-C`== T &
											 	`mouse-RS` == T) %>%
	pull(Genus) %>% 
	print()

print("Mouse C genera")

core_genus %>% filter(`human-C`== F &
											 	`human-RS`== F &
											 	`mouse-C`== T &
											 	`mouse-RS` == F) %>%
	pull(Genus) %>% 
	print()

print("Mouse RS genera")

core_genus %>% filter(`human-C`== F &
											 	`human-RS`== F &
											 	`mouse-C`== F &
											 	`mouse-RS` == T) %>%
	pull(Genus) %>% 
	print()

# family ------------------------------------------------------------------
# convert to relative abundance
starch_ra <- transform_sample_counts(starch_family, fun=function(x) x/sum(x))

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
# image parameters
path <- "../40_core-microbiome/venn/"


# venn diagram 
	
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
names(core_group) <- group_names

core_gender <- list(f_core, m_core)
names(core_gender) <- gender_names

core_group_gender <- list(human_f_core, human_m_core, mouse_f_core, mouse_m_core)
names(core_group_gender) <- group_gender_names

ggvenn(core_group, 
			 fill_color = pal, 
			 stroke_size = 0.5, 
			 set_name_size = 4,
			 show_percentage = TRUE)

ggsave(paste0(path,"group",name,"family",ext), 
			 units = "px", width = w, height = h, 
			 create.dir = TRUE)

ggvenn(core_gender, 
			 fill_color = pal3, 
			 stroke_size = 0.5, 
			 set_name_size = 4,
			 show_percentage = TRUE)

ggsave(paste0(path,"gender",name,"family",ext), units = "px", width = w, height = h)

ggvenn(core_group_gender, 
			 fill_color = pal3, 
			 stroke_size = 0.5, 
			 set_name_size = 4,
			 show_percentage = TRUE)

ggsave(paste0(path,"group_gender",name,"family",ext), 
			 units = "px", width = w, height = h)

core_family <- list_to_data_frame(core_group)
tax_table_family <- tax_table(starch_family) %>% 
	as.data.frame() %>% 
	rownames_to_column(var = "key") %>% 
	select(c(key, Family))

core_family <- inner_join(tax_table_family, core_family)

write.csv(core_family, "../40_core-microbiome/shared_family.csv")

print("Shared families between all groups")

core_family %>% filter(`human-C`== TRUE &
												`human-RS`== TRUE &
												`mouse-C`== TRUE &
												`mouse-RS` == TRUE) %>%
	pull(Family) %>% 
	print()

print("Shared families between all except mouse-RS")

core_family %>% filter(`human-C`== TRUE &
												`human-RS`== TRUE &
												`mouse-C`== TRUE &
												`mouse-RS` == FALSE) %>%
	pull(Family) %>% 
	print()

print("Human specific families")

core_family %>% filter(`human-C`== TRUE &
												`human-RS`== TRUE &
												`mouse-C`== FALSE &
												`mouse-RS` == FALSE) %>%
	pull(Family) %>% 
	print()

print("Human C families")

core_family %>% filter(`human-C`== TRUE &
											 	`human-RS`== FALSE &
											 	`mouse-C`== FALSE &
											 	`mouse-RS` == FALSE) %>%
	pull(Family) %>% 
	print()

print("Human RS families")

core_family %>% filter(`human-C`== F &
											 	`human-RS`== T &
											 	`mouse-C`== FALSE &
											 	`mouse-RS` == FALSE) %>%
	pull(Family) %>% 
	print()

print("Mouse specific families")

core_family %>% filter(`human-C`== F &
												`human-RS`== F &
												`mouse-C`== T &
												`mouse-RS` == T) %>%
	pull(Family) %>% 
	print()

print("Mouse C families")

core_family %>% filter(`human-C`== F &
											 	`human-RS`== F &
											 	`mouse-C`== T &
											 	`mouse-RS` == F) %>%
	pull(Family) %>% 
	print()

print("Mouse RS families")

core_family %>% filter(`human-C`== F &
											 	`human-RS`== F &
											 	`mouse-C`== F &
											 	`mouse-RS` == T) %>%
	pull(Family) %>% 
	print()

