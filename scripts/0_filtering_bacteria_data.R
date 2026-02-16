#****************************#
#### CAMARGUE - BACTERIES ####
#****************************#

#Packages----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"), update = FALSE)
install.packages("microViz",
 repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos")))
install.packages("paletteer")


library(tidyverse)
library(phyloseq)
# library(microbiome)
library(microViz)
library(paletteer)
# library(ggrepel)

# library(vegan)

#Silva digest----
#source("scripts/silva_work.R") 

silva_levels <- "data_ref/taxa2ranks_ssu_138.1.table"

silva_levels %>%
  read_tsv(show_col_types = FALSE) -> taxa2ranks

good_taxonomic_levels <- c("domain", "phylum", "class",
                           "order", "family", "genus")
#*******************----
#Prepare data for analyses----
#**Load data----
illumina_data = "data/Camargue_diversity_MiSeq_16S_341F_785R_20240417_405_samples.OTU.filtered.cleaved.nosubstringOTUs.mumu.table2"

illumina_data %>%
  read_tsv(show_col_types = FALSE, progress = FALSE) %>%
  colnames()

illumina_data %>%
  read_tsv(show_col_types = FALSE, progress = FALSE) %>%
  filter(str_detect(taxonomy, "yanobacter")) %>% 
  view()

#**Count organelles----
illumina_data %>%
  read_tsv(show_col_types = FALSE, progress = FALSE) %>%
  mutate(organelle = str_detect(taxonomy, "Chloroplast")|str_detect(taxonomy, "Mitochondria")) %>% 
  group_by(organelle) %>% 
  summarise(OTU_number = length(total), read_count = sum(total))

#read count
4822435/(4412547+4822435) #52% des séquences proviennent d'organelles
#OTU count
1341/(1341+29074) #4.4% des OTU sont organelles

#**Filter out organelles and stuff----
#voir la distribution des No_hit
illumina_data %>%
  read_tsv(na = c("", "NA", "*"), 
           show_col_types = FALSE, 
           progress = FALSE) %>% #30,415 OTUs
  select(OTU, taxonomy, identity, length) %>% 
  filter(!str_detect(taxonomy, "Bacteria"))   #1,034 OTUS non assignées aux bactéries

#Filtrer les données
illumina_data %>%
  read_tsv(na = c("", "NA", "*"), 
           show_col_types = FALSE, 
           progress = FALSE) %>% #30,415 OTUs
  filter(str_detect(taxonomy, "Bacteria")) %>% #29,381 OTUs
  filter(!str_detect(taxonomy, "Chloroplast")) %>% 
  filter(!str_detect(taxonomy, "Mitochondria")) %>%
  filter(!str_detect(taxonomy, "endosymbionts")) %>%
  filter(length > 300) %>% 
  filter(identity > 80) %>% 
  select(-total, -amplicon, -cloud, -chimera, -spread, -length, -quality, -identity, -abundance, -sequence, -references) -> filtered

##**Prepare metadata----
filtered %>%
  select(-OTU, -taxonomy) %>% 
  colnames() -> all_samples

filtered %>%
  select(-OTU, -taxonomy) %>%
  select(contains("-CAM-")) %>% 
  colnames() -> only_environmental_samples

filtered %>%
  select(-OTU, -taxonomy) %>%
  select(contains("Zymo")) %>% 
  colnames() -> zymo_samples

filtered %>%
  select(-OTU, -taxonomy) %>%
  select(contains(c("Tgen","Tneg"))) %>% 
  colnames() -> neg_samples

data.frame(
  sample_ID = all_samples, 
  type = case_when(
    str_detect(all_samples, "-CAM-") ~ "ENV",
    str_detect(all_samples, "Zymo") ~ "ZYMO",
    .default = "NEG"),
  year = if_else(str_detect(all_samples, "-CAM-"), 
                 str_sub(all_samples, 1, 2), 
                 NA), 
  grid = if_else(str_detect(all_samples, "-CAM-"), 
                 str_sub(all_samples, 1, 9), 
                 NA), 
  quadrat = if_else(str_detect(all_samples, "-CAM-"), 
                    str_sub(all_samples, 1, 11), 
                    NA)
) %>% 
  column_to_rownames("sample_ID") %>% 
  mutate(sample_ID = rownames(.)) -> my_mtd

read_delim("data/Metadata_grid_CAM.csv")

#**Stats descriptives----
columns_to_keep <- c("OTU","taxonomy",only_environmental_samples)

filtered %>% 
  select(all_of(columns_to_keep)) %>% 
  replace(is.na(.),0) %>% 
  mutate(.before = 1,total = rowSums(.[,3:ncol(.)])) %>%  
  filter(total > 0) %>% 
  select(-total) %>% 
  dim()
#Sur les 377 échantillons des grilles, il y a 27 100 OTU bactériennes non organelles

filtered %>% 
  select(all_of(columns_to_keep)) %>% 
  replace(is.na(.),0) %>% 
  mutate(.before = 1,total = rowSums(.[,3:ncol(.)])) %>%  
  pull(total) %>% 
  sum()
#Sur les 377 échantillons des grilles, il y a 3 843 518 reads bactériens non organelles

#*******************----
#Building a phyloseq object based on OTUs----
#**OTU----
columns_to_keep <- c("OTU", all_samples)

filtered %>% 
  select(all_of(columns_to_keep)) %>% 
  replace(is.na(.),0) %>% 
  mutate(.before = 1,
         total = rowSums(.[,colnames(.)!="OTU"])) %>% 
  # Remove OTU absent in all samples and controls (don't remove any here)
  filter(total>0) %>% 
  select(-total) %>% 
  column_to_rownames(var = "OTU") -> my_otu

#**TAX à l'ancêtre commun----

# Cleans and reshapes taxonomic data by mapping OTUs to their respective ranks
# and filling missing levels with the lowest known classification (prefixed with "Unk_").
# It ultimately converts the data from a long format into a structured wide 
# table for downstream analysis.

filtered %>% 
  select(OTU,taxonomy) %>% 
  filter(OTU %in% rownames(my_otu)) %>% # Remove absent OTU (useless here)
  separate_longer_delim(col = taxonomy, delim = "|") %>% 

  left_join(x = .,
            y = taxa2ranks,
            by = c("taxonomy" = "taxonomic_path")) %>%
  # Filter out unknwown taxonomic levels
  filter(taxonomic_rank %in% good_taxonomic_levels) %>% 
  # relevel taxonomic_rank
  mutate(taxonomic_rank = fct_relevel(taxonomic_rank, good_taxonomic_levels),
         # More practical naming
         taxonomy = str_replace(taxonomy, "Peptostreptococcales-Tissierellales", "Peptostreptococcaceae"), 
         taxonomy = str_replace(taxonomy, "Methylobacterium-Methylorubrum", "Methylobacterium"),
         taxonomy = str_replace(taxonomy, "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium","Allorhizobium")) %>%
  group_by(OTU) %>%
  # Remove duplicated family rank (OTU 2128 and 25741)
  filter(!duplicated(taxonomic_rank)) %>% 
  # Lower taxonomy level known per OTU
  mutate(lower_kwn_tax = last(taxonomy, order_by = taxonomic_rank )) %>%
  ungroup() %>% 
  # Add missing rank (Unkwown OTU taxonomy) before pivot longer
  # (function is faster outside group_by)
  complete(OTU, taxonomic_rank) %>%  
  group_by(OTU) %>% 
  # fill added NA rows
  fill(lower_kwn_tax, .direction = "downup") %>%
  mutate(taxonomy = case_when(is.na(taxonomy) ~ paste0("Unk_", lower_kwn_tax),
                              .default = taxonomy 
                           )) %>% 
  ungroup() %>%
  # Pivot wider with all unkonw rank rename -> Unk_lower-rank-known
  select(-lower_kwn_tax) %>% 
  pivot_wider(names_from = taxonomic_rank, values_from = taxonomy) %>%
  mutate(OTU = as.character(OTU)) %>%
  column_to_rownames("OTU")-> my_tax

#**Phyloseq----
sample_data(my_mtd) -> MTD
otu_table(my_otu,taxa_are_rows = T) -> OTU
tax_table(as.matrix(my_tax)) -> TAX

sample_names(MTD)
taxa_names(OTU)
taxa_names(TAX)
ps = phyloseq(OTU, TAX, MTD)

#*******************----
##CONTROLS----
prune_samples(sample_data(ps)$type == "NEG", ps) %>% 
  tax_fix() %>% 
  comp_barplot(tax_level = "class", n=30, tax_transform_for_plot = "identity")+
  theme_minimal() +
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.justification = "top",
        legend.box = "vertical",
        legend.key.size = unit(0.5, "cm"),
        axis.text.x = element_text(size = 8, angle=90, vjust=0.4, hjust=1.2),
        legend.text = element_text(size = 11))

##MOCKS----
#Qu'y a-t-il dans les mocks ? 
prune_samples(sample_data(ps)$type == "ZYMO", ps) %>% 
  tax_fix() %>% 
  comp_barplot(tax_level = "family", n=15)+
  theme_minimal() +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.justification = "top",
        legend.box = "vertical",
        legend.key.size = unit(0.5, "cm"),
        axis.text.x = element_text(size = 8, angle=90, vjust=0.4, hjust=1.2),
        legend.text = element_text(size = 11))

#**Theoretical mock----
read_delim("data_ref/mock_zymo_expected.csv") %>% 
  rename_with(., tolower) %>% 
  rename(frq = rel_ab) %>% 
  mutate(species = str_replace(species, pattern = " ", replacement = "_")) -> mock_expected

#**Genus barplot----
zymo_samples -> my_ordered_samples

my_ordered_samples %>% 
  str_sub(., end = -12) -> my_labels

mock_expected %>%
  select( genus, frq, frq_log) %>%
  pivot_longer(cols = c(frq,frq_log), names_to = "sample_ID") %>%
  mutate(sample_ID = str_c("Expected_",sample_ID)) %>%
  dplyr::rename(frq = "value") -> mock_exp

my_otu %>% 
  select(all_of(zymo_samples)) %>% 
  filter(rowSums(.) > 0) %>% 
  mutate(.before = 1, OTU = rownames(.)) %>%
  pivot_longer(cols = c(2:ncol(.)), names_to = "sample_ID", values_to = "reads") %>% 
  filter(reads > 0) -> my_df_ZYMO

my_tax %>%
  filter(rownames(.) %in% my_df_ZYMO$OTU) %>% 
  mutate(OTU = rownames(.))-> my_tax_ZYMO

my_df_ZYMO %>% 
  left_join(., my_tax_ZYMO, by = "OTU") -> my_df_ZYMO 
my_df_ZYMO$sample_ID %>% unique

my_df_ZYMO %>% 
  group_by(sample_ID, genus) %>% 
  summarise(read_count = sum(reads)) %>%
  group_by(sample_ID) %>% 
  mutate(frq = read_count/sum(read_count)) %>%
  ungroup() %>%
 # mutate(genus = ifelse(frq < 0.0001, "Others", genus)) %>% 
  bind_rows(., mock_exp) %>% 
  ggplot(., aes(x = sample_ID, y = frq, fill = genus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.justification = "top",
        legend.box = "vertical",
        legend.key.size = unit(0.5, "cm"),
        axis.text.x = element_text(angle=90, vjust=0.4, hjust=1.2),
        legend.text = element_text(size = 11)) +
  guides(fill = guide_legend(ncol = 1))+
  labs(x="Runs", y = "Relative abundances", fill="Genera") +  
  #scale_x_discrete(limits = c(my_ordered_samples, "Mock"), labels = c(my_labels, "Expected"))+
  scale_fill_paletteer_d("ggsci::default_igv")

my_df_ZYMO %>% 
  group_by(sample_ID, genus) %>% 
  summarise(read_count = sum(reads)) %>%
  group_by(sample_ID) %>% 
  mutate(frq = read_count/sum(read_count)) %>%
  ungroup() %>%
  # mutate(genus = ifelse(frq < 0.0001, "Others", genus)) %>% 
  bind_rows(., mock_exp) %>%
  filter(str_detect(sample_ID, "log")) %>% 
  ggplot(., aes(x = sample_ID, y = frq, fill = genus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.justification = "top",
        legend.box = "vertical",
        legend.key.size = unit(0.5, "cm"),
        axis.text.x = element_text(angle=90, vjust=0.4, hjust=1.2),
        legend.text = element_text(size = 11)) +
  guides(fill = guide_legend(ncol = 1))+
  labs(x="Runs", y = "Relative abundances", fill="Genera") +  
  #scale_x_discrete(limits = c(my_ordered_samples, "Mock"), labels = c(my_labels, "Expected"))+
  scale_fill_paletteer_d("ggsci::default_igv") +
  scale_y_log10()

my_df_ZYMO %>% 
  group_by(sample_ID, genus) %>% 
  summarise(read_count = sum(reads)) %>%
  group_by(sample_ID) %>% 
  mutate(frq = read_count/sum(read_count)) %>%
  ungroup() %>%
  # mutate(genus = ifelse(frq < 0.0001, "Others", genus)) %>% 
  bind_rows(., mock_exp) %>%
  filter(str_detect(sample_ID, "log")) %>% 
  ggplot(., aes(x = sample_ID, y = frq, fill = genus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.justification = "top",
        legend.box = "vertical",
        legend.key.size = unit(0.5, "cm"),
        axis.text.x = element_text(angle=90, vjust=0.4, hjust=1.2),
        legend.text = element_text(size = 11)) +
  guides(fill = guide_legend(ncol = 1))+
  labs(x="Runs", y = "Relative abundances", fill="Genera") +  
  #scale_x_discrete(limits = c(my_ordered_samples, "Mock"), labels = c(my_labels, "Expected"))+
  scale_fill_paletteer_d("ggsci::default_igv") +
  scale_y_log10()
#Sans les unidentified
my_df_ZYMO %>%
  filter(!str_detect(genus, "Unk")) %>% 
  group_by(sample_ID, genus) %>% 
  summarise(read_count = sum(reads)) %>%
  group_by(sample_ID) %>% 
  mutate(frq = read_count/sum(read_count)) %>%
  ungroup() %>%
  mutate(genus = ifelse(frq < 0.001, "Others", genus)) %>% 
  bind_rows(., mock_exp) %>% 
  ggplot(., aes(x = sample_ID, y = frq, fill = genus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.justification = "top",
        legend.box = "vertical",
        legend.key.size = unit(0.5, "cm"),
        axis.text.x = element_text(angle=90, vjust=0.4, hjust=1.2),
        legend.text = element_text(size = 11)) +
  guides(fill = guide_legend(ncol = 1))+
  labs(x="Runs", y = "Relative abundances", fill="Genera") +  
  scale_x_discrete(limits = c(my_ordered_samples, "Mock"), labels = c(my_labels, "Expected"))+
  scale_fill_paletteer_d("ggsci::default_igv")

#je récupère les especes inattendues (qui ne sont pas Unknown)
my_df_ZYMO %>% 
  filter(!genus %in% mock_exp$genus) %>% 
  filter(!is.na(genus)) %>% 
  filter(!str_detect(genus, "Unk")) -> false_positives

false_positives %>% 
  view()

#The false positive most frequent in a mock is below 0.0005 in its sample and below 0.00006 in total ab.
my_df_ZYMO$genus %>% unique
my_df_ZYMO %>% 
  filter(genus %in% mock_exp$genus) %>% 
  filter(!str_detect(genus, "Unk")) %>% 
  group_by(sample_ID, genus) %>% 
  summarise(reads = sum(reads)) -> true_positives

true_positives %>% 
  view() 

true_positives %>% 
  group_by(sample_ID) %>% 
  mutate(local_rel_ab = reads/sum(reads)) %>% 
  arrange(local_rel_ab) %>% 
  view()

#**Family barplot----
zymo_samples -> my_ordered_samples

my_ordered_samples %>% 
  str_sub(., end = -12) -> my_labels

mock_expected %>%
  mutate(sample_ID = "Mock") %>%
  group_by(sample_ID, family, frq) %>% 
  summarise(frq = sum(frq)) -> mock_exp

my_otu %>%
  select(all_of(zymo_samples)) %>% 
  filter(rowSums(.) > 0) %>% 
  mutate(.before = 1, OTU = rownames(.)) %>% 
  pivot_longer(cols = c(2:ncol(.)), names_to = "sample_ID", values_to = "reads") %>% 
  filter(reads > 0) -> my_df_ZYMO

my_tax %>%
  filter(rownames(.) %in% my_df_ZYMO$OTU) %>% 
  mutate(OTU = rownames(.))-> my_tax_ZYMO

my_df_ZYMO %>% 
  left_join(., my_tax_ZYMO, by = "OTU") -> my_df_ZYMO 

my_df_ZYMO %>% 
  group_by(sample_ID, family) %>% 
  summarise(read_count = sum(reads)) %>%
  group_by(sample_ID) %>% 
  mutate(frq = read_count/sum(read_count)) %>%
  ungroup() %>%
  mutate(family = ifelse(frq < 0.001, "Others", family)) %>% 
  bind_rows(., mock_exp) %>% 
  ggplot(., aes(x = sample_ID, y = frq, fill = family)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.justification = "top",
        legend.box = "vertical",
        legend.key.size = unit(0.5, "cm"),
        axis.text.x = element_text(angle=90, vjust=0.4, hjust=1.2),
        legend.text = element_text(size = 11)) +
  guides(fill = guide_legend(ncol = 1))+
  labs(x="Runs", y = "Relative abundances", fill="Families") +  
  scale_x_discrete(limits = c(my_ordered_samples, "Mock"), labels = c(my_labels, "Expected"))+
  scale_fill_paletteer_d("ggsci::default_igv")

#Sans les unidentified
my_df_ZYMO %>%
  filter(!str_detect(family, "Unk")) %>% 
  group_by(sample_ID, family) %>% 
  summarise(read_count = sum(reads)) %>%
  group_by(sample_ID) %>% 
  mutate(frq = read_count/sum(read_count)) %>%
  ungroup() %>%
  mutate(family = ifelse(frq < 0.001, "Others", family)) %>% 
  bind_rows(., mock_exp) %>% 
  ggplot(., aes(x = sample_ID, y = frq, fill = family)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.justification = "top",
        legend.box = "vertical",
        legend.key.size = unit(0.5, "cm"),
        axis.text.x = element_text(angle=90, vjust=0.4, hjust=1.2),
        legend.text = element_text(size = 11)) +
  guides(fill = guide_legend(ncol = 1))+
  labs(x="Runs", y = "Relative abundances", fill="Families") +  
  scale_x_discrete(limits = c(my_ordered_samples, "Mock"), labels = c(my_labels, "Expected"))+
  scale_fill_paletteer_d("ggsci::default_igv")

#je récupère les especes inattendues, et je vire les undidentified
my_df_ZYMO %>% 
  filter(!family %in% mock_exp$family) %>% 
  filter(!is.na(family)) %>% 
  filter(!str_detect(family, "Unk")) -> false_positives

false_positives %>% 
  view()

#The false positive most frequent in a mock is below 0.0005 in its sample and below 0.00006 in total ab.

my_df_ZYMO %>% 
  filter(family %in% mock_exp$family) %>% 
  filter(!str_detect(family, "Unk")) %>%
  group_by(sample_ID, family) %>% 
  summarise(reads = sum(reads)) -> true_positives

true_positives %>% 
  view() 

true_positives %>% 
  group_by(sample_ID) %>% 
  mutate(local_rel_ab = reads/sum(reads)) %>% 
  arrange(local_rel_ab) %>% 
  view()

#*******************----
#MICRODECON----
install.packages("remotes")
remotes::install_github("donaldtmcknight/microDecon")
library("microDecon")

neg <- my_otu[,neg_samples]
env <- my_otu[, only_environmental_samples]

data.frame(OTU = rownames(my_otu)) %>% 
  cbind(., neg) %>% 
  cbind.data.frame(., env) -> df_decon

decontaminated <- decon(data = df_decon,numb.blanks=length(colnames(neg)),numb.ind=length(colnames(env)), taxa=F)
save(decontaminated, file = "outputs/decontaminated.Rdata")

load("outputs/decontaminated.Rdata")

#**OTUs à retirer----
decontaminated$OTUs.removed %>% 
  pull(OTU)-> otus_to_remove
length(otus_to_remove) #il y en a 213

# my_tax %>% 
#   filter(OTU %in% otus_to_remove) %>% view()

#**OTUs seulement présentes dans les puits négatifs----
my_otu %>% 
  mutate(.before = 1, total_env = rowSums(.[,only_environmental_samples])) %>%
  filter(total_env==0) %>% 
  rownames(.)-> OTU_absent_from_env_samples

length(OTU_absent_from_env_samples) #Il y en a 51

#**OTUs dont les reads sont corrigés----
decontaminated$reads.removed %>%
  replace(.,is.na(.),0) %>%
  filter(!OTU %in% otus_to_remove) %>% 
  filter(!OTU %in% OTU_absent_from_env_samples) %>% 
  pull(OTU) -> partially_removed 
length(partially_removed) #Il y en a 121

#**Phyloseq object pour faire graph OTU conta----

otu_table(my_otu,taxa_are_rows = TRUE) -> OTU.decon 
tax_table(as.matrix(my_tax)) -> TAX.decon
sample_data(my_mtd) -> MTD

ps_decon = phyloseq(OTU.decon, TAX.decon, MTD)

#**Visualisation of decontamination----
ps.pa <- transform_sample_counts(ps_decon, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$type == "NEG", ps.pa) 
ps.pa.pos <- prune_samples(!(sample_data(ps.pa)$type == "NEG"), ps.pa)

df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), 
                    pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=case_when(
                      taxa_names(ps.pa) %in% OTU_absent_from_env_samples ~ "Only neg",
                      taxa_names(ps.pa) %in% otus_to_remove ~ "Full contaminant",
                      taxa_names(ps.pa) %in% partially_removed ~ "Partial contaminant", 
                      TRUE ~ "True OTU"
                    ))

ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + 
  geom_point(position = "jitter") +
  xlab("Prevalence (Negative Controls)") + 
  ylab("Prevalence (True Samples)")  + 
  ggtitle("Prevalence of conserved and removed OTUs following microDecon")

#Graphique avec les nombres de reads 
ps.decon.neg <- prune_samples(sample_data(ps_decon)$type == "NEG", ps_decon) 
ps.decon.pos <- prune_samples(!(sample_data(ps_decon)$type == "NEG"), ps_decon)

df.pa1 <- data.frame(
  OTU = taxa_names(ps.pa),
  pa.pos1 = log10(taxa_sums(ps.decon.pos) + 1),
  pa.neg1 = log10(taxa_sums(ps.decon.neg) + 1),
  contaminant = case_when(
    taxa_names(ps.pa) %in% OTU_absent_from_env_samples ~ "Only neg",
    taxa_names(ps.pa) %in% otus_to_remove ~ "Full contaminant",
    taxa_names(ps.pa) %in% partially_removed ~ "Partial contaminant",
    TRUE ~ "True OTU"
  )
)

ggplot(data=df.pa1, aes(x=pa.neg1, y=pa.pos1, color=contaminant)) + 
  geom_point(position = "jitter")+
  #geom_label_repel(aes(label = OTU,  color = contaminant), size = 3)+
  xlab("Number of reads in negative Controls (Log scale))") + 
  ylab("Number of reads in true Samples (Log-scale)")  + 
  ggtitle("Decontamination following microDecon")

my_tax %>% 
  filter(OTU == "5")

#*******************----
#Objets finaux----
#**Table d'OTU décontaminée----
decontaminated$decon.table %>% 
  mutate(total = rowSums(.[,3:ncol(.)])) %>%
  filter(total > 0) %>%
  select(-total,-Mean.blank) %>% 
  column_to_rownames(var = "OTU") -> my_cleaned_otu_table
#cleaned OTU table est bien débarrassée des otus uniquement présentes dans les neg et des otu à enlever complètement (26886 OTU restantes)

my_cleaned_otu_table %>% 
  summarise(across(everything(), sum)) %>% 
  t() %>% as.data.frame() %>% 
  pull(1)-> total_reads_per_sample

hist(total_reads_per_sample)
summary(total_reads_per_sample)

#clearNames(ps.pa, ps.pa.neg, ps.pa.pos, ps.decon, ps.decon.neg, ps.decon.pos)

#**Taxonomie à l'ancetre commun sur la table d'OTU décontaminée----
my_tax %>%
  rownames_to_column("OTU") %>%
  filter(OTU %in% rownames(my_cleaned_otu_table)) %>%
  column_to_rownames("OTU")-> my_cleaned_tax_table

#**Objets phyloseq ----
#uniquement les puits environnementaux 

my_mtd %>% 
  filter(type=="ENV") -> my_mtd_env 
sample_data(my_mtd_env) -> MTD
otu_table(my_cleaned_otu_table,taxa_are_rows = T) -> OTU
tax_table(as.matrix(my_cleaned_tax_table)) -> TAX

sample_names(MTD)

ps_bacteria = phyloseq(OTU, TAX, MTD)

#*******************----
#Qualité du séquençage----
#**Profondeurs de séquençage----
#On ajoute les nombres de reads totaux de chaque puit dans le tableau des métadonnées
my_cleaned_otu_table %>% 
  summarise(across(everything(), sum)) %>% 
  t() %>% 
  data.frame() %>% 
  rownames_to_column("sample_ID") %>% 
  rename(total_reads = 2) %>%
  left_join(my_mtd_env, ., by = "sample_ID") -> my_mtd_env

summary(my_mtd_env$total_reads)

#Plein de graphiques A TRIER
ggplot(my_mtd_env, aes(log10(total_reads), fill=grid))+
  geom_histogram()+
  theme_minimal()+
  labs(title="Camargue",
       x ="Total reads", y = "Frequency")

# ggsave("res/Distribution_of_bacterial_reads.pdf")

my_mtd_env <- my_mtd_env %>% arrange(total_reads)
sample_order<-my_mtd_env$sample_ID
label_order<-my_mtd_env$quadrat

ggplot(my_mtd_env, aes(x=sample_ID, y=total_reads, col=grid))+
  geom_point(size=2)+
  theme_minimal()+
  labs(title="Camargue",
       x ="Bacterial reads", y = "Frequency")+
  theme(legend.position="none",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4, size = 4))+
  scale_x_discrete(limits=sample_order, labels=label_order)

my_mtd_env %>%
  group_by(grid) %>% 
  mutate(grid_mean = mean(total_reads)) %>% 
  arrange(grid_mean) -> my_mtd_env

sample_order<-my_mtd_env$sample_ID
label_order<-my_mtd_env$grid

ggplot(my_mtd_env, aes(x=grid, y=total_reads, col=grid))+
  geom_boxplot(alpha = 0.5)+
  geom_point(size=1)+
  geom_point(aes(x=grid, y = grid_mean))+
  theme_minimal()+
  labs(title="Camargue",
       x ="Bacterial reads", y = "Frequency")+
  theme(legend.position="none",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4, size = 4))+
  scale_x_discrete(limits=label_order)

# ggsave("res/RUN4_Distribution_of_total_reads.pdf")

ps_bacteria %>% 
  prune_samples(sample_sums(.) >= 1000, .) %>% 
  tax_fix() %>% 
  comp_barplot(tax_level = "family", n=15)+
  theme_minimal() +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.justification = "top",
        legend.box = "vertical",
        legend.key.size = unit(0.5, "cm"),
        axis.text.x = element_text(size = 8, angle=90, vjust=0.4, hjust=1.2),
        legend.text = element_text(size = 11))

save(ps_bacteria, file = "outputs/ps_bacteria.Rdata")
