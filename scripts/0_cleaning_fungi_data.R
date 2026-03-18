#Packages----
library(tidyverse)
library(phyloseq)
library(microbiome)
library(microViz)


#Load data----
input_table <- "data/fungi/Camargue_diversity_MiSeq_ITS2_ITS86F_ITS4r_20240220_818_samples.OTU.filtered.cleaved.nosubstringOTUs.mumu.table2"
unite_taxonomy_table <- "data/fungi/unite_ITS86F_20230725_modified_taxonomy.table"
taxonomic_levels <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

#Taxonomy----

##Deal with multiple references----
#Seuil de longueur
LSEUIL = 0

###Functions
. %>%
  read_tsv(show_col_types = FALSE) %>%
  filter(length>=LSEUIL) %>% 
  select(OTU, references) %>%
  separate_rows(references, sep = ",") %>%
  left_join(reference_to_taxonomy, by = "references") %>%
  separate(taxonomy, into = taxonomic_levels, sep = "[|]") -> prepare_table

. %>%
  add_count(OTU, kingdom, name = "n_kingdom") %>%
  add_count(OTU, kingdom, phylum, name = "n_phylum") %>%
  add_count(OTU, kingdom, phylum, class, name = "n_class") %>%
  add_count(OTU, kingdom, phylum, class, order, name = "n_order") %>%
  add_count(OTU, kingdom, phylum, class, order, family, name = "n_family") %>%
  add_count(OTU, kingdom, phylum, class, order, family, genus, name = "n_genus") %>%
  add_count(OTU, kingdom, phylum, class, order, family, genus, species,
            name = "n_species") -> count_taxonomic_assignments_per_OTU

. %>%
  group_by(OTU) %>%
  mutate(number_of_references = n(),
         kingdom_frequency = 100 * n_kingdom / number_of_references,
         phylum_frequency = 100 * n_phylum / number_of_references,
         class_frequency = 100 * n_class / number_of_references,
         order_frequency = 100 * n_order / number_of_references,
         family_frequency = 100 * n_family / number_of_references,
         genus_frequency = 100 * n_genus / number_of_references,
         phylum_frequency = 100 * n_phylum / number_of_references,
         species_frequency = 100 * n_species / number_of_references) %>%
  ungroup() -> compute_taxonomic_assignment_frequencies

. %>%
  select(! starts_with("n_")) -> clean_up

input_table %>%
  read_tsv(show_col_types = FALSE)  %>% 
  filter(str_detect(taxonomy, "k__Fungi"))-> tmp2
  
###Main
#Seuil de certitude sur assignations multiples
P_SEUIL = 75

#Taxonomie correspondant à chaque référence
unite_taxonomy_table %>%
  read_tsv(show_col_types = FALSE) -> reference_to_taxonomy

#Comptage du nb de refs par OTU et de la fréquence de chaque taxo
input_table %>%
  prepare_table %>%
  count_taxonomic_assignments_per_OTU %>% 
  compute_taxonomic_assignment_frequencies %>%
  clean_up -> assignments

assignments %>% 
  #filter(number_of_references>1) %>% 
  mutate(kingdom = paste(kingdom, kingdom_frequency, sep="%"), 
         phylum = paste(phylum, phylum_frequency, sep="%"), 
         class = paste(class, class_frequency, sep="%"), 
         order = paste(order, order_frequency, sep="%"), 
         family = paste(family, family_frequency, sep="%"), 
         genus = paste(genus, genus_frequency, sep="%"), 
         species = paste(species, species_frequency, sep="%")) %>% 
  select(-references, -contains("frequency"), -number_of_references) %>% 
  pivot_longer(., cols = c(kingdom, phylum, class, order, family, genus, species), names_to = "tax_level", values_to = "tax_name") %>% 
  mutate(tax_name = if_else(str_detect(tax_name, "NA%100"), tax_name, str_sub(tax_name, 4, nchar(tax_name)))) %>% 
  separate_wider_delim(cols = "tax_name", delim = "%", names = c("taxa", "freq")) %>% 
  mutate(freq = as.numeric(freq)) %>% 
  filter(freq >= P_SEUIL)  %>% 
  select(-freq) %>% 
  distinct() %>% 
  pivot_wider(names_from = "tax_level", values_from = "taxa", values_fill = NA) -> new_assignment_table

##Cleaning taxa names----
#Avec ce script les chemins avec kingdom = NA 
new_assignment_table %>%
  mutate(kingdom = if_else(str_detect(kingdom, "Incertae_sedis"), NA, kingdom), 
         phylum = if_else(str_detect(phylum, "Incertae_sedis"), NA, phylum),
         class = if_else(str_detect(class, "Incertae_sedis"), NA, class), 
         order = if_else(str_detect(order, "Incertae_sedis"), NA, order), 
         family = if_else(str_detect(family, "Incertae_sedis"), NA, family), 
         genus = if_else(str_detect(genus, "Incertae_sedis"), NA, genus), 
         species = if_else(str_detect(species, "Incertae_sedis")|str_ends(species, "_sp"), NA, species)) %>%
  mutate(phylum = if_else(is.na(kingdom), NA, phylum),
         class = if_else(is.na(kingdom)|is.na(phylum), NA, class), 
         order = if_else(is.na(kingdom)|is.na(phylum)|is.na(class), NA, order), 
         family = if_else(is.na(kingdom)|is.na(phylum)|is.na(class)|is.na(order), NA, family), 
         genus = if_else(is.na(kingdom)|is.na(phylum)|is.na(class)|is.na(order)|is.na(family), NA, genus), 
         species = if_else(is.na(kingdom)|is.na(phylum)|is.na(class)|is.na(order)|is.na(family)|is.na(genus), NA, species), NA, species) %>% 
  mutate(phylum = if_else(is.na(phylum), paste0("unidentified_", kingdom), phylum),
         class = if_else(is.na(class), paste0("unidentified_", phylum), class), 
         order = if_else(is.na(order), paste0("unidentified_", class), order), 
         family = if_else(is.na(family), paste0("unidentified_", order), family), 
         genus = if_else(is.na(genus), paste0("unidentified_", family), genus), 
         species = if_else(is.na(species), paste0("unidentified_", genus), species)) %>%  
  mutate(phylum = str_replace_all(phylum, "unidentified_unidentified_","unidentified_"), 
         class = str_replace_all(class, "unidentified_unidentified_","unidentified_"),
         order = str_replace_all(order, "unidentified_unidentified_","unidentified_"),
         family = str_replace_all(family, "unidentified_unidentified_","unidentified_"),
         genus = str_replace_all(genus, "unidentified_unidentified_","unidentified_"),
         species = str_replace_all(species, "unidentified_unidentified_","unidentified_")) %>% 
  mutate(phylum = str_replace_all(phylum, "unidentified_unidentified_","unidentified_"), 
         class = str_replace_all(class, "unidentified_unidentified_","unidentified_"), 
         order = str_replace_all(order, "unidentified_unidentified_","unidentified_"), 
         family = str_replace_all(family, "unidentified_unidentified_","unidentified_"), 
         genus = str_replace_all(genus, "unidentified_unidentified_","unidentified_"), 
         species = str_replace_all(species, "unidentified_unidentified_","unidentified_")) %>% 
  mutate(phylum = str_replace_all(phylum, "unidentified_unidentified_","unidentified_"), 
         class = str_replace_all(class, "unidentified_unidentified_","unidentified_"), 
         order = str_replace_all(order, "unidentified_unidentified_","unidentified_"), 
         family = str_replace_all(family, "unidentified_unidentified_","unidentified_"), 
         genus = str_replace_all(genus, "unidentified_unidentified_","unidentified_"), 
         species = str_replace_all(species, "unidentified_unidentified_","unidentified_")) -> my_tax

my_tax %>% 
  column_to_rownames("OTU") %>% 
  select(-"NA") -> my_tax

my_tax %>% 
  filter(! is.na(kingdom)) -> my_tax

#OTU table----
#Seuil d'identité
I_SEUIL = 0.80

read_tsv(input_table, show_col_types = FALSE) %>% 
  filter(identity >= I_SEUIL) %>% 
  select(-c("total", "cloud", "amplicon", "length", "abundance", "chimera", "spread", "quality", "sequence", "identity", "taxonomy", "references")) %>% 
  column_to_rownames("OTU") -> my_otu
  #filter(!rownames(.) %in% to_delete) 

#Metadata----
as_tibble(colnames(my_otu)) %>% 
  rename(full_sample_ID = value) %>% 
  separate_wider_delim(cols = full_sample_ID, delim = "-Pl", cols_remove = F, names = c("sample_ID", "plate")) %>% 
  separate_wider_delim(cols = plate, delim = "-R", names = c("plate", "run")) %>% 
  separate_wider_delim(cols = run, delim = "_S", names = c("run", "library")) %>% 
  select(-library) %>% 
  #filter(str_starts(sample_ID, "T")) %>% 
  mutate(type = case_when(str_detect(sample_ID, "Tneg") ~"NEG",
                          str_detect(sample_ID, "Zymo") ~"ZYMO",
                          str_detect(sample_ID, "genseq") ~"NEG",
                          .default = "ENV"))  %>%
  mutate(row_ID = full_sample_ID) %>% 
  column_to_rownames("row_ID")  -> my_mtd

#PHYLOSEQ----
library(phyloseq)

sample_data(my_mtd) -> MTD

my_tax %>% 
  filter(kingdom == "Fungi") %>% 
  filter(rownames(.) %in% rownames(my_otu)) -> my_tax

my_otu %>% 
  filter(rownames(.) %in% rownames(my_tax)) -> my_otu

# otu_table(my_otu, taxa_are_rows = T) -> OTU
# 
# tax_table(as.matrix(my_tax)) -> TAX
# 
# ps <- phyloseq(OTU, TAX, MTD)
# 
# taxa_sums(ps) %>% 
#   log10() %>% 
#   hist()
# 
# taxa_sums(ps) %>% 
#   summary()
# 
# sample_sums(ps) %>% hist()
# 
# my_otu %>% 
#   colSums() %>% 
#   as.data.frame() %>% 
#   rename(sample_read_count = ".") %>% 
#   filter(!sample_read_count==0)

#*******************----
#MicroDecon----
# install.packages("remotes")
# remotes::install_github("donaldtmcknight/microDecon")
library("microDecon")

my_mtd %>% 
  filter(type=="NEG") %>% 
  pull(full_sample_ID) -> neg_samples

my_mtd %>% 
  filter(type=="ENV") %>% 
  pull(full_sample_ID) -> env_samples

neg <- my_otu[,neg_samples]
env <- my_otu[, env_samples]

data.frame(OTU = rownames(my_otu)) %>% 
  cbind(., neg) %>% 
  cbind.data.frame(., env) -> df_decon

decontaminated <- decon(data = df_decon,numb.blanks=length(colnames(neg)),numb.ind=length(colnames(env)), taxa=F) 

#**OTUs à retirer----
decontaminated$OTUs.removed %>% 
  pull(OTU)-> otus_to_remove
length(otus_to_remove) #il y en a 23

my_tax %>%
  filter(rownames(.) %in% otus_to_remove) %>% view()

#**OTUs seulement présentes dans les puits négatifs----
my_otu %>% 
  mutate(.before = 1, total_env = rowSums(.[,env_samples])) %>%
  filter(total_env==0) %>% 
  rownames(.)-> OTU_absent_from_env_samples

length(OTU_absent_from_env_samples) #Il y en a 4

my_tax %>%
  filter(rownames(.) %in% OTU_absent_from_env_samples) %>% view()

#**OTUs dont les reads sont corrigés----
decontaminated$reads.removed %>%
  replace(.,is.na(.),0) %>%
  filter(!OTU %in% otus_to_remove) %>% 
  filter(!OTU %in% OTU_absent_from_env_samples) %>% 
  pull(OTU) -> partially_removed 
length(partially_removed) #Il y en a 21

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
                      TRUE ~ "True OTU"))

ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + 
  geom_point(position = "jitter") +
  xlab("Prevalence (Negative Controls)") + 
  ylab("Prevalence (True Samples)")  + 
  ggtitle("Conserved and removed OTUs following microDecon")

ggsave("decontamination_fungi.pdf", width = 15, units = "cm", scale =1.618)

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

library(ggrepel)

ggplot(data=df.pa1, aes(x=pa.neg1, y=pa.pos1, color=contaminant)) + 
  geom_point(position = "jitter")+
  geom_label_repel(aes(label = OTU,  color = contaminant), size = 3)+
  xlab("Number of reads in negative Controls (Log scale))") + 
  ylab("Number of reads in true Samples (Log-scale)")  + 
  ggtitle("Decontamination following microDecon")

ggsave("decontamination_2.pdf", width = 15, units = "cm", scale =1.618)

#*******************----
#Objets finaux----
#**Table d'OTU décontaminée----
decontaminated$decon.table %>% 
  mutate(total = rowSums(.[,3:ncol(.)])) %>%
  filter(total > 0) %>%
  select(-total,-Mean.blank) %>% 
  column_to_rownames(var = "OTU") -> my_cleaned_otu_table

my_cleaned_otu_table %>% 
  summarise(across(everything(), sum)) %>% 
  t() %>% as.data.frame() %>% 
  pull(1) -> total_reads_per_sample

hist(total_reads_per_sample)
summary(total_reads_per_sample)

my_tax %>% 
  filter(rownames(.) %in% rownames(my_cleaned_otu_table)) -> my_cleaned_tax_table

#clearNames(ps.pa, ps.pa.neg, ps.pa.pos, ps.decon, ps.decon.neg, ps.decon.pos)

my_mtd %>% 
  filter(type=="ENV") -> my_mtd_env 
rownames(my_mtd_env) <- my_mtd_env$sample_ID
sample_data(my_mtd_env) -> MTD

my_cleaned_otu_table %>%  
  rename_with(~str_sub(., start = 1, end = 13)) %>% 
  otu_table(.,taxa_are_rows = T) -> OTU

tax_table(as.matrix(my_cleaned_tax_table)) -> TAX

sample_names(MTD)

ps_fungi = phyloseq(OTU, TAX, MTD)

#TRANFERT pour MATHIS----
ps_fungi %>% 
  tax_fix() %>% 
  comp_barplot(tax_level = "class", n=30)+
  theme_minimal() +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.justification = "top",
        legend.box = "vertical",
        legend.key.size = unit(0.5, "cm"),
        axis.text.x = element_text(size = 8, angle=90, vjust=0.4, hjust=1.2),
        legend.text = element_text(size = 11))


save(ps_fungi, file = "outputs/clean_data_to_analise/ps_fungi.Rdata")

library(vegan)

rarecurve(fungi_OTU, step = 50, sample = 1000, label = T, tidy=T) %>% 
  ggplot(., aes(x = Sample, y = Species, group = Site))+
  geom_line()+
  xlim(0,1000)
  
  