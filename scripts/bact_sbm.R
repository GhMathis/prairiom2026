library(tidyverse)
library(sbm)
library(phyloseq)
library(RColorBrewer)
library(alluvial)
library(GGally)
library(vegan)

main_theme = theme_minimal()+
  theme(line = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.ticks =  element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", size=22),
        axis.text.y = element_text(colour = "black", size=22),
        legend.title = element_text(colour = "black", size=20,
                                    hjust =0.5),
        legend.text = element_text(colour = "black", size=18),
        axis.title= element_text(size=28),
        strip.text = element_text(colour = "black", size=15, face ="italic"))

plot.OTUcumSums <-function(x){plot(1:length(colSums(x)), log10(sort(colSums(x))+1))}

# Load data ----
load(file = "outputs/clean_data_to_analise/ps_bacteria.Rdata")

read.table("outputs/clean_data_to_analise/Metadata_grid_CAM.txt", header = T) -> metadata_grid


metadata_grid %>% 
  dplyr::select(log_biomass, Vegetation)%>% 
  ggpairs()
# N read qaudrat ----

ps_bacteria %>%
  otu_table() %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "quadrat_code") %>%
  group_by(quadrat_code) %>%
  summarise(tot_reads= sum(value)) %>%
  mutate(grid_code = str_sub(quadrat_code,1,9)) %>% 
  ggplot() +
  facet_wrap(~grid_code, scale="free_x") +
  geom_point(aes(quadrat_code, tot_reads), cex = 3) +
  scale_y_log10() +
  theme_minimal() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Bacteria community ---

ps_bacteria %>%
  otu_table() %>%
  t() %>%
  plot.OTUcumSums

ps_bacteria_grid <- merge_samples(ps_bacteria, "grid")

ps_bacteria_grid %>%
  otu_table() %>% 
  t() %>%
  as.data.frame() -> bacteria_grid_df

ps_bacteria_grid %>% 
  otu_table() %>% 
  colSums() %>%
  sort(decreasing = T)
ps_bacteria_grid %>%
  otu_table() %>%
  t() %>%  
  plot.OTUcumSums

## Find sur-abundante bacteria -----

### 10 max abundances per grid -----

bacteria_grid_df %>%
  rownames_to_column("OTU") %>%
  pivot_longer(-OTU) %>%
  # 1. Calculer le total par échantillon (name)
  group_by(name) %>%
  mutate(tot_read = sum(value)) %>%
  # 2. IMPORTANT : Dégrouper pour que la réorganisation soit globale
  ungroup() %>%
  # 3. Réordonner le facteur 'name' selon 'tot_read'
  mutate(name = fct_reorder(name, tot_read)) %>%
  # 4. Reprendre le groupement si vous voulez les 10 meilleurs OTU par échantillon
  group_by(name) %>%
  slice_max(order_by = value, n = 10) -> highest_read_per_grid
  
highest_read_per_grid %>%
  ggplot() +
  geom_label(aes(x = name, y = value,label = OTU, col = tot_read)) +
  # Optionnel : inverser les axes si vous avez beaucoup de noms
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 10 OTUs par échantillon (ordonnés par reads totaux)",
       x = "Échantillon (Name)",
       y = "Abondance (Value)",
       color = "Total Reads")

### 10 max abundances per quadra 21_CAM_10 -----

bacteria_grid_df %>%
  rownames_to_column("OTU") %>%
  pivot_longer(-OTU) %>%
  filter(str_detect(name, "21-CAM-10")) %>%
  # 1. Calculer le total par échantillon (name)
  group_by(name) %>%
  mutate(tot_read = sum(value)) %>%
  # 2. IMPORTANT : Dégrouper pour que la réorganisation soit globale
  ungroup() %>%
  # 3. Réordonner le facteur 'name' selon 'tot_read'
  mutate(name = fct_reorder(name, tot_read)) %>%
  # 4. Reprendre le groupement si vous voulez les 10 meilleurs OTU par échantillon
  group_by(name) %>%
  slice_max(order_by = value, n = 10) -> highest_read_grid_21_10

highest_read_grid_21_10 %>%
  ggplot() +
  geom_label(aes(x = name, y = value,label = OTU, col = tot_read)) +
  # Optionnel : inverser les axes si vous avez beaucoup de noms
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 10 OTUs par quadara (ordonnés par reads totaux)",
       x = "Quadara",
       y = "Abondance (Value)",
       color = "Total Reads")

highest_read_per_grid %>%
  arrange(desc(value))

ps_bacteria %>%
  tax_table() %>%
  as.data.frame() %>%
  rownames_to_column("OTU") %>%
  filter(OTU %in% c(9, 14, 17, 29, 35))

illumina_data = "data/bacteria/Camargue_diversity_MiSeq_16S_341F_785R_20240417_405_samples.OTU.filtered.cleaved.nosubstringOTUs.mumu.table2"
illumina_data %>%
  read_tsv(show_col_types = FALSE, progress = FALSE) %>%
  filter(OTU %in% c(9, 14, 17, 29, 35)) %>%
  select(OTU, amplicon)
  
illumina_data %>%
  read_tsv(show_col_types = FALSE, progress = FALSE) %>%
  filter(OTU %in% c(9, 14, 17, 29, 35)) %>%
  select(OTU, sequence)

illumina_data %>%
  read_tsv(show_col_types = FALSE, progress = FALSE) %>%
  filter(str_detect(taxonomy,"Pantoea" ))


bacteria_grid_df %>%
  pivot_longer(everything(), names_to = "sample_ID", values_to = "reads") %>%
  group_by(sample_ID) %>%
  summarise(reads = sum(reads)) %>%
  mutate(sample_ID = fct_reorder(sample_ID, reads)) %>%
  ggplot() +
  geom_point(aes(sample_ID, reads)) 


### Rarefation -----

set.seed(4892763)
bacteria_grid_df %>%
  rownames_to_column("OTU") %>%
  pivot_longer(-OTU, names_to = "sample_ID", values_to = "reads") %>%
  pivot_wider(names_from = "sample_ID", values_from = "reads") %>%
  column_to_rownames("OTU") %>% 
  t() %>%
  vegan::rrarefy(x=.,sample = min(rowSums(.))) %>%
  t() %>%
  as.data.frame() %>%
  filter(rowSums(.)>0)-> bacteria_grid_rare

str(bacteria_grid_rare)
dim(bacteria_grid_rare)
dim(bacteria_grid_df)



bacteria_grid_df %>%
  t() %>%
  rarecurve(x=., step = 500, sample = min(rowSums(.))) -> bacteria_grid_rarecurve

# SBM ----

## Data visualisation ----

bacteria_grid_rare %>%
  rownames_to_column("OTU") %>%
  pivot_longer(-OTU, names_to = "grid_code", values_to = "reads") %>%
  ggplot() +
  facet_wrap(~grid_code) +
  geom_histogram(aes(log1p(reads)))

# bacteria_grid_rare %>%
#   as.matrix() %>%
#   log1p() %>%
#   plotMyMatrix()

bacteria_grid_rare %>%
  rownames_to_column("OTU") %>%
  pivot_longer(-OTU, names_to = "grid_code", values_to = "reads") %>%
  filter(reads > 10) %>%
  pivot_wider(names_from = "grid_code", values_from = "reads", values_fill = 0) %>%
  column_to_rownames("OTU") -> bacteria_grid_rare_filter

dim(bacteria_grid_rare)
bacteria_grid_rare_filter %>%
  as.matrix() %>%
  t() %>%
  log1p %>%
  plotMyMatrix()

bacteria_grid_rare_filter %>%
  as.matrix() %>%
  t() %>%
  log1p %>%
  hist

bacteria_grid_rare_filter %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("OTU") %>%
  pivot_longer(-OTU, names_to = "grid_code", values_to = "reads") %>%
  filter(reads > 0) %>%
  pull(reads) %>%
  hist

bacteria_grid_rare_filter %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("OTU") %>%
  pivot_longer(-OTU, names_to = "grid_code", values_to = "reads") %>%
  filter(reads > 0) %>%
  pull(reads) %>%
  log() %>%
  hist

bacteria_grid_rare %>%
  as.data.frame() %>%
  rownames_to_column("OTU") %>%
  pivot_longer(-OTU, names_to = "grid_code", values_to = "reads") %>%
  filter(reads > 0) %>%
  pull(reads) %>%
  log() %>%
  hist

## Matrice subsampling ----

### Sub-sampling function ----
subsample.equally <- function(x, n_subsample=5){
  # All position which are not 0
  pos <- which(x > 0)
  
  # Incase non 0 are less represented than n_subsample
  n_to_draw <- min(n_subsample, length(pos))
  
  # Sub sample across all non 0 abundance
  kept <- sample(pos, n_to_draw)
  
  # Put subsampled OTU/ASC n reads in a vector with the correct order 
  replace(numeric(length(x)), kept, x[kept])
}
subsample.equally.rowise <- function(x, n_subsample=5){
  S <- rownames(x)
  kept <- sample(S,size = n_subsample,replace = F)
  
  x[kept,]
}

set.seed(1234)
bacteria_grid_subsample <- bacteria_grid_df %>%
  mutate(across(everything(),
                ~ subsample.equally(.x, n = 500)))  %>%
  filter(rowSums(.)>0)

bacteria_grid_subsample %>% colSums()

### Compuation ----

library(furrr)

#### Random subsample of abundance matrix and SBM ----

if(file.exists("outputs/heavy_computation_save/subsample_sbm_bact.Rdata")){
  
  load(file = "outputs/heavy_computation_save/subsample_list_bact.Rdata")
  load(file = "outputs/heavy_computation_save/subsample_sbm_bact.Rdata")

  }else{
  
  n_sim = 1000
  set.seed(1234)

  purrr::map(1:n_sim, \(x) subsample.equally.rowise(bacteria_grid_df, n_subsample = 500)) -> subsample_list
  save(subsample_list, file = "outputs/heavy_computation_save/subsample_list_bact.Rdata")
  
  plan(multisession, workers = 7)
  subsample_list %>%
    furrr::future_map(\(x) x%>% 
                        as.matrix() %>%
                        log1p %>% 
                        estimateBipartiteSBM(netMat= ., model = "poisson", dimLabels = c(row = "ASV", col = "grid_code"),
                                                                                         estimOptions = list(plot = FALSE,
                                                                                                             verbosity = 0)),
                      .options = furrr_options(scheduling = 2),
                      .progress = T
  ) -> subsample_sbm_list

save(subsample_sbm_list, file = "outputs/heavy_computation_save/subsample_sbm_bact.Rdata")
  }

#### Random subsample of occurence matrix and SBM ----

if(file.exists("outputs/heavy_computation_save/subsample_sbm_bact_binar.Rdata")){
  
  load(file = "outputs/heavy_computation_save/subsample_sbm_bact.Rdata")
  
  }else{
  
  plan(multisession, workers = 7)
  subsample_list %>%
    furrr::future_map(\(x) x %>% 
                        mutate(across(where(is.numeric), ~ as.integer(.x > 0))) %>%
                        as.matrix() %>%
                        estimateBipartiteSBM(netMat= ., model = "bernoulli", 
                                             dimLabels = c(row = "ASV", col = "grid_code"),
                                             estimOptions = list(plot = FALSE,
                                                                 verbosity = 0)),
                      .options = furrr_options(scheduling = 2),
                      .progress = T
    ) -> subsample_sbm_binar_list
  
  save(subsample_sbm_binar_list, file = "outputs/heavy_computation_save/subsample_sbm_bact_binar.Rdata")
  }

subsample_list[[1]] %>% 
  mutate(across(where(is.numeric), ~ as.integer(.x > 0))) %>%
  as.matrix() %>%
  estimateBipartiteSBM(netMat= ., model = "bernoulli", dimLabels = c(row = "ASV", col = "grid_code")) -> test_sbm_binary

plot(test_sbm_binary)
test_sbm_binary$connectParam

subsample_sbm_binar_list[[1]] %>% plot
subsample_sbm_binar_list[[2]] %>% plot
subsample_sbm_binar_list[[3]] %>% plot

subsample_sbm_binar_list[[1]]$memberships
subsample_sbm_binar_list[[2]]$memberships
subsample_sbm_binar_list[[3]]$memberships

subsample_sbm_binar_list[[1]]$connectParam
subsample_sbm_binar_list[[2]]$connectParam
subsample_sbm_binar_list[[3]]$connectParam

#### Agglomération at genus level ----

ps_bacteria_genus_grid <- ps_bacteria_grid %>%
  tax_glom(taxrank = "genus")

ps_bacteria_genus_grid %>%
  tax_table() %>%
  head

ps_bacteria_genus_grid %>%
  otu_table() %>%
  t() %>%
  rowSums() %>%
  sort

ps_bacteria_genus_grid %>%
  otu_table() %>%
  t() %>%
  as.vector() %>%
  hist

ps_bacteria_genus_grid %>%
  otu_table() %>%
  log1p %>% 
  plotMyMatrix()

ps_bacteria_genus_grid %>% 
  otu_table() %>%
  t() %>%
  as.data.frame() %>%
  mutate(across(where(is.numeric), ~ as.integer(.x > 0))) %>%
  as.matrix() %>%
  estimateBipartiteSBM(netMat= ., model = "bernoulli", dimLabels = c(row = "ASV", col = "grid_code")) -> smb_bin_genus

plot(smb_bin_genus)
smb_bin_genus$memberships$grid_code %>% table
smb_bin_genus$connectParam
smb_bin_genus$memberships$ASV %in% c(9,10) -> rare_genus

ps_bacteria_genus_grid %>%
  tax_table() %>%
  as.data.frame() %>%
  rownames_to_column("tax_id") %>%
  filter(rare_genus)
ps_bacteria_genus_grid %>%
  otu_table() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("tax_id") %>%
  filter(rare_genus) %>%
  column_to_rownames("tax_id") %>%
  rowSums() %>% sort
ps_bacteria_genus_grid %>%
  tax_table() %>%
  as.data.frame() %>%
  rownames_to_column("tax_id") %>%
  filter(tax_id == 312)
ps_bacteria_genus_grid %>%
  otu_table() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("tax_id") %>%
  filter(tax_id == 75)

ps_bacteria_genus_grid %>%
  otu_table() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("tax_id") %>%
  filter(rare_genus) %>%
  column_to_rownames("tax_id") %>%
  as.matrix() %>%
  log1p %>% plotMyMatrix()

ps_bacteria_genus_grid %>%
  otu_table() %>%
  t() %>%
  colSums()
#### Multiplex and abundance are way too long to run

# ps_bacteria_genus_grid %>%
#   otu_table() %>%
#   t() %>%
#   estimateBipartiteSBM(netMat= ., model = "poisson", dimLabels = c(row = "ASV", col = "grid_code")) -> smb_poisson_genus
# 
# ps_bacteria_genus_grid %>% 
#   otu_table() %>%
#   t() %>%
#   as.data.frame() %>%
#   mutate(across(where(is.numeric), ~ as.integer(.x > 0))) %>%
#   as.matrix() %>%
#   defineSBM(model = "bernoulli", dimLabels = c(row = "ASV", col = "grid_code")) -> net_bin
# 
# ps_bacteria_genus_grid %>%
#   otu_table() %>%
#   t() %>%
#   defineSBM(model = "poisson", dimLabels = c(row = "ASV", col = "grid_code")) -> net_pois
# (net_bin)
# plotMyMultiplexMatrix(list(net_bin, net_pois))
# MultiplexFitIndep <- estimateMultiplexSBM(list(net_bin, net_pois), dependent = F)
