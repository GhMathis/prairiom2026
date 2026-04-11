# Load pacakges -----

library(tidyverse)
library(sbm)
library(phyloseq)
library(iNEXT.3D)

# SBM at grid level ----

## Bacteria -----

load(file = "outputs/clean_data_to_analise/ps_bacteria.Rdata")

ps_bacteria_grid <- merge_samples(ps_bacteria, "grid") %>%
  prune_samples(!str_detect(sample_names(.), "21-CAM-1[45]"), .)


ps_bacteria_grid %>%
  otu_table() %>% 
  t() %>%
  as.data.frame() %>%
  rename_with(~str_replace_all(., "-","_")) -> bacteria_grid_matrix

bacteria_grid_matrix %>%
  rownames_to_column("OTU") %>%
  pivot_longer(-OTU, names_to = "grid_code", values_to = "reads") %>%
  mutate(reads = floor(reads/10)) %>%
  filter(reads > 0) %>%
  pivot_wider(names_from = "grid_code", values_from = "reads", values_fill = 0) %>%
  column_to_rownames("OTU") -> bacteria_grid_truncated_matrix

reference_sample_sort <- colnames(bacteria_grid_truncated_matrix)

dim(bacteria_grid_matrix)
dim(bacteria_grid_truncated_matrix)
## Fungi -----

load(file = "outputs/clean_data_to_analise/ps_fungi.Rdata")

ps_fungi_grid <- merge_samples(ps_fungi, "grid_code") %>%
  prune_samples(!str_detect(sample_names(.), "21_CAM_1[45]"), .)

ps_fungi_grid %>%
  otu_table() %>% 
  t() %>%
  as.data.frame() -> fungi_grid_matrix

fungi_grid_matrix <- fungi_grid_matrix[,reference_sample_sort]

fungi_grid_matrix %>%
  rownames_to_column("OTU") %>%
  pivot_longer(-OTU, names_to = "grid_code", values_to = "reads") %>%
  mutate(reads = floor(reads/10)) %>%
  filter(reads > 0) %>%
  pivot_wider(names_from = "grid_code", values_from = "reads", values_fill = 0) %>%
  column_to_rownames("OTU") -> fungi_grid_truncated_matrix

dim(fungi_grid_matrix)
dim(fungi_grid_truncated_matrix)

## Virus -----

read.table(file = "outputs/clean_data_to_analise/OTU_virus_CAM.txt", header = T) -> virus_quadra_matrix
virus_quadra_matrix %>%
  mutate(grid_code = str_sub(Host_code,1,9)) %>%
  filter(!str_detect(grid_code, "21_CAM_1[45]")) %>%
  select(-Host_code) %>%
  group_by(grid_code) %>%
  summarise(across(where(is.numeric), ~sum(.x))) %>%
  column_to_rownames("grid_code") %>%
  as.matrix() -> virus_grid_matrix

virus_grid_matrix <- virus_grid_matrix[reference_sample_sort,]

dim(virus_grid_matrix)

## Plant -----

read.table(file = "outputs/clean_data_to_analise/abund_plant_grid.txt", header = T) %>%
  column_to_rownames("Grid_code") %>%
  as.matrix() -> plant_grid_matrix

plant_grid_matrix <- plant_grid_matrix[reference_sample_sort,]

dim(plant_grid_matrix)






## Metadata -----

read.table("outputs/clean_data_to_analise/Metadata_grid_CAM.txt", header = T) -> metadata_grid

dim(metadata_grid)
# Rarefaction at reference coverage ----

## Compute sample coverage ----

# Bacteria
iNext_object_OTU_bact <- bacteria_grid_truncated_matrix %>%
  iNEXT.3D::iNEXT3D( q = c(0), datatype = "abundance", nboot = 0)

iNext_object_OTU_bact$TDiNextEst$coverage_based %>% 
  filter(Order.q == 0 & Method == "Observed") %>%
  filter(SC == min(SC, na.rm = TRUE)) %>%
  pull(SC) -> min_coverage_OTU_bact

bacteria_grid_truncated_matrix %>%
  iNEXT.3D::estimate3D(level = min_coverage_OTU_bact, q = 0, nboot = 0) %>%
  select(Assemblage,m) -> required_sizes_OTU_bact

# Fungi
iNext_object_OTU_fungi <- fungi_grid_truncated_matrix %>%
  iNEXT.3D::iNEXT3D( q = c(0), datatype = "abundance", nboot = 0)

iNext_object_OTU_fungi$TDiNextEst$coverage_based %>% 
  filter(Order.q == 0 & Method == "Observed") %>%
  filter(SC == min(SC, na.rm = TRUE)) %>%
  pull(SC) -> min_coverage_OTU_fungi

fungi_grid_truncated_matrix %>%
  iNEXT.3D::estimate3D(level = min_coverage_OTU_fungi, q = 0, nboot = 0) %>%
  select(Assemblage,m) -> required_sizes_OTU_fungi

# Virus
rowSums(virus_grid_matrix)
rbind(n = rep(9,ncol(t(virus_grid_matrix))),t(virus_grid_matrix)) %>% 
  as.data.frame() -> virus_grid_matrix_Inext 
str(virus_grid_matrix_Inext)

# library(iNEXT)
# dim(virus_grid_matrix_Inext)
# iNEXT::iNEXT(virus_grid_matrix_Inext, q=c(0,1), datatype="incidence_freq", se=T) -> iNext_object_OTU_virus
# 
# iNext_object_OTU_virus$TDiNextEst$coverage_based %>% filter(Method == "Observed")
# iNext_object_OTU_virus$TDiNextEst$coverage_based$Assemblage %>% unique %>% length()
# iNext_object_OTU_virus <- virus_quadra_matrix %>%
#   mutate(grid_code = str_sub(Host_code, 1,9)) %>%
#   pivot_longer(-c(Host_code, grid_code)) %>%
#   group_by(grid_code) %>%
#   filter(sum(value)>5) %>%
#   ungroup() %>%
#   pivot_wider(names_from = "name", values_from = "value") %>%
#   arrange(Host_code) %>%
#   select(-c(Host_code,grid_code)) %>%
#   as.matrix() %>%
#   t() %>%
#   iNEXT.3D::iNEXT3D(q = c(0), datatype = "incidence_raw", nT = rep(9,27), nboot = 0)

# iNext_object_OTU_virus$TDiNextEst$coverage_based %>% 
#   filter(Order.q == 0 & Method == "Observed") %>%
#   filter(SC == min(SC, na.rm = TRUE)) %>%
#   pull(SC) -> min_coverage_OTU_virus
# virus_grid_matrix %>%
#   t() %>%
#   iNEXT.3D::estimate3D(level = min_coverage_OTU_fungi, q = 0, nboot = 0) %>%
#   select(Assemblage,m) -> required_sizes_OTU_fungi

required_sizes_OTU_fungi %>%
  left_join(required_sizes_OTU_bact, by = "Assemblage") %>%
  ggplot() +
  geom_point(aes(m.x,m.y)) +
  geom_smooth(aes(m.x,m.y), method = "glm")

## Rarefaction to sample coverage reference  ----

# Function
coverage_rarefy_matrix <- function(otu_matrix, 
                                   required_sizes,
                                   drop_lowcoverage = FALSE){
  
  # Calculer le coverage de chaque échantillon
  sample_coverages <- apply(otu_matrix, 2, function(sample_vector) {
    sample_vector <- sample_vector[sample_vector > 0]
    n <- sum(sample_vector)
    f1 <- sum(sample_vector == 1)
    f2 <- sum(sample_vector == 2)
    # 
    # if (f1 > 0 && f2 > 0) {
    #   1 - (f1 / n) * (((n - 1) * f1) / ((n - 1) * f1 + 2 * f2))
    # } else if (f1 > 0 && f2 == 0) {
    #   1 - (f1 / n) * (((n - 1) * f1) / ((n - 1) * f1 + 2))
    # } else {
    #   1
    # }
  })
  
  # Raréfier chaque échantillon
  rarefied_matrix <- sapply(colnames(otu_matrix), function(sample_name) {
    sample_vector <- otu_matrix[, sample_name]
    target_size <- required_sizes %>%
      filter(Assemblage == sample_name) %>%
      pull(m)
    
    # Créer un vecteur répété de chaque OTU selon son abondance
    otu_pool <- rep(rownames(otu_matrix), times = sample_vector)
    
    # Échantillonner sans remise
    sampled_otus <- sample(otu_pool, size = target_size, replace = FALSE)
    
    # Compter les abondances
    rarefied <- table(factor(sampled_otus, levels = rownames(otu_matrix)))
    
    return(as.numeric(rarefied))
  })
  
  rownames(rarefied_matrix) <- rownames(otu_matrix)
  rarefied_matrix %>%
    as.data.frame() %>%
    filter(rowSums(across(where(is.numeric)))!=0) %>%
    as.matrix() -> rarefied_matrix
  
  return(rarefied_matrix)
}

# Computation

n_sim = 100
set.seed(2365)
1:n_sim %>%
  map(\(x) bacteria_grid_truncated_matrix %>%
        coverage_rarefy_matrix(required_sizes = required_sizes_OTU_bact)) -> coverage_OTU_bact_list

coverage_OTU_bact_list %>%
  map(\(mat) {
    # convert to tibble with OTU and sample as identifiers
    mat %>%
      as.data.frame() %>%
      rownames_to_column("OTU") %>%
      pivot_longer(-OTU, names_to = "sample", values_to = "count") 
  }) %>% 
  reduce(full_join, by = c("OTU", "sample")) %>%   # union of all OTUs x samples
  # sum across all count columns (one per simulation)
  mutate(total = rowSums(across(starts_with("count")), na.rm = TRUE)) %>%
  select(OTU, sample, total) %>%
  pivot_wider(names_from = "sample", values_from = "total") %>%
  column_to_rownames("OTU") %>%
  as.matrix() -> test2

test3<- test2[rowSums(round(test2/n_sim)) >3,]



n_sim = 100
set.seed(2365)
1:n_sim %>%
  map(\(x) fungi_grid_truncated_matrix %>%
        coverage_rarefy_matrix(required_sizes = required_sizes_OTU_fungi)) -> coverage_OTU_fungi_list
coverage_OTU_fungi_list %>%
  map(\(mat) {
    # convert to tibble with OTU and sample as identifiers
    mat %>%
      as.data.frame() %>%
      rownames_to_column("OTU") %>%
      pivot_longer(-OTU, names_to = "sample", values_to = "count") 
  }) %>% 
  reduce(full_join, by = c("OTU", "sample")) %>%   # union of all OTUs x samples
  # sum across all count columns (one per simulation)
  mutate(total = rowSums(across(starts_with("count")), na.rm = TRUE)) %>%
  select(OTU, sample, total) %>%
  pivot_wider(names_from = "sample", values_from = "total") %>%
  column_to_rownames("OTU") %>%
  as.matrix() -> test2_fungi

test3_fungi <- test2_fungi[rowSums(round(test2_fungi/n_sim)) >3,]

colnames(coverage_OTU_bact_list[[1]]) == colnames(coverage_OTU_fungi_list[[1]])
colnames(coverage_OTU_fungi_list[[1]]) == rownames(virus_grid_matrix)
rownames(virus_grid_matrix) == rownames(plant_grid_matrix)

### Bact SBMs ----

# sbm_Bact_bin <- estimateBipartiteSBM(+(t(test3)>0), 
#                                  model = "bernoulli",
#                                  dimLabels = c("grid", "bact"))
# save(sbm_Bact_bin, file = "outputs/heavy_computation_save/sbm_Bact_bin.Rdata")
load(file = "outputs/heavy_computation_save/sbm_Bact_bin.Rdata")
source("scripts/plot_multipartite_sbm.R")
plot_bipartite_sbm(sbm_Bact_bin)

# Doesn't converge, to many groups
# sbm_Bact_poison <- estimateBipartiteSBM(t(test3), 
#                                         model = "poisson",
#                                         dimLabels = c("grid", "bact"))

# sbm_Bact_gaussian <- estimateBipartiteSBM(log1p(t(test3)), 
#                                           model = "gaussian",
#                                           dimLabels = c("grid", "bact"))
# save(sbm_Bact_gaussian, file = "outputs/heavy_computation_save/sbm_Bact_gaussian.Rdata")
load(file = "outputs/heavy_computation_save/sbm_Bact_gaussian.Rdata")

plot_bipartite_sbm(sbm_Bact_gaussian)

### Fungi SBMs ----
str(test3_fungi)

# sbm_Fungi_bin <- estimateBipartiteSBM(+(t(test3_fungi)>0),
#                                  model = "bernoulli",
#                                  dimLabels = c("grid", "fungi"))
# save(sbm_Fungi_bin, file = "outputs/heavy_computation_save/sbm_Fungi_bin.Rdata")
load(file = "outputs/heavy_computation_save/sbm_Fungi_bin.Rdata")
plot_bipartite_sbm(sbm_Fungi_bin)



sbm_Fungi_gaussian <- estimateBipartiteSBM(log1p(t(test3_fungi)),
                       model = "gaussian",
                       dimLabels = c("grid", "fungi"))
save(sbm_Fungi_gaussian, file = "outputs/heavy_computation_save/sbm_Fungi_gaussian.Rdata")
load(file = "outputs/heavy_computation_save/sbm_Fungi_gaussian.Rdata")
plot_bipartite_sbm(sbm_Fungi_gaussian)

sbm_Virus <- estimateBipartiteSBM(virus_grid_matrix,
                       model = "poisson",
                       dimLabels = c("grid", "virus"))

sbm_Plant <- estimateBipartiteSBM(log1p(plant_grid_matrix),
                       model = "gaussian",
                       dimLabels = c("grid", "plant"))

gg_bact<- plot_bipartite_sbm(sbm_Bact_gaussian)
gg_fungi<- plot_bipartite_sbm(sbm_Fungi_gaussian)
gg_virus <- plot_bipartite_sbm(sbm_Virus)
gg_plant <- plot_bipartite_sbm(sbm_Plant)

library(patchwork)

list(gg_bact,
     gg_fungi,
     gg_virus,
     gg_plant) %>%
  wrap_plots(nrow = 1)

gg_bact_expect <- plot_bipartite_sbm(sbm_Bact_gaussian, mat = sbm_Bact_gaussian$expectation)
gg_fungi_expect <- plot_bipartite_sbm(sbm_Fungi_gaussian, mat = sbm_Fungi_gaussian$expectation)
gg_virus_expect <- plot_bipartite_sbm(sbm_Virus, mat = sbm_Virus$expectation)
gg_plant_expect <- plot_bipartite_sbm(sbm_Plant, mat = sbm_Plant$expectation)

list(gg_bact_expect,
     gg_fungi_expect,
     gg_virus_expect,
     gg_plant_expect) %>%
  wrap_plots(nrow = 1)

plot_bipartite_sbm(sbm_Bact_gaussian, mat = sbm_Bact_gaussian$expectation) 

# graph sbms -----

all(rownames(sbm_Bact_gaussian$networkData) == rownames(sbm_Fungi_gaussian$networkData))
all(rownames(sbm_Fungi_gaussian$networkData) == rownames(sbm_Virus$networkData))
all(rownames(sbm_Virus$networkData) == rownames(sbm_Plant$networkData))

B_multi_regne <- as.data.frame(
  table(
    rownames(sbm_Bact_gaussian$networkData),
    sbm_Bact_gaussian$memberships$grid,
    sbm_Fungi_gaussian$memberships$grid,
    sbm_Virus$memberships$grid,
    sbm_Plant$memberships$grid
  )
)

colnames(B_multi_regne) <- c("grid_code", "Bact", "Fungi", "Virus","Plant", "Freq")


w_filter <- which(B_multi_regne$Freq != 0)
B_multi_regne <- B_multi_regne[w_filter,]

library(alluvial)

alluvial(B_multi_regne[, c(2, 3, 4, 5)],
         freq = B_multi_regne$Freq,
         alpha = 0.7,
         cex = 1.5,
         cex.axis = 1.5,
         xw = 0.2)
# Extract block per species -----


sbm_Bact_gaussian$memberships$grid -> row_cluster
unique(sbsbm_Bact_gaussian$memberships$bact) -> col_cluster
sbm_Bact_gaussian -> bsbm

str(sbm_Bact_gaussian)
sbm_Bact_gaussian$networkData %>% str
sbm_Bact_gaussian$memberships %>% str
sbm_Bact_gaussian$expectation %>%
  as.data.frame() %>%
  rownames_to_column("grid") %>%
  pivot_longer(-grid,names_to = "ASV", values_to = "log_expected_reads") %>%
  mutate(ASV = as.character(ASV)) %>%
  left_join(data.frame(cluster = sbm_Bact_gaussian$memberships$bact,
                       ASV = colnames(sbm_Bact_gaussian$networkData)), by = "ASV") %>% 
  group_by(grid, cluster) %>%
  summarise(log_expected_reads = mean(log_expected_reads)) %>%
  pivot_wider(names_from = "cluster", values_from = "log_expected_reads") %>% 
  column_to_rownames("grid") %>%
  as.matrix %>% str
plot_bipartite_sbm(sbm_Bact_gaussian, mat = sbm_Bact_gaussian$expectation)

plot_bipartite_sbm_avg(sbm_Bact_gaussian) -> sbm_Bact_average
sbm_Bact_average$plot
plot_bipartite_sbm_avg(sbm_Fungi_gaussian) -> sbm_Fungi_average
sbm_Fungi_average$plot

# SBM at quadra level ----

## Bacteria ----

ps_bacteria_quadra <- ps_bacteria %>%
  prune_samples(!str_detect(sample_names(.), "21-CAM-1[45]"), .)

ps_bacteria_quadra %>%
  sample_data %>%
  pull(grid) %>% unique %>% length

ps_bacteria_quadra %>%
  otu_table() %>% 
  as.data.frame() %>% 
  rename_with(~str_replace_all(., "-","_")) %>%
  rename_with(~substring(.,1,13))-> bacteria_quadra_matrix

plot(1:dim(bacteria_quadra_matrix)[2],
     sort(log10(colSums(bacteria_quadra_matrix))))

table(sort(log10(colSums(bacteria_quadra_matrix)))> 2.5)
sort(colSums(bacteria_quadra_matrix))
dim(bacteria_quadra_matrix)

bacteria_quadra_matrix %>%
  rownames_to_column("OTU") %>%
  pivot_longer(-OTU, names_to = "quadra_code", values_to = "reads") %>%
  group_by(quadra_code) %>%
  filter(log10(sum(reads))> 2.5) %>%
  ungroup() %>%
  pivot_wider(names_from = "quadra_code", values_from = "reads", values_fill = 0) %>%
  column_to_rownames("OTU") -> bacteria_quadra__filter_matrix
dim(bacteria_quadra__filter_matrix)
plot(1:dim(bacteria_quadra__filter_matrix)[2],
     sort(log10(colSums(bacteria_quadra__filter_matrix))))
sort(colSums(bacteria_quadra__filter_matrix))

bacteria_quadra__filter_matrix %>%
  rownames_to_column("OTU") %>%
  pivot_longer(-OTU, names_to = "quadra_code", values_to = "reads") %>%
  mutate(reads = floor(reads/10)) %>%
  filter(reads > 0) %>%
  pivot_wider(names_from = "quadra_code", values_from = "reads", values_fill = 0) %>%
  column_to_rownames("OTU") -> bacteria_quadra_truncated_matrix

dim(bacteria_quadra_matrix)
dim(bacteria_quadra_truncated_matrix)


plot(1:dim(bacteria_quadra_truncated_matrix)[2],
     sort(log10(colSums(bacteria_quadra_truncated_matrix))))

sort(colSums(bacteria_quadra_matrix))
sort(colSums(bacteria_quadra__filter_matrix))
sort(colSums(bacteria_quadra_truncated_matrix))


## Fungi ----
ps_fungi_quadra %>%
  sample_data %>%
  pull( grid_code) %>% unique %>% length
ps_fungi_quadra <- ps_fungi %>%
  prune_samples(!str_detect(sample_names(.), "21-CAM-1[45]"), .)

ps_fungi_quadra %>%
  otu_table() %>% 
  as.data.frame() %>% 
  rename_with(~str_replace_all(., "-","_")) -> fungi_quadra_matrix

fungi_quadra_matrix %>%
  rownames_to_column("OTU") %>%
  pivot_longer(-OTU, names_to = "grid_code", values_to = "reads") %>%
  mutate(reads = floor(reads/5)) %>%
  filter(reads > 0) %>%

  pivot_wider(names_from = "grid_code", values_from = "reads", values_fill = 0) %>%
  column_to_rownames("OTU") -> fungi_quadra_truncated_matrix

dim(fungi_quadra_matrix)
dim(fungi_quadra_truncated_matrix)



## Compute sample coverage ----

iNext_object_bact_quadra <- bacteria_quadra_truncated_matrix %>%
  iNEXT.3D::iNEXT3D( q = c(0), datatype = "abundance", nboot = 0)

iNext_object_bact_quadra$TDiNextEst$coverage_based %>% 
  filter(Order.q == 0 & Method == "Observed") %>%
  arrange(SC) 
iNext_object_bact_quadra$TDiNextEst$coverage_based %>% 
  filter(Order.q == 0 & Method == "Observed") %>%
  arrange(SC) %>%
  filter(SC == min(SC, na.rm = TRUE)) %>%
  pull(SC) -> min_coverage_bact_coverage

bacteria_quadra_truncated_matrix %>%
  iNEXT.3D::estimate3D(level = min_coverage_bact_coverage, q = 0, nboot = 0) %>%
  select(Assemblage,m) -> required_sizes_bact_quadra

### Fungi ----
iNext_object_OTU_fungi <- fungi_grid_truncated_matrix %>%
  iNEXT.3D::iNEXT3D( q = c(0), datatype = "abundance", nboot = 0)

iNext_object_OTU_fungi$TDiNextEst$coverage_based %>% 
  filter(Order.q == 0 & Method == "Observed") %>%
  filter(SC == min(SC, na.rm = TRUE)) %>%
  pull(SC) -> min_coverage_OTU_fungi

fungi_grid_truncated_matrix %>%
  iNEXT.3D::estimate3D(level = min_coverage_OTU_fungi, q = 0, nboot = 0) %>%
  select(Assemblage,m) -> required_sizes_OTU_fungi

