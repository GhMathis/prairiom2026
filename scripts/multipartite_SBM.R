
# Load pacakges and data -----

library(tidyverse)
library(sbm)
library(phyloseq)
library(iNEXT.3D)

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
# Rarefation at reference coverage ----

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

library(iNEXT)
dim(virus_grid_matrix_Inext)
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
    
    if (f1 > 0 && f2 > 0) {
      1 - (f1 / n) * (((n - 1) * f1) / ((n - 1) * f1 + 2 * f2))
    } else if (f1 > 0 && f2 == 0) {
      1 - (f1 / n) * (((n - 1) * f1) / ((n - 1) * f1 + 2))
    } else {
      1
    }
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

n_sim = 100
set.seed(2365)
1:n_sim %>%
  map(\(x) fungi_grid_truncated_matrix %>%
        coverage_rarefy_matrix(required_sizes = required_sizes_OTU_fungi)) -> coverage_OTU_fungi_list

# Multipartite SBM -----

colnames(coverage_OTU_bact_list[[1]]) == colnames(coverage_OTU_fungi_list[[1]])
colnames(coverage_OTU_fungi_list[[1]]) == rownames(virus_grid_matrix)
rownames(virus_grid_matrix) == rownames(plant_grid_matrix)

GridBact <- defineSBM(log1p(t(coverage_OTU_bact_list[[1]])), 
                         model = "gaussian",
                         type = "bipartite", 
                         directed = F, dimLabels = c("grid", "bact"))

GridFungi <- defineSBM(log1p(t(coverage_OTU_fungi_list[[1]])),
                      model = "gaussian",
                      type = "bipartite", 
                      directed = F, dimLabels = c("grid", "fungi"))
                      
GridVirus <- defineSBM(virus_grid_matrix,
                       model = "poisson",
                       type = "bipartite", 
                       directed = F, dimLabels = c("grid", "virus"))

GridPlant <- defineSBM(log1p(plant_grid_matrix),
                       model = "gaussian",
                       type = "bipartite", 
                       directed = F, dimLabels = c("grid", "plant"))

plotMyMultipartiteMatrix(list(GridBact, GridFungi, GridVirus, GridPlant))

estimOptions = list(initBM = FALSE, nbCores = 2)
listSBM <- list(GridBact, GridFungi, GridPlant, GridVirus)

myMSBM <- estimateMultipartiteSBM(listSBM, estimOptions)                       

source("scripts/plot_multipartite_sbm.R")

plot_multipartite_sbm(myMSBM)

myMSBM$blockProp
myMSBM$connectParam

GridBact_bin <- defineSBM(+(t(coverage_OTU_bact_list[[1]])>0), 
                      model = "bernoulli",
                      type = "bipartite", 
                      directed = F, dimLabels = c("grid", "bact"))

GridFungi_bin <- defineSBM(+(t(coverage_OTU_fungi_list[[1]])>0),
                       model = "bernoulli",
                       type = "bipartite", 
                       directed = F, dimLabels = c("grid", "fungi"))

GridVirus<- defineSBM(virus_grid_matrix,
                           model = "bernoulli",
                           type = "bipartite", 
                           directed = F, dimLabels = c("grid", "virus"))

GridVirus_bin <- defineSBM(+(virus_grid_matrix>0),
                       model = "bernoulli",
                       type = "bipartite", 
                       directed = F, dimLabels = c("grid", "virus"))

GridPlant_bin <- defineSBM(+(plant_grid_matrix>0),
                       model = "bernoulli",
                       type = "bipartite", 
                       directed = F, dimLabels = c("grid", "plant"))

plotMyMultipartiteMatrix(list(GridBact_bin , GridFungi_bin , GridVirus_bin , GridPlant_bin ))

listSBM_bin <- list(GridBact_bin, GridFungi_bin, GridPlant_bin)
if(file.exists( file = "outputs/heavy_computation_save/Multipartie_sbm_no_virus.Rdata")){
 load(file = "outputs/heavy_computation_save/Multipartie_sbm_no_virus.Rdata")
} else{
  myMSBM_bin <- estimateMultipartiteSBM(listSBM_bin, estimOptions) 
  save(myMSBM_bin, file = "outputs/heavy_computation_save/Multipartie_sbm_no_virus.Rdata")
}

listSBM_bin2 <- list(GridBact_bin, GridFungi_bin, GridPlant_bin,GridVirus)
if(file.exists( file = "outputs/heavy_computation_save/Multipartie_sbm_virus.Rdata")){
  load(file = "outputs/heavy_computation_save/Multipartie_sbm_virus.Rdata")
} else{
  myMSBM_bin2 <- estimateMultipartiteSBM(listSBM_bin2, estimOptions) 
  save(myMSBM_bin2, file = "outputs/heavy_computation_save/Multipartie_sbm_virus.Rdata")
}

listSBM_bin3 <- list(GridBact_bin, GridFungi_bin, GridPlant_bin,GridVirus_bin)
if(file.exists( file = "outputs/heavy_computation_save/Multipartie_sbm_virus_bin.Rdata")){
  load(file = "outputs/heavy_computation_save/Multipartie_sbm_virus_bin.Rdata")
} else{
  myMSBM_bin3 <- estimateMultipartiteSBM(listSBM_bin3, estimOptions) 
  save(myMSBM_bin3,file =  "outputs/heavy_computation_save/Multipartie_sbm_virus_bin.Rdata")
}

# Equivalent coverage for all, except viruses ----

# Fungi

fungi_grid_truncated_matrix %>%
  iNEXT.3D::estimate3D(level = min_coverage_OTU_bact, q = 0, nboot = 0) %>%
  select(Assemblage,m) -> required_sizes_OTU_fungi2

n_sim = 100
set.seed(2365)
1:n_sim %>%
  map(\(x) fungi_grid_truncated_matrix %>%
        coverage_rarefy_matrix(required_sizes = required_sizes_OTU_fungi2)) -> coverage_OTU_fungi_list2

GridFungi_bactcoverage_bin <- defineSBM(+(t(coverage_OTU_fungi_list2[[1]])>0),
                                        model = "bernoulli",
                                        type = "bipartite", 
                                        directed = F, dimLabels = c("grid", "fungi"))

listSBM_bin4 <- list(GridBact_bin, GridFungi_bactcoverage_bin, GridPlant_bin,GridVirus_bin)
if(file.exists( file = "outputs/heavy_computation_save/Multipartie_sbm_equi_covearge_bin.Rdata")){
  load( file = "outputs/heavy_computation_save/Multipartie_sbm_equi_covearge_bin.Rdata")
} else{
  
  myMSBM_equi_covearge_bin <- estimateMultipartiteSBM(listSBM_bin4, estimOptions) 
  save(myMSBM_equi_covearge_bin, file = "outputs/heavy_computation_save/Multipartie_sbm_equi_covearge_bin.Rdata")
}

# Plot all multipartite sbm tests ----
source("scripts/plot_multipartite_sbm.R")

plot_multipartite_sbm(myMSBM_bin)
plot_multipartite_sbm(myMSBM_bin2)
plot_multipartite_sbm(myMSBM_bin3)
plot_multipartite_sbm(myMSBM_equi_covearge_bin)

# Paralelisation multi sbm ----
# Pass ALL matrices as explicit arguments — no hidden globals
define.multiSBM <- function(mat_bact, mat_fungi, mat_virus, mat_plant) {
  
  GridBact  <- defineSBM(t(mat_bact),  model = "bernoulli", type = "bipartite",
                         directed = FALSE, dimLabels = c("grid", "bact"))
  GridFungi <- defineSBM(t(mat_fungi), model = "bernoulli", type = "bipartite",
                         directed = FALSE, dimLabels = c("grid", "fungi"))
  GridVirus <- defineSBM(mat_virus,    model = "bernoulli", type = "bipartite",
                         directed = FALSE, dimLabels = c("grid", "virus"))
  GridPlant <- defineSBM(mat_plant,    model = "bernoulli", type = "bipartite",
                         directed = FALSE, dimLabels = c("grid", "plant"))
  
  estimOptions <- list(initBM = FALSE, nbCores = 1)
  estimateMultipartiteSBM(list(GridBact, GridFungi, GridVirus, GridPlant), estimOptions)
}

library(furrr)

if(file.exists( file = "outputs/heavy_computation_save/multipartite_sbm_list.Rdata")){
  load(file = "outputs/heavy_computation_save/multipartite_sbm_list.Rdata")
} else{
  # Use multicore (fork-based) on Linux/Mac — much safer with C++ backends
  # On Windows, keep multisession but limit workers and increase timeout
  # plan(multicore, workers = 10)          # Linux/Mac
  plan(multisession, workers = 10)     # Windows
  options(future.globals.maxSize = 2 * 1024^3)  # 2 GB per worker
  function_time <- system.time({
    future_map2(
      coverage_OTU_bact_list,
      coverage_OTU_fungi_list,
      .f = ~ define.multiSBM(
        mat_bact  = +(.x > 0),
        mat_fungi = +(.y > 0),
        mat_virus = +(virus_grid_matrix > 0),   # now explicitly passed → furrr exports it
        mat_plant = +(plant_grid_matrix > 0)
      ),
      .options = furrr_options(
        globals = c("virus_grid_matrix", "plant_grid_matrix",
                    "define.multiSBM"),
        packages = c("sbm"),             # make sure sbm is loaded in each worker
        scheduling = Inf,
        seed     = TRUE
      )
    ) -> multipartite_sbm_list
  }
  )
  save(multipartite_sbm_list,file = "outputs/heavy_computation_save/multipartite_sbm_list.Rdata")
}

str(multipartite_sbm_list)
x = 1
map(1:n_sim, \(x) data.frame(subsample = multipartite_sbm_list[[x]]$memberships$grid,
                             grid_code = colnames(subsample_list[[1]]),
                             sim = x)
) %>%
  bind_rows() %>%
  pivot_wider(names_from = sim, names_glue = "simulation{sim}", values_from = "subsample") %>%
  column_to_rownames("grid_code") -> cluster_df
n_echantillons <- nrow(cluster_df)
grid_cluster_mat <- matrix(0, nrow = n_echantillons, ncol = n_echantillons)

for (sim in 1:n_sim) {
  grid_cluster_mat <- grid_cluster_mat + outer(cluster_df[[sim]], cluster_df[[sim]], "==")
}


