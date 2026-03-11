# Load data and packages ----

library(tidyverse)
library(viridis)
library(alluvial)
library(sbm)
library(pals)
library(RColorBrewer)
library(aricode)
library(cowplot)
library(vegan)
library(FactoMineR)


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



read.table("data/data_clean/OTU_virus_CAM.txt", header = T)%>%
  filter(!str_detect(Host_code, "22_CAM_15...."))%>%
  mutate(Grid_code = str_extract(Host_code, ".._CAM_.."))%>%
  group_by(Grid_code)%>%
  summarise(across(-c(Host_code), sum))%>%
  column_to_rownames("Grid_code")%>%
  t() -> otu_virus_grid

hist(otu_virus_grid)

read.table("data/data_clean/Metadata_grid_CAM.txt", header = T) -> metadata_grid

otu_virus_grid%>%
  as.matrix()%>%
  t()%>%
  plotMyMatrix(dimLabels = c("Grids", "Virus"), plotOptions = list(rowNames = T,colNames = F))+
  theme(element_text(angle = 90, vjust = 1, hjust = 1, size = 0.1, face="italic"),
        strip.text.x = element_text(colour = "gray", size=15, face ="italic"),
        strip.text.y  = element_text(colour = "gray", size=15, face ="italic"),
        strip.background = element_rect(fill="white"),
        panel.border = element_rect(color = "black",  fill = NA))

# SBM virus ----
sbm_virus <-otu_virus_grid%>%
  as.matrix() %>%
  estimateBipartiteSBM(
    model = 'poisson', 
    dimLabels = c(row = "Virus", col = "Grids"))

if(!any(colnames(metadata_grid) %in% "sbm_virus")){ #add the column only if it is not present in the df (avoid duplicated column)

  
  metadata_grid%>%
    left_join(data.frame(sbm_virus = sbm_virus$memberships$Grids,
                         Grid_code = colnames(otu_virus_grid)), by ="Grid_code")->metadata_grid
  metadata_grid %>% write.table("data/data_clean/Metadata_grid_CAM.txt")
} 

## aluvial plot ----
B <-
  as.data.frame(
    table(
      
      metadata_grid$Grid_code,
      metadata_grid$clust_smb_grid_plant ,
      metadata_grid$clust_land_use,
      metadata_grid$clust_landscp,
      metadata_grid$sbm_virus,
      metadata_grid$sbm_soil)
  )


colnames(B) =  c( "Grid_code","sbm_plant", "land_use",  "landscp", "sbm_virus","smb_soil","Freq")


w   <- which(B$Freq != 0)
B <- B[w, ]

data.frame(sbm_plant = unique(B$sbm_plant),
           col_plant = brewer.pal(7, "Set1")) -> col_cluster_plant
data.frame(sbm_virus = unique(B$sbm_virus),
           col_virus = brewer.pal(4, "Set2")) -> col_cluster_virus

B%>%
  left_join(col_cluster_virus, by = "sbm_virus")%>%
  left_join(col_cluster_plant, by = "sbm_plant")-> B

B%>%
  mutate(landscp =  case_when(landscp == 1 ~"(1) Cultivés\n  dominant+\n Anthropisés",
                              landscp == 2 ~"(2) Naturels\n  non Humide\n dominant",
                              landscp == 3 ~"(3) Mixtes",
                              landscp == 4 ~"(4) Naturels\n  Humide\n dominant")) -> B

B$sbm_plant <- factor(B$sbm_plant, levels = c("5", "3", "2","1", "4","6","7"))
B$sbm_virus <- factor(B$sbm_virus, levels = c("2", "1", "3", "4"))

B$landscp <- factor(B$landscp, levels = c( "(1) Cultivés\n  dominant+\n Anthropisés",
  "(3) Mixtes", "(2) Naturels\n  non Humide\n dominant","(4) Naturels\n  Humide\n dominant"))
colnames(B)[ c(4,2,5)] = c("Paysage", "Végétale", "Virale")
alluvial(B[, c(4,5,2)],
         freq = B$Freq,
         col = B$col_virus,
         alpha = 0.7,
         cex = 1.2,
         cex.axis = 1.5,
         xw = 0.2)

alluvial(B[, c(5,3)],
         freq = B$Freq,
         col = B$col_virus,
         alpha = 0.7,
         xw = 0.2)

alluvial(B[, c(5,4)],
         freq = B$Freq,
         col = B$col_virus,
         alpha = 0.7,
         xw = 0.2)

alluvial(B[, c(5,6)],
         freq = B$Freq,
         col = B$col_virus,
         alpha = 0.7,
         xw = 0.2)
## NMI ----
nmi_land_use = NMI(c1 =metadata_grid$sbm_virus ,c2 = metadata_grid$clust_land_use)
nmi_landscape = NMI(c1 =metadata_grid$sbm_virus ,c2 = metadata_grid$clust_landscp)
nmi_soil = NMI(c1 =metadata_grid$sbm_virus ,c2 = metadata_grid$sbm_soil)
nmi_plant = NMI(c1 =metadata_grid$sbm_virus ,c2 = metadata_grid$clust_smb_grid_plant)
nmi_land_use 
nmi_landscape
nmi_soil
nmi_plant 

## Randomization ----

sapply(1:10000,function(x) sample(metadata_grid$sbm_virus)) -> rowperm_nmi
data.frame(row_land_use  = apply(rowperm_nmi,2 ,function(x) NMI(c1 = x ,c2 = metadata_grid$clust_land_use)) ,
           row_soil = apply(rowperm_nmi,2, function(x) NMI(c1 = x ,c2 = metadata_grid$sbm_soil)) ,
           row_landscape = apply(rowperm_nmi, 2, function(x) NMI(c1 = x ,c2 = metadata_grid$clust_landscp)),
           row_plant = apply(rowperm_nmi, 2, function(x) NMI(c1 = x ,c2 = metadata_grid$clust_smb_grid_plant))
)-> nmi_rand_df
nmi_rand_df%>%
  pivot_longer(everything(), names_to = "rand_type", values_to = "nmi_rand_values")%>%
  group_by(rand_type)%>%
  left_join(data.frame(nmi_values = c(nmi_land_use, nmi_landscape, nmi_soil, nmi_plant),
                       rand_type = c("row_land_use", "row_landscape",  "row_soil", "row_plant")),
            by = "rand_type")-> nmi_rand_df

ggplot(nmi_rand_df, scale = "free")+
  facet_wrap(~rand_type)+#,
  #labeller = as_labeller(setNames(plot_names, sort(unique(nmi_rand_df$rand_type)))))+
  geom_histogram(aes(nmi_rand_values),bins = 100)+
  geom_vline(aes(xintercept = nmi_values),col ="red", linetype =2,linewidth =1.5)+
  labs(x ="Values of NMI", y= "Count")+
  main_theme+
  theme_classic()
nmi_rand_df%>%
  group_by(rand_type)%>%
  summarise(mean(nmi_rand_values))
# CCA ----
Y = metadata_grid$sbm_virus%>%as.factor()
X = metadata_grid%>%dplyr::select(clust_landscp,clust_land_use, sbm_soil, clust_smb_grid_plant)%>%
  mutate_all(as.factor)
full_cca_virus = cca(tab.disjonctif(Y)~sbm_soil +clust_land_use+ clust_landscp+clust_smb_grid_plant, data = X)

mod_varpart_virus_chi = varpart(tab.disjonctif(as.factor(Y)), ~clust_smb_grid_plant,
                                 ~clust_land_use, ~clust_landscp, data=  X, permutations = how(nperm = 10000), chisquare=T)

bg_colors <- c("#00735c", "#8c52ff", "#c0344a")
bg_colors_alpha <- adjustcolor(bg_colors, alpha.f = 1) 

plot(mod_varpart_virus_chi,Xnames = c(list("Plant \n community profiles"),
                                        list("Agricultural \n practice profiles"),
                                      list("Land use \n profiles")),
     bg = bg_colors_alpha, lwd = 2, cex = 1.5)

source("code/functions_modified.R")

test_env_ana = analysis.function.alt(incid = otu_virus_grid%>%t(), memb_obs = as.factor(metadata_grid$sbm_virus),
                                                   var1 = as.factor(metadata_grid$clust_land_use), 
                                                   var2 = as.factor(metadata_grid$clust_landscp),
                                                   var3 = as.factor(metadata_grid$clust_smb_grid_plant), 
                                                   NPERM = 1000,
                                                   configs = NULL #permutation_plant_quadclust%>%as.data.frame(),
                                                   )
select_cca_virus = cca(tab.disjonctif(Y)~clust_land_use+ clust_landscp+clust_smb_grid_plant, data = X)
test_env_ana$R2_proxy = test_env_ana$chi2/select_cca_virus$tot.chi
test_env_ana

## Clustert profils ----
virus_id = sbm_virus$memberships[1]%>%
  as.data.frame()%>%
  rownames_to_column("id")%>%
  mutate(id = as.numeric(id))%>%
  arrange(Virus)
grid_id = sbm_virus$memberships[2]%>%
  as.data.frame()%>%
  rownames_to_column("id")%>%
  mutate(id = as.numeric(id))%>%
  arrange(Grids)
otu_virus_grid[virus_id$id,grid_id$id]%>%
  t()%>%
  plotMyMatrix(dimLabels = c("Grilles", "Virus"), plotOptions = list(rowNames = T,colNames = T))+
  geom_hline(yintercept = abs(which(!duplicated(grid_id$Grids))[-1]-nrow(grid_id)-1)+0.5, col ="red", linetype = 5)+
  geom_vline(xintercept = which(!duplicated(virus_id$Virus))[-1]-0.5, col = "red")+
  
  theme(element_text(angle = 45, vjust = 1, hjust = 1, size = 0.1, face="italic"),
        strip.text.x = element_text(colour = "gray", size=20, face ="italic"),
        strip.text.y  = element_text(colour = "gray", size=20, face ="italic"),
        axis.text.x =element_text( size = 8, face="italic"),
        axis.text.y = element_text(size = 10,face="italic"),
        strip.background = element_rect(fill="white"),
        panel.border = element_rect(color = "black",  fill = NA))

str(otu_virus_grid)
reshape2::melt(otu_virus_grid)%>%
  dplyr::rename(Virus = "Var1", Grid_code = "Var2" )%>%
  full_join(metadata_grid, by = "Grid_code")-> virus_profil


virus_profil%>%
  dplyr::select(Virus, value,sbm_virus )%>%
  filter(value>0)%>%
  
  group_by(sbm_virus, Virus)%>%
  summarise(n = sum(value))%>%
  group_by(sbm_virus)%>%
  add_tally(n)%>%
  ungroup()%>%
  group_by(sbm_virus,Virus)%>%
  arrange(desc(n))%>%
  group_by(sbm_virus)%>%
  slice_head( n = 10)%>%
  mutate(prop = n/nn) -> virussp_df

virussp_df$Virus -> best_virus_cluster
virus_id$virus_names = rownames(otu_virus_grid[virus_id$id, grid_id$id])
virus_id2 = virus_id%>%
  filter(virus_names %in% best_virus_cluster)
abund_virus_grid_fomated = as.matrix(otu_virus_grid)
colnames(abund_virus_grid_fomated) = as.character(1:41)
rownames(abund_virus_grid_fomated) = str_replace_all(rownames(abund_virus_grid_fomated),"[.]", " ")
rownames(abund_virus_grid_fomated) = str_replace_all(rownames(abund_virus_grid_fomated),"sp", "sp. ")

abund_virus_grid_fomated[virus_id2$id,grid_id$id]%>%
  plotMyMatrix(dimLabels = c("Virus", "Grilles"), plotOptions = list(rowNames = T,colNames = T))+
  geom_hline(yintercept = abs(which(!duplicated(virus_id2$Virus))[-1]-nrow(virus_id2)-1)+0.5, col ="red", linetype = 5)+
  geom_vline(xintercept = which(!duplicated(grid_id$Grids))[-1]-0.5, col = "red")+
  
  theme(element_text(angle = 45, vjust = 1, hjust = 1, size = 0.1, face="italic"),
        strip.text.x = element_text(colour = "gray", size=20, face ="italic"),
        strip.text.y  = element_text(colour = "gray", size=20, face ="italic"),
        axis.text.x =element_text( size = 8, face="italic"),
        axis.text.y = element_text(size = 10,face="italic"),
        strip.background = element_rect(fill="white"),
        panel.border = element_rect(color = "black",  fill = NA))

virussp_df %>%
  mutate(Virus =  str_replace_all(Virus,"[.]", " "),
         Virus = str_replace_all(Virus,"sp", "sp. "),
         prop = round(prop,3)*100,
         Virus =paste(Virus," (", prop, "%)"))%>%
  group_by(sbm_virus)%>%
  mutate(id_row = 1:10)%>%
  ungroup()%>%
  dplyr::select(sbm_virus, Virus, id_row)%>%
  pivot_wider(values_from = Virus, names_from =id_row )%>%
  as.data.frame()%>%
  t()%>%
  as_tibble()%>%
  tail(10)%>%
  knitr::kable()

  
