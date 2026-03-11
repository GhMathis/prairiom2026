# Load data and packages ----

library(tidyverse)
library(viridis)
library(alluvial)
library(sbm)
library(pals)
library(RColorBrewer)
library(aricode)
library(cowplot)
library(reshape2)


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



read.table("data/data_clean/Metadata_grid_CAM.txt", header = T)%>%
  filter(!str_detect(Grid_code, "22_CAM_15")) -> metadata_grid

read.table("data/data_clean/Metadata_quadra_CAM.txt", header = T)%>%
  filter(!str_detect(Host_code, "22_CAM_15...."))-> metadata_quadra

read.table("data/data_clean/abund_plant_grid.txt", header = T)%>%
  filter(!str_detect(Grid_code, "22_CAM_15"))%>%
  column_to_rownames("Grid_code")%>%
  as.matrix()-> abund_plant_grid

read.table("data/data_clean/abund_plant_grid.txt", header = T)%>%
  filter(!str_detect(Grid_code, "22_CAM_15"))%>%
  mutate(across(where(is.numeric), ~ifelse(. != 0,1,0 )))%>%
  column_to_rownames("Grid_code")%>%
  as.matrix() -> bin_plant_grid
summary(rowSums(bin_plant_grid))
hist(log1p(abund_plant_grid))
hist(abund_plant_grid)

# abund plant
abund_plant_grid%>%
  as.matrix()%>%
  plotMyMatrix(dimLabels = c("Grids", "Plants"), plotOptions = list(rowNames = T,colNames = F))+
  theme(element_text(angle = 90, vjust = 1, hjust = 1, size = 0.1, face="italic"),
        strip.text.x = element_text(colour = "gray", size=15, face ="italic"),
        strip.text.y  = element_text(colour = "gray", size=15, face ="italic"),
        strip.background = element_rect(fill="white"),
        panel.border = element_rect(color = "black",  fill = NA))

# SBM Plants ----

if(file.exists("outputs/SBM_abund_grid.Rdata")){
  load(file="outputs/SBM_abund_grid.Rdata")
  
}else{
  SBM_abund_grid <-
    as.matrix(abund_plant_grid) %>%
    estimateBipartiteSBM(
      model = 'poisson', 
      dimLabels = c(row = "Grids", col = "Plants"))
  save(SBM_abund_grid, file="outputs/SBM_abund_grid.Rdata")
}


plot(SBM_abund_grid)


SBM_abund_grid$expectation
plant_id = SBM_abund_grid$memberships[1]%>%
  as.data.frame()%>%
  rownames_to_column("id")%>%
  mutate(id = as.numeric(id))%>%
  arrange(Plants)
grid_id = SBM_abund_grid$memberships[2]%>%
  as.data.frame()%>%
  rownames_to_column("id")%>%
  mutate(id = as.numeric(id))%>%
  arrange(Grids)


abund_plant_grid[grid_id$id,plant_id$id]%>%
  log1p()%>%
  as.matrix()%>%
  plotMyMatrix(dimLabels = c("Grilles", "Plantes"), plotOptions = list(rowNames = T,colNames = T))+
  geom_hline(yintercept = abs(which(!duplicated(grid_id$Grids))[-1]-nrow(grid_id)-1)+0.5, col ="red")+
  geom_vline(xintercept = which(!duplicated(plant_id$Plants))[-1]-0.5, col = "red", linetype = 5)+
  theme(element_text(angle = 90, vjust = 1, hjust = 1, size = 0.1, face="italic"),
        strip.text.x = element_text(colour = "gray", size=15, face ="italic"),
        strip.text.y  = element_text(colour = "gray", size=15, face ="italic"),
        strip.background = element_rect(fill="white"),
        panel.border = element_rect(color = "black",  fill = NA),
        axis.text.x = element_blank(),  # Retire les étiquettes des ticks de l'axe x
        axis.text.y = element_blank())

## Clustert profils ----

data.frame(Grid_code = rownames(abund_plant_grid),
           clust_smb_grid_plant = as.factor(SBM_abund_grid$memberships$Grids))%>%
  left_join(metadata_grid%>%dplyr::select(Grid_code, clust_land_use,
                                                  clust_landscp,
                                                  sbm_soil), by = "Grid_code" ) ->alluvial_site_df

reshape2::melt(t(abund_plant_grid)[plant_id$id,grid_id$id])%>%
  dplyr::rename(Plants = "Var1", Grid_code = "Var2" )%>%
  full_join(alluvial_site_df, by = "Grid_code")-> plant_profil

plant_profil%>%
  dplyr::select(Plants, value, clust_smb_grid_plant)%>%
  filter(value>0)%>%
  
  group_by(clust_smb_grid_plant, Plants)%>%
  summarise(n = sum(value))%>%
  group_by(clust_smb_grid_plant)%>%
  add_tally(n)%>%
  ungroup()%>%
  group_by(clust_smb_grid_plant,Plants)%>%
  arrange(desc(n))%>%
  group_by(clust_smb_grid_plant)%>%
  slice_head( n = 10)%>%
  mutate(prop = n/nn) -> plantsp_df

plantsp_df$Plants -> best_plant_cluster
plant_id$Plants_names = colnames(abund_plant_grid[grid_id$id, plant_id$id])
plant_id2 = plant_id%>%
  filter(Plants_names %in% best_plant_cluster)
abund_plant_grid_fomated = t(as.matrix(abund_plant_grid))
colnames(abund_plant_grid_fomated) = as.character(1:41)
rownames(abund_plant_grid_fomated) = str_replace_all(rownames(abund_plant_grid_fomated),"[.]", " ")
rownames(abund_plant_grid_fomated) = str_replace_all(rownames(abund_plant_grid_fomated),"sp", "sp. ")
  
abund_plant_grid_fomated[plant_id2$id,grid_id$id]%>%
  log1p()%>%
  plotMyMatrix(dimLabels = c("Plantes", "Grilles"), plotOptions = list(rowNames = T,colNames = T))+
  geom_hline(yintercept = abs(which(!duplicated(plant_id2$Plants))[-1]-nrow(plant_id2)-1)+0.5, col ="red", linetype = 5)+
  geom_vline(xintercept = which(!duplicated(grid_id$Grids))[-1]-0.5, col = "red")+
  
  theme(element_text(angle = 45, vjust = 1, hjust = 1, size = 0.1, face="italic"),
        strip.text.x = element_text(colour = "gray", size=20, face ="italic"),
        strip.text.y  = element_text(colour = "gray", size=20, face ="italic"),
        axis.text.x =element_text( size = 8, face="italic"),
        axis.text.y = element_text(size = 9,face="italic"),
        strip.background = element_rect(fill="white"),
        panel.border = element_rect(color = "black",  fill = NA))

plantsp_df %>%
  mutate(Plants =  str_replace_all(Plants,"[.]", " "),
         Plants = str_replace_all(Plants,"sp", "sp. "),
         prop = round(prop,3)*100,
         Plants =paste(Plants," (", prop, "%)"))%>%
  group_by(clust_smb_grid_plant)%>%
  mutate(id_row = 1:10)%>%
  ungroup()%>%
  dplyr::select(clust_smb_grid_plant, Plants, id_row)%>%
  pivot_wider(values_from = Plants, names_from =id_row )%>%
  as.data.frame()%>%
  t()%>%
  as_tibble()%>%
  tail(10)%>%
  knitr::kable()
## Aluvial plot ----

if(!any(colnames(metadata_grid) %in% "clust_smb_grid_plant")){ #add the column only if it is not present in the df (avoid duplicated column)
  metadata_grid%>%
    left_join(alluvial_site_df%>%dplyr::select(Grid_code, clust_smb_grid_plant), by ="Grid_code")%>%
    write.table("data/data_clean/Metadata_grid_CAM.txt")
}
  

B <-
  as.data.frame(
    table(
 
      alluvial_site_df$Grid_code,
      alluvial_site_df$clust_smb_grid_plant,
      alluvial_site_df$clust_land_use,
      alluvial_site_df$clust_landscp,
      alluvial_site_df$sbm_soil,
      metadata_grid$Ecosystem)
  )


colnames(B) =  c( "Grid_code","sbm_plant", "land_use", "landscp", "sbm_soil", "Ecosystem","Freq")

w   <- which(B$Freq != 0)
B <- B[w, ]

data.frame(sbm_plant = unique(B$sbm_plant),
           col_plant_sbm = brewer.set1(7)) -> col_cluster_sbm_plant

B%>%
  left_join(col_cluster_sbm_plant, by = "sbm_plant") -> B

B%>%
  mutate(sbm_soil =  case_when(sbm_soil == 1 ~"(1) Salés+\n peu Fertiles",
                               sbm_soil == 2 ~"(2) Salés+\n Fertiles",
                               sbm_soil == 3 ~"(3) peu Salés+\n Fertiles",
                               sbm_soil == 4 ~"(4) peu Salés+\n peu Fertiles"),
         land_use = case_when(land_use==1 ~"(1) Pâture",
                      land_use==2 ~"(2) Rien",
                      land_use==3 ~"(3) Fauche", 
                      land_use==4 ~"(4) Pâture+\n Fauche")) -> B


alluvial(B[, c(2,6)],
         freq = B$Freq,
         alpha = 0.7,
         col = B$col_plant_sbm,
         xw = 0.2)

B$sbm_plant <- factor(B$sbm_plant, levels = c("1","2","3", "6","5", "7", "4"))
alluvial(B[, c(2,3)],
         freq = B$Freq,
         col =  B$col_plant_sbm,
         alpha = 0.7,
         xw = 0.2)

B$sbm_plant <- factor(B$sbm_plant, levels = c("1","2","3", "4","5", "6", "7"))
alluvial(B[, c(2,4)],
         freq = B$Freq,
         col = B$col_plant_sbm,
         alpha = 0.7,
         xw = 0.2)

B$sbm_plant <- factor(B$sbm_plant, levels = c("3", "6", "2","1","5",  "7", "4"))
B$land_use <- factor(B$land_use, levels = c("(1) Pâture","(2) Rien","(3) Fauche", "(4) Pâture+\n Fauche"))
B$sbm_soil<- factor(B$sbm_soil, levels =c(  "(2) Salés+\n Fertiles","(3) peu Salés+\n Fertiles", "(1) Salés+\n peu Fertiles",
                                             "(4) peu Salés+\n peu Fertiles"))

colnames(B)[ c(2,3,5)] = c("Végétale","Agricole", "Sol")
alluvial(B[, c(3,2,5)],
         freq = B$Freq,
         col = B$col_plant_sbm,
         alpha = 0.7,
         cex = 1.2,
         cex.axis = 1.5,
         xw = 0.2)

## NMI ----

nmi_land_use = NMI(c1 =alluvial_site_df$clust_smb_grid_plant ,c2 = alluvial_site_df$clust_land_use)
nmi_soil = NMI(c1 =alluvial_site_df$clust_smb_grid_plant ,c2 = alluvial_site_df$sbm_soil)
nmi_landscape = NMI(c1 =alluvial_site_df$clust_smb_grid_plant ,c2 = alluvial_site_df$clust_landscp)
nmi_land_use
nmi_soil
nmi_landscape

## NMI randomization ----
library(bipartite)
library(doParallel)
sapply(1:1000,function(x) sample(alluvial_site_df$clust_smb_grid_plant)) -> rowperm_nmi
trials = 1000


cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)
if(file.exists("outputs/permutation_plant.Rdata")){
  load(file = "outputs/permutation_plant.Rdata")
}else{
  permutation_plant <- foreach(icount(trials),
                               .packages = c("tidyverse", "sbm","bipartite")) %dopar% {
                                 
                                 abund_plant_grid%>%
                                   swap.web(N = 1 ) -> temp
                                 
                                 temp[[1]]%>%
                                   as.data.frame()%>%
                                   set_names(colnames(abund_plant_grid))%>%
                                   mutate(grid_code = rownames(abund_plant_grid))%>%
                                   column_to_rownames("grid_code")%>%
                                   as.matrix() ->temp
                                 
                                 estimateBipartiteSBM(temp,
                                                      model = 'poisson', 
                                                      dimLabels = c(row = "Grid", col = "Plant"),
                                                      estimOptions  = list(plot = F,exploreMax = 11),
                                 ) -> SBM_abund_plant_rand
                                 
                                 
                                 list(SBM_abund_plant_rand$memberships$Grid, SBM_abund_plant_rand$memberships$Plant)
                               }
  stopCluster(cl)
  save(permutation_plant ,file = "outputs/permutation_plant.Rdata")
  
}


permutation_plant_gridclust = foreach(i = icount(trials), .combine = cbind) %do% {
    permutation_plant[[i]][[1]]}
data.frame(row_land_use  = apply(rowperm_nmi,2 ,function(x) NMI(c1 = x ,c2 = alluvial_site_df$clust_land_use)) ,
           row_soil = apply(rowperm_nmi,2, function(x) NMI(c1 = x ,c2 = alluvial_site_df$sbm_soil)) ,
           row_landscape = apply(rowperm_nmi, 2, function(x) NMI(c1 = x ,c2 = alluvial_site_df$clust_landscp)),
           ntw_land_use = apply(permutation_plant_gridclust, 2, function(x) NMI(c1 = x ,c2 = alluvial_site_df$clust_land_use)) ,
           ntw_soil = apply(permutation_plant_gridclust, 2, function(x) NMI(c1 = x,c2 = alluvial_site_df$sbm_soil)) ,
           ntw_landscape = apply(permutation_plant_gridclust, 2, function(x) NMI(c1 = x,c2 = alluvial_site_df$clust_landscp)))-> nmi_rand_df

nmi_rand_df%>%
  pivot_longer(everything(), names_to = "rand_type", values_to = "nmi_rand_values")%>%
  group_by(rand_type)%>%
  left_join(data.frame(nmi_values = c(nmi_land_use,nmi_land_use, nmi_landscape, nmi_landscape, nmi_soil, nmi_soil),
                       rand_type = c("row_land_use","ntw_land_use", "row_landscape", "ntw_landscape", "row_soil", "ntw_soil")),
            by = "rand_type")-> nmi_rand_df

plot_names <- c("Edge x land use", "Edge x landscape","Edge x soil", "Row x land use", "Row x landscape", "Row x soil"  )

ggplot(nmi_rand_df)+
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

#Multivar
library(vegan)
library(FactoMineR)
library(doParallel)

#Graph
library(igraph)
library(incidentally)

source("code/functions_modified.R")


str(alluvial_site_df)
Y = alluvial_site_df$clust_smb_grid_plant%>%as.factor()
X = alluvial_site_df%>%dplyr::select(clust_landscp,clust_land_use, sbm_soil)%>%
  mutate_all(as.factor)

abund_plant_grid

## Variance partitioning with created function ----
test_env_ana = analysis.function.alt(incid = abund_plant_grid%>%t(),
                                     memb_obs = Y,
                                     var1 = X$sbm_soil, 
                                    
                                     var3 =X$clust_land_use,
                                     var2 = X$clust_landscp,
                                     
                                     NPERM = 1000,
                                     configs = NULL)

full_cca = cca(tab.disjonctif(Y)~sbm_soil + clust_landscp+clust_land_use, data = X)

# Some ouputs from cca
predict(full_cca)%>%round(2)
plot(full_cca)
test_env_ana$R2_proxy = test_env_ana$chi2/full_cca$tot.chi

test_env_ana%>%
  as.data.frame()%>%
  mutate(across(where(is.numeric),~round(.,3)))%>%
  knitr::kable()
anova.cca(full_cca,by = "margin")

## Variance partitioning ----
mod_varpart_chi = varpart(tab.disjonctif(as.factor(Y)),~sbm_soil,~clust_land_use, ~clust_landscp,
                          data=  X, chisquare=T, permutations = how(nperm = 10000))

plot(mod_varpart_chi,Xnames = c(list("Sol"), list("Usage agri."), list("Paysage")),
     bg = c("hotpink","skyblue","forestgreen"),
     lwd = 2, cex = 1.5)

mod2_varpart_chi = varpart(tab.disjonctif(as.factor(Y)),~sbm_soil,~clust_land_use,
                          data=  X, chisquare=T)
plot(mod2_varpart_chi,Xnames = c(list("Soil"), list("Land use")),
     bg = c("hotpink","skyblue"))
summary(mod_varpart_chi)

## list of plant cluster ----

reshape2::melt(abund_plant_grid[grid_id$id,plant_id$id])%>%
  dplyr::rename(Grid_code = "Var1", Plants = "Var2" )%>%
  full_join(alluvial_site_df, by = "Grid_code")-> Grid_profil

Grid_profil%>%
  dplyr::select(Grid_code, Plants, value, clust_smb_grid_plant)%>%
  group_by(clust_smb_grid_plant ,Plants)%>%
  summarise(grp_size = n(),
            n = sum(value)) %>%
  filter(n!=0)%>%
  # mutate(freq = n / grp_size)%>%
  slice_max(n,n=50) -> temp

for(i in 1:7){
  temp%>%filter(clust_smb_grid_plant  == i)%>%print(n = 20)%>%knitr::kable()
}

## map of plant cluster ----
library(sf)
library(ggspatial) 

st_read("data/shapefiles/crop_shapefile.shp") -> soil_occup_crop

metadata_grid%>%
  dplyr::select(X,Y,Grid_code)%>%
  right_join(alluvial_site_df, by ="Grid_code") %>%
  st_as_sf(coords = c("X", "Y"), crs = 4326) %>%
  st_transform(grid_pos_metadata, crs = st_crs(soil_occup_crop))%>%
  mutate(X = as.numeric(st_coordinates(geometry)[, 1]),
         Y =  as.numeric(st_coordinates(geometry)[, 2]))-> metadata_grid_grp

metadata_grid_grp$clust_smb_grid_plant

ggplot(soil_occup_crop)+
  geom_sf(col = "gray", fill = "white")+
  geom_point(data = metadata_grid_grp, aes(X,Y, fill =clust_smb_grid_plant),
             cex = 6, alpha = 0.7, shape = 21 )+
  #geom_text_repel(data = metadata_grid,aes(X,Y, label = 1:41),  cex = 3.5,
  #                box.padding = 0.5)+
  labs(fill = "Classification communautés végétales")+
  scale_fill_brewer(palette = "Pastel1") +
  theme_void()+
  theme(legend.title = element_text(colour = "black", size=20,
                                    hjust =0.5),
        legend.text = element_text(colour = "black", size=16))

