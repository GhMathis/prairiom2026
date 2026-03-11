# Packages ----
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(vegan)
library(car)
library(GGally)
library(corrplot)
library(aricode)
library(alluvial)
library(RColorBrewer)
library(sbm)

# Load data and setup function ----
read.table("data/data_clean/Metadata_grid_CAM.txt", header  = T)-> metadata_grid
str(metadata_grid)
main_theme = theme_minimal()+
  theme(line = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.ticks =  element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", size=22,angle = ),
        axis.text.y = element_text(colour = "black", size=22),
        legend.title = element_text(colour = "black", size=20,
                                    hjust =0.5),
        legend.text = element_text(colour = "black", size=18),
        axis.title= element_text(size=28),
        strip.text = element_text(colour = "black", size=15, face ="italic"))


zscore_normalization <- function(data) {
  # Step 1: Rank the data
  ranked_data <- rank(data)
  n <- length(data)
  
  # Step 2: Calculate quantile positions
  quantile_positions <- (ranked_data) / (n + 1)
  
  # Step 3: Transform quantiles to z-scores
  mean_std_normal <- qnorm(quantile_positions)
  
  return(mean_std_normal)
}



metadata_grid%>%
  dplyr::select(pHwater, MO, Phos, K, Mg, Ca, Na, N, C, 
                CN, clay, SiltF, SiltC, SandF, SandC, Cl, Res, Cond,
                dep_oxy_num , Pature, Fauche, natural_landscape,
                cultivated, artificial, non_emitting, wetland, Grid_code) %>%
  #remove_rownames()%>%
  #column_to_rownames("Grid_code")%>%
  mutate(across(c(Fauche, Pature), as.factor),
         across(c(pHwater, MO, Phos, K, Mg, Ca, Na, N, C, 
                  CN, clay, SiltF, SiltC, SandF, SandC, Cl, Res, Cond,
                  dep_oxy_num ),  ~decostand(., method = "standardize")[,1]),
         across(c(natural_landscape, cultivated, artificial,
                  non_emitting, wetland),~logit(.)))%>%
  filter(!str_detect(Grid_code, "22_CAM_15")) -> metadata_grid_standard # warning messages from the logit => normal don't worry
#metadata_grid_standard%>%ggpairs

# Soil cluster ----

# z-score on ranks
metadata_grid%>%
  dplyr::select(pHwater, MO, Phos, K, Mg, Ca, Na, N, C, 
         CN, clay, SiltF, SiltC, SandF, SandC, Cl, Res, Cond,
         dep_oxy_num,Grid_code)%>%
  filter(!str_detect(Grid_code, "22_CAM_15"))%>%
  remove_rownames()%>%
  column_to_rownames("Grid_code")-> soil_only
soil_only%>%
  mutate_all(~zscore_normalization(.))%>%t() -> soil_standard_rank

# z-score on values
metadata_grid_standard%>%
  dplyr::select(pHwater, MO, Phos, K, Mg, Ca, Na, N, C, 
         CN, clay, SiltF, SiltC, SandF, SandC, Cl, Res, Cond,
         dep_oxy_num)%>%t() -> soil_standard_zscore 

#### SBM soil
par(mfrow =c(1,1))
soil_standard_rank%>%
  as.matrix() %>%
  estimateBipartiteSBM(
    model = 'gaussian', 
    dimLabels = c(row = "Soil", col = "Grid")) -> smb_soil

sbm_soil = as.factor(smb_soil$memberships$Grid) ##### Grids cluster for soil

plot(smb_soil)
rownames(soil_standard_zscore)
soil_standard_zscore%>%
  as.data.frame()%>%
  rownames_to_column("soil")%>%
  
  mutate(soil = c("pH", "MO", "Phos", "K", "Mg", "Ca", "Na",  "N", "C", "CN", "Argile",
           "Limon F.", "Limon G.", "Sable F.", "Sable G.",
           "Cl", "Résistivité", "Conductivité", "Profondeur nappe"))%>%
  column_to_rownames("soil") -> soil_standard_zscore

#standadized matrix, not arrange
soil_standard_zscore%>%
  as.matrix()%>%
  plotMyMatrix(dimLabels = c("Soil", "Grid"), plotOptions = list(rowNames = T,colNames = T))+
  theme(
    legend.position="none",
    axis.text.x =element_text(angle = 45, vjust = 1, hjust = 1, size = 8, face="italic"),
    axis.text.y = element_text(size = 8,face="italic"))+
  theme(strip.text.x = element_text(colour = "gray", size=15, face ="italic"),
        strip.text.y  = element_text(colour = "gray", size=15, face ="italic"),
        strip.background = element_rect(fill="white"),
        panel.border = element_rect(color = "black",  fill = NA))

grid_id = smb_soil$memberships[2]%>%
  as.data.frame()%>%
  rownames_to_column("id")%>%
  mutate(id = as.numeric(id))%>%
  arrange(Grid)

soil_id = smb_soil$memberships[1]%>%
  as.data.frame()%>%
  rownames_to_column("id")%>%
  mutate(id = as.numeric(id))%>%
  arrange(Soil)


soil_standard_rank[soil_id$id,grid_id$id]%>%
  as.matrix()%>%
  plotMyMatrix(dimLabels = c("Soil", "Grid"), plotOptions = list(rowNames = T,colNames = T))+
  theme(
    legend.position="none",
    axis.text.x =element_text(angle = 45, vjust = 1, hjust = 1, size = 6, face="italic"),
    axis.text.y = element_text(size = 5,face="italic"))+
  geom_vline(xintercept = which(!duplicated(grid_id$Grid))[-1]-0.5, col ="red")+
  geom_hline(yintercept = abs(which(!duplicated(soil_id$Soil))[-1]-nrow(soil_id)-1)+0.5, col = "red", linetype = 5)+
  theme(strip.text.x = element_text(colour = "gray", size=15, face ="italic"),
        strip.text.y  = element_text(colour = "gray", size=15, face ="italic"),
        strip.background = element_rect(fill="white"),
        panel.border = element_rect(color = "black",  fill = NA))


soil_standard_zscore[soil_id$id,grid_id$id]%>%
  as.matrix()%>%
  plotMyMatrix(dimLabels = c("Varaibles physico-chimiques", "Grilles"), plotOptions = list(rowNames = T,colNames = T))+
  theme(
    legend.position="none",
    axis.text.x =element_text(angle = 45, vjust = 1, hjust = 1, size = 8, face="italic"),
    axis.text.y = element_text(size = 8,face="italic"))+
  geom_vline(xintercept = which(!duplicated(grid_id$Grid))[-1]-0.5, col ="red")+
  geom_hline(yintercept = abs(which(!duplicated(soil_id$Soil))[-1]-nrow(soil_id)-1)+0.5, col = "red", linetype = 5)+
  annotate("text", x =  c(1,which(!duplicated(grid_id$Grid))[-1]),y = 1.2, 
           label= c("1","2","3","4"),colour = "white", size =8)+
  theme(strip.text.x = element_text(colour = "gray", size=15, face ="italic"),
        strip.text.y  = element_text(colour = "gray", size=15, face ="italic"),
        strip.background = element_rect(fill="white"),
        panel.border = element_rect(color = "black",  fill = NA))


# Landscape cluster ----

metadata_grid_standard%>%
  dplyr::select(natural_landscape, cultivated, artificial, non_emitting, wetland) -> landscp_only

landscp_only%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column("landscp")%>%
  mutate(landscp = c("Naturels secs", "Cultivés", "Anthropisés", 
                     "Non emmeteurs", "Naturels humides"))%>%
  column_to_rownames("landscp")%>%
  t()%>%
  as.matrix()%>%
  plotMyMatrix(dimLabels = c("Grid", "Occupation des sols"),
               plotOptions = list(rowNames = T,colNames = T))+
  theme(
    legend.position="none",
    axis.text.x =element_text(angle = 45, vjust = 1, hjust = 1, size = 8, face="italic"),
    axis.text.y = element_text(size = 8,face="italic"))+
  theme(strip.text.x = element_text(colour = "gray", size=15, face ="italic"),
        strip.text.y  = element_text(colour = "gray", size=15, face ="italic"),
        strip.background = element_rect(fill="white"),
        panel.border = element_rect(color = "black",  fill = NA))

ggpairs(landscp_only)

## clustering
landscp_only_named = landscp_only
colnames(landscp_only_named) = c("Antropized", "Cultivated", "Natural non-wet",
                                 "No propagules emitter", "Natural\n wet")
PCA_landscape = PCA(landscp_only_named,scale.unit = F)
HCPC_landscape =PCA_landscape%>%
  HCPC(-1)
PCA_landscape$call

clust_landscp = HCPC_landscape$data.clust$clust

fviz_pca_biplot(PCA_landscape,col.ind = as.factor(HCPC_landscape$data.clust$clust))+
  scale_colour_manual(values = c(brewer.pal(4, "Set1")))+
  labs(col = "Land use\n classification", shape =   "Land use\n classification")+
  main_theme

## landscape groups profiles
metadata_grid %>%
  filter(!str_detect(Grid_code, "22_CAM_15")) %>%
  dplyr::select(natural_landscape, cultivated, artificial, non_emitting, wetland, Grid_code, clust_landscp) %>%
  mutate(clust_landscp = as.factor(clust_landscp)) %>%
  pivot_longer(-c(Grid_code, clust_landscp), values_to = "landscp_val", names_to = "landscp_name") %>%
  ggplot() +
  facet_wrap(~clust_landscp, labeller = labeller(clust_landscp = c(
    `1` = "Cultivated \n dominant",
    `2` = "Natural non-wet \n dominant",
    `3` = "Mixed  land use",
    `4` = "Natural wet \n dominant"
  ))) +
  geom_col(aes(x = Grid_code, y = landscp_val, fill = landscp_name)) +
  scale_fill_manual(values = c("#d73027", "#66c2a5", "#fee08b", "#cfd8dc", "#4575b4"),
                    labels = c("Antropized", "Cultivated", "Natural non-wet",
                               "No propagules emitter", "Natural wet")) +
  main_theme +
  theme(axis.text.x = element_blank()) +
  labs(x = "Grids", y = "Percentage of Land-use", fill = "Land-use types")

# land_use (agricultural usages)----

metadata_grid_standard%>%
  dplyr::select(Fauche,Pature) -> land_use_only

land_use_only%>%
  as.matrix()%>%
  plotMyMatrix(dimLabels = c("Grid", "Usages agricoles"),
               plotOptions = list(rowNames = T,colNames = T))+
  theme(
    legend.position="none",
    axis.text.x =element_text(angle = 45, vjust = 1, hjust = 1, size = 8, face="italic"),
    axis.text.y = element_text(size = 8,face="italic"))+
  theme(strip.text.x = element_text(colour = "gray", size=15, face ="italic"),
        strip.text.y  = element_text(colour = "gray", size=15, face ="italic"),
        strip.background = element_rect(fill="white"),
        panel.border = element_rect(color = "black",  fill = NA))

land_use_only%>%
  mutate(clust_land_use = case_when(Pature == 0 & Fauche ==0 ~2 ,
                                    Pature == 0 & Fauche ==1 ~3 ,
                                    Pature == 1 & Fauche ==0 ~1 ,
                                    Pature == 1 & Fauche ==1 ~4 )) -> land_use_only

clust_land_use = land_use_only$clust_land_use

land_use_only%>%
  dplyr::select(Fauche,Pature)%>%
  table()

# Soil sbm ----

NMI(sbm_soil,clust_landscp)
NMI(sbm_soil,clust_land_use)
NMI(clust_landscp,clust_land_use)

freq_df<-
  as.data.frame(
    table(
      sbm_soil,
      clust_landscp,
      clust_land_use)
  )


colnames(freq_df) = c("sbm_soil",  "landscp", "land_use","Freq")
str(freq_df)
w   <- which(freq_df$Freq != 0)
freq_df <- freq_df[w, ]
col_clust = brewer.pal(4, "Set1")
alluvial(freq_df[, c(1, 2, 3)],
         freq = freq_df$Freq,
         col = case_when(freq_df$sbm_soil == 1 ~ col_clust[1],
                         freq_df$sbm_soil == 2 ~ col_clust[2],
                         freq_df$sbm_soil == 3 ~ col_clust[3],
                         .default = col_clust[4]))

if(any(colnames(metadata_grid)%in%c("clust_landscp", "clust_land_use", "sbm_soil"))){#add columns only if they are not present in the df (avoid duplicated columns)
  cat("Metadata already contain classifications")
}else{
  metadata_grid%>%
    mutate(clust_landscp = clust_landscp,
           clust_land_use = clust_land_use,
           sbm_soil = sbm_soil)-> metadata_grid
    write.table(metadata_grid, file = "data/data_clean/Metadata_grid_CAM.txt")
}
metadata_grid%>%
  mutate(clust_landscp = as.factor(clust_landscp),
         clust_land_use = as.factor(clust_land_use),
         sbm_soil = as.factor(sbm_soil)) -> metadata_grid

# map of groups ----
st_read("data/shapefiles/crop_shapefile.shp") -> soil_occup_crop
metadata_grid%>%
  dplyr::select(X,Y,Grid_code)%>%
  st_as_sf(coords = c("X", "Y"), crs = 4326) %>%
  st_transform(grid_pos_metadata, crs = st_crs(soil_occup_crop))%>%
  mutate(X = as.numeric(st_coordinates(geometry)[, 1]),
         Y =  as.numeric(st_coordinates(geometry)[, 2]))-> metadata_grid_grp

## Soil ----
ggplot(soil_occup_crop)+
  geom_sf(col = "gray", fill = "white")+
  geom_point(data = metadata_grid_grp, aes(X,Y, fill =sbm_soil),
             cex = 6, alpha = 0.7, shape = 21 )+
  #geom_text_repel(data = metadata_grid,aes(X,Y, label = 1:41),  cex = 3.5,
  #                box.padding = 0.5)+
  labs(fill = "Classification types de Sols")+
  scale_fill_brewer(palette = "Pastel1",
                    labels = c("Salés/Argileux +\n Nappe surface +  Fertiles",
                                                  "Salés/Argileux +\n  Nappe profonde +  peu Fertiles",
                                                  "peu Salés/Sableux +\n  Nappe profonde +  Fertiles", 
                                                  "peu Salés/Sableux +\n  Nappe profonde +  peu Fertiles")) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl", width_hint = 0.5,text_cex = 1) +
  theme_void()+
  theme(legend.title = element_text(colour = "black", size=20,
                                    hjust =0.5),
        legend.text = element_text(colour = "black", size=16))

## Landscape ----
ggplot(soil_occup_crop)+
  geom_sf(col = "gray", fill = "white")+
  geom_point(data = metadata_grid_grp, aes(X,Y, fill =clust_landscp),
             cex = 6, alpha = 0.7, shape = 21 )+
  #geom_text_repel(data = metadata_grid,aes(X,Y, label = 1:41),  cex = 3.5,
  #                box.padding = 0.5)+
  labs(fill = "Land use classification")+
  scale_fill_brewer(palette = "Pastel1",
                    labels = c("Cultivated dominant",
                               "Natural non-wet dominant",
                               "Mixed land use", 
                               "Natural wet dominant")) +
  theme_void()+
  theme(legend.title = element_text(colour = "black", size=20,
                                    hjust =0.5),
        legend.text = element_text(colour = "black", size=16))
c("Antropized", "Cultivated", "Natural non-wet",
  "No propagules emitter", "Natural wet")
## Agricultural usages ----
ggplot(soil_occup_crop)+
  geom_sf(col = "gray", fill = "white")+
  geom_point(data = metadata_grid_grp, aes(X,Y, fill =as.factor(clust_land_use)),
             cex = 6, alpha = 0.7, shape = 21 )+
  #geom_text_repel(data = metadata_grid,aes(X,Y, label = 1:41),  cex = 3.5,
  #                box.padding = 0.5)+
  
  labs(fill = "Classification Usages agricoles")+
  scale_fill_brewer(palette = "Pastel1",
                    labels = c("Pâture",
                               "Rien",
                               "Fauche", 
                               "Pâture + Fauche")) +
  theme_void()+
  theme(legend.title = element_text(colour = "black", size=20,
                                    hjust =0.5),
        legend.text = element_text(colour = "black", size=16))
