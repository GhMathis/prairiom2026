library(tidyverse)
library(readxl)
library(vegan)

main_theme = theme_bw()+
  theme(line = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.ticks =  element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", size=15),
        axis.text.y = element_text(colour = "black", size=15),
        axis.title.x = element_text(colour = "black", size=18),
        axis.title.y = element_text(colour = "black", size=18),
        legend.title = element_text(colour = "black", size=20),
        legend.title.align=0.5,
        legend.text = element_text(colour = "black", size=18),
        axis.title=element_text(size=10),
        strip.background = element_rect(fill="cornsilk"))

# OTU Camargue ----

## PLANT ----
read_xlsx("data/OTU_plant.xlsx")%>%
  filter(str_detect( Host_code, "CAM"),
         !str_detect(Host_code, "21_CAM_1[45]")) %>%
  replace(is.na(.),0)%>%
  select(-Richness)%>%
  select(Host_code, where(~is.numeric(.x) && sum(.x) != 0 )) -> otu_plant_cam

write.table(otu_plant_cam, "data/data_clean/OTU_plant_CAM.txt",row.names = F)

otu_plant_cam%>%
  mutate( Grid_code = str_extract(Host_code, ".._CAM_.."))%>%
  group_by(Grid_code) %>%
  mutate(across(where(is.numeric), ~.*2))%>%
  summarise_if(is.numeric, sum)-> abund_plant_grid 


write.table(abund_plant_grid, "data/data_clean/abund_plant_grid.txt",row.names = F)

## VIRUS ----
read.delim2("data/S1_Viral_OTU.txt",header = T)%>%
  rename(Host_code = "X")%>%
  filter(str_detect( Host_code, "CAM"),
         !str_detect(Host_code, "21_CAM_1[45]")) %>%
  select(Host_code, where(~is.numeric(.x) && sum(.x) != 0 ))-> otu_virus_cam

write.table(otu_virus_cam, "data/data_clean/OTU_virus_CAM.txt",row.names = F)

## Fungi ----

load("data/fungi_for_mathis.Rdata")

fungi_OTU%>%
  mutate(sp = fungi_TAX$species)%>%
  group_by(sp)%>%
  summarise(across(where(is.numeric), sum))%>%
  column_to_rownames("sp")-> fungi_OTU_sp

## unidentified proportion of fungi
fungi_OTU_sp%>%
  rownames_to_column("sp")%>%
  pivot_longer(-sp)%>%
  group_by(sp)%>%
  summarise(abond = sum(value))%>%
  mutate(identified = case_when(str_detect(sp,"unidentified") ~ 0,
                                !str_detect(sp,"unidentified") ~1))%>%
  ggplot()+
  facet_wrap(~identified)+
  geom_histogram(aes(log10(abond+1)))

## remove unidentified at species level
fungi_OTU_sp%>%
  rownames_to_column("sp")%>%
  filter(!str_detect(sp,"unidentified"))%>%
  column_to_rownames("sp") %>%
  rownames_to_column("sp")%>%
  pivot_longer(-sp)%>%
  mutate(name = str_replace_all(name, "-","_"))%>%
  filter(!str_detect(name, "22_CAM_15...."),
         !str_detect(name, "21_CAM_1[45]...."))%>%
  pivot_wider(names_from = "name",values_from ="value", values_fill = 0)%>%
  column_to_rownames("sp") -> fungi_OTU_sp_filter

sort(apply(fungi_OTU_sp, 1, sum), decreasing = T)
sort(apply(fungi_OTU_sp_filter, 1, sum), decreasing = T)
sort(apply(fungi_OTU_sp_filter, 2, sum))

nrow(fungi_OTU_sp_filter) # n species

### Filter and rarefaction ----

fungi_OTU_sp_filter%>%
  rownames_to_column("sp")%>%
  pivot_longer(-sp)%>%
  group_by(name)%>%
  summarise(value_sum = sum(value))%>%
  mutate(grid_code = str_extract(name, ".._CAM_.."))%>%
  arrange(desc(value_sum)) -> tmp2

tmp2$name -> grid_list

ggplot(tmp2)+
  geom_point(aes(name, log10(value_sum+1), col = grid_code))+
  scale_x_discrete(limits = grid_list)

## filter
fungi_OTU_sp_filter%>%
  rownames_to_column("sp")%>%
  pivot_longer(-sp)%>%
  group_by(name)%>%
  mutate(tot_reads = sum(value),
         grid_code = str_extract(name, ".._CAM_.."))%>%
  filter(tot_reads > 10^3.5) %>%
  ungroup()%>%
  dplyr::select(-tot_reads)%>%
  group_by(grid_code, sp)%>%
  summarise(value_sum = sum(value))%>%
  group_by(sp)%>%
  ungroup()%>%
  pivot_wider(names_from = "grid_code",values_from ="value_sum", values_fill = 0)%>%
  column_to_rownames("sp")%>%
  as.data.frame()-> fungi_OTU_sp_grid_filter2

if(!file.exists("data/data_clean/OTU_fungi.txt")){
  write.table(fungi_OTU_sp_grid_filter2, file = "data/data_clean/OTU_fungi.txt")
}
str(fungi_OTU_sp_grid_filter2)
hist(log1p(apply(fungi_OTU_sp_grid_filter2, 1, sum)))  

# otu fungi table for sbm
fungi_OTU_sp_filter%>%
  rownames_to_column("sp")%>%
  pivot_longer(-sp)%>%
  group_by(name)%>%
  mutate(tot_reads = sum(value),
         grid_code = str_extract(name, ".._CAM_.."))%>%
  filter(tot_reads > 10^3.5) %>%
  ungroup()%>%
  dplyr::select(-tot_reads)%>%
  group_by(grid_code, sp)%>%
  summarise(value_sum = sum(value))%>%
  group_by(sp)%>%
  filter(sum(value_sum)>10)%>%
  ungroup()%>%
  pivot_wider(names_from = "grid_code",values_from ="value_sum", values_fill = 0)%>%
  column_to_rownames("sp")%>%
  as.matrix()-> fungi_OTU_sp_grid_filter3

nrow(fungi_OTU_sp_grid_filter3)

if(!file.exists("data/data_clean/OTU_fungi_for_sbm.txt")){
  write.table(fungi_OTU_sp_grid_filter3, file = "data/data_clean/OTU_fungi_for_sbm.txt")
}

# Metadata ----
## Metadata Quadra ----
read_xlsx("data/Metadata_Sample.xlsx")%>%
  filter(str_detect(Locality, "Arles"))%>%
  select(-c(Rocks, Additional_information))-> metadata_quad_cam
str(metadata_quad_cam)

## add information to quadra metadata

read.table("data/data_clean/OTU_plant_CAM.txt",header = T)%>%
  reframe(Host_code = Host_code,
         year = as.factor(substring(Host_code,1,2)),
         grid = as.factor(substring(Host_code,8,9)),
         quadra = as.factor(substring(Host_code,10,11)),
         total_cover = rowSums(across(where(is.numeric))),
         Plant_richness = rowSums(across(where(is.numeric), ~.x != 0))
         
  )%>%
  select(Plant_richness,year,grid,quadra, total_cover,Host_code)%>%
  left_join(metadata_quad_cam, by ="Host_code")%>%
  mutate(across(c(Vegetation, Litter, Soil),~as.numeric(.x)),
         Soil = case_when(str_detect(Host_code, "20_CAM_01")~ 0,
         .default = Soil),
         Free_space = Litter +Soil)%>%
  select(-Host_community)-> metadata_quad_cam

write.table(metadata_quad_cam, "data/data_clean/Metadata_quadra_CAM.txt")

## Metadata grid (joint all metadata at grid level) ----


read_xlsx("data/Habitats_EDGG_patures.xlsx")%>%
  select(Grid_code = Sampling_n, Pature,Fauche)-> paturage


read_xlsx("data/donnee_sol.xlsx",col_names = F, skip = 3)%>%
  rename_with(~c("Grid_code", "Num", "pHwater", "lime_tot", "MO", "Phos", "K", "Mg",
                 "Ca", "Na", "N","C","CN","clay",  "SiltF", "SiltC","SandF", "SandC",
                 "Cl", "Res", "Cond"))%>%
  filter(str_detect(Grid_code,"CAM" ))%>%
  mutate(Grid_code = str_extract(Grid_code,".._CAM_.." ))%>%
  left_join(y = read_xlsx("data/distance_fer_sol_EDGG_CAM.xlsx"), by = join_by(Grid_code == EDGG))%>%
  mutate(dep_oxy_num = as.numeric(str_remove(Profondeur, "[>]")),
         depth_oxy = cut(dep_oxy_num,
                         breaks = c(0, 9, 19, 29, 39, Inf),
                         labels = c("0-9", "10-19", "20-29", "30-39", ">40"),
                         include.lowest = TRUE)) -> soil_data

read_xlsx("data/Metadata_grid.xlsx")%>%
  filter(Locality == "Arles") %>%
  mutate(
    X = as.numeric(str_split_fixed(GPS_localisation, pattern = ",", n = 2)[, 2]),
    Y = as.numeric(str_split_fixed(GPS_localisation, pattern = ",", n = 2)[, 1])
  )%>%
  select(-GPS_localisation)%>%
  left_join(soil_data, by= "Grid_code")%>%
  left_join(paturage, by = "Grid_code")-> metadata_grid_cam

## Diversity ----
## Virus----
### Lire et formater les données ----

read.table("data/data_clean/OTU_virus_CAM.txt",header = T)%>%
  mutate( Grid_code = str_extract(Host_code, ".._CAM_.."))%>%
  group_by(Grid_code) %>%
  mutate(across(where(is.numeric), ~ifelse(. != 0,1,0 ))) %>%
  summarise_if(is.numeric, sum)%>%
  pivot_longer(-Grid_code, names_to = "Virus", values_to = "occur")%>%
  group_by(Grid_code)%>%
  filter(sum(occur)!=0)%>%
  pivot_wider(values_from = occur, names_from = Grid_code )%>%
  column_to_rownames("Virus") -> OUT_virus_grid_binary_df

rbind(n = rep(9,ncol(OUT_virus_grid_binary_df)),OUT_virus_grid_binary_df) -> inext_formatted_data

## Estimations de diversité ----
if(file.exists("outputs/inext_estimates_of_viral_diversity.Rdata")){
  load("outputs/inext_estimates_of_viral_diversity.Rdata")
  
}else{
  z <- iNEXT(inext_formatted_data, q=c(0,1), datatype="incidence_freq", se=TRUE)
  z
  save(inext_formatted_data, z, file = "outputs/inext_estimates_of_viral_diversity.Rdata")
}

z$iNextEst$size_based %>% 
  filter(Order.q==0) %>% 
  mutate(qD_obs = if_else(Method == "Observed", qD, NA),
         qD_rar = if_else(Method == "Extrapolation", NA, qD),
         qD_ext = if_else(Method == "Rarefaction", NA, qD)) %>% 
  ggplot(., aes(x = t, y = qD, group = Assemblage))+
  geom_point(aes(x = t, y = qD_obs)) +
  geom_line(aes(x = t, y = qD_rar, group = Assemblage), linewidth = 0.5) +
  geom_line(aes(x = t, y = qD_ext, group = Assemblage), linetype = 2, linewidth = 0.5)+
  geom_ribbon(aes(ymin=qD.LCL, ymax=qD.UCL), alpha = 0.05, linetype = 0)+
  ylab("Virus richness")+
  xlab("Sample size")+
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    legend.title = element_text(size=12, face="bold"),
    legend.text = element_text(size=12, face="italic")
  )

z$iNextEst$size_based %>% 
  filter(Order.q==0 & Method == "Observed") %>% 
  dplyr::select(Assemblage, qD, qD.LCL, qD.UCL) %>% 
  rename(Grid_code = Assemblage, 
         virus_richness = qD, 
         v_rich_LCL = qD.LCL, 
         v_rich_UCL = qD.UCL)-> virus_richness

z$iNextEst$size_based %>% 
  filter(Method == "Observed" & Order.q == 1) %>% 
  dplyr::select(Assemblage, qD, qD.LCL, qD.UCL) %>% 
  rename(Grid_code = Assemblage, 
         virus_H = qD, 
         v_H_LCL = qD.LCL, 
         v_H_UCL = qD.UCL) -> virus_H

z$AsyEst %>% 
  filter(Diversity == "Species richness") %>% 
  ggplot(., aes(x=Observed, y = Estimator))+
  geom_point()+
  main_theme

## Plants----
### Lire et formater les données ----
read.table("data/data_clean/OTU_plant_CAM.txt", header = T)%>%
  mutate( Grid_code = str_extract(Host_code, ".._CAM_.."))%>%
  group_by(Grid_code) %>%
  mutate(across(where(is.numeric), ~ifelse(. != 0,1,0 ))) %>%
  summarise_if(is.numeric, sum)%>%
  pivot_longer(-c(Grid_code), names_to = "Plant", values_to = "occur")%>%
  pivot_wider(values_from = occur, names_from = Grid_code )%>%
  column_to_rownames("Plant")-> OUT_plant_grid_binary_df

rbind(n = rep(9,42),OUT_plant_grid_binary_df) -> OUT_plant_grid_binary_df

read.table("data/data_clean/abund_plant_grid.txt", header = T)%>%
  pivot_longer(-c(Grid_code), names_to = "Plant", values_to = "occur")%>%
  pivot_wider(values_from = occur, names_from = Grid_code )%>%
  column_to_rownames("Plant")-> OUT_plant_grid_abund
rbind(n = rep(1800,length(colnames(OUT_plant_grid_abund))),OUT_plant_grid_abund) -> plant_inext_formatted_data

### Estimations de diversité ----
if(file.exists("outputs/plant_diversity_estimates_with_inext.Rdata")){
  load(file = "outputs/plant_diversity_estimates_with_inext.Rdata")
  
}else{
  
  y <- iNEXT(plant_inext_formatted_data, q=c(0,1), datatype="incidence_freq", se=T)
  y
  save(plant_inext_formatted_data, y, file = "outputs/plant_diversity_estimates_with_inext.Rdata")
  
}

y$iNextEst$size_based %>% 
  filter(Order.q==0) %>%
  mutate(qD_obs = if_else(Method == "Observed", qD, NA),
         qD_rar = if_else(Method == "Extrapolation", NA, qD),
         qD_ext = if_else(Method == "Rarefaction", NA, qD)) %>% 
  ggplot(., aes(x = t, y = qD, group = Assemblage))+
  geom_point(aes(x = t, y = qD_obs)) +
  geom_line(aes(x = t, y = qD_rar, group = Assemblage), linewidth = 0.5) +
  geom_line(aes(x = t, y = qD_ext, group = Assemblage), linetype = 2, linewidth = 0.5)+
  geom_ribbon(aes(ymin=qD.LCL, ymax=qD.UCL), alpha = 0.05, linetype = 0)+
  ylab("Plant richness")+
  xlab("Sample size")+
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    legend.title = element_text(size=12, face="bold"),
    legend.text = element_text(size=12, face="italic")
  )

y$iNextEst$size_based %>% 
  filter(Method == "Observed" & Order.q == 0) %>% 
  dplyr::select(Assemblage, qD, qD.LCL, qD.UCL) %>% 
  rename(Grid_code = Assemblage, 
         plant_richness = qD, 
         pl_rich_LCL = qD.LCL, 
         pl_rich_UCL = qD.UCL) -> plant_richness


y$iNextEst$size_based %>% 
  filter(Method == "Observed" & Order.q == 1) %>% 
  dplyr::select(Assemblage, qD, qD.LCL, qD.UCL) %>% 
  rename(Grid_code = Assemblage, 
         plant_H = qD, 
         pl_H_LCL = qD.LCL, 
         pl_H_UCL = qD.UCL) -> plant_H


## Fungi----
##read and transform data ----
fungi_OTU_sp_grid_filter2%>%
  as.matrix()-> OTU_fungi
fungi_richness= round(as.data.frame(rarefy(t(OTU_fungi), 10^3.5)))%>%
  rownames_to_column("Grid_code")
summary(fungi_richness) 
colnames(fungi_richness)[2] = "fungi_richness"

 
# merge and save metadata grid ----

if(!file.exists("data/data_clean/Metadata_grid_CAM.txt")){
  metadata_quad_cam%>%
    group_by(Grid_code) %>% 
    summarise(Biomass = sum(Biomass),
              Vegetation = 2*sum(as.numeric(Vegetation))/1800) -> biom_df


  
  metadata_grid_cam%>%
    left_join( plant_richness, by = "Grid_code")%>%
    left_join( plant_H, by = "Grid_code")%>%
    left_join( virus_richness, by = "Grid_code")%>%
    left_join( virus_H, by = "Grid_code")%>%
    left_join( biom_df, by = "Grid_code")%>%
    left_join( fungi_richness, by = "Grid_code")%>%
    
  mutate(virus_richness = replace(virus_richness, is.na(virus_richness), 0),
         virus_H = replace(virus_H, is.na(virus_H), 0),
         
         plant_richness2 = plant_richness^2, 
         p_richness_class5 = cut(plant_richness, breaks = seq(0,30, 5)),
         p_richness_class2 = cut(plant_richness, breaks = seq(0,30, 2)),
         log_virus_richness = log(virus_richness+1),
         log_plant_richness = log(plant_richness+1),
         log_plant_richness2 = (log(plant_richness+1))^2,
         log_plant_H = log(plant_H+1),
         log_virus_H = log(virus_H+1),
         log_biomass = log(Biomass+1),
         log_plant_equit = log_plant_H-log_plant_richness,
         plant_equit = plant_H/plant_richness) -> metadata_grid_cam

  str(metadata_grid_cam)

  write.table(metadata_grid_cam, "data/data_clean/Metadata_grid_CAM.txt")
}