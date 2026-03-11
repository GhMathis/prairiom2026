# Packages ----
# General
library(RColorBrewer)
library(tidyverse)
library(readxl)

# Multivar
library(GGally)
library(tidyverse)
library(corrplot)
library(vegan)
library(car)

# Geographical
library(sp)
library(sf)
library(siland)
library(terra)
library(ggspatial) # 

# theme for ggplot graph
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

# Load data ----
read.table("data/data_clean/Metadata_Grid_CAM.txt", header = T) %>%
  dplyr::select(X, Y,Grid_code, plant_richness,log_plant_richness, virus_richness,log_virus_richness,fungi_richness)%>%
  data.frame()-> metadata_grid


st_read("data/shapefiles/sql_statement_d551600.shp")  %>% 
  #rename(lib = lib4_16)%>%
  # mutate(lib1_16 = str_replace_all(lib1_16, "[,/ ']","_"),
  #        lib1_16 = str_replace_all(lib1_16, "[.]","_"),
  #        lib1_16 = str_replace_all(lib1_16, "[êéè]","e"),
  #        lib1_16 = str_remove(lib1_16, "[*]$"),
  #        lib1_16 = str_remove(lib1_16, "_$"))%>%
  mutate(lib = case_when(
    lib4_16 %in% c("Bâti diffus", "Bâti individuel dense", "Tissu urbain continu", "Bâti individuel lâche",            
                   "Bâti isolé", "Décharge",  "Bâti collectif",  "Terrain vague en zone urbaine",     
                   "Équipement sportif et de loisirs", "Bâti léger ou informel", "Zone d'activité et équipement","Extraction de matériaux",           
                   "Chantier", "Place", "Jardins familiaux", "Espace vert urbain",                
                   "Zone portuaire", "Bâti individuel dans parc paysager", "Cimetière" ) ~ "artificial",
    
    lib4_16 %in% c("Marais ouvert", "Feuillu", "Formation arbustive et arborée semi-ouverte", 
                   "Conifère",  "Formation arbustive et arborée fermée","Plage",
                   "Etang et/ou lagune", "Canal", "Forêt mélangée","Ripisylve") ~ "non_emitting",
    
    lib4_16 %in% c("Riz", "Luzerne", "Prairie temporaire", "Blé", "Tournesol", "Verger, oliveraie", "Friche récente",    
                   "Culture maraichère", "Sorgho, soja", "Vignoble", "Colza", "Maïs" ) ~ "cultivated",
    
    lib4_16 %in% c("Prairie naturelle","Coussoul", "Dune embryonnaire",                         
                   "Dune végétalisée", "Dune à végétation arbustive") ~ "natural_landscape",
    
    lib4_16 %in% c("Roselière", "Sansouire basse", "Marais ouvert", "Sansouire haute",                  
                   "Jonchaie", "Autre marais à végétation émergée", "Marais à marisque", "Sol nu",                           
                   "Lagune de pré-concentration", "Table saunante", "Friche salicole récente" ) ~ "wetland",
    
    
    .default = lib4_16
  ))%>%
  dplyr::select(lib)%>%
  st_transform( 2154)%>%
  
  mutate(area =  st_area(.), #area of each polygone
         n =1,
         lib2 = as.factor(lib))%>%
  pivot_wider(names_from = lib2, values_from = n, values_fill = 0)-> soil_occup

# Arrange the map ----
## Change positions of grids to Lambert projection (for siland package) ----
metadata_grid%>%
  st_as_sf(coords = c("X", "Y"), crs = 4326)%>%
  st_transform(grid_pos_metadata, crs = st_crs(soil_occup))%>%
  mutate(X = as.numeric(st_coordinates(geometry)[, 1]),
         Y =  as.numeric(st_coordinates(geometry)[, 2]))%>%
  st_join(.,
          soil_occup%>%
            dplyr::select(lib,area)
  )%>%
  distinct(Grid_code, .keep_all = TRUE)%>% # if duplication append in the st_joit above... just in case
  as.data.frame(.)-> metadata_grid

## setup limits and crop ----
metadata_grid%>%
  #dplyr::select(X,Y)%>%
  st_as_sf(coords = c("X", "Y"), crs = st_crs(soil_occup))%>%
  st_buffer(4000)%>%
  st_union()%>%
  st_cast(to = "POLYGON") -> limit_map

soil_occup_crop = st_crop(x = soil_occup, y =limit_map)
if(!file.exists("data/shapefiles/crop_shapefile.shp")){
  st_write(soil_occup_crop, "data/shapefiles/crop_shapefile.shp")
  }

ggplot(limit_map)+
  geom_sf()



ggplot() +
  geom_sf(data= soil_occup,col = "gray", fill = "white") +
  geom_rect(aes(xmin = 825589.5, xmax = 842084.4, ymin = 6253261.7, ymax = 6281143.1), 
            fill = NA, col = "red", linewidth = 1) +
  theme_void() +
  theme(legend.position = "none")

col_landscape =brewer.pal(5, "Set2")[c(4,2,1,3,5)]
ggplot(soil_occup_crop)+
  geom_sf(col = "gray",alpha = 0.6,aes(fill = lib) )+
  geom_point(data = metadata_grid,aes(X,Y),fill = "black",  cex = 5, alpha = 0.3, shape = 21 )+
  #geom_text_repel(data = metadata_grid,aes(X,Y, label = 1:41),  cex = 3.5,
  #                box.padding = 0.5)+
  scale_fill_manual(values = col_landscape,
                    labels = c("Antropisés", "Cultivés", "Naturels non humide",
                               "Non emmeteur de propagules", "Naturel humide")) +
  labs(colour = "Grid cluster", fill ="Types d'occupation des sols", x =NULL, y = NULL)+
  annotation_north_arrow(location = "tr", which_north = "true", 
                         style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl", width_hint = 0.5,text_cex = 1) +
  main_theme+
  theme(line = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))



## Area per class cover on the overall landscape ----
as.data.frame(soil_occup_crop)%>%
  dplyr::select(lib, area)%>%
  group_by(lib)%>%
  summarise(sum_area = sum(area, na.rm = T))%>%
  mutate(perc_area = as.numeric(sum_area/sum(sum_area)*100))%>%
  arrange(desc(perc_area)) -> area_per_class

area_per_class

area_per_class%>%
  #filter(perc_area >2.5)%>%
  pull(lib) -> cover_names 

# Compute buffers sizes estimations with siland ----


(fmla <- as.formula(paste("log_plant_richness  ~1 +", paste(cover_names, collapse= "+"))))
if(file.exists("outputs/siland_mod_plants.Rdata")){
  load("outputs/siland_mod_plants.Rdata")
}else{
  mod1pois = siland(fmla, land = soil_occup_crop, init = c(100, 100, 100, 100, 100), data = metadata_grid, wd = 10, family = "gaussian")
  save(mod1pois, file ="outputs/siland_mod_plants.Rdata")
}


summary(mod1pois)


likres = siland.lik(mod1pois,land = soil_occup_crop,data = metadata_grid, varnames = cover_names)
likres+
  main_theme
siland.quantile(mod1pois)

plotsiland.sif(mod1pois)+
  main_theme

# Compute buffer and percentage around points ----


buffer.around.point = function(pt, geo_data, layer_name, buffer_size){
  pt = as.matrix(pt)
  if(nrow(pt) == 2){
    pt = t(pt)
  }
  pt_vect = terra::vect(pt, type="points", atts=NULL, crs=terra::crs(geo_data))
  buffer_vec <- terra::buffer(pt_vect, buffer_size)
  crop_suface = terra::crop(geo_data, buffer_vec)
  
  names(crop_suface)[names(crop_suface) == layer_name] = "focal_layer"
  
  sufaces_class = data.frame( class= crop_suface$focal_layer, surface = expanse(crop_suface, unit="m", transform=TRUE))
  sufaces_per_class = sufaces_class%>%
    group_by(class)%>%
    summarise(sum_area = sum(surface, na.rm = T))%>%
    mutate(perc_area = sum_area/sum(sum_area)*100,
           X = pt[1],
           Y = pt[2])%>%
    arrange(desc(perc_area))
  
  return(list(sufaces_per_class = sufaces_per_class, crop_suface = st_as_sf(crop_suface) ))
}

## Compute area percentage within a buffer ----
buffer_size = 300 # buffer radius
croped_surfaces = apply(data.frame(metadata_grid$X, metadata_grid$Y), 1, function(x)
  buffer.around.point(pt = x, geo_data = vect(soil_occup_crop), layer_name = "lib", buffer_size = buffer_size))

# croped_surfaces contain percentage area dataframes and also a sp objects of 
# each buffer (if needed for graphical representation of the landscape for exemple)

## Extract and group percentage area dataframes into one dataframe ----

area_per_buffer = croped_surfaces[[1]][[1]]
area_per_buffer$Grid_code =metadata_grid$Grid_code[1]
for (i in 2:length(croped_surfaces)){
  temp = croped_surfaces[[i]][[1]]
  temp$Grid_code = metadata_grid$Grid_code[i]
  area_per_buffer = rbind(area_per_buffer ,temp)
}

ggplot(area_per_buffer)+
  geom_point(aes(x=class, perc_area),cex = 2)+
  facet_wrap(~Grid_code)

area_per_buffer$buffer_size = buffer_size
area_per_buffer%>%
  select(class, perc_area, Grid_code,buffer_size)%>%
  pivot_wider(names_from = class, values_from = perc_area, values_fill = 0) -> area_per_buffer_wide 

#add columns only if they are not present in the df (avoid duplicated columns)
if(any(colnames( read.table("data/data_clean/Metadata_grid_CAM.txt", header = T))%in%c("wetland", "natural_landscape", "non_emitting", "cultivated", "artificial"))){
  cat("Metadata already contain landscapes covers values")
}else{
  read.table("data/data_clean/Metadata_grid_CAM.txt", header = T)%>%
    left_join(area_per_buffer_wide, by = "Grid_code")%>%
    write.table(file = "data/data_clean/Metadata_grid_CAM.txt")
}

