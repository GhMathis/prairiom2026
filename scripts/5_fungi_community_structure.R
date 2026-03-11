library(tidyverse)
library(sbm)
library(RColorBrewer)
library(alluvial)
library(GGally)

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



# GLM richness ----
read.table("data/data_clean/Metadata_grid_CAM.txt", header = T) -> metadata_grid
read.table("data/data_clean/OTU_fungi_for_sbm.txt", check.names = FALSE) -> fungi_otu_sbm

metadata_grid %>% 
  dplyr::select(log_virus_richness, log_plant_richness, log_plant_H,
                log_plant_equit, log_biomass, Vegetation,
                fungi_richness)%>% 
  ggpairs()


# SBM ----


if(file.exists("outputs/sbm_fungi_grid.Rdata")){
  load("outputs/sbm_fungi_grid.Rdata")
  
}else{
  fungi_otu_sbm%>%
    log1p()%>%
    t()%>%
    plotMyMatrix(dimLabels = c("Grids", "Fungi"), plotOptions = list(rowNames = T,colNames = F))+
    theme(element_text(angle = 90, vjust = 1, hjust = 1, size = 0.1, face="italic"),
          strip.text.x = element_text(colour = "gray", size=15, face ="italic"),
          strip.text.y  = element_text(colour = "gray", size=15, face ="italic"),
          strip.background = element_rect(fill="white"),
          panel.border = element_rect(color = "black",  fill = NA))
  
  
  fungi_otu_sbm%>%
    log1p()%>%
    estimateBipartiteSBM(
      model = 'poisson', 
      dimLabels = c(row = "Fungi", col = "Grid"),
      estimOptions = list(nbCores = 8,plot = F)) -> sbm_fungi_grid
  save(sbm_fungi_grid, file = "outputs/sbm_fungi_grid.Rdata")
}


fungi_id = sbm_fungi_grid$memberships[1]%>%
  as.data.frame()%>%
  rownames_to_column("id")%>%
  mutate(id = as.numeric(id))%>%
  arrange(Fungi)
grid_id = sbm_fungi_grid$memberships[2]%>%
  as.data.frame()%>%
  rownames_to_column("id")%>%
  mutate(id = as.numeric(id))%>%
  arrange(Grid)

fungi_otu_sbm[fungi_id$id,grid_id$id]%>%
  as.matrix()%>%
  log1p()%>%
  t()%>%
  plotMyMatrix(dimLabels = c("Grid", "Fungi"), plotOptions = list(rowNames = T,colNames = F))+
  geom_hline(yintercept = abs(which(!duplicated(grid_id$Grid))[-1]-nrow(grid_id)-1)+0.5, col ="red")+
  geom_vline(xintercept = which(!duplicated(fungi_id$Fungi))[-1]-0.5, col = "red", linetype = 5)+
  theme(element_text(angle = 90, vjust = 1, hjust = 1, size = 0.1, face="italic"),
        strip.text.x = element_text(colour = "gray", size=15, face ="italic"),
        strip.text.y  = element_text(colour = "gray", size=15, face ="italic"),
        strip.background = element_rect(fill="white"),
        panel.border = element_rect(color = "black",  fill = NA))


str(metadata_grid)
if (!any(colnames(metadata_grid) %in% c("sbm_fungi"))){ #add the column only if it is not present in the df (avoid duplicated column)
  metadata_grid%>%
    left_join(data.frame(sbm_fungi = sbm_fungi_grid$memberships$Grid,
                       Grid_code = str_replace_all(colnames(fungi_otu_sbm),
                                              "-","_")), by = "Grid_code") -> metadata_grid
  write.table(metadata_grid, file = "data/data_clean/Metadata_grid_CAM.txt")
  }

## aluvial plot ----
B <-
  as.data.frame(
    table(
      
      metadata_grid$Grid_code,
      metadata_grid$clust_smb_grid_plant ,
      metadata_grid$clust_land_use,
      metadata_grid$clust_landscp,
      metadata_grid$sbm_fungi,
      metadata_grid$sbm_soil)
  )


colnames(B) =  c( "Grid_code","sbm_plant", "land_use",  "landscp", "sbm_fungi","smb_soil","Freq")


w   <- which(B$Freq != 0)
B <- B[w, ]

data.frame(sbm_plant = unique(B$sbm_plant),
           col_plant = brewer.pal(7, "Set1")) -> col_cluster_plant
data.frame(sbm_fungi = unique(B$sbm_fungi),
           col_fungi = brewer.pal(8, "Set3")) -> col_cluster_fungi

B%>%
  left_join(col_cluster_fungi, by = "sbm_fungi")%>%
  left_join(col_cluster_plant, by = "sbm_plant")-> B

B%>%
  mutate(landscp =  case_when(landscp == 1 ~"(1) Cultivés dominant+\n Anthropisés",
                              landscp == 2 ~"(2) Naturels non Humide\n dominant",
                              landscp == 3 ~"(3) Naturels dominants+\n Non Emmeteurs",
                              landscp == 4 ~"(4) Naturels Humide\n dominant")) -> B
B%>%
  mutate(land_use =  case_when(land_use == 1 ~"(1) Pâture",
                               land_use == 2 ~"(2) Rien",
                               land_use == 3 ~"(3) Fauche",
                               land_use == 4 ~"(4) Pâture +\n  Fauche")) -> B
B$sbm_plant <- factor(B$sbm_plant, levels = c( "7", "1",  "5", "2", "4", "3", "6"))
B$sbm_fungi <- factor(B$sbm_fungi, levels = c( "1", "2",  "4", "5", "6", "7", "8","3"))
colnames(B)[ c(3,2,5)] = c("Agricole","Végétale", "Fongique")

alluvial(B[, c(3,5,2)],
         freq = B$Freq,
         col = B$col_fungi,
         alpha = 0.7,
         cex = 1.5,
         cex.axis = 1.5,
         xw = 0.2)

alluvial(B[, c(5,3)],
         freq = B$Freq,
         col = B$col_fungi,
         alpha = 0.7,
         xw = 0.2)

alluvial(B[, c(5,4)],
         freq = B$Freq,
         col = B$col_fungi,
         alpha = 0.7,
         xw = 0.2)

alluvial(B[, c(5,6)],
         freq = B$Freq,
         col = B$col_fungi,
         alpha = 0.7,
         xw = 0.2)

alluvial(B[, c(3,5,2)],
         freq = B$Freq,
         col = B$col_fungi,
         alpha = 0.7,
         xw = 0.2)

alluvial(B[, c(2,3)],
         freq = B$Freq,
         col = B$col_fungi,
         alpha = 0.7,
         xw = 0.2)

## NMI ----
library(aricode)
nmi_land_use = NMI(c1 =metadata_grid$sbm_fungi ,c2 = metadata_grid$clust_land_use)
nmi_landscape = NMI(c1 =metadata_grid$sbm_fungi ,c2 = metadata_grid$clust_landscp)
nmi_soil = NMI(c1 =metadata_grid$sbm_fungi ,c2 = metadata_grid$sbm_soil)
nmi_plant = NMI(c1 =metadata_grid$sbm_fungi ,c2 = metadata_grid$clust_smb_grid_plant)
nmi_land_use
nmi_landscape
nmi_soil
nmi_plant

## Randomization ----

sapply(1:10000,function(x) sample(metadata_grid$sbm_fungi)) -> rowperm_nmi
data.frame(row_land_use  = apply(rowperm_nmi,2 ,function(x) NMI(c1 = x ,c2 = metadata_grid$clust_land_use)) ,
           row_soil = apply(rowperm_nmi,2, function(x) NMI(c1 = x ,c2 = metadata_grid$sbm_soil)) ,
           row_landscape = apply(rowperm_nmi, 2, function(x) NMI(c1 = x ,c2 = metadata_grid$clust_landscp)),
           row_plant = apply(rowperm_nmi, 2, function(x) NMI(c1 = x ,c2 = metadata_grid$clust_smb_grid_plant))
)-> nmi_rand_metadata_grid

nmi_rand_metadata_grid%>%
  pivot_longer(everything(), names_to = "rand_type", values_to = "nmi_rand_values")%>%
  group_by(rand_type)%>%
  left_join(data.frame(nmi_values = c(nmi_land_use, nmi_landscape, nmi_soil, nmi_plant),
                       rand_type = c("row_land_use", "row_landscape",  "row_soil", "row_plant")),
            by = "rand_type")-> nmi_rand_metadata_grid

ggplot(nmi_rand_metadata_grid)+
  facet_wrap(~rand_type)+#,
  #labeller = as_labeller(setNames(plot_names, sort(unique(nmi_rand_metadata_grid$rand_type)))))+
  geom_histogram(aes(nmi_rand_values),bins = 100)+
  geom_vline(aes(xintercept = nmi_values),col ="red", linetype =2,linewidth =1.5)+
  labs(x ="Values of NMI", y= "Count")+
  main_theme+
  theme_classic()

data.frame(nmi_values = c(nmi_land_use,nmi_plant, nmi_soil, nmi_landscape
                          ),
           nmi_names = factor(c("Usages agri.","Communauté\n végétale", "Sol" ,"Paysage"),
                                levels= c("Communauté\n végétale","Usages agri.", "Sol" ,"Paysage")))%>%
  ggplot()+
    geom_point(aes(nmi_names, nmi_values), cex = 5)+
  labs(x = 'Classification', "Valeurs de congruence (NMI)")+
  main_theme
  

## CCA ----
library(vegan)
library(FactoMineR)
Y = metadata_grid$sbm_fungi%>%as.factor()
X = metadata_grid%>%dplyr::select(clust_landscp,clust_land_use, sbm_soil, clust_smb_grid_plant)%>%
  mutate_all(as.factor)
full_cca_virus = cca(tab.disjonctif(Y)~sbm_soil +clust_land_use+clust_smb_grid_plant, data = X)
full_cca_virus = cca(tab.disjonctif(Y)~sbm_soil+
                       clust_land_use+
                       clust_smb_grid_plant+
                       clust_landscp, data = X)

vif.cca(full_cca_virus) # to much colinearity

mod_varpart_champi_chi3 = varpart(tab.disjonctif(as.factor(Y)),
                                 ~clust_smb_grid_plant,
                                 ~clust_land_use, 
                                 ~clust_landscp, data=  X, chisquare=T,
                                 permutations = how(nperm = 10000))
plot(mod_varpart_champi_chi3,Xnames = c(list("Végétale"), list("Usage agri."), list("Paysage")),
     bg = c( "#CB6CE6","#CB6CE6","#C0344A"), lwd = 2, cex = 1.5, alpha = 170)

##
metadata_grid%>%
  mutate(sbm_fungi = as.factor(sbm_fungi_grid$memberships[2]$Grid))%>%
  group_by(sbm_fungi)%>%
  summarise(median(fungi_richness),
            n(),max(fungi_richness), min(fungi_richness))%>%t()%>%knitr::kable()

fungi_id = sbm_fungi_grid$memberships[1]%>%
  as.data.frame()%>%
  rownames_to_column("id")%>%
  mutate(id = as.numeric(id))%>%
  arrange(Fungi)
grid_id = sbm_fungi_grid$memberships[2]%>%
  as.data.frame()%>%
  rownames_to_column("id")%>%
  mutate(id = as.numeric(id))%>%
  arrange(Grid)
str(fungi_otu_sbm)
reshape2::melt(as.matrix(fungi_otu_sbm))%>%
  dplyr::rename(Fungi = "Var1", Grid_code = "Var2")%>%
  full_join(metadata_grid, by = "Grid_code")-> fungi_profil

fungi_profil%>%
  dplyr::select(Fungi, value,sbm_fungi )%>%
  filter(value>0)%>%
  
  group_by(sbm_fungi, Fungi)%>%
  summarise(n = sum(value))%>%
  group_by(sbm_fungi)%>%
  add_tally(n)%>%
  ungroup()%>%
  group_by(sbm_fungi,Fungi)%>%
  arrange(desc(n))%>%
  group_by(sbm_fungi)%>%
  slice_head(n= 10)%>%
  mutate(prop = n/nn) -> fungisp_metadata_grid

fungisp_metadata_grid$Fungi -> best_fungi_cluster
fungi_id$fungi_names = rownames(fungi_otu_sbm[fungi_id$id, grid_id$id])
fungi_id2 = fungi_id%>%
  filter(fungi_names %in% best_fungi_cluster)
abund_fungi_grid_fomated = as.matrix(fungi_otu_sbm)
colnames(abund_fungi_grid_fomated) = as.character(1:41)
rownames(abund_fungi_grid_fomated) = str_replace_all(rownames(abund_fungi_grid_fomated),"[_]", " ")


abund_fungi_grid_fomated[fungi_id2$id, grid_id$id]%>%
  log1p()%>%
  plotMyMatrix(dimLabels = c("Virus", "Grilles"), plotOptions = list(rowNames = T,colNames = T))+
  geom_hline(yintercept = abs(which(!duplicated(fungi_id2$Fungi))[-1]-nrow(fungi_id2)-1)+0.5, col ="red", linetype = 5)+
  geom_vline(xintercept = which(!duplicated(grid_id$Grid))[-1]-0.5, col = "red")+
  
  theme(element_text(angle = 45, vjust = 1, hjust = 1, size = 0.1, face="italic"),
        strip.text.x = element_text(colour = "gray", size=20, face ="italic"),
        strip.text.y  = element_text(colour = "gray", size=20, face ="italic"),
        axis.text.x =element_text( size = 8, face="italic"),
        axis.text.y = element_text(size = 10,face="italic"),
        strip.background = element_rect(fill="white"),
        panel.border = element_rect(color = "black",  fill = NA))

fungisp_metadata_grid %>%
  mutate(Fungi =  str_replace_all(Fungi,"[_]", " "),
         prop = round(prop,3)*100,
         Fungi =paste(Fungi," (", prop, "%)"))%>%
  group_by(sbm_fungi)%>%
  mutate(id_row = 1:10)%>%
  ungroup()%>%
  dplyr::select(sbm_fungi, Fungi, id_row)%>%
  pivot_wider(values_from = Fungi, names_from =id_row )%>%
  as.data.frame()%>%
  t()%>%
  as_tibble()%>%
  tail(10)%>%
  knitr::kable()
