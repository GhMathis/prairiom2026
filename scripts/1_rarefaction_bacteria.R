library(tidyverse)
library(phyloseq)
library(iNEXT)

#Courbes de raréfaction----

load(file = "outputs/ps_bacteria.Rdata")
ps_bacteria -> ps
ps_grid_merge <- merge_samples(ps, "grid")

ps %>%
  otu_table() %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "sample_ID", values_to = "reads") %>%
  group_by(sample_ID) %>%
  summarise(reads = sum(reads)) %>%
  mutate(sample_ID = fct_reorder(sample_ID, reads)) %>%
  ggplot() +
  geom_point(aes(sample_ID, reads)) 

# Coverage estimation

ps %>%
  otu_table() %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "sample_ID", values_to = "reads") %>%
  group_by(sample_ID) %>% 
  mutate(sample_size = sum(reads)) %>% 
  ungroup() %>%
  filter(reads <= 2, reads!=0) %>%
  group_by(sample_ID, reads) %>% 
  summarize(count = n(),
            sample_size = unique(sample_size)) %>% 
  ungroup %>%
  pivot_wider(names_from =reads , values_from = count) %>%
  mutate(coverage = 1-(`1`/sample_size) * ((sample_size-1)*`1`)/((sample_size-1)*`1` + 2*`2`),
         max_coverage_extrapolation = coverage+(exp(-2*`2`/`1`))*(1-coverage)) -> coverage_reads_df

coverage_reads_df %>%
  mutate(sample_ID = reorder(sample_ID, coverage)) %>%
  ggplot() +
  geom_point(aes(sample_ID, coverage))
coverage_reads_df %>%
  mutate(sample_ID = reorder(sample_ID, coverage)) %>%
  ggplot() +
  geom_point(aes(sample_ID, sample_size, colour = coverage)) +
  scale_y_log10()

coverage_reads_df %>%
  mutate(sample_ID = fct_reorder(sample_ID, sample_size)) %>%
  ggplot() +
  geom_point(aes(sample_ID, sample_size, colour = coverage)) +
  scale_y_log10()

coverage_reads_df %>%
  ggplot() +
  geom_point(aes(sample_size,coverage )) +
  geom_hline(aes(yintercept= min(max_coverage_extrapolation,na.rm = T)), linetype = 2, linewidth = 1.5, color = "orange")+
  geom_text(aes(x = 1,y = 1,label = round(min(max_coverage_extrapolation,na.rm = T),2))) +
  scale_x_log10()
  
coverage_reads_df %>%
  filter(coverage >0.60)%>%
    ggplot() +
    geom_point(aes(sample_size,coverage )) +
    geom_hline(aes(yintercept= min(max_coverage_extrapolation,na.rm = T)), linetype = 2, linewidth = 1.5, color = "orange") +
    geom_text(aes(x = 1,y = 1,label = round(min(max_coverage_extrapolation,na.rm = T),2))) +
    scale_x_log10()

ps %>% 
  prune_samples(sample_sums(.) >= 1000, .) %>% 
  otu_table(.) %>% 
  as.data.frame(.) -> tab
tab %>% colSums %>% sort

out <- iNEXT(tab, q=c(0, 1, 2), datatype="abundance", size = c(1000, 2000, 10000, 20000, 50000))
save(out, file="outputs/diversities.Rdata")

load("outputs/diversities.Rdata")

out$iNextEst$size_based %>% 
  #filter(Assemblage %in% my_mtd_env$sampli_ID[1:72]) %>%
  filter(Order.q == 0) %>% 
  ggplot(., aes(x = m, y = qD, col = Assemblage))+
  geom_line()+
  xlab("Number of reads") +
  ylab("Number of ASV") +
  theme_minimal()+
  theme(
    legend.position = "none", 
    axis.text.x = element_text(angle=90,vjust=0.4,hjust=1.2))

out$iNextEst$size_based %>% 
  filter(Order.q == 1) %>%
  filter(m<=50000) %>% 
  ggplot(., aes(x = m, y = qD, col = Assemblage))+
  geom_line()+
  xlab("Number of reads") +
  ylab("Number of ASV") +
  theme_minimal()+
  theme(
    legend.position = "none", 
    axis.text.x = element_text(angle=90,vjust=0.4,hjust=1.2))

out$iNextEst$coverage_based %>% 
  filter(Order.q == 1) %>%
  filter(m<=50000) %>% 
  ggplot()+
  geom_path(aes(x =  SC, y = qD, col = Method,group = Assemblage))+
  xlab("Coverage") +
  ylab("Number of ASV") +
  theme_minimal()+
  theme(
    legend.position = "none", 
    axis.text.x = element_text(angle=90,vjust=0.4,hjust=1.2)) +
  scale_y_log10()+
  scale_x_log10()

out$iNextEst$size_based %>% 
  filter(Order.q == 1) %>%
  filter(m<=50000) %>% 
  ggplot(., aes(x = m, y = qD, col = Assemblage))+
  geom_line()+
  xlab("Number of reads") +
  ylab("Number of ASV") +
  theme_minimal()+
  theme(
    legend.position = "none", 
    axis.text.x = element_text(angle=90,vjust=0.4,hjust=1.2))

out$iNextEst$size_based %>% 
  filter(Order.q == 1) %>%
  left_join(., my_mtd_env, by = c("Assemblage" = "sample_ID")) %>% 
  ggplot(., aes(x = m, y = qD, col = quadrat))+
  geom_line()+
  facet_wrap(~grid, scales = "free_y")+
  xlab("Number of reads") +
  ylab("Number of ASV") +
  theme_minimal()+
  theme(
    legend.position = "none", 
    axis.text.x = element_text(angle=90,vjust=0.4,hjust=1.2))

# out_small <- iNEXT(tab, q=c(0, 1), datatype="abundance", size = seq(1000, 2000, by=200), nboot = 25)
# 
# save(out_small, file="outputs/diversities_small.Rdata")

load("outputs/diversities_small.Rdata")

out_small$iNextEst$coverage_based %>% 
  filter(Order.q == 0) %>%

  ggplot()+
  geom_path(aes(x =  SC, y = qD, col = Method,group = Assemblage))+
  xlab("Coverage") +
  ylab("Number of ASV") +
  theme_minimal()+
  theme(
    legend.position = "bottom", 
    axis.text.x = element_text(angle=90,vjust=0.4,hjust=1.2)) +
  scale_y_log10()+
  scale_x_log10()

out_small$iNextEst$coverage_based %>%
  filter(Method == "Observed") %>%
  pull(SC) %>%
  min() -> coverage_min

richesse_standardisee <- estimateD(
  x = tab,
  datatype = "abundance",
  base = "coverage",           # Standardisation par coverage
  level = coverage_min,        # Niveau de coverage cible
  conf = 0.95
)
#*******************----
#BARPLOT - Composition bactérienne par échantillon----
#Une façon de faire basée sur les abondances relatives

ps %>% 
  prune_samples(sample_sums(.) >= 1000, .) %>% 
  otu_table(.) %>% 
  as.data.frame(.) %>% 
  filter(rowSums(.) > 0) %>% 
  mutate(.before = 1, OTU = rownames(.)) %>%
  pivot_longer(cols = c(2:ncol(.)), names_to = "sample_ID", values_to = "reads") %>% 
  filter(reads > 0) -> my_df

my_tax %>%
  filter(rownames(.) %in% my_df$OTU) %>% 
  mutate(OTU = rownames(.))-> my_tax_df

my_df %>% 
  left_join(., my_tax_df, by = "OTU") -> my_df 

my_df %>% 
  count(family,sample_ID, wt = reads, name = "total") %>% 
  group_by(sample_ID) %>% 
  mutate(frq = total/sum(total)) %>% 
  ungroup() -> df_with_NA

# Définition de la palette de couleurs
library(paletteer)
ncol = 1
colors <- paletteer_d("ggsci::default_igv")[seq(ncol, ncol + 36, 1)]

df_with_NA %>%
  group_by(family) %>% 
  mutate(family = ifelse(all(frq < 0.05), "Autres", family)) %>% 
  ungroup() %>%
  left_join(., my_mtd_env, by = "sample_ID" )%>% 
  ggplot(., aes(x=sample_ID, y = frq, fill = family)) +
  geom_bar(stat = "identity")+
  theme_minimal()+
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.justification = "top",
        legend.box = "vertical",
        legend.key.size = unit(0.5, "cm"),
        axis.text.x = element_text(angle=90, vjust=0.4, hjust=1.2, size = 0),
        legend.text = element_text(size = 11))+
  guides(fill = guide_legend(ncol = 1))+
  labs(x="Echantillons environnementaux", y = "Abondances relatives", fill="Familles")+
  facet_wrap(~grid, scales="free_x", ncol = 11)+
  scale_fill_manual(values = colors)

# Abondance des familles par échantillon
my_df  %>% 
  count(family,sample_ID, wt = reads, name = "total") %>% 
  group_by(sample_ID) %>% 
  mutate(frq = total/sum(total)) %>% 
  ungroup() %>% 
  group_by(family) %>%
  mutate(family = ifelse(all(frq < 0.05), "Autres", family)) %>% 
  ungroup() %>% 
  group_by(family, sample_ID) %>% 
  summarise(frq = sum(frq)) %>% 
  ungroup() %>% 
  mutate(frq = frq*100) %>% 
  # filter(family == "Enterobacteriaceae") %>%
  pivot_wider(
    id_cols = 1,
    names_from = sample_ID,
    values_from = frq,
    values_fill = NA
  ) %>%
  view()


#*******************----
#NMDS----
#Avec matrice présence/absence
#Env/sex
ps_rare = rarefy_even_depth(ps, sample.size = 2000)
ps_bin = transform_sample_counts(ps_rare, function(x) 1*(x>0))

# sample_data(ps_bin) %>% View()
# 
# otu_table(ps) %>% 
#   as.data.frame() %>% 
#   summarise(across(everything(), sum)) %>% 
#   select_if(~ sum(.x) < 70000) %>% View()

library(paletteer)
ncol = 1
colors <- paletteer_d("ggsci::default_igv")[seq(ncol, ncol + 42, 1)]
library(ggrepel)

jaccard <- ordinate(ps_bin, "NMDS", "jaccard") #stress trop haut

plot_ordination(ps_bin,
                jaccard,
                type = "samples",
                color = "grid") +
  theme_bw() +
  geom_point(size = 1) +
  theme(legend.position = "none",
        legend.text = element_text(size = 10)) +
  theme(legend.text = element_text(face = "italic")) + 
  scale_color_manual(name = "Grid", values = colors) +
  geom_text_repel(aes(label = grid, fontface = 'italic'), size = 2,max.overlaps = 50)

ggsave("results/NMDS.svg", width=15, height = 10, units = "cm")

#NMDS abondance
bray <- ordinate(ps_rare, "NMDS", "bray") #stress élevé

plot_ordination(ps_rare,
                bray,
                type = "samples",
                color = "grid") +
  theme_bw() +
  geom_point(size = 2)+
  theme(legend.position = "none") +
  theme(legend.text = element_text(face = "italic"))+
  scale_color_manual(name = "Grid", values = colors) +
  geom_text_repel(aes(label = grid, fontface = 'italic'), size = 2,max.overlaps = 30)

#*******************----
#Analyse de diversité----
load("results/diversities.Rdata")

div_df<-out$AsyEst[out$AsyEst$Diversity=="Shannon diversity",]

div_df<-left_join(div_df, my_mtd_env, by = c("Assemblage" = "sample_ID"))

#div alpha par site 
div_df %>% 
  group_by(grid) %>% 
  mutate(.before = 4, "mean_ASV_richness" = mean(Estimator)) %>% 
  distinct(mean_ASV_richness, .keep_all = TRUE) -> div_df

#Box plot env
div_df %>% 
  arrange(mean_ASV_richness) %>% 
  pull(grid) -> sample_order

ggplot(div_df, aes(x = grid, y = Estimator, fill = grid, col = grid)) +
  geom_jitter(alpha = 0.5, width = 0.25) + 
  geom_point(aes(x = grid, y = mean_ASV_richness), size = 3, shape = 1) + 
  theme_minimal() +
  theme(legend.position = "none") +
  theme(
    axis.text.x = element_text(size = 4, angle = 90),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(face = "italic")
  ) +
  scale_x_discrete(name = "Grilles", limits = sample_order) +
  scale_y_continuous(name = "Diversité alpla") +
  scale_fill_manual(name = "Grids", values = colors) +
  scale_color_manual(name = "Grids", values = colors)


#*******************----
#CLUSTERING FROM READ COUNTS----
library("sbm")

#il y a 19000 ASV même après raréfaction, je dois agglomérer avant de faire des calculs de SBM
hist(log10(taxa_sums(ps_rare)+1))

#set.seed(1)
#ps_rare = rarefy_even_depth(ps, sample.size = 24000)

ps_family = tax_glom(ps, taxrank="family")
summary(taxa_sums(ps_family))
hist(log10(taxa_sums(ps_family)+1))

summary(sample_sums(ps_family))
hist(log10(sample_sums(ps_family)+1))
tax_table(ps_family) %>% view()

ps_rare = rarefy_even_depth(ps_family, sample.size = 2000)


#Représentation graphique de la matrice
Log_FB <- as.matrix(log10(t(otu_table(ps_rare))+1))
plotMyMatrix(Log_FB, dimLabels = c(row = "Samples", col = "ASV"))

Log_BF <- log10(otu_table(ps_rare)+1)
plotMyMatrix(Log_BF, dimLabels = c(row = "ASV", col = "Samples"), plotOptions = list(colNames = T))+
  theme(axis.text.x = element_text(size = 4))

#Estimation du LBM
LBM_log_reads_gaussian <-
  as.matrix(Log_FB) %>%
  estimateBipartiteSBM(
    model = 'gaussian')

save(LBM_log_reads_gaussian, ps_rare, file = "results/LBM_log_reads_gaussian2000.Rdata")
# ggsave("LBM_log_reads_gaussian.svg")

load("results/LBM_log_reads_gaussian.Rdata")

memb_spl_obs <- LBM_log_reads_gaussian$memberships$row
memb_bac_obs <- LBM_log_reads_gaussian$memberships$col

LBM_log_reads_gaussian$storedModels %>% arrange(ICL)

# ggsave("results/weigthed_incid_mat.svg")

# LBM_log_reads_poisson <-
#   as.matrix(Log_FB) %>%
#   estimateBipartiteSBM(
#     model = 'poisson')

# save(LBM_log_reads_poisson, ps_rare, file = "results_fruit/LBM_log_reads_poisson24000.Rdata")
# load("results_fruit/LBM_log_reads_poisson24000.Rdata")
# 
# memb_spl_POISSON <- LBM_log_reads_poisson$memberships$row
# memb_bac_POISSON <- LBM_log_reads_poisson$memberships$col

remotes::install_github("GrossSBM/sbm")
library(sbm)
LBM_log_reads_ZIgaussian <-
  as.matrix(Log_FB) %>%
  estimateBipartiteSBM(
    model = 'ZIgaussian')

#représentation du clustering obtenu----
plotOptions <- list(rowNames = T,colNames = TRUE,legend = TRUE,legend.position = 'right',line.width = 1)
plotOptions$nodeNames <-T

plot(LBM_log_reads_gaussian, dimLabels = list(row = 'fungis', col= 'tree'), plotOptions= plotOptions)+
  theme(
    legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 8, face="italic"),
    axis.text.y = element_text(size = 6))

#
df=as.data.frame(t(Log_FB))
df$taxon=rownames(df)
df_long<-pivot_longer(df, cols=names(df)[-ncol(df)], names_to =
                        "sample_ID", values_to = "reads")
#df_long<-left_join(df_long,mtd, by="sample_ID")

bac_list<-colnames(Log_FB)
sample_list<-rownames(Log_FB)

LBM_log_reads_gaussian$nbBlocks #how many clusters
Qspl<-LBM_log_reads_gaussian$nbBlocks["Samples"][[1]]
Qasv<-LBM_log_reads_gaussian$nbBlocks["ASV"][[1]]

LBM_log_reads_gaussian$memberships #cluster composition
lapply(1:Qspl,function(q){sample_list[LBM_log_reads_gaussian$memberships$Samples
                                      == q]})
lapply(1:Qasv,function(q){bac_list[LBM_log_reads_gaussian$memberships$ASV == q]})

df_tax_cluster=data.frame(taxon=bac_list,tax_cluster=LBM_log_reads_gaussian$memberships$ASV)
df_long<-left_join(df_long, df_tax_cluster, by="taxon")

df_long %>% arrange(tax_cluster) -> df_long
correct_taxon_order <- unique(df_long$taxon)

df_spl_cluster=data.frame(sample_ID=sample_list,spl_cluster=LBM_log_reads_gaussian$memberships$Samples)
df_long<-left_join(df_long, df_spl_cluster, by="sample_ID")

df_long %>% arrange(spl_cluster) -> df_long
correct_sample_order <- unique(df_long$sample_ID)

ggplot(df_long,aes(x=taxon, y=sample_ID, fill=reads))+
  geom_tile()+
  theme_minimal()+
  theme(legend.position="none",
        axis.text.y = element_text(vjust = 0, hjust = 1, size = 6),
        axis.text.x = element_text(angle = 90, size = 6))+
  scale_fill_gradient(low="white", high="Black")+
  scale_x_discrete(limits = correct_taxon_order, label =
                     my_tax$class[match(my_tax$OTU, colnames(FB_bin))])+
  scale_y_discrete(limits = correct_sample_order)+
  geom_hline(yintercept=0.5+cumsum(table(LBM_log_reads_gaussian$memberships$Samples))[-length(table(LBM_log_reads_gaussian$memberships$Samples))],
             color = "red")+
  geom_vline(xintercept=0.5+cumsum(table(LBM_log_reads_gaussian$memberships$ASV))[-length(table(LBM_log_reads_gaussian$memberships$ASV))],
             color = "red")+
  ylab("Echantillons")+
  xlab("Bactéries")

# #--------- Null model 
# NPERM=10000
# NPERM=10
# 
# randomize.WEDD <- function(bipmat) {
#   s <- sum(bipmat)
#   m1 <- apply(bipmat, 1, sum)
#   m2 <- apply(bipmat, 2, sum)
#   apply(m1 %o% m2 / s, c(1, 2), function(x)
#     rnorm(1, mean = x, sd = sd(m1 %o% m2 / s)))
# }
# 
# many.random.WEDD <- function(bipmat, nperm = NPERM) {
#   res <- array(0, dim = c(dim(bipmat)[1], dim(bipmat)[2], nperm))
#   bidule <- lapply(1:nperm, function(x)
#     randomize.WEDD(bipmat))
#   for (i in 1:nperm) {
#     res[, , i] <- bidule[[i]]
#   }
#   res
# }
# 
# ###production of NPERM random matrices
# random.mat <- many.random.WEDD(Log_FB)
# 
# # save(random.mat, ps_rare, file = "results_fruit/random.mat_10000.Rdata")
# load("results_fruit/random.mat_10000.Rdata")
# 
# dim(random.mat)
# random.mat
# 
# ###LBM classification on these NPERM random matrices
# clust_sim = sapply(1:NPERM, function(i) {
#   Log_FB_sim <- random.mat[, , i]
#   LBM_sim <-
#     as.matrix(Log_FB_sim) %>%
#     estimateBipartiteSBM(model = 'gaussian')
#   
#   memb_spl_sim <- LBM_sim$memberships$row
# })
# 
# save(clust_sim, file = "results_fruit/sim_weighted_NPERM10000_log.Rdata")
# load("results_fruit/sim_weighted_NPERM10.Rdata")
# 

set.seed(1)
ps_rare = rarefy_even_depth(ps, sample.size = 24000)
hist(log10(taxa_sums(ps_rare)+1))

otu_table(ps_rare) %>% View()

#Représentation graphique de la matrice
Log_FB <- log10(t(otu_table(ps_rare)+1))
Log_FB

# ## --------------------------------- variables
library(parallel)

set.seed(1)
nbCores <- 1
n_iterations <- 10000

sum_reads <- sum(Log_FB)
sum_reads_per_OTU <- apply(Log_FB, 1, sum)
sum_reads_per_samples <- apply(Log_FB, 2, sum)
nrow_log_FB <- nrow(Log_FB)

sum_reads_per_OTU %o% sum_reads_per_samples / sum_reads %>% str()

## --------- Null model ------

randomize.WEDD <- function(OTU_matrix_t) {
  apply(sum_reads_per_OTU %o% sum_reads_per_samples / sum_reads, 
        c(1, 2), function(x)
          rnorm(1, mean = x, 
                sd = sd(sum_reads_per_OTU %o% sum_reads_per_samples / sum_reads)))
}

simulate_membership <- function(OTU_matrix_t) {
  sbm::estimateBipartiteSBM(netMat = OTU_matrix_t,
                            model = "gaussian",
                            estimOptions = list(
                              verbosity = 0,
                              plot = FALSE,
                              nbCores = nbCores)) -> LBM_sim
  
  LBM_sim$memberships$row
}

randomize_and_simulate <- function(n) {
  #argument n non utilisé mais obligatoire pour mclapply
  Log_FB %>%
    randomize.WEDD(.) %>%
    simulate_membership(.)
}

nrow_log_FB <- nrow(FB_bin)

system.time(
  parallel::mclapply(1:n_iterations,
                     randomize_and_simulate,
                     mc.cores = (parallel::detectCores()/2)-1) %>% 
    unlist(.) %>%
    matrix(., nrow = nrow_log_FB , byrow = TRUE) -> clust
)
rm(sum_reads, sum_reads_per_OTU, sum_reads_per_samples)
