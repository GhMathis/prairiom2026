library(tidyverse)
library(phyloseq)
library(iNEXT)

#Courbes de raréfaction----

load(file = "outputs/ps_bacteria.Rdata")
ps_bacteria -> ps

ps %>% 
  otu_table() %>%
  as.data.frame() %>%
  rownames_to_column("ASV") %>%
  pivot_longer(-ASV, names_to = "grid_id", values_to = "reads") %>%
  mutate(grid_id = str_extract(grid_id, "..-CAM-.."),
         grid_id = str_replace_all(grid_id, "-", "_")) %>%
  group_by(ASV, grid_id) %>%
  summarise(reads = sum(reads)) -> otu_table_grid_long

otu_table_grid_long %>%
  filter(reads != 0) %>%
  ggplot() +
  geom_histogram(aes(reads), bins = 1000) +
  scale_y_log10()

otu_table_grid_long %>%
  filter(reads != 0) %>%
  ggplot() +
  facet_wrap(~grid_id) +
  geom_histogram(aes(reads), bins = 500) +
  scale_y_log10() +
  scale_x_log10()

otu_table_grid_long %>%
  pivot_wider(names_from = grid_id, values_from = reads) %>%
  column_to_rownames("ASV")-> otu_table_grid
str(otu_table_grid)

ps %>% sample_data() %>% head
ps %>% sample_data() %>% head
ps_grid_merge %>% tax_table() %>% str
otu_table_grid_long %>%
  group_by(grid_id   ) %>%
  filter(reads != 0) %>%
  summarise(reads = sum(reads),
            n_OTU = n()) %>%
  mutate(grid_id = fct_reorder(grid_id, reads)) %>%
  ggplot() +
  geom_point(aes(grid_id, reads, col = log10(n_OTU))) +
  scale_y_log10()

otu_table_grid_long %>%
  group_by(grid_id   ) %>%
  filter(reads != 0) %>%
  summarise(reads = sum(reads),
            n_OTU = n()) %>%
  mutate(grid_id = fct_reorder(grid_id, reads)) %>%
  ggplot() +
  geom_point(aes(n_OTU, reads)) 

# Coverage estimation

otu_table_grid_long %>%
  group_by(grid_id) %>% 
  mutate(sample_size = sum(reads)) %>% 
  ungroup() %>%
  filter(reads <= 2, reads!=0) %>%
  group_by(grid_id, reads) %>% 
  summarize(count = n(),
            sample_size = unique(sample_size)) %>% 
  ungroup %>%
  pivot_wider(names_from =reads , values_from = count) %>%
  mutate(coverage = 1-(`1`/sample_size) * ((sample_size-1)*`1`)/((sample_size-1)*`1` + 2*`2`),
         max_coverage_extrapolation = coverage+(exp(-2*`2`/`1`))*(1-coverage)) -> coverage_reads_df

coverage_reads_df 
coverage_reads_df %>%
  mutate(grid_id = reorder(grid_id, coverage)) %>%
  ggplot() +
  geom_point(aes(grid_id, coverage)) +
  geom_hline(aes(yintercept= min(max_coverage_extrapolation)), linetype = 2, linewidth = 1.5, color = "orange")
coverage_reads_df %>%
  mutate(grid_id = reorder(grid_id, coverage)) %>%
  ggplot() +
  geom_point(aes(grid_id, sample_size, colour = coverage)) +
  scale_y_log10()

coverage_reads_df %>%
  mutate(grid_id = fct_reorder(grid_id, sample_size)) %>%
  ggplot() +
  geom_point(aes(grid_id, sample_size, colour = coverage)) +
  geom_hline(aes(yintercept= min(sample_size)*2), linetype = 2, linewidth = 1.5, color = "orange")+
  scale_y_log10()


otu_table_grid %>% colSums %>% sort %>% min *2 %>% round-> max_extrapolation

out <- iNEXT(otu_table_grid, q = 0, datatype="abundance", size = seq(round(max_extrapolation/2), max_extrapolation, length.out = 10 ))

out$iNextEst$coverage_based %>% 
  filter(Order.q == 0) %>%
  
  ggplot()+
  geom_path(aes(x =  SC, y = qD, col = Method,group = Assemblage))+
  geom_point(aes(x =  SC, y = qD, col = Method,group = Assemblage))+
  geom_vline(xintercept = min(coverage_reads_df$max_coverage_extrapolation))+
  xlab("Coverage") +
  ylab("Number of ASV") +
  theme(
    legend.position = "bottom", 
    axis.text.x = element_text(angle=90,vjust=0.4,hjust=1.2)) +
  scale_y_log10()+
  scale_x_log10()

min(coverage_reads_df$max_coverage_extrapolation) -> coverage_min

richesse_standardisee <- estimateD(
  x = otu_table_grid,
  datatype = "abundance",
  base = "coverage",           # Standardisation par coverage
  level = coverage_min,        # Niveau de coverage cible
  conf = 0.95
)
