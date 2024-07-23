# EDA proteinortho
# cds show good performance for three transcriptomes
dir <- "~/Documents/GitHub/Hybryd_Haliotis/Results/05.Meta_assembly/Protein_orth/"

f <- list.files(path = dir, pattern = "tsv_numeric.txt", full.names = T)

library(tidyverse)

df <- read_tsv(f[1])


recode_to <- structure(c("Hybryd (RF)", "Red (RR)", "Green (FF)"), 
  names = c("RF_cross.cds", "RR_cross.cds", "FF_cross.cds"))

UPSETDF <- df %>%
  pivot_longer(-FAM_ID) %>%
  dplyr::mutate(name = dplyr::recode_factor(name, !!!recode_to)) %>%
  filter(value >0) %>%
  group_by(FAM_ID) %>%
  summarise(across(name, .fns = list), n = n()) %>% arrange(desc(n))


library(ggupset)

UPSETDF %>%
  ggplot(aes(x = name)) +
  geom_bar(position = position_dodge(width = 1), 
    color = "black", linewidth = 0, fill = "gray80") +
  # geom_segment(aes(xend = CONTRAST_DE, yend = n, y = Inf, color = SIGN), linewidth = 1) +
  geom_text(stat='count', aes(label = after_stat(count)), 
    position = position_dodge(width = 1), vjust = -0.2, family = "GillSans", size = 5) +
  ggupset::scale_x_upset(order_by = "degree", reverse = F) +
  theme_bw(base_family = "GillSans", base_size = 18) +
  ggupset::theme_combmatrix(combmatrix.panel.point.color.fill = "black",
    combmatrix.panel.line.size = 0, 
    combmatrix.panel.point.size = 5,
    base_family = "GillSans") +
  labs(x = '', y = 'Number of orthologs')

m <- UPSETDF %>% filter(n == 3) %>% left_join(df, by = "FAM_ID") %>% select_if(is.double)
# m <- df %>% select_if(is.double)

sample_cor <- cor(m, method='spearman', use='pairwise.complete.obs')

h <- heatmap(sample_cor, col = cm.colors(12), keep.dendro = T)

hc_samples <- as.hclust(h$Colv)

plot(hc_samples)

hc_order <- hc_samples$labels[h$colInd]

sample_cor %>% 
  as_tibble(rownames = 'Strain') %>%
  pivot_longer(cols = colnames(sample_cor), values_to = 'cor') %>%
  # left_join(.colData) %>% 
  mutate(Strain = factor(Strain, levels = rev(hc_order))) -> sample_cor_long

library(ggh4x)

sample_cor_long %>%
  ggplot(aes(x = Strain, y = name, fill = cor)) +  
  geom_tile(linewidth = 0.2) +
  ggsci::scale_fill_material(name = "", "blue-grey") +
  ggsci::scale_color_material(name = "", "blue-grey") +
  scale_x_discrete(position = 'bottom') +
  ggh4x::scale_y_dendrogram(hclust = hc_samples, position = "left", labels = NULL) +
  guides(y.sec = guide_axis_manual(labels = hc_order, label_size = 5, label_family = "GillSans")) +
  ggh4x::scale_x_dendrogram(hclust = hc_samples, position = "top", labels = NULL) +
  theme_bw(base_size = 7, base_family = "GillSans")
