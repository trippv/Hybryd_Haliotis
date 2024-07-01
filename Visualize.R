library(tidyverse, help, pos = 2, lib.loc = NULL)

dir <- "/Users/cigom/Documents/GitHub/Hybryd_Haliotis/Report"

f1 <- list.files(dir, pattern = "02.Assembly-alignment.tsv",  full.names = T)
f2 <- list.files(dir, pattern = "BUSCO.tsv", full.names = T)

# /Users/cigom/Documents/GitHub/Hybryd_Haliotis/Report/03.quantification

f3 <- list.files(dir, pattern = "FeatureCount", full.names = T) 


cols <- read_tsv(f1) %>%  select(starts_with("PE")) %>% names()

df1 <- read_tsv(f1) %>% select(-Treads, -`% Alignment`)

df1 <- df1 %>% 
  mutate(Method = ifelse(Method %in% "Hisat to Reference", "Reference-guide", Method)) %>%
  pivot_longer(cols = all_of(cols), names_to = "Category", values_to =  "Reads") 

df1 <- df1 %>% mutate(Sample = sapply(strsplit(Sample, "_"), `[`, 1))

df1 <- df1 %>% mutate(Category = factor(Category, levels = cols))

df1 <- df1 %>%
  group_by(Sample, Method) %>%
  mutate(pct_align = Reads / sum(Reads))


category <- "PE neither mate aligned"

p <- df1 %>%
  mutate(Method = stringr::str_to_title(Method)) %>%
  mutate(Method = factor(Method, levels = c("Trinity", "Rnaspades", "Reference-Guide"))) %>%
  mutate(Align = ifelse(Category != category, "Transcriptome coverage", "Unalignment")) %>%
  # group_by(Sample, Method, Align) %>%
  # summarise(pct_align = sum(pct_align)) %>%
  filter(Align == "Transcriptome coverage") %>%
  ggplot(aes(x = Method, y = pct_align, fill = Category)) +
  geom_col() +
  facet_grid(Sample ~ Align, scales = "free_x") +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  labs(x = "Assembly method", y = "% Alignment") +
  coord_flip() +
  scale_fill_grey("") +
  # scale_fill_manual("Assembly method", values = c("black", "grey89")) +
  guides(fill=guide_legend(nrow = 5)) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "bottom", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank())

# p

ggsave(p, filename = 'alignment-methods.png', path = dir, width = 4, height = 7, device = png, dpi = 300)



df2 <- read_tsv(f2) %>%
  mutate(my_species = gsub("transcripts_|Trinity_|rnaspades_", "", my_species)) %>%
  mutate(my_species = stringr::str_to_title(my_species)) %>%
  mutate(Method = stringr::str_to_title(Method)) %>%
  mutate(Method = ifelse(Method %in% "Refbased", "Reference-Guide", Method)) %>%
  filter(Method != "Denovo")

col <- c("#ED4647", "#EFE252", "#3A93E5", "#5BB5E7")

#names <- c("Complete", "Duplicated", "Fragmented", "Missing")

names <- c("S", "D", "F", "M")

labels <- c("Complete (C) and single-copy (S)",
            "Complete (C) and duplicated (D)",
            "Fragmented (F)  ",
            "Missing (M)")

col <- structure(col, names = rev(names))


my_sp_lev <- c("Mollusca","Metazoa","Eukaryota","Mammalia","Bacteria")

figure <- df2 %>%
  mutate(Method = factor(Method, levels = c("Trinity", "Rnaspades", "Reference-Guide"))) %>%
  # filter(my_species != "Bacteria") %>%
  mutate(facet = "Completeness") %>%
  mutate(category = factor(category, levels = rev(names))) %>%
  mutate(my_species = factor(my_species, levels = my_sp_lev)) %>%
  ggplot(aes(x = Method, y = my_percentage, fill = category, group = Method)) +
  # facet_grid(my_species ~ ., scales = "free") +
  facet_grid(my_species~ facet) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.75) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(y = "% BUSCOs", x = "", caption = "Transcriptome assembly") +
  coord_flip() +
  scale_fill_manual("", values = col, labels = rev(labels)) +
  guides(fill=guide_legend(nrow = 4)) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "bottom", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank())

figure

ggsave(figure, filename = 'BUSCO-methods.png', path = dir, width = 5, height = 6, device = png, dpi = 300)



f0 <- list.files(dir, pattern = "01.Contig-metrics.tsv", full.names = T)

cols <- read_tsv(f0) %>%  select(!starts_with("Metric")) %>% names()

df0 <- read_tsv(f0)

n_seqs <- c(76658,74226,766989,398352)
n_orfs <- c(56684,55862,92573,81604)

df0 <- rbind(df0, c("orfs", as.numeric(n_orfs/n_seqs)))

df0 <- df0 %>% 
  pivot_longer(cols = all_of(cols), names_to = "Method", values_to =  "Value") %>%
  mutate(Method = stringr::str_to_title(Method)) %>%
  filter(Method != "Denovo") %>%
  mutate(Method = factor(Method, levels = c("Trinity", "Rnaspades", "Reference-Guide"))) %>%
  mutate(Value = as.numeric(Value))

which_metrics <- c("n seqs", "n90","n50", "largest","orfs")

df0 %>%
  # mutate(facet = "Contig Metrics") %>%
  filter(Metric %in% which_metrics) %>%
  mutate(Metric = factor(Metric, levels = which_metrics)) %>%
  ggplot(aes(x = Method, y = Value, group = Method)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.75) +
  # ggh4x::facet_grid2(Metric ~., scales = "free_y") +
  # facet_grid(Metric ~., scales = "free_y") +
  facet_wrap(Metric ~ ., ncol = 1, scales = "free") +
  coord_flip() +
  # scale_fill_manual("", values = col, labels = rev(labels)) +
  guides(fill=guide_legend(nrow = 4)) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "bottom", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()) -> p0

p0

ggsave(p0, filename = 'metrics-methods.png', path = dir, width = 3.5, height = 6, device = png, dpi = 300)



df0 %>%
  filter(Metric %in% c("n under 200")) %>%
  # mutate(Metric = factor(Metric, levels = which_metrics)) %>%
  ggplot(aes(x = Method, y = Value, group = Method)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.75) +
  # ggh4x::facet_grid2(Metric ~., scales = "free_y") +
  # facet_grid(Metric ~., scales = "free_y") +
  facet_wrap(Metric ~ ., ncol = 1, scales = "free") +
  coord_flip() +
  # scale_fill_manual("", values = col, labels = rev(labels)) +
  guides(fill=guide_legend(nrow = 4)) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "bottom", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank())

for(i in rev(c(1:length(levels(my_species))))){
    detailed_values <- my_values[my_species==my_species[my_species==levels(my_species)[i]]]
    total_buscos <- sum(detailed_values)
    figure <- figure + 
    annotate("text", label=paste("C:", detailed_values[1] + detailed_values[2], " [S:", detailed_values[1], ", D:", detailed_values[2], "], F:", detailed_values[3], ", M:", detailed_values[4], ", n:", total_buscos, sep=""), 
             y=3, x = i, size = labsize*4*my_size_ratio, colour = "black", hjust=0, family=my_family)
    
    
    


}
