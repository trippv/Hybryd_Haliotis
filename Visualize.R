library(tidyverse, help, pos = 2, lib.loc = NULL)

dir <- "/Users/cigom/Documents/GitHub/Hybryd_Haliotis/Report"

# f1 <- list.files(dir, pattern = "02.Assembly-alignment.tsv",  full.names = T)
f2 <- list.files(dir, pattern = "BUSCO.tsv", full.names = T)

df2 <- read_tsv(f2) %>%
  mutate(my_species = gsub("transcripts_", "", my_species)) %>%
  mutate(my_species = stringr::str_to_title(my_species)) %>%
  mutate(Method = stringr::str_to_title(Method))

col <- c("#ED4647", "#EFE252", "#3A93E5", "#5BB5E7")

#names <- c("Complete", "Duplicated", "Fragmented", "Missing")

names <- c("S", "D", "F", "M")

labels <- c("Complete (C) and single-copy (S)",
            "Complete (C) and duplicated (D)",
            "Fragmented (F)  ",
            "Missing (M)")

col <- structure(col, names = rev(names))

figure <- df2 %>%
  filter(my_species != "Bacteria") %>%
  mutate(category = factor(category, levels = rev(names))) %>%
  ggplot(aes(x = Method, y = my_percentage, fill = category, group = Method)) +
  facet_grid(my_species ~ ., scales = "free") +
  geom_col(position = position_stack(reverse = TRUE), width = 0.75) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(y = "% BUSCOs", x = "", caption = "Transcriptome assembly") +
  coord_flip() +
  scale_fill_manual("", values = col, labels = rev(labels)) +
  theme_grey(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top")

figure


for(i in rev(c(1:length(levels(my_species))))){
    detailed_values <- my_values[my_species==my_species[my_species==levels(my_species)[i]]]
    total_buscos <- sum(detailed_values)
    figure <- figure + 
    annotate("text", label=paste("C:", detailed_values[1] + detailed_values[2], " [S:", detailed_values[1], ", D:", detailed_values[2], "], F:", detailed_values[3], ", M:", detailed_values[4], ", n:", total_buscos, sep=""), 
             y=3, x = i, size = labsize*4*my_size_ratio, colour = "black", hjust=0, family=my_family)
  }