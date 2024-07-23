# calculate nx distribution


rm(list = ls())

if(!is.null(dev.list())) dev.off()

library(tidyverse, help, pos = 2, lib.loc = NULL)


require(Biostrings)
require(dplyr)
require(ggplot2)

# options(base.)

dir <- "/Users/cigom/Documents/GitHub/Hybryd_Haliotis/Report"


f <- list.files(dir, "fasta", full.names = T)

contig_Nx <- function(f) {
  
  DNA <- Biostrings::readDNAStringSet(f)
  
  contig_width <- sort(Biostrings::width(DNA), decreasing = T)
  
  # contig_df <- data.frame(width = contig_width, Assembly = basename(f))
  
  return(contig_width)
}


metrics_df <- function(f) {
  
  width <- contig_Nx(f)
  
  n_seqs <- length(width)
  
  v <- seq(0.1,1, by = 0.1)
  
  Lx <- function(x) { sum(cumsum(width) < (sum(width) * x)) + 1 }
  
  # add the number of sequences per Nx (field n_seqs)
  
  # Lx <- function(x) { sum(cumsum(width) < (sum(width) * x)) + 1 }
  
  l <- unlist( lapply(v, Lx))
  
  metrics_df <- data.frame(x = paste0("N", v*100), n = width[l], l = l,
    n_seqs = n_seqs-l, Assembly = basename(f))

  
  return(metrics_df)
  
}


df <- lapply(f, metrics_df)

df <- do.call(rbind, df)

df <- mutate(df, x = factor(x, levels = unique(df$x)))

unique(df$Assembly)

recode_to <- structure(c("Reference-Guide","Rnaspades","Trinity"), names = unique(df$Assembly))

df <- mutate(df, Assembly = dplyr::recode_factor(Assembly, !!!recode_to))

# rnsps <- mean(contig_Nx(f[[1]]))
# trnt <- mean(contig_Nx(f[[2]]))

p <- ggplot(df, aes(x = x, y = n, group = Assembly, color = Assembly)) +
  geom_vline(xintercept = "N50", linetype="dashed", alpha=0.5) +
  geom_point() +
  ggplot2::geom_path() +
  labs(x = "Nx", y = "Contig length", color = "Assembly method") +
  scale_color_grey("") +
  # scale_fill_manual("Assembly method", values = c("black", "grey89")) +
  guides(color=guide_legend(nrow = 1)) +
  theme_bw(base_size = 11, base_family = "GillSans") +
  theme(legend.position = "top", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()) 
# annotate("text", y = rnsps, x = "N50", angle = 90, label = "label")

# p

ggsave(p, filename = 'Nx-methods.png', path = dir, width = 4, height = 3, device = png, dpi = 300)

# calculat pca for every count-matrix, 
# plot variance of explanation barplot per method
# add metadata

dir <-  "/Users/cigom/Documents/GitHub/Hybryd_Haliotis/Report/03.quantification/"

f <- list.files(dir, pattern = "_count.txt", full.names = T) 

# cols <- read_tsv(f, skip = 1) %>% select(contains(".sorted.bam")) %>% names()

PCA <- function(f) {
  
  df <- read_tsv(f, skip = 1) 
  
  m <- df %>% select(contains(".sorted.bam")) %>% as("matrix")
  
  data <- DESeq2::vst(round(m) + 1)

  PCA <- prcomp(t(data), center = T, scale. = FALSE)
  
  PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], Method = basename(f))
  
  # percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
  
  percentVar <- round(PCA$sdev/sum(PCA$sdev)*100,1)

  PCAvar <- data.frame(
                      Eigenvalues = PCA$sdev,
                      percentVar = percentVar,
                      Varcum = cumsum(percentVar),
                      Method = basename(f))
  
  return(PCAdf)
}

PCAdf <- lapply(f, PCA)
PCAdf <- do.call(rbind, PCAdf)

recode_to <- structure(c("Reference-Guide","Rnaspades","Trinity"), names = unique(PCAdf$Method))

PCAdf <- mutate(PCAdf, Method = dplyr::recode_factor(Method, !!!recode_to))

PCAdf %>%
  group_by(Method) %>%
  #dplyr::as_tibble(rownames = "ID")
  mutate(Dim = row_number()) %>%
  ggplot(., aes(y = percentVar, x = as.factor(Dim), fill = Method, color = Method)) +
  geom_col(position = position_dodge2()) +
  geom_line(aes(y = Varcum, group = Method)) +
  geom_point(aes(y = Varcum)) +
  # labs(x = "Component Number", y = "Eigenvalue") +
  labs(x = "Principal component", y = "Fraction variance explained (%)") +
  scale_fill_grey("") +
  scale_color_grey("") +
  guides(color=guide_legend(nrow = 1)) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()) -> p 

ggsave(p, filename = 'Eigenevalues.png', path = dir, width = 4, height = 3, device = png, dpi = 300)

PCAdf <- lapply(f, PCA)

PCAdf <- do.call(rbind, PCAdf)

PCAdf %>%
  mutate(LIBRARY_ID = rownames(.)) %>%
  ggplot(., aes(PC1, PC2, color = Method)) +
  scale_color_grey("") +
  facet_grid(~ Method) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  geom_point(size = 7, alpha = 0.7) 
  # ggforce::geom_mark_ellipse(aes(group = as.factor(hpf), label = as.factor(hpf)),
  #   fill = 'grey', colour = NA) +


cor_heights <- function(f) {
  
  df <- read_tsv(f, skip = 1) 
  
  m <- df %>% select(contains(".sorted.bam")) %>% as("matrix")
  
  data <- DESeq2::vst(round(m) + 1) # vst if cols > 10
  
  sample_cor <- cor(data, method='pearson', use='pairwise.complete.obs')
  
  h <- heatmap(sample_cor, keep.dendro = T)
  
  hc_samples <- as.hclust(h$Colv)
  
  out <- data.frame(height = hc_samples$height, 
                    # labels = hc_samples$labels, 
                    # order = hc_samples$order,
                    Method = basename(f))
  
}

cor_df <- lapply(f, cor_heights)
cor_df <- do.call(rbind, cor_df)

recode_to <- structure(c("Reference-Guide","Rnaspades","Trinity"), names = unique(cor_df$Method))

# cor_df <- mutate(cor_df, Method = dplyr::recode_factor(Method, !!!recode_to))
cor_df <- mutate(cor_df, Method = dplyr::recode_factor(Method, !!!rev(recode_to)))

# The height of the branch points indicates how similar or different they are from each other: the greater the height, the greater the difference.

cor_df %>% 
  group_by(Method) %>%
  mutate(row_number = row_number()) %>%
  ggplot(aes(x = height, y = Method)) + 
  geom_path() +
  geom_point(shape = 21, size = 1.5, stroke = 1.2, fill = "white", color = "black") +
  #1.4 Align labels on the top or bottom edge
  #Use hjust or vjust to justify the text neatly:
  #hjust = 0 for left-align
  #hjust = 0.5 for center
  #hjust = 1 for right-align
  ggrepel::geom_text_repel(
    aes(label = row_number),
    nudge_y      = 0.05,
    direction    = "x",
    angle        = 0,
    vjust        = 0,
    hjust        = 1,
    segment.size = 0.2,
    family = "GillSans"
  ) +
  guides(color=guide_legend(nrow = 1)) +
  theme_bw(base_size = 14, base_family = "GillSans") +
  theme(legend.position = "top", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()) -> p


ggsave(p, filename = 'heights.png', path = dir, width = 7, height = 3, device = png, dpi = 300)

# from refdata

read_FtrCnt <- function(f) {
  
  .m <- read_tsv(f, skip = 1) 
  
  m <- .m %>% select(contains(".sorted.bam")) %>% as("matrix")
  
  colNames <- colnames(m)
  
  colnames(m) <- gsub("../|.sorted.bam", "", colNames)
  
  # m <- DESeq2::vst(round(m)+1) # vst if cols > 10
  
  m <- select(.m, !contains(".sorted.bam")) %>% cbind(m)
  
  m <- as_tibble(m)
  
  return(m)
}

df <- read_FtrCnt(f[1]) 

df %>% count(Chr)

samNames <- df %>% select(starts_with("C")) %>% select_if(is_double) %>% colnames()

df_summ <- df %>% select(contains(c("Chr","Geneid", samNames))) %>%
  pivot_longer(cols = all_of(samNames)) %>%
  mutate(name = substr(name, 1,3)) %>%
  group_by(name, Geneid) %>% summarise(value = sum(value)) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  ungroup() 

# Prevalence <- df_summ %>% select_if(is.double)

recode_to <- structure(c("Hybryd", "Red", "Green"), names = c("CRF","CRR","CFF"))

UPSETDF <- df_summ %>%
  pivot_longer(-Geneid) %>%
  dplyr::mutate(name = dplyr::recode_factor(name, !!!recode_to)) %>%
  filter(value >0) %>%
  group_by(Geneid) %>%
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
  labs(x = '', y = 'Genes number')


df %>%
  mutate(Chr = strsplit(Chr, ";")) %>%
  unnest(Chr) %>% 
  # group_by(Geneid) %>% # turn on to calculate per gene 
  count(Chr) %>% 
  mutate(Genome = gsub("^HiC_scaffold_[0-9]+", "HiC_scaffold", as.character(Chr))) %>%
  mutate(Genome = gsub("^JALGQA010000.[0-9]+", "JALGQA010000", as.character(Genome))) %>%
  ungroup() %>%
  count(Genome)
  

as_tibble(df)

recode_to <- structure(c("Reference-Guide","Rnaspades","Trinity"), 
  names = unique(cor_df$Method))

refbase <- vsn::meanSdPlot(read_FtrCnt(f[1]), plot = F)
rnasp <- vsn::meanSdPlot(read_FtrCnt(f[2]), plot = F)
trin <- vsn::meanSdPlot(read_FtrCnt(f[3]), plot = F)


rbind(data.frame(py = refbase$sd, px = refbase$rank, col = "Reference-Guide"),
  data.frame(py = rnasp$sd, px = rnasp$rank, col = "Rnaspades"),
  data.frame(py = trin$sd, px = trin$rank, col = "Trinity")) %>%
  ggplot(aes(px, py, color = col)) +
  labs(x = "Ranks", y = "sd", color = "") +
  geom_line(orientation = NA, position = position_identity(), size = 2) +
  theme_bw(base_family = "GillSans", base_size = 20) +
  theme(legend.position = "top") +
  scale_fill_grey("") +
  scale_color_grey("") 

# Prevalence

apply(m, 1, function(x) sum(x > 0)) %>% table()

prevalence <- function(f) {
  
  m <- read_tsv(f, skip = 1) 
  
  m <- m %>% select(contains(".sorted.bam")) %>% as("matrix")
  
  colnames(m) <- gsub("../|.sorted.bam", "", colnames(m))
  
  prevelancedf <- apply(m, 1, function(x) sum(x > 0))
  
  data.frame(Prevalence = prevelancedf, 
    TotalAbundance = rowSums(m)) %>% # mean_se
    as_tibble(rownames = "GeneRank") %>%
    arrange(desc(TotalAbundance)) -> prevelancedf
  
  prevelancedf <- prevelancedf %>% count(Prevalence) %>% 
    mutate(Method = basename(f), pct = n/sum(n))
  
  return(prevelancedf)
}

prevelancedf <- lapply(f, prevalence)

prevelancedf <- do.call(rbind, prevelancedf)

recode_to <- structure(c("Reference-Guide","Rnaspades","Trinity"), names = unique(prevelancedf$Method))

prevelancedf <- mutate(prevelancedf, Method = dplyr::recode_factor(Method, !!!recode_to))


prevelancedf %>% group_by(Method) %>% summarise(sum(n))

prevelancedf %>% 
  ggplot(aes(x = as.factor(Prevalence), y = pct, fill = Method)) + 
  geom_col(position = position_dodge2()) +
  theme_classic(base_family = "GillSans") + 
  scale_y_continuous("Frac. of genes", labels = scales::percent) +
  labs(x = "Frequency in samples") +
  scale_fill_grey("") +
  scale_color_grey("") +
  theme_bw(base_family = "GillSans", base_size = 20) +
  theme(legend.position = "top")
  # ylim(0,60000)

# z_scores <- function(x) {(x-mean(x))/sd(x)}

# data <- apply(data, 1, z_scores)

# calculate cor matrix 

m <- read_FtrCnt(f[1])

sample_cor <- cor(m, method='pearson', use='pairwise.complete.obs')

h <- heatmap(sample_cor, col = cm.colors(12), keep.dendro = T)

hc_samples <- as.hclust(h$Colv)

hc_samples$height

plot(hc_samples)
hc_order <- hc_samples$labels[h$colInd]


sample_cor %>% 
  as_tibble(rownames = 'LIBRARY_ID') %>%
  pivot_longer(cols = colnames(sample_cor), values_to = 'cor') %>%
  # left_join(.colData) %>% 
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = rev(hc_order))) -> sample_cor_long

library(ggh4x)

sample_cor_long %>%
  ggplot(aes(x = LIBRARY_ID, y = name, fill = cor)) +  
  geom_tile(linewidth = 0.2) +
  ggsci::scale_fill_material(name = "", "blue-grey") +
  ggsci::scale_color_material(name = "", "blue-grey") +
  scale_x_discrete(position = 'bottom') +
  ggh4x::scale_y_dendrogram(hclust = hc_samples, position = "left", labels = NULL) +
  guides(y.sec = guide_axis_manual(labels = hc_order, label_size = 5, label_family = "GillSans")) +
  # ggh4x::scale_x_dendrogram(hclust = hc_samples, position = "top", labels = NULL) +
  theme_bw(base_size = 7, base_family = "GillSans")

# make expression weighted gene length

len <- df$Length; # Biostrings::width(DNA)
expr <- df %>% select(contains(".sorted.bam")) %>% rowSums()

sum_expr_n_len <- len * expr

sum_expr <- sum(expr)
gene_len <- sum_expr_n_len / sum_expr;

ExN50df <- data.frame(lenth = len, max_expr_over_samples = expr, 
  sum_expr_over_samples = sum_expr_n_len)

