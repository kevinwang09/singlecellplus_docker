## ----load qc pkg, warning=FALSE, message=FALSE---------------------------
library(DropletUtils)
library(dplyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(stringr)
library(forcats)

theme_set(theme_classic(16))


## ------------------------------------------------------------------------
liver = read.csv(file = "data/GSE87795_counts.csv", 
                 header = TRUE,
                 row.names = 1)


## ------------------------------------------------------------------------
dim(liver)
liver[1:5, 1:5]


## ------------------------------------------------------------------------
stage = str_split(colnames(liver), "_") %>% 
  sapply("[[", 1)

table(stage)


## ------------------------------------------------------------------------
liver = liver[rowSums(liver) != 0, ]
dim(liver)


## ---- fig.height=6, fig.width=8------------------------------------------
barcode = DropletUtils::barcodeRanks(liver)
barcode_data = as.data.frame(barcode)

barcode_points = data.frame(
  type = c("inflection", "knee"),
  value = c(barcode@metadata$inflection, barcode@metadata$knee))


ggplot(data = barcode_data, aes(x = rank, y = total)) +
  geom_point() +
  geom_hline(data = barcode_points, 
             aes(yintercept = value, 
                 colour = type), linetype = 2) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(title = "Waterfall plot of read counts (log)", 
       x = "Log Rank", 
       y = "Log counts")


## ---- fig.height=6, fig.width=8------------------------------------------
cell_plotdf_full = data.frame(
  cell_name = colnames(liver),
  stage = stage, 
  library = colSums(liver),
  num_genes = colSums(liver != 0)
) 

head(cell_plotdf_full)

cell_plotdf = cell_plotdf_full %>% 
 dplyr::filter(stage %in% c("E11.5", "E12.5", "E13.5"))


cell_plotdf_mean = cell_plotdf %>% 
  dplyr::group_by(stage) %>% 
  dplyr::summarise(library_mean = mean(library), 
                   num_genes_mean = mean(num_genes))


cell_plotdf_mean

ggplot(cell_plotdf, 
       aes(x = library/1e6, colour = stage, fill = stage)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", bins = 30) + 
  geom_density(alpha = 0.4) + 
  geom_vline(data = cell_plotdf_mean, 
             aes(xintercept = library_mean/1e6), 
             colour = "blue", linetype = "dashed", size = 1.5) +
  # scale_x_continuous(breaks = 1:6) +
  labs(title = "Histogram of library size", 
       x = "library size (in millions)", 
       y = "density") + 
  facet_wrap(~ stage, ncol = 3) 


## ---- fig.height=6, fig.width=8, message=FALSE, warning=FALSE------------
ggplot(cell_plotdf, 
       aes(x = num_genes/1000, colour = stage, fill = stage)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", bins = 30) + 
  geom_density(alpha = 0.4) + 
  geom_vline(data = cell_plotdf_mean,
             aes(xintercept = num_genes_mean/1000), 
             colour = "blue", linetype = "dashed", size = 1.5) +
  labs(title = "Histogram of number of unique genes expressed", 
       x = "Number of unique genes expressed (in thousands)", 
       y = "density") + 
  facet_wrap(~ stage, ncol = 3) 


## ---- message = FALSE, warning = FALSE-----------------------------------
mito_genes <- c("ENSMUSG00000064336", "ENSMUSG00000064337", "ENSMUSG00000064338", "ENSMUSG00000064339", "ENSMUSG00000064340", "ENSMUSG00000064341", "ENSMUSG00000064342", "ENSMUSG00000064343", "ENSMUSG00000064344", "ENSMUSG00000064345", "ENSMUSG00000064346", "ENSMUSG00000064347", "ENSMUSG00000064348", "ENSMUSG00000064349", "ENSMUSG00000064350", "ENSMUSG00000064351", "ENSMUSG00000064352", "ENSMUSG00000064353", "ENSMUSG00000064354", "ENSMUSG00000064355", "ENSMUSG00000064356", "ENSMUSG00000064357", "ENSMUSG00000064358", "ENSMUSG00000064359", "ENSMUSG00000064360", "ENSMUSG00000064361", "ENSMUSG00000065947", "ENSMUSG00000064363", "ENSMUSG00000064364", "ENSMUSG00000064365", "ENSMUSG00000064366", "ENSMUSG00000064367", "ENSMUSG00000064368", "ENSMUSG00000064369", "ENSMUSG00000064370", "ENSMUSG00000064371", "ENSMUSG00000064372", "ENSMUSG00000096105")


liver_mito = liver[mito_genes, ]

mito_plotdf_full = data.frame(
  cell_name = colnames(liver_mito),
  stage = stage, 
  mito_percent = colSums(liver_mito)/colSums(liver) * 100
) 

mito_plotdf = mito_plotdf_full %>% dplyr::filter(stage %in% c("E11.5", "E12.5", "E13.5"))

mito_plotdf_mean = mito_plotdf %>% 
  group_by(stage) %>% 
  dplyr::summarise(mito_percent_mean = mean(mito_percent))

mito_plotdf_mean

g1 = ggplot(mito_plotdf, 
       aes(x = mito_percent, colour = stage, fill = stage)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", bins = 30) + 
  geom_density(alpha = 0.4) + 
  geom_vline(data = mito_plotdf_mean,
             aes(xintercept = mito_percent_mean), 
             colour = "blue", linetype = "dashed", size = 1.5) +
  labs(x = "Mitochondrial percentage", 
       y = "Density") + 
  facet_wrap(~ stage, nrow = 3) 


g2 = ggplot(mito_plotdf, 
            aes(x = stage, 
                y = mito_percent, 
                fill = stage)) +
  geom_violin(alpha = 0.8, draw_quantiles = 0.5) + 
  labs(y = "Mitochondrial percentage")

ggpubr::ggarrange(g1, g2, ncol=  2, common.legend = TRUE, legend = "right")


## ----  fig.height=12, fig.width=8----------------------------------------
## We will go to each column of the matrix, and divide the column by its sum 
## and convert this into a percentage. 
## The result is a matrix with each column that sum up to 100%. 
liver_percent = sweep(liver, 2, STATS = colSums(liver), FUN = "/") * 100

liver_percent %>% colSums %>% head

liver_percent_meanEachGene = rowMeans(liver_percent)

top50Genes = liver_percent_meanEachGene %>% sort(decreasing = TRUE) %>% head(50) %>% names


liver_percent_plotdf = liver_percent[top50Genes, ] %>% 
  as.matrix %>% 
  reshape2::melt(varnames = c("gene_name", "cell_name"),
                 value.name = "percent") %>% 
  dplyr::mutate(gene_name = fct_reorder(gene_name, percent, mean))

liver_percent_plotdf %>% 
  ggplot(aes(x = gene_name, y = percent, colour = gene_name)) +
  geom_violin() +
  geom_jitter(alpha = 0.2) +
  coord_flip() +
  labs(x = "percent of library", 
       y = "genes",
       title = "Top 50 genes for all cell stages") +
  theme(legend.position = "none")


## ----  fig.height=12, fig.width=8----------------------------------------
liver_percent_split = split.data.frame(x = t(liver_percent), f = stage) %>% lapply(t)

liver_percent_split_meanEachGene = liver_percent_split %>% lapply(rowMeans)


e115_top50 = liver_percent_split_meanEachGene$E11.5 %>% sort(decreasing = TRUE) %>% head(50) %>% names

e115_liver_percent = liver_percent[e115_top50, stage == "E11.5"] %>% 
  as.matrix %>% 
  reshape2::melt(varnames = c("gene_name", "cell_name"),
                 value.name = "percent") %>% 
  dplyr::mutate(gene_name = fct_reorder(gene_name, percent, mean))

e115_liver_percent %>% 
  ggplot(aes(x = gene_name, y = percent, colour = gene_name)) +
  geom_violin() +
  geom_jitter(alpha = 0.2) +
  coord_flip() +
  labs(x = "genes", 
       y = "percent of library",
       title = "Top 50 expressed genes for E11.5 cells") +
  theme(legend.position = "none")


## ------------------------------------------------------------------------
table(barcode_data$total >= barcode@metadata$inflection)
liver_new = liver[, barcode_data$total >= barcode@metadata$inflection]
dim(liver_new)


## ---- eval = FALSE-------------------------------------------------------
## dim(liver)
## liver_new2 = liver[, barcode_data$total >= barcode@metadata$inflection & mito_plotdf_full$mito_percent <= 10]
## dim(liver_new2)


## ------------------------------------------------------------------------
sessionInfo()

