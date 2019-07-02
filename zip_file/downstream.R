## ----load downstream pkg,  warning=FALSE, message=FALSE------------------
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(dplyr)
  library(edgeR)
  library(scdney)
  library(mclust)
  library(Rtsne)
  library(parallel)
  library(cluster)
  library(ggplot2)
  library(MAST)
  library(viridis)
  library(ggpubr)
  library(plyr)
  library(monocle)
})

theme_set(theme_classic(16))



## ------------------------------------------------------------------------
sce_scMerge = readRDS("data/liver_scMerge.rds")
## We will subset Su et al. and Yang et al. datasets.
ids = colData(sce_scMerge)$batch %in% c("GSE87795", "GSE90047")

subset_data = sce_scMerge[,ids]

lab = colData(subset_data)$cellTypes

nCs = length(table(lab))

mat = SummarizedExperiment::assay(subset_data, "scMerge")


## ---- fig.width=12, fig.height=6, message=FALSE, echo = FALSE------------
cell_data = colData(subset_data) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("cell_name") %>% 
  dplyr::mutate(stage = stringr::str_sub(cell_name, 1, 3)) %>% 
  dplyr::group_by(cellTypes, stage, batch) %>%
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  tidyr::complete(cellTypes, stage, batch, fill = list(n = 0))

cell_data %>% 
  ggplot(aes(x = stage, y = cellTypes, fill = n, label = n)) +
  geom_tile() +
  geom_text() +
  facet_wrap(~batch) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  labs(title = "Number of cells split by batch, celltypes and stage")


## ---- warning = FALSE, eval = FALSE--------------------------------------
## ## For demonstration purpose, we will run k = 6 (which is actually the number of cell types in our dataset)
## simlr_result_k6 = scClust(mat, 6, similarity = "pearson", method = "simlr", seed = 1, cores.ratio = 0, geneFilter = 0)
## 
## load("data/simlr.results.RData")


## ---- eval = FALSE-------------------------------------------------------
## ## We will NOT run for various `k` to save time. Instead, we will load pre-computed results for `k` between 3 to 8
## 
## ## This is an easy way to run `scClust` for k = 3, 4, 5, 6, 7, 8.
## all_k = 3:8
## simlr_results = sapply(as.character(all_k), function(k) {
##   scClust(mat, as.numeric(k), similarity = "pearson", method = "simlr", seed = 1, cores.ratio = 0, geneFilter = 0)
## }, USE.NAMES = TRUE, simplify = FALSE)


## ---- include = FALSE----------------------------------------------------
## To make the knitting faster, we will load the precomputed results.
load("data/simlr.results.RData")
simlr_result_k6 = simlr_results$`6`


## ---- fig.height=6, fig.width=8------------------------------------------
# Find total WSS from all cluster outputs
all_wss = sapply(simlr_results, function(result) {
  sum(result$y$withinss)
}, USE.NAMES = TRUE, simplify = TRUE)

plot_data = data.frame(
  k = as.integer(names(all_wss)),
  total_wss = all_wss
)

ggplot(plot_data, 
       aes(x = k, 
           y = total_wss)) +
  geom_point(size = 3) +
  stat_smooth(method = loess, col = "red", 
              method.args = list(degree = 1), se = FALSE) +
  labs(title = "Compare Total WSS for each k",
       y = "Total WSS")


## ------------------------------------------------------------------------
## To run scClust with euclidean distance, uncommnet the following lines.
## simlr_result_eucl_k6 = scClust(mat, 6, similarity = "euclidean", method = "simlr", seed = 1, cores.ratio = 0, geneFilter = 0)

## for convenience, we will load our pre-computed result
load("data/simlr_result_eucl_k6.RData")


## ---- fig.width=10, fig.height=5-----------------------------------------
# create tsne object
set.seed(123)
tsne_result = Rtsne(t(mat), check_duplicates = FALSE)
#################################################
tmp_lab = as.numeric(factor(lab))
pear_cluster = plyr::mapvalues(
  simlr_result_k6$y$cluster,
  from = c(1,2,3,4,5,6),
  to = c(2,3,1,4,6,5)
)
eucl_cluster = plyr::mapvalues(
  simlr_result_eucl_k6$y$cluster,
  from = c(1,2,3,4,5,6),
  to = c(6,2,4,3,1,5)
)
#################################################

plot_data = data.frame(
  tsne1 = rep(tsne_result$Y[,1], 3),
  tsne2 = rep(tsne_result$Y[,2], 3),
  cluster = factor(c(tmp_lab, pear_cluster, eucl_cluster)),
  label = rep(c("Truth", "Pearson", "Euclidean"), each = length(lab)))


ggplot(plot_data, aes(x = tsne1, y = tsne2, colour = cluster)) +
  geom_point(size = 2) +
  labs(title = "t-SNE plot") +
  facet_grid(~label) + 
  theme(legend.position = "none")


## ---- fig.width=10, fig.height=5-----------------------------------------
plot_data2 = data.frame(
  Truth = lab,
  computed_cluster = as.factor(c(pear_cluster, eucl_cluster)),
  label = rep(c("Pearson", "Euclidean"), each = length(lab))
)


plot_data2 %>% 
  dplyr::group_by(Truth, computed_cluster, label) %>%
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  tidyr::complete(Truth, computed_cluster, label, fill = list(n = 0)) %>% 
  ggplot(aes(x = computed_cluster, 
             y = Truth,
             fill = n, label = n)) + 
  geom_tile() +
  geom_text() +
  facet_wrap(~label) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  labs(title = "")


## ---- fig.height=5, fig.width=10-----------------------------------------
# ARI
ari = c(mclust::adjustedRandIndex(lab, simlr_result_eucl_k6$y$cluster),
        mclust::adjustedRandIndex(lab, simlr_result_k6$y$cluster))

# NMI
nmi = c(igraph::compare(as.numeric(factor(lab)), 
                        simlr_result_eucl_k6$y$cluster, method = "nmi"), 
        igraph::compare(as.numeric(factor(lab)), 
                        simlr_result_k6$y$cluster, method = "nmi"))

plot_data = data.frame(
  dist = rep(c("Euclidean", "Pearson"), 2),
  value = c(ari, nmi),
  eval = rep(c("ARI", "NMI"), each = 2)
)

ggplot(plot_data, aes(x = dist, y = value, fill = dist)) + 
  geom_bar(stat="identity") + 
  facet_grid(col = vars(eval)) +
  labs(x = "Similarity metrics", 
       y = "Evalution score", 
       title = "Affect of similarity metrics in scRNA-seq data") +
  theme(legend.position = "none")


## ---- eval = FALSE-------------------------------------------------------
## marker_cluster4 = findMarker(mat = mat,
##                              cluster = simlr_result_k6$y$cluster,
##                              cluster_id = 4)


## ---- include = FALSE----------------------------------------------------
marker_cluster4 = readRDS("data/marker_cluster4.rds")


## ------------------------------------------------------------------------
marker_cluster4[1:10, ]


## ---- message=FALSE, warning=FALSE---------------------------------------
tsne_plotdf = data.frame(
  tsne1 = tsne_result$Y[, 1],
  tsne2 = tsne_result$Y[, 2]) %>% 
  dplyr::mutate(
    cluster = as.factor(simlr_result_k6$y$cluster),
    Gys2 = mat["Gys2", ])


ggplot(data = tsne_plotdf, aes(x = tsne1, y = tsne2, colour = Gys2) ) +
  geom_point(alpha = 0.5) +
  scale_color_viridis() +
  labs(col = "Gys2 expression", x = "tsne1", y = "tsne2")


## ------------------------------------------------------------------------
marker_cluster3 = readRDS("data/marker_cluster3.rds")


## ---- message=FALSE, warning=FALSE---------------------------------------
tsne_plotdf = tsne_plotdf %>% 
  dplyr::mutate(Erich5 = mat["Erich5",])

ggplot(data = tsne_plotdf, 
       mapping = aes(x = tsne1, y = tsne2, colour = Erich5)) +
  geom_point(alpha = 0.5) + 
  scale_color_viridis() +
  labs(col="Erich5 expression")


## ---- fig.height=6, fig.width=12-----------------------------------------
tsne_plotdf = tsne_plotdf %>% 
  dplyr::mutate(Hnf4a = mat["Hnf4a",])


fig1 = ggplot(data = tsne_plotdf, 
              mapping = aes(x = tsne1, y = tsne2)) + 
  geom_point(aes(color = ifelse(cluster == 4, 'Yellow', 'Purple')), alpha = 0.5) +
  scale_colour_viridis_d() + 
  labs(x = "", 
       y = "", 
       title = "Cluster 4") +
  theme(legend.position = "none") 


fig2 = ggplot(data = tsne_plotdf, 
              mapping = aes(x = tsne1, y = tsne2, colour = Hnf4a) ) +
  geom_point(alpha = 0.5) +
  scale_color_viridis() +
  labs(x = "", 
       y = "",
       title = "Hnf4a expression pattern")

ggarrange(fig1,fig2, ncol= 2, nrow = 1)


## ---- fig.height=6, fig.width=12-----------------------------------------
tsne_plotdf = tsne_plotdf %>% 
  dplyr::mutate(Epcam = mat["Epcam",])


fig1 = ggplot(data = tsne_plotdf, mapping = aes(x = tsne1, y = tsne2) ) + 
  geom_point(aes(color = ifelse(cluster == 3, 'Yellow', 'Purple')), alpha = 0.5) +
  scale_colour_viridis_d() + 
  labs(x = "", 
       y = "", 
       title = "Cluster 3") +
  theme(legend.position = "none") 


fig2 = ggplot(data = tsne_plotdf, 
               mapping = aes(x = tsne1, y = tsne2, colour = Epcam) ) +
  geom_point(alpha=0.5) +
  scale_color_viridis() +
  labs(x = "", 
       y = "",
       title = "Epcam expression pattern")

ggarrange(fig1,fig2, ncol= 2, nrow = 1)


## ---- fig.width = 10, fig.height=10--------------------------------------
plot_data = data.frame(table(lab)) %>% 
  dplyr::mutate(lab = reorder(lab, Freq))

ggplot(plot_data, 
       aes(x = lab, 
           y = Freq, 
           fill = lab)) +
  geom_bar(stat = "identity") +
  labs(x = "Cell types", 
       y = "Frequency", 
       title = "Composition of cell types") +
  theme(legend.position = "none")


## ---- message = FALSE, warning = FALSE-----------------------------------
## Subsetting data to "hepatoblast/hepatocyte"
monocle_data = subset_data[,colData(subset_data)$cellTypes %in% c("hepatoblast/hepatocyte")]
## Add a "stage" column to the colData of the monocle_data
colData(monocle_data)$stage = stringr::str_sub(colnames(monocle_data), 1, 3)
table(colData(monocle_data)$stage)
## monocle needs a rowData (data about each gene)
rowData(monocle_data) = DataFrame(gene_short_name = rownames(monocle_data))
monocle_data

## moncole requires a `CellDataSet` object to run. 
## You can convert monocle_data into a `CellDataSet` object using the scran package.
monocle_CellDataSet = scran::convertTo(
  monocle_data,
  type = "monocle",
  assay.type = "scMerge",
  col.fields = c("cellTypes", "stage", "batch"),
  row.fields = c("gene_short_name")) %>%
  estimateSizeFactors()

## Performing differential gene test using "stage". 
diff_test_res = differentialGeneTest(
  monocle_CellDataSet, fullModelFormulaStr = "~stage")

## We will select the top genes to be used for clustering and 
## calculate dispersion (variability) parameters before constructing the trajectory
ordering_genes = row.names(subset(diff_test_res, qval < 0.00001))
length(ordering_genes)
monocle_CellDataSet = setOrderingFilter(monocle_CellDataSet, ordering_genes)
monocle_CellDataSet = estimateDispersions(monocle_CellDataSet) %>% suppressWarnings()
# plot_ordering_genes(monocle_CellDataSet)

## Construcing the trajectory
monocle_CellDataSet = reduceDimension(monocle_CellDataSet, 
                                      max_components = 2,
                                      method = 'DDRTree')
monocle_CellDataSet = orderCells(monocle_CellDataSet)
plot_cell_trajectory(monocle_CellDataSet, color_by = "stage")
# monocle::plot_cell_clusters(monocle_CellDataSet)
plot_genes_in_pseudotime(monocle_CellDataSet[c("Sarm1", "Gm38388"),], color_by = "stage")


## ------------------------------------------------------------------------
sessionInfo()

