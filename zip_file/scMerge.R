## ----load scMerge pkg,  warning=FALSE, message=FALSE---------------------
library(scMerge)
library(scater)
library(dplyr)
library(ggpubr)
library(forcats)
library(dplyr)
library(tidyr)

theme_set(theme_classic(16))


## ------------------------------------------------------------------------
su = readRDS("data/sce_GSE87795.rds")
yang = readRDS("data/sce_GSE90047.rds")


## ------------------------------------------------------------------------
sce_list = list(
  su = su, 
  yang = yang
)

sce_list

sce_combine = scMerge::sce_cbind(sce_list = sce_list, 
                                 method = "union", 
                                 colData_names = c("cellTypes", "stage"),
                                 batch_names = c("Su", "Yang"))

sce_combine


## ------------------------------------------------------------------------
table(
  colData(sce_combine)$cellTypes, 
  colData(sce_combine)$batch
)

colData(sce_combine)$cellTypes = colData(sce_combine)$cellTypes %>% 
  forcats::fct_recode(Hepatoblast = "hepatoblast/hepatocyte") %>% 
  droplevels()

sce_combine = sce_combine[rowSums(SingleCellExperiment::counts(sce_combine)) != 0,
                          colSums(SingleCellExperiment::counts(sce_combine)) != 0]

table(
  colData(sce_combine)$cellTypes, 
  colData(sce_combine)$batch
)


## ---- fig.width=12, fig.height=6, message=FALSE, echo = FALSE------------
cell_data = colData(sce_combine) %>% 
  as.data.frame() %>% 
  dplyr::group_by(cellTypes, stage, batch) %>%
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  tidyr::complete(cellTypes, stage, batch, fill = list(n = 0)) %>% 
  dplyr::mutate(cellTypes = as.character(cellTypes))

cell_data %>% 
  ggplot(aes(x = stage, y = cellTypes, fill = n, label = n)) +
  geom_tile() +
  geom_text() +
  facet_wrap(~batch) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  labs(title = "Number of cells split by batch, celltypes and stage")


## ---- fig.width=12, fig.height=6, message=FALSE--------------------------
set.seed(1234)
tsne_logcounts_cellTypes = scater::plotTSNE(sce_combine, 
                                          colour_by = "cellTypes",
                                          run_args = list(exprs_values = "logcounts")) +
  scale_fill_brewer(palette = "Set1")


set.seed(1234)
tsne_logcounts_batch = scater::plotTSNE(sce_combine, 
                                      colour_by = "batch",
                                      run_args = list(exprs_values = "logcounts")) +
  scale_fill_brewer(palette = "Dark2")

ggpubr::ggarrange(tsne_logcounts_cellTypes, tsne_logcounts_batch, ncol = 2, nrow = 1)


## ---- fig.width=5, fig.height=5------------------------------------------
data("segList_ensemblGeneID", package = "scMerge")

scMerge_supervised = scMerge(
  sce_combine = sce_combine,
  ctl = which(rownames(sce_combine) %in% segList_ensemblGeneID$mouse$mouse_scSEG),
  cell_type = sce_combine$cellTypes,
  replicate_prop = 1,
  assay_name = "scMerge_supervised",
  verbose = TRUE)


## ---- fig.width=12, fig.height=6, message=FALSE--------------------------
scMerge_supervised

set.seed(1234)
tsne_scMerge_supervised_cellTypes = scater::plotTSNE(scMerge_supervised, 
                                                   colour_by = "cellTypes",
                                                   run_args = list(exprs_values = "scMerge_supervised")) +
  scale_fill_brewer(palette = "Set1")

set.seed(1234)
tsne_scMerge_supervised_batch = scater::plotTSNE(scMerge_supervised, 
                                               colour_by = "batch",
                                               run_args = list(exprs_values = "scMerge_supervised")) +
  scale_fill_brewer(palette = "Dark2")

ggpubr::ggarrange(tsne_scMerge_supervised_cellTypes, tsne_scMerge_supervised_batch, ncol = 2, nrow = 1)


## ---- fig.width=5, fig.height=5------------------------------------------
scMerge_unsupervised = scMerge(
  sce_combine = sce_combine,
  ctl = which(rownames(sce_combine) %in% segList_ensemblGeneID$mouse$mouse_scSEG),
  kmeansK = c(6, 2),
  replicate_prop = 1,
  assay_name = "scMerge_unsupervised",
  verbose = TRUE)


## ---- fig.width=12, fig.height=6, message=FALSE--------------------------
scMerge_unsupervised

set.seed(1234)
tsne_scMerge_unsupervised_cellTypes = scater::plotTSNE(scMerge_unsupervised, 
                                                     colour_by = "cellTypes",
                                                     run_args = list(exprs_values = "scMerge_unsupervised")) +
  scale_fill_brewer(palette = "Set1")

set.seed(1234)
tsne_scMerge_unsupervised_batch = scater::plotTSNE(scMerge_unsupervised, 
                                                 colour_by = "batch",
                                                 run_args = list(exprs_values = "scMerge_unsupervised")) +
  scale_fill_brewer(palette = "Dark2")

ggpubr::ggarrange(tsne_scMerge_unsupervised_cellTypes, tsne_scMerge_unsupervised_batch, ncol = 2, nrow = 1)


## ---- fig.width=5, fig.height=5------------------------------------------
scMerge_semisupervised = scMerge(
  sce_combine = sce_combine,
  ctl = which(rownames(sce_combine) %in% segList_ensemblGeneID$mouse$mouse_scSEG),
  kmeansK = c(6, 2),
  replicate_prop = 1,
  WV = sce_combine$stage,
  WV_marker = c("ENSMUSG00000045394","ENSMUSG00000054932","ENSMUSG00000045394"),
  assay_name = "scMerge_semisupervised",
  verbose = TRUE)


## ---- fig.width=12, fig.height=6, message=FALSE--------------------------
scMerge_semisupervised


set.seed(1234)
tsne_scMerge_semisupervised_cellTypes = scater::plotTSNE(scMerge_semisupervised, 
                                                       colour_by = "cellTypes",
                                                       run_args = list(exprs_values = "scMerge_semisupervised")) +
  scale_fill_brewer(palette = "Set1")


set.seed(1234)
tsne_scMerge_semisupervised_batch = scater::plotTSNE(scMerge_semisupervised, 
                                                   colour_by = "batch",
                                                   run_args = list(exprs_values = "scMerge_semisupervised")) +
  scale_fill_brewer(palette = "Dark2")

ggpubr::ggarrange(tsne_scMerge_semisupervised_cellTypes, tsne_scMerge_semisupervised_batch, ncol = 2, nrow = 1)


## ---- fig.width=12, fig.height=6, message=FALSE--------------------------
set.seed(1234)
tsne_scMerge_unsupervised_hept = scater::plotTSNE(scMerge_unsupervised[, scMerge_unsupervised$cellTypes == "Hepatoblast"], 
                                                       colour_by = "stage",
                                                       run_args = list(exprs_values = "scMerge_unsupervised")) +
  scale_fill_brewer(palette = "Set1")

tsne_scMerge_semisupervised_hept = scater::plotTSNE(scMerge_semisupervised[, scMerge_semisupervised$cellTypes == "Hepatoblast"], 
                                                       colour_by = "stage",
                                                       run_args = list(exprs_values = "scMerge_semisupervised")) +
  scale_fill_brewer(palette = "Set1")

ggpubr::ggarrange(tsne_scMerge_unsupervised_hept, tsne_scMerge_semisupervised_hept, ncol = 2, nrow = 1)


## ---- eval = FALSE, fig.width=5, fig.height=5----------------------------
## scMerge_fast <- scMerge(
##   sce_combine = sce_combine,
##   ctl = which(rownames(sce_combine) %in% segList_ensemblGeneID$mouse$mouse_scSEG),
##   cell_type = sce_combine$cellTypes,
##   replicate_prop = 1,
##   assay_name = "scMerge_supervised",
##   verbose = TRUE,
##   fast_svd = TRUE,
##   rsvd_prop = 0.1)


## ---- eval = FALSE, fig.width=12, fig.height=6, message=FALSE------------
## scMerge_fast
## 
## set.seed(1234)
## tsne_scMerge_fast_cellTypes = scater::plotTSNE(scMerge_fast,
##                                                colour_by = "cellTypes",
##                                                run_args = list(exprs_values = "scMerge_fast")) +
##   scale_fill_brewer(palette = "Set1")
## 
## 
## set.seed(1234)
## tsne_scMerge_fast_batch = scater::plotTSNE(scMerge_fast,
##                                            colour_by = "batch",
##                                            run_args = list(exprs_values = "scMerge_fast")) +
##   scale_fill_brewer(palette = "Dark2")
## 
## ggpubr::ggarrange(tsne_scMerge_fast_cellTypes, tsne_scMerge_fast_batch, ncol = 2, nrow = 1)


## ------------------------------------------------------------------------
sessionInfo()

