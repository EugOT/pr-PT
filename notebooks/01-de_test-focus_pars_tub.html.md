---
title: "Differential expression analysis of Hypothalamus datasets from Kim DW et al., 2020 and Romanov et al., 2020 with focus on Pars Tuberalis"
author: "Evgenii O. Tretiakov"
date: "2025-01-30"
format:
  elsevier-html:
    toc: true
    df-print: paged
    code-fold: true
    fig-width: 9
    fig-height: 12
    fig-format: retina
    fig-responsive: true
    fig-dpi: 120
  elsevier-pdf:
    colorlinks: true
    fontsize: 12pt
execute:
  keep-md: true
  echo: true
  error: false
  message: false
  warning: false
  debug: false
knitr:
  opts_chunk:
    autodep: true
---


::: {.cell .hidden}

```{.r .cell-code .hidden}
#| label: setup
#| include: false
DOCNAME <- "pars-tuberalis-analysis"
NOW <- Sys.time()

# Time chunks during knitting
knitr::knit_hooks$set(timeit = function(before) {
  if (before) {
    print(paste("Start:", Sys.time()))
    NOW <<- Sys.time()
  } else {
    print(paste("Stop:", Sys.time()))
    print(Sys.time() - NOW)
  }
})

knitr::knit_hooks$set(debug = function(before, options, envir) {
  if (!before) {
    message(
      paste(names(envir), as.list(envir),
        sep = " = ", collapse = "\n"
      )
    )
  }
})

knitr::opts_chunk$set(
  cache          = FALSE,
  dev            = c("png", "pdf"),
  timeit         = TRUE
)
```
:::





## Load data and setup parameters





::: {.cell}

```{.r .cell-code .hidden}
#| label: libraries
#| cache: false
# Load tidyverse infrastructure packages
suppressPackageStartupMessages({
  library(future)
  library(here)
  library(tidyverse)
  library(magrittr)
  library(stringr)
  library(skimr)
  library(RColorBrewer)
  library(viridis)
})


# Load packages for scRNA-seq analysis and visualisation
suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(patchwork)
  library(ggstatsplot)
  library(Seurat)
  library(SeuratWrappers)
  library(scCustomize)
})
```
:::





### Set paths





::: {.cell}

```{.r .cell-code .hidden}
#| label: paths
src_dir <- here("code")
data_dir <- here("data")
output_dir <- here("output")
plots_dir <- here(output_dir, "figures/")
tables_dir <- here(output_dir, "tables/")
```
:::





### Load helper functions and gene-sets





::: {.cell}

```{.r .cell-code .hidden}
#| label: source
#| cache: false
source(here(src_dir, "genes.R"))
source(here(src_dir, "functions.R"))
```
:::





### Set fixed variables





::: {.cell}

```{.r .cell-code .hidden}
#| label: params-computation
#| cache: false
# set seed
reseed <- 42
set.seed(seed = reseed)

# Parameters for parallel execution
n_cores <- parallelly::availableCores() / 2 - 1
plan("multisession", workers = n_cores)
options(
  future.globals.maxSize = 100000 * 1024^2,
  future.rng.onMisuse = "ignore"
)
plan()
```

::: {.cell-output .cell-output-stdout}

```
multisession:
- args: function (..., workers = c(cgroups2.cpu.max = 6), envir = parent.frame())
- tweaked: TRUE
- call: plan("multisession", workers = n_cores)
```


:::

```{.r .cell-code .hidden}
#| label: params-computation
#| cache: false
# ggplot2 theme
theme_set(theme_minimal(base_size = 12))
```
:::

::: {.cell}

```{.r .cell-code .hidden}
#| label: params
bioproject <- "PRJNA547712"
project <- "kim2020_Hypoth-dev"
cb_fpr <- 0.001
low_cutoff_gene <- 500
high_cutoff_gene <- NULL
high_cutoff_gene <- 5000
low_cutoff_umis <- NULL
low_cutoff_umis <- -Inf
high_cutoff_umis <- 25000
high_cutoff_pc_mt <- 15
high_cutoff_pc_ribo <- 20
high_cutoff_pc_hb <- 0.1
high_cutoff_doublet_score <- 0.33
high_cutoff_complexity <- 0.85
connectivity_model <- "min_tree"
k <- 10
metric <- "euclidean"
signature <- 100
```
:::





## Load data from Kim DW et al 2020 (bioproject PRJNA547712)





::: {.cell}

```{.r .cell-code .hidden}
#| label: convert-to-seurat
srt.kim <- schard::h5ad2seurat(here(
  data_dir,
  "kim2020_combined.h5ad"
), use.raw = TRUE)

X_umap <- srt.kim@meta.data |>
  select(X, Y, Z) |>
  as.matrix()
colnames(X_umap) <- c("UMAP_1", "UMAP_2", "UMAP_3")
rownames(X_umap) <- colnames(srt.kim)
srt.kim[["umap"]] <- CreateDimReducObject(embeddings = X_umap, key = "umap_", assay = DefaultAssay(srt.kim))
srt.kim$Age %<>% forcats::fct(levels = c(
  "E10", "E11", "E12", "E13", "E14",
  "E15", "E16", "E17", "E18", "P0",
  "P2", "P4", "P8", "P10", "P14", "P23", "P45"
))

Idents(srt.kim) <- "Age"
srt.kim <- Store_Palette_Seurat(seurat_object = srt.kim, palette = rev(brewer.pal(n = 11, name = "Spectral")), palette_name = "expr_Colour_Pal")
```

::: {.cell-output .cell-output-stderr .hidden}

```
Seurat Object now contains the following items in @misc slot:
ℹ 'expr_Colour_Pal'
```


:::
:::

::: {.cell}

```{.r .cell-code .hidden}
#| label: load-seurat
print(srt.kim)
```

::: {.cell-output .cell-output-stdout}

```
An object of class Seurat 
27998 features across 128006 samples within 1 assay 
Active assay: RNA (27998 features, 0 variable features)
 2 layers present: counts, data
 1 dimensional reduction calculated: umap
```


:::

```{.r .cell-code .hidden}
#| label: load-seurat
srt.romanov.pub <- readRDS("/data/1_heteroAstrocytes/PRJNA548917/old/oldCCA_nae_srt.rds") # Load the Seurat object deposited to GEO by Romanov et al. (2020)
srt.romanov.pub <- UpdateSeuratObject(srt.romanov.pub)
```

::: {.cell-output .cell-output-stderr .hidden}

```
Validating object structure
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Updating object slots
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Ensuring keys are in the proper structure
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Updating matrix keys for DimReduc 'pca'
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Updating matrix keys for DimReduc 'tsne'
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Updating matrix keys for DimReduc 'umap'
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: Assay RNA changing from Assay to Assay
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: DimReduc pca changing from DimReduc to DimReduc
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: DimReduc tsne changing from DimReduc to DimReduc
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: DimReduc umap changing from DimReduc to DimReduc
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Ensuring keys are in the proper structure
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Ensuring feature names don't have underscores or pipes
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Updating slots in RNA
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Updating slots in pca
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Updating slots in tsne
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Setting tsne DimReduc to global
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Updating slots in umap
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Setting umap DimReduc to global
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Validating object structure for Assay 'RNA'
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Validating object structure for DimReduc 'pca'
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Validating object structure for DimReduc 'tsne'
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Validating object structure for DimReduc 'umap'
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Object representation is consistent with the most current Seurat version
```


:::

```{.r .cell-code .hidden}
#| label: load-seurat
srt.romanov.pub$Age <-
  Cells(srt.romanov.pub) |>
  str_split(pattern = ":", simplify = T) %>%
  .[, 1] %>%
  str_split_fixed(pattern = "_", n = 3) %>%
  .[, 3]
srt.romanov.pub$Age %<>% forcats::fct(levels = c(
  "E10", "E11", "E12", "E13", "E14",
  "E15", "E16", "E17", "E18", "P0",
  "P2", "3P2", "P4", "P8", "1P10",
  "P10", "P14", "P23", "P45"
))
srt.romanov.pub$Age %<>% fct_collapse(
  P2 = c("P2", "3P2"),
  P10 = c("1P10", "P10")
)
srt.romanov.pub <- Store_Palette_Seurat(
  seurat_object = srt.romanov.pub,
  palette = read_lines(here(data_dir, "colours_wtree.tsv")), palette_name = "wtree_Colour_Pal"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
Seurat Object now contains the following items in @misc slot:
ℹ 'vst.out' and 'wtree_Colour_Pal'
```


:::

```{.r .cell-code .hidden}
#| label: load-seurat
names(srt.romanov.pub@misc$wtree_Colour_Pal) <- 1:45
print(srt.romanov.pub)
```

::: {.cell-output .cell-output-stdout}

```
An object of class Seurat 
24340 features across 51199 samples within 1 assay 
Active assay: RNA (24340 features, 3500 variable features)
 3 layers present: counts, data, scale.data
 3 dimensional reductions calculated: pca, tsne, umap
```


:::

```{.r .cell-code .hidden}
#| label: load-seurat
glimpse(srt.romanov.pub@meta.data)
```

::: {.cell-output .cell-output-stdout}

```
Rows: 51,199
Columns: 20
$ nGene            <int> 1652, 782, 447, 1706, 1106, 894, 727, 734, 669, 617, …
$ nUMI             <dbl> 2787, 1090, 544, 2709, 1817, 1220, 995, 1036, 920, 86…
$ orig.ident       <fct> Hypothalamus, Hypothalamus, Hypothalamus, Hypothalamu…
$ res.0.2          <chr> "23", "23", "23", "23", "23", "23", "23", "23", "23",…
$ res.0.4          <chr> "34", "34", "34", "34", "34", "34", "34", "34", "34",…
$ res.0.8          <chr> "42", "42", "42", "42", "42", "42", "42", "42", "42",…
$ res.1.2          <chr> "47", "47", "47", "47", "47", "47", "47", "47", "47",…
$ res.2            <chr> "54", "54", "54", "54", "54", "54", "54", "54", "54",…
$ tree.ident       <int> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
$ pro_Inter        <chr> "41", "41", "41", "41", "41", "41", "41", "41", "41",…
$ pro_Enter        <chr> "41", "41", "41", "41", "41", "41", "41", "41", "41",…
$ tree_final       <fct> 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 1…
$ subtree          <fct> 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 4…
$ prim_walktrap    <fct> 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 3…
$ umi_per_gene     <dbl> 1.687046, 1.393862, 1.217002, 1.587925, 1.642857, 1.3…
$ log_umi_per_gene <dbl> 0.22712693, 0.14421974, 0.08529138, 0.20082998, 0.215…
$ nCount_RNA       <dbl> 2787, 1090, 544, 2709, 1817, 1220, 995, 1036, 920, 86…
$ nFeature_RNA     <int> 1652, 782, 447, 1706, 1106, 894, 727, 734, 669, 617, …
$ wtree            <fct> 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 3…
$ Age              <fct> P23, P2, P2, P2, P2, P2, P2, P2, P2, P2, P2, P2, P2, …
```


:::

```{.r .cell-code .hidden}
#| label: load-seurat
table(Idents(srt.romanov.pub))
```

::: {.cell-output .cell-output-stdout}

```

    1     2     3     4     5     6     7     8     9    10    11    12    13 
 2344  8146   395   402  3234   712   552   374   259   952 13727  1615   765 
   14    15    16    17    18    19    20    21    22    23    24    25    26 
  832  1244   792   590   808  2486  1683   628  1039  1750   292   394   547 
   27    28    29    30    31    32    33    34    35    36    37    38    39 
  391   407   507    93    81   402   143   701   222   353   324    73    78 
   40    41    42    43    44    45 
  328   190    73    37   179    55 
```


:::

```{.r .cell-code .hidden}
#| label: load-seurat
srt.romanov.pt <- subset(srt.romanov.pub, idents = c("38", "42", "45"))
print(srt.romanov.pt)
```

::: {.cell-output .cell-output-stdout}

```
An object of class Seurat 
24340 features across 201 samples within 1 assay 
Active assay: RNA (24340 features, 3500 variable features)
 3 layers present: counts, data, scale.data
 3 dimensional reductions calculated: pca, tsne, umap
```


:::
:::

::: {#cell-fig-feature-romanov2020 .cell}

```{.r .cell-code .hidden}
#| label: fig-feature-romanov2020
#| fig-cap: "Feature plot of selected genes in hypothalamus across different developmental stages in the Romanov et al. (2020) dataset (original UMAP embedding). Cells are colored by expression level.  Note the distinct localization patterns of each gene."
#| fig-width: 24
#| fig-height: 36
#| fig-dpi: 90
FeaturePlot(
  srt.romanov.pub,
  features = c(
    "Tshb", "Cck", "Pitx1",
    "Eya1", "Eya2", "Eya3", "Eya4",
    "Sox2", "Hlf", "Tshr",
    "Cckar", "Cckbr", "Gpr173",
    "Foxl2", "Lhx3", "Lhx4", "Pit1", "Gata2"
  ),
  label = F,
  blend = F,
  order = TRUE,
  pt.size = 1.2,
  raster.dpi = c(512, 512),
  alpha = 0.5,
  split.by = "Age"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: The following requested variables were not found: Pit1
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Foxl2"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Lhx3"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Foxl2"
All cells have the same value (0) of "Foxl2"
All cells have the same value (0) of "Foxl2"
All cells have the same value (0) of "Foxl2"
```


:::

::: {.cell-output-display}
![Feature plot of selected genes in hypothalamus across different developmental stages in the Romanov et al. (2020) dataset (original UMAP embedding). Cells are colored by expression level.  Note the distinct localization patterns of each gene.](01-de_test-focus_pars_tub_files/figure-html/fig-feature-romanov2020-1.png){#fig-feature-romanov2020 width=2160}
:::
:::

::: {#cell-fig-feature-pt-romanov2020 .cell}

```{.r .cell-code .hidden}
#| label: fig-feature-pt-romanov2020
#| fig.cap: "Feature plot of selected genes in hypothalamus across different developmental stages in the subset of the pars tuberalis clusters from the Romanov et al. (2020) dataset (original UMAP embedding). Cells are colored by a gene expression level.  Note the distinct localization patterns of each gene."
#| fig-width: 24
#| fig-height: 36
FeaturePlot(
  srt.romanov.pt,
  features = c(
    "Tshb", "Cck", "Pitx1",
    "Eya1", "Eya2", "Eya3", "Eya4",
    "Sox2", "Hlf", "Tshr",
    "Cckar", "Cckbr", "Gpr173",
    "Foxl2", "Lhx3", "Lhx4", "Pit1", "Gata2"
  ),
  label = T,
  blend = F,
  order = TRUE,
  pt.size = 1.2,
  raster.dpi = c(512, 512),
  alpha = 0.5,
  split.by = "Age"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: The following requested variables were not found: Pit1
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Tshb"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Cck"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Pitx1"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Eya1"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Eya3"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Hlf"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Tshr"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Cckar"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Cckbr"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Gpr173"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Foxl2"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Lhx3"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Lhx4"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Gata2"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Tshb"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Cck"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Pitx1"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Eya1"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Eya2"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Eya3"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Eya4"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Sox2"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Tshr"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Cckar"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Cckbr"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Gpr173"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Foxl2"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Lhx3"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Lhx4"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Gata2"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Cckar"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Cckbr"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Foxl2"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Eya2"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Eya3"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Foxl2"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Lhx4"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Cckar"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Cckbr"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Foxl2"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Lhx4"
```


:::

::: {.cell-output-display}
![Feature plot of selected genes in hypothalamus across different developmental stages in the subset of the pars tuberalis clusters from the Romanov et al. (2020) dataset (original UMAP embedding). Cells are colored by a gene expression level.  Note the distinct localization patterns of each gene.](01-de_test-focus_pars_tub_files/figure-html/fig-feature-pt-romanov2020-1.png){#fig-feature-pt-romanov2020 width=2304}
:::
:::

::: {.cell}

```{.r .cell-code .hidden}
#| label: norm-scale-matrix
# Using gene name patterns
srt.kim <- Add_Mito_Ribo(object = srt.kim, species = "Mouse", ensembl_ids = FALSE)
```

::: {.cell-output .cell-output-stderr .hidden}

```
Adding Percent Mitochondrial genes for Mouse using gene symbol pattern: "^mt-".
Adding Percent Ribosomal genes for Mouse using gene symbol pattern: "^Rp[sl]".
Adding Percent Mito+Ribo by adding Mito & Ribo percentages.
```


:::

```{.r .cell-code .hidden}
#| label: norm-scale-matrix
srt.kim$umi_per_gene <-
  (srt.kim$nCount_RNA/srt.kim$nFeature_RNA)
srt.kim$log_umi_per_gene <-
  log10(srt.kim$umi_per_gene)
srt.kim <- NormalizeData(srt.kim)
srt.kim <- FindVariableFeatures(srt.kim, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(srt.kim)
srt.kim <- ScaleData(srt.kim, features = all.genes, vars.to.regress = c("log_umi_per_gene", "percent_mito_ribo"))
```

::: {.cell-output .cell-output-stderr .hidden}

```
Regressing out log_umi_per_gene, percent_mito_ribo
Centering and scaling data matrix
```


:::

```{.r .cell-code .hidden}
#| label: norm-scale-matrix
# srt.kim <- ScaleData(srt.kim)
```
:::

::: {.cell}

```{.r .cell-code .hidden}
#| label: transfer-annotations
hypoth.anchors <- FindTransferAnchors(
  reference = srt.romanov.pub,
  query = srt.kim,
  dims = 1:30,
  npcs = 30,
  reference.reduction = "pca",
  reduction = "rpca",
  k.anchor = 15,
  k.filter = 30,
  k.score = 50,
  max.features = 500,
  nn.method = "annoy",
  n.trees = 200
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
Performing PCA on the provided query using 3500 features as input.
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Centering and scaling data matrix
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Projecting new data onto SVD
Projecting new data onto SVD
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Projecting cell embeddings
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Finding neighborhoods
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Finding anchors
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
	Found 6393 anchors
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Filtering anchors
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
	Retained 397 anchors
```


:::

```{.r .cell-code .hidden}
#| label: transfer-annotations
predictions <- TransferData(anchorset = hypoth.anchors, refdata = srt.romanov.pub$wtree, dims = 1:30) # TODO: check if I can randomly subset the reference data to handle the imbalance of clusters
```

::: {.cell-output .cell-output-stderr .hidden}

```
Finding integration vectors
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Finding integration vector weights
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Predicting cell labels
```


:::

```{.r .cell-code .hidden}
#| label: transfer-annotations
srt.kim <- AddMetaData(srt.kim, metadata = predictions)
table(srt.kim$predicted.id)
```

::: {.cell-output .cell-output-stdout}

```

    1    12    13    15    19     2    20    21    22    23    26    27    29 
 1504  2751  7692  5423 76464  1508  2952  6239   264   141    20  3053  1733 
    3    32    38     4    44    45     5     6     7     8 
 1580     3   264   124    38 10653  1897  3291   161   251 
```


:::
:::

::: {.cell}

```{.r .cell-code .hidden}
#| label: transfer-umap
srt.romanov.pub <- RunUMAP(srt.romanov.pub, dims = 1:30, reduction = "pca", return.model = TRUE)
```

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
This message will be shown once per session
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
UMAP will return its model
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
10:55:01 UMAP embedding parameters a = 0.9922 b = 1.112
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
10:55:01 Read 51199 rows and found 30 numeric columns
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
10:55:01 Using Annoy for neighbor search, n_neighbors = 30
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
10:55:01 Building Annoy index with metric = cosine, n_trees = 50
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
0%   10   20   30   40   50   60   70   80   90   100%
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
[----|----|----|----|----|----|----|----|----|----|
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
**************************************************|
10:55:10 Writing NN index file to temp file /tmp/Rtmpx1NRAQ/file268d1a7104f
10:55:10 Searching Annoy index using 6 threads, search_k = 3000
10:55:15 Annoy recall = 100%
10:55:16 Commencing smooth kNN distance calibration using 6 threads with target n_neighbors = 30
10:55:18 Initializing from normalized Laplacian + noise (using RSpectra)
10:55:19 Commencing optimization for 200 epochs, with 2530018 positive edges
10:55:41 Optimization finished
```


:::

```{.r .cell-code .hidden}
#| label: transfer-umap
srt.kim <- IntegrateEmbeddings(
  anchorset = hypoth.anchors, reference = srt.romanov.pub, query = srt.kim,
  new.reduction.name = "ref.pca"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
Requested to reuse weights matrix, but no weights found. Computing new weights.
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: Layer counts isn't present in the assay object; returning NULL
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: Layer counts isn't present in the assay object; returning NULL
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```

Integrating dataset 2 with reference dataset
Finding integration vectors
Finding integration vector weights
Integrating data
```


:::

```{.r .cell-code .hidden}
#| label: transfer-umap
srt.kim <- ProjectUMAP(
  query = srt.kim, query.reduction = "ref.pca", reference = srt.romanov.pub,
  reference.reduction = "pca", reduction.model = "umap"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
Computing nearest neighbors
Running UMAP projection
10:59:38 Read 128006 rows
10:59:38 Processing block 1 of 1
10:59:38 Commencing smooth kNN distance calibration using 6 threads with target n_neighbors = 30
10:59:38 Initializing by weighted average of neighbor coordinates using 6 threads
10:59:39 Commencing optimization for 67 epochs, with 3840180 positive edges
11:00:06 Finished
```


:::
:::

::: {#cell-fig-reference-umap-transfered .cell}

```{.r .cell-code .hidden}
#| label: fig-reference-umap-transfered
#| fig-cap: "UMAP plot of the Romanov et al. (2020) dataset with reference annotations and transferred labels on the the Kim DW et al. (2020) dataset. Cells are colored by published *wtree* clusters' labels identified as pars tuberalis. Note the successful embedding of query dataset according to the reference dataset and slight change of the reference embedding model in the process."
Idents(srt.romanov.pub) <- "wtree"
p1 <- Cluster_Highlight_Plot(
  seurat_object = srt.romanov.pub,
  reduction = "umap",
  label = TRUE,
  label.size = 3,
  repel = TRUE,
  cluster_name = c("38", "42", "45"),
  highlight_color = srt.romanov.pub@misc$wtree_Colour_Pal[c("38", "42", "45")]
) +
  NoLegend() +
  ggtitle("Reference annotations")

Idents(srt.kim) <- "predicted.id"
p2 <- Cluster_Highlight_Plot(
  seurat_object = srt.kim,
  reduction = "ref.umap",
  label = TRUE,
  label.size = 3,
  repel = TRUE,
  cluster_name = c("38", "42", "45"),
  highlight_color = srt.romanov.pub@misc$wtree_Colour_Pal[c("38", "42", "45")]
) +
  NoLegend() +
  ggtitle("Query transferred labels")
```

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: The following `cluster_name` were omitted as they were not found the
active.ident slot: 42
```


:::

```{.r .cell-code .hidden}
#| label: fig-reference-umap-transfered
#| fig-cap: "UMAP plot of the Romanov et al. (2020) dataset with reference annotations and transferred labels on the the Kim DW et al. (2020) dataset. Cells are colored by published *wtree* clusters' labels identified as pars tuberalis. Note the successful embedding of query dataset according to the reference dataset and slight change of the reference embedding model in the process."
p1 + p2
```

::: {.cell-output-display}
![UMAP plot of the Romanov et al. (2020) dataset with reference annotations and transferred labels on the the Kim DW et al. (2020) dataset. Cells are colored by published *wtree* clusters' labels identified as pars tuberalis. Note the successful embedding of query dataset according to the reference dataset and slight change of the reference embedding model in the process.](01-de_test-focus_pars_tub_files/figure-html/fig-reference-umap-transfered-1.png){#fig-reference-umap-transfered width=672}
:::
:::

::: {#cell-fig-feature-romanov2020-integrated .cell}

```{.r .cell-code .hidden}
#| label: fig-feature-romanov2020-integrated
#| fig.cap: "Feature plot of selected genes in hypothalamus across different developmental stages in the Romanov et al. (2020) dataset (new integrated UMAP embedding). Cells are colored by expression level.  Note the distinct localization patterns of each gene."
#| fig-width: 24
#| fig-height: 36
#| fig-dpi: 90
FeaturePlot(
  srt.romanov.pub,
  features = c(
    "Tshb", "Cck", "Pitx1",
    "Eya1", "Eya2", "Eya3", "Eya4",
    "Sox2", "Hlf", "Tshr",
    "Cckar", "Cckbr", "Gpr173",
    "Foxl2", "Lhx3", "Lhx4", "Pit1", "Gata2"
  ),
  label = F,
  blend = F,
  order = TRUE,
  pt.size = 1.2,
  raster.dpi = c(512, 512),
  alpha = 0.5,
  split.by = "Age"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: The following requested variables were not found: Pit1
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Foxl2"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Lhx3"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Foxl2"
All cells have the same value (0) of "Foxl2"
All cells have the same value (0) of "Foxl2"
All cells have the same value (0) of "Foxl2"
```


:::

::: {.cell-output-display}
![Feature plot of selected genes in hypothalamus across different developmental stages in the Romanov et al. (2020) dataset (new integrated UMAP embedding). Cells are colored by expression level.  Note the distinct localization patterns of each gene.](01-de_test-focus_pars_tub_files/figure-html/fig-feature-romanov2020-integrated-1.png){#fig-feature-romanov2020-integrated width=2160}
:::
:::





# Calculate and plot chi2 test of independence between Sox2 and Tshr expression in hypothalamus across different developmental stages





::: {.cell}

```{.r .cell-code .hidden}
#| label: get-goi-sox2-tshr
sbs_mtx <- GetAssayData(object = srt.kim, layer = "counts", assay = "RNA")[c("Sox2", "Tshr"), ] %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  select(Sox2, Tshr) %>%
  dplyr::bind_cols(srt.kim@meta.data) %>%
  select(Age, Sox2, Tshr) %>%
  mutate(
    Sox2_pos = Sox2 > 0,
    Tshr_pos = Tshr > 0
  )

sbs_mtx %>% skimr::skim()
```

::: {.cell-output-display}

Table: Data summary

|                         |           |
|:------------------------|:----------|
|Name                     |Piped data |
|Number of rows           |128006     |
|Number of columns        |5          |
|_______________________  |           |
|Column type frequency:   |           |
|factor                   |1          |
|logical                  |2          |
|numeric                  |2          |
|________________________ |           |
|Group variables          |None       |


**Variable type: factor**

|skim_variable | n_missing| complete_rate|ordered | n_unique|top_counts                                    |
|:-------------|---------:|-------------:|:-------|--------:|:---------------------------------------------|
|Age           |         0|             1|FALSE   |       12|P45: 17025, E16: 16781, E14: 16543, P8: 13677 |


**Variable type: logical**

|skim_variable | n_missing| complete_rate| mean|count                   |
|:-------------|---------:|-------------:|----:|:-----------------------|
|Sox2_pos      |         0|             1| 0.14|FAL: 110343, TRU: 17663 |
|Tshr_pos      |         0|             1| 0.00|FAL: 127714, TRU: 292   |


**Variable type: numeric**

|skim_variable | n_missing| complete_rate| mean|   sd| p0| p25| p50| p75| p100|hist  |
|:-------------|---------:|-------------:|----:|----:|--:|---:|---:|---:|----:|:-----|
|Sox2          |         0|             1| 0.27| 1.21|  0|   0|   0|   0|   82|▇▁▁▁▁ |
|Tshr          |         0|             1| 0.00| 0.06|  0|   0|   0|   0|    5|▇▁▁▁▁ |


:::
:::

::: {#cell-fig-sox2-tshr-stats .cell}

```{.r .cell-code .hidden}
#| label: fig-sox2-tshr-stats
#| fig-cap: "Proportion of Sox2 and Tshr expression in hypothalamus across different developmental stages. Cells are colored by Sox2 and Tshr expression status. Note the distinct Sox2 and Tshr expression patterns across different developmental stages."
#| fig-width: 8
#| fig-height: 24
write_csv(sbs_mtx, here(tables_dir, "Sox2-Tshr-expression-status-between-Ages-on-evaluation-datasets.csv"))


# plot
grouped_ggpiestats(
  data = sbs_mtx,
  x = Tshr_pos,
  y = Sox2_pos,
  grouping.var = Age,
  perc.k = 1,
  package = "ggsci",
  palette = "category10_d3",
  title.text = "Sox2 specification of Tshr-positive hypothalamic development",
  caption.text = "Asterisks denote results from proportion tests; \n***: p < 0.001, ns: non-significant",
  plotgrid.args = list(nrow = 8)
)
```

::: {.cell-output-display}
![Proportion of Sox2 and Tshr expression in hypothalamus across different developmental stages. Cells are colored by Sox2 and Tshr expression status. Note the distinct Sox2 and Tshr expression patterns across different developmental stages.](01-de_test-focus_pars_tub_files/figure-html/fig-sox2-tshr-stats-1.png){#fig-sox2-tshr-stats width=768}
:::
:::





## Calculate and plot hexagonal cells representation in hypothalamus across different developmental stages with meta information





::: {.cell}

```{.r .cell-code .hidden}
srt.kim
```

::: {.cell-output .cell-output-stdout}

```
An object of class Seurat 
27998 features across 128006 samples within 1 assay 
Active assay: RNA (27998 features, 5000 variable features)
 3 layers present: counts, data, scale.data
 3 dimensional reductions calculated: umap, ref.pca, ref.umap
```


:::
:::

::: {#cell-fig-hexbin-umap-kim2020 .cell}

```{.r .cell-code .hidden}
#| label: fig-hexbin-umap-kim2020
#| fig-cap: "Hexagonal binning plot of the Kim DW et al. (2020) dataset in the hypothalamus of all developmental stages. Cells are colored by cell count."
library(hexbin)
# Extract UMAP coordinates
umap_coords <- Embeddings(srt.kim, reduction = "ref.umap")

# Create hexbin object
hb <- hexbin(umap_coords[, 1], umap_coords[, 2], xbins = 64)

# Create a data frame for plotting
hex_data <- data.frame(
  x = hcell2xy(hb)$x,
  y = hcell2xy(hb)$y,
  count = hb@count
)

# Create the plot
ggplot(hex_data, aes(x = x, y = y, fill = count)) +
  geom_hex(stat = "identity") +
  scale_fill_gradientn(colors = ggsci::pal_material("amber")(9)) +
  labs(x = "UMAP_1", y = "UMAP_2", fill = "Cell\nCount") +
  coord_fixed()
```

::: {.cell-output-display}
![Hexagonal binning plot of the Kim DW et al. (2020) dataset in the hypothalamus of all developmental stages. Cells are colored by cell count.](01-de_test-focus_pars_tub_files/figure-html/fig-hexbin-umap-kim2020-1.png){#fig-hexbin-umap-kim2020 width=672}
:::
:::

::: {#cell-fig-feature-kim2020 .cell}

```{.r .cell-code .hidden}
#| label: fig-feature-kim2020
#| fig-cap: "Feature plot of selected genes in hypothalamus across different developmental stages in the Kim DW et al. (2020) dataset (integrated UMAP embedding). Cells are colored by expression level.  Note the distinct localization patterns of each gene."
#| fig-width: 32
#| fig-height: 36
FeaturePlot(
  srt.kim,
  features = c(
    "Tshb", "Cck", "Pitx1",
    "Eya1", "Eya2", "Eya3", "Eya4",
    "Sox2", "Hlf", "Tshr",
    "Cckar", "Cckbr", "Gpr173",
    "Foxl2", "Lhx3", "Lhx4", "Pit1", "Gata2"
  ),
  reduction = "ref.umap",
  label = F,
  blend = F,
  order = TRUE,
  pt.size = 1.2,
  raster.dpi = c(512, 512),
  alpha = 0.5,
  split.by = "Age"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: The following requested variables were not found: Pit1
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Tshb"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Cck"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Cckar"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Cckbr"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Foxl2"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Tshb"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Cckbr"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Foxl2"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Pitx1"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Foxl2"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Lhx4"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Cckbr"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Foxl2"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Lhx3"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Lhx4"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Tshr"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Foxl2"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Lhx3"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Tshb"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Pitx1"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Foxl2"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Lhx3"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Pitx1"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Foxl2"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Lhx3"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Pitx1"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Foxl2"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Lhx3"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Lhx4"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Foxl2"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Lhx4"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Foxl2"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Lhx4"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Foxl2"
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: All cells have the same value (0) of "Lhx4"
```


:::

::: {.cell-output-display}
![Feature plot of selected genes in hypothalamus across different developmental stages in the Kim DW et al. (2020) dataset (integrated UMAP embedding). Cells are colored by expression level.  Note the distinct localization patterns of each gene.](01-de_test-focus_pars_tub_files/figure-html/fig-feature-kim2020-1.png){#fig-feature-kim2020 width=3072}
:::
:::

::: {#cell-fig-prediction-scores-kim2020 .cell}

```{.r .cell-code .hidden}
#| label: fig-prediction-scores-kim2020
#| fig-cap: "Feature plot of prediction scores in hypothalamus of all developmental stages in the Kim DW et al. (2020) dataset (integrated UMAP embedding). Cells are colored by prediction score of pars tuberalis identified clusters. Note the distinct localization patterns of prediction score."
#| fig-width: 21
#| fig-height: 5
Idents(srt.kim) <- "predicted.id"
FeaturePlot(
  srt.kim,
  features = c(
    "prediction.score.38",
    "prediction.score.42",
    "prediction.score.45"
  ),
  reduction = "ref.umap",
  label = T,
  repel = T,
  blend = F,
  order = TRUE,
  pt.size = 4,
  raster.dpi = c(512, 512),
  alpha = 0.8,
  max.cutoff = "q90",
  ncol = 3
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
Rasterizing points since number of points exceeds 100,000.
To disable this behavior set `raster=FALSE`
Rasterizing points since number of points exceeds 100,000.
To disable this behavior set `raster=FALSE`
Rasterizing points since number of points exceeds 100,000.
To disable this behavior set `raster=FALSE`
```


:::

::: {.cell-output-display}
![Feature plot of prediction scores in hypothalamus of all developmental stages in the Kim DW et al. (2020) dataset (integrated UMAP embedding). Cells are colored by prediction score of pars tuberalis identified clusters. Note the distinct localization patterns of prediction score.](01-de_test-focus_pars_tub_files/figure-html/fig-prediction-scores-kim2020-1.png){#fig-prediction-scores-kim2020 width=2016}
:::
:::

::: {#cell-fig-prediction-scores-split-kim2020 .cell}

```{.r .cell-code .hidden}
#| label: fig-prediction-scores-split-kim2020
#| fig-cap: "Feature plot of prediction scores in hypothalamus across different developmental stages in the Kim DW et al. (2020) dataset (integrated UMAP embedding). Cells are colored by prediction score of pars tuberalis identified clusters. Note the distinct localization patterns of prediction score."
#| fig-width: 32
#| fig-height: 6.353
FeaturePlot(
  srt.kim,
  features = c(
    "prediction.score.38",
    "prediction.score.42",
    "prediction.score.45"
  ),
  reduction = "ref.umap",
  label = F,
  blend = F,
  order = TRUE,
  pt.size = 1.2,
  raster.dpi = c(512, 512),
  alpha = 0.8,
  max.cutoff = "q90",
  split.by = "Age"
)
```

::: {.cell-output-display}
![Feature plot of prediction scores in hypothalamus across different developmental stages in the Kim DW et al. (2020) dataset (integrated UMAP embedding). Cells are colored by prediction score of pars tuberalis identified clusters. Note the distinct localization patterns of prediction score.](01-de_test-focus_pars_tub_files/figure-html/fig-prediction-scores-split-kim2020-1.png){#fig-prediction-scores-split-kim2020 width=3072}
:::
:::

::: {#cell-fig-prediction-id-kim2020 .cell}

```{.r .cell-code .hidden}
#| label: fig-prediction-id-kim2020
#| fig-cap: "Feature plot of predicted clusters in hypothalamus of all developmental stages in the Kim DW et al. (2020) dataset (integrated UMAP embedding). Cells are colored by predicted clusters. Note the distinct localization patterns of clusters 38 and 45."
#| fig-width: 7
#| fig-height: 5
DimPlot(
  srt.kim,
  group.by = c("predicted.id"),
  reduction = "ref.umap",
  label = T,
  repel = T,
  pt.size = 4,
  raster.dpi = c(300, 300),
  alpha = 0.5
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
Rasterizing points since number of points exceeds 100,000.
To disable this behavior set `raster=FALSE`
```


:::

::: {.cell-output-display}
![Feature plot of predicted clusters in hypothalamus of all developmental stages in the Kim DW et al. (2020) dataset (integrated UMAP embedding). Cells are colored by predicted clusters. Note the distinct localization patterns of clusters 38 and 45.](01-de_test-focus_pars_tub_files/figure-html/fig-prediction-id-kim2020-1.png){#fig-prediction-id-kim2020 width=672}
:::
:::

::: {#cell-fig-prediction-split-kim2020 .cell}

```{.r .cell-code .hidden}
#| label: fig-prediction-split-kim2020
#| fig-cap: "Feature plot of projected Pars Tuberalis clusters in hypothalamus across different developmental stages in the Kim DW et al. (2020) dataset (integrated UMAP embedding). Cells are colored by original colour code of wtree clusters. Note the distinct localization patterns of clusters 38 and 45."
#| fig-width: 14
#| fig-height: 5
p38kim <- Cluster_Highlight_Plot(
  seurat_object = srt.kim,
  reduction = "ref.umap",
  label = TRUE,
  label.size = 1.2,
  repel = TRUE,
  # split.by = "Age",
  cluster_name = "38",
  highlight_color = srt.romanov.pub@misc$wtree_Colour_Pal["38"]
) + NoLegend()

# Not found
# p42kim <- Cluster_Highlight_Plot(
#   seurat_object = srt.kim,
#   reduction = "ref.umap",
#   label = TRUE,
#   label.size = 1.2,
#   repel = TRUE,
#   # split.by = "Age",
#   cluster_name = "42",
#   highlight_color = srt.romanov.pub@misc$wtree_Colour_Pal["42"]
# ) + NoLegend()

p45kim <- Cluster_Highlight_Plot(
  seurat_object = srt.kim,
  reduction = "ref.umap",
  label = TRUE,
  label.size = 1.2,
  repel = TRUE,
  # split.by = "Age",
  cluster_name = "45",
  highlight_color = srt.romanov.pub@misc$wtree_Colour_Pal["45"]
) + NoLegend()

# (p38kim | p42kim | p45kim)
p38kim | p45kim
```

::: {.cell-output-display}
![Feature plot of projected Pars Tuberalis clusters in hypothalamus across different developmental stages in the Kim DW et al. (2020) dataset (integrated UMAP embedding). Cells are colored by original colour code of wtree clusters. Note the distinct localization patterns of clusters 38 and 45.](01-de_test-focus_pars_tub_files/figure-html/fig-prediction-split-kim2020-1.png){#fig-prediction-split-kim2020 width=1344}
:::
:::





## Quantify and plot representation of gene interactions with feature expression Spearman’s correlation

### Kim DW et al. 2020 dataset





::: {#cell-fig-violin-gene-interactions-kim2020 .cell}

```{.r .cell-code .hidden}
#| label: fig-violin-gene-interactions-kim2020
#| fig-cap: Gene expression of projected Pars Tuberalis clusters in hypothalamus across different developmental stages in the Kim DW et al. (2020) dataset.
#| fig-width: 12
#| fig-height: 12
#| warning: false
cells_to_check <- WhichCells(srt.kim, idents = c("38", "42", "45") %>% .[. %in% levels(srt.kim)])
dat_kim <- srt.kim@assays$RNA@counts[, cells_to_check]
genes_to_check <- c(
  "Tshb",
  "Cck",
  "Pitx1",
  "Eya1",
  "Eya2",
  "Eya3",
  "Eya4",
  "Sox2",
  "Hlf",
  "Tshr",
  "Cckar",
  "Cckbr",
  "Gpr173",
  "Foxl2",
  "Lhx3",
  "Lhx4",
  "Pit1",
  "Gata2"
) %>%
  .[. %in% (srt.kim@assays$RNA@counts[, cells_to_check] |> rowSums() %>% .[. > 3] %>% names())]

dat_kim <- dat_kim[genes_to_check, cells_to_check] %>% as.data.frame() %>% t()

dat_kim <- dat_kim |> as.data.frame() |> mutate(Age = srt.kim[["Age"]][rownames(dat_kim), ])

# 1. Complete gene pairs - Ensure all combinations are present
gene_pairs <- combn(
  genes_to_check, 2
) |>
  as_tibble() |>
  array_tree(margin = 2)
```

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: The `x` argument of `as_tibble.matrix()` must have unique column names if
`.name_repair` is omitted as of tibble 2.0.0.
ℹ Using compatibility `.name_repair`.
```


:::

```{.r .cell-code .hidden}
#| label: fig-violin-gene-interactions-kim2020
#| fig-cap: Gene expression of projected Pars Tuberalis clusters in hypothalamus across different developmental stages in the Kim DW et al. (2020) dataset.
#| fig-width: 12
#| fig-height: 12
#| warning: false
ages <- c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P8", "P45")


# Violin plots of individual gene expression split by cluster
Idents(srt.kim) <- "Age" # Reset identity back to Age

Stacked_VlnPlot(subset(srt.kim, cells = cells_to_check), features = genes_to_check, split.by = "predicted.id", colors_use = srt.romanov.pub@misc$wtree_Colour_Pal[c("38", "42", "45")])
```

::: {.cell-output .cell-output-stderr .hidden}

```
The default behaviour of split.by has changed.
Separate violin plots are now plotted side-by-side.
To restore the old behaviour of a single split violin,
set split.plot = TRUE.
      
This message will be shown once per session.
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
```


:::

::: {.cell-output-display}
![Gene expression of projected Pars Tuberalis clusters in hypothalamus across different developmental stages in the Kim DW et al. (2020) dataset.](01-de_test-focus_pars_tub_files/figure-html/fig-violin-gene-interactions-kim2020-1.png){#fig-violin-gene-interactions-kim2020 width=1152}
:::
:::

::: {.cell .hidden}

```{.r .cell-code .hidden}
#| label: setup-gene-correlation-analysis
#| include: false

# Helper functions
calculate_correlations <- function(seurat_obj, gene_pairs, group_var, cells) {
  map_dfr(gene_pairs, function(genes) {
    data <- seurat_obj |>
      subset(cells = cells) |>
      GetAssayData(layer = "data") |>
      t() |>
      as.data.frame() |>
      select(all_of(genes))

    groups <- seurat_obj[[group_var]][cells, , drop = FALSE]

    # Calculate correlations for each group
    map_dfr(unique(groups[[1]]), function(group) {
      group_cells <- groups[[1]] == group
      group_data <- data[group_cells, ]

      # Spearman correlation test
      test_result <- pspearman::spearman.test(
        group_data[[genes[1]]],
        group_data[[genes[2]]]
      ) |>
        broom::tidy()

      tibble(
        group = group,
        gene1 = genes[1],
        gene2 = genes[2],
        correlation = test_result$estimate,
        p.value = test_result$p.value,
        n = sum(group_cells)
      )
    })
  })
}

generate_correlation_plots <- function(seurat_obj, genes, group_var, cells, age = NULL) {
  plot_data <- seurat_obj |>
    subset(cells = cells) |>
    GetAssayData(layer = "data") |>
    t() |>
    as.data.frame() |>
    dplyr::select(genes) |>
    mutate(group = as_factor(seurat_obj[[group_var]][colnames(seurat_obj), 1]))

  grouped_ggscatterstats(
    data = plot_data,
    x = !!genes[1],  # Force evaluation of gene names as symbols
    y = !!genes[2],  # Force evaluation of gene names as symbols
    grouping.var = group,
    type = "nonparametric",
    xlab = genes[1],
    ylab = genes[2]
  )
}
```
:::

::: {.cell}

````{.python .cell-code .hidden}
# | label: generate-kim2020-correlation-plots
# | output: false
# | echo: false

import pandas as pd

# Access R variables in Python using py$
gene_pairs = r.gene_pairs
ages = r.ages
dat_kim = r.dat_kim

# Convert gene_pairs and data to a Pandas DataFrame
gene_pairs = pd.DataFrame(gene_pairs).T
dat_kim = pd.DataFrame(dat_kim)


def generate_plot_chunk(genes, age):
    # Filter the data for the current age
    age_data = dat_kim[dat_kim["Age"] == age]

    # Check if both genes have at least two unique values in the age subset
    gene1_values = age_data[genes[0]].unique()
    gene2_values = age_data[genes[1]].unique()

    if len(gene1_values) >= 2 and len(gene2_values) >= 2:
        # Generate the plot chunk
        chunk_text = f"""```{{r}}
#| label: fig-kim2020-correlation-{genes[0]}-{genes[1]}-{age}
#| fig-cap: "Correlation between {genes[0]} and {genes[1]} expression at {age} in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "{age}", cells = cells_to_check),
    c("{genes[0]}", "{genes[1]}"),
    "predicted.id",
    cells_to_check,
    age = "{age}"
)
```"""
        return chunk_text
    else:
        # Skip this gene pair for this age if one or both genes don't have enough values
        return ""


# Generate the Quarto markdown content
qmd_content = "#### Gene Correlation Plots\n\n"

for _, row in gene_pairs.iterrows():
    for age in ages:
        qmd_content += generate_plot_chunk([row[0], row[1]], age) + "\n\n"

# Write the content to a .qmd file
with open("gene_correlation_plots_kim2020.qmd", "w") as f:
    f.write(qmd_content)
````

::: {.cell-output .cell-output-stdout}

```
166347
```


:::
:::


#### Gene Correlation Plots

















::: {#cell-fig-kim2020-correlation-Tshb-Cck-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshb-Cck-P8
#| fig-cap: "Correlation between Tshb and Cck expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Tshb", "Cck"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: Using an external vector in selections was deprecated in tidyselect 1.1.0.
ℹ Please use `all_of()` or `any_of()` instead.
  # Was:
  data %>% select(genes)

  # Now:
  data %>% select(all_of(genes))

See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
Registered S3 method overwritten by 'ggside':
  method from   
  +.gg   ggplot2
```


:::

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Cck expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshb-Cck-P8-1.png){#fig-kim2020-correlation-Tshb-Cck-P8 width=1152}
:::
:::



















::: {#cell-fig-kim2020-correlation-Tshb-Pitx1-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshb-Pitx1-P8
#| fig-cap: "Correlation between Tshb and Pitx1 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Tshb", "Pitx1"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Pitx1 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshb-Pitx1-P8-1.png){#fig-kim2020-correlation-Tshb-Pitx1-P8 width=1152}
:::
:::



















::: {#cell-fig-kim2020-correlation-Tshb-Eya1-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshb-Eya1-P8
#| fig-cap: "Correlation between Tshb and Eya1 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Tshb", "Eya1"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Eya1 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshb-Eya1-P8-1.png){#fig-kim2020-correlation-Tshb-Eya1-P8 width=1152}
:::
:::



















::: {#cell-fig-kim2020-correlation-Tshb-Eya2-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshb-Eya2-P8
#| fig-cap: "Correlation between Tshb and Eya2 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Tshb", "Eya2"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Eya2 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshb-Eya2-P8-1.png){#fig-kim2020-correlation-Tshb-Eya2-P8 width=1152}
:::
:::



















::: {#cell-fig-kim2020-correlation-Tshb-Eya3-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshb-Eya3-P8
#| fig-cap: "Correlation between Tshb and Eya3 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Tshb", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Eya3 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshb-Eya3-P8-1.png){#fig-kim2020-correlation-Tshb-Eya3-P8 width=1152}
:::
:::



















::: {#cell-fig-kim2020-correlation-Tshb-Eya4-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshb-Eya4-P8
#| fig-cap: "Correlation between Tshb and Eya4 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Tshb", "Eya4"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Eya4 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshb-Eya4-P8-1.png){#fig-kim2020-correlation-Tshb-Eya4-P8 width=1152}
:::
:::



















::: {#cell-fig-kim2020-correlation-Tshb-Sox2-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshb-Sox2-P8
#| fig-cap: "Correlation between Tshb and Sox2 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Tshb", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Sox2 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshb-Sox2-P8-1.png){#fig-kim2020-correlation-Tshb-Sox2-P8 width=1152}
:::
:::



















::: {#cell-fig-kim2020-correlation-Tshb-Hlf-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshb-Hlf-P8
#| fig-cap: "Correlation between Tshb and Hlf expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Tshb", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Hlf expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshb-Hlf-P8-1.png){#fig-kim2020-correlation-Tshb-Hlf-P8 width=1152}
:::
:::



















::: {#cell-fig-kim2020-correlation-Tshb-Tshr-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshb-Tshr-P8
#| fig-cap: "Correlation between Tshb and Tshr expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Tshb", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Tshr expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshb-Tshr-P8-1.png){#fig-kim2020-correlation-Tshb-Tshr-P8 width=1152}
:::
:::







































::: {#cell-fig-kim2020-correlation-Tshb-Cckbr-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshb-Cckbr-P8
#| fig-cap: "Correlation between Tshb and Cckbr expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Tshb", "Cckbr"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Cckbr expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshb-Cckbr-P8-1.png){#fig-kim2020-correlation-Tshb-Cckbr-P8 width=1152}
:::
:::



















::: {#cell-fig-kim2020-correlation-Tshb-Gpr173-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshb-Gpr173-P8
#| fig-cap: "Correlation between Tshb and Gpr173 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Tshb", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Gpr173 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshb-Gpr173-P8-1.png){#fig-kim2020-correlation-Tshb-Gpr173-P8 width=1152}
:::
:::



















::: {#cell-fig-kim2020-correlation-Tshb-Lhx3-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshb-Lhx3-P8
#| fig-cap: "Correlation between Tshb and Lhx3 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Tshb", "Lhx3"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Lhx3 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshb-Lhx3-P8-1.png){#fig-kim2020-correlation-Tshb-Lhx3-P8 width=1152}
:::
:::







































::: {#cell-fig-kim2020-correlation-Tshb-Gata2-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshb-Gata2-P8
#| fig-cap: "Correlation between Tshb and Gata2 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Tshb", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Gata2 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshb-Gata2-P8-1.png){#fig-kim2020-correlation-Tshb-Gata2-P8 width=1152}
:::
:::





::: {#cell-fig-kim2020-correlation-Cck-Pitx1-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Pitx1-E11
#| fig-cap: "Correlation between Cck and Pitx1 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Cck", "Pitx1"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Pitx1 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Pitx1-E11-1.png){#fig-kim2020-correlation-Cck-Pitx1-E11 width=1152}
:::
:::













::: {#cell-fig-kim2020-correlation-Cck-Pitx1-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Pitx1-P8
#| fig-cap: "Correlation between Cck and Pitx1 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Cck", "Pitx1"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Pitx1 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Pitx1-P8-1.png){#fig-kim2020-correlation-Cck-Pitx1-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Pitx1-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Pitx1-P45
#| fig-cap: "Correlation between Cck and Pitx1 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Cck", "Pitx1"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Pitx1 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Pitx1-P45-1.png){#fig-kim2020-correlation-Cck-Pitx1-P45 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Cck-Eya1-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya1-E11
#| fig-cap: "Correlation between Cck and Eya1 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Cck", "Eya1"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya1 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya1-E11-1.png){#fig-kim2020-correlation-Cck-Eya1-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Eya1-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya1-E12
#| fig-cap: "Correlation between Cck and Eya1 expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Cck", "Eya1"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya1 expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya1-E12-1.png){#fig-kim2020-correlation-Cck-Eya1-E12 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Eya1-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya1-E13
#| fig-cap: "Correlation between Cck and Eya1 expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Cck", "Eya1"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya1 expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya1-E13-1.png){#fig-kim2020-correlation-Cck-Eya1-E13 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Eya1-E14 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya1-E14
#| fig-cap: "Correlation between Cck and Eya1 expression at E14 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E14", cells = cells_to_check),
    c("Cck", "Eya1"),
    "predicted.id",
    cells_to_check,
    age = "E14"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya1 expression at E14 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya1-E14-1.png){#fig-kim2020-correlation-Cck-Eya1-E14 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Eya1-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya1-E15
#| fig-cap: "Correlation between Cck and Eya1 expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Cck", "Eya1"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya1 expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya1-E15-1.png){#fig-kim2020-correlation-Cck-Eya1-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Eya1-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya1-E16
#| fig-cap: "Correlation between Cck and Eya1 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Cck", "Eya1"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya1 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya1-E16-1.png){#fig-kim2020-correlation-Cck-Eya1-E16 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Eya1-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya1-E18
#| fig-cap: "Correlation between Cck and Eya1 expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Cck", "Eya1"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya1 expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya1-E18-1.png){#fig-kim2020-correlation-Cck-Eya1-E18 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Eya1-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya1-P8
#| fig-cap: "Correlation between Cck and Eya1 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Cck", "Eya1"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya1 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya1-P8-1.png){#fig-kim2020-correlation-Cck-Eya1-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Eya1-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya1-P45
#| fig-cap: "Correlation between Cck and Eya1 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Cck", "Eya1"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya1 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya1-P45-1.png){#fig-kim2020-correlation-Cck-Eya1-P45 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Cck-Eya2-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya2-E11
#| fig-cap: "Correlation between Cck and Eya2 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Cck", "Eya2"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya2 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya2-E11-1.png){#fig-kim2020-correlation-Cck-Eya2-E11 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Cck-Eya2-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya2-E13
#| fig-cap: "Correlation between Cck and Eya2 expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Cck", "Eya2"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya2 expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya2-E13-1.png){#fig-kim2020-correlation-Cck-Eya2-E13 width=1152}
:::
:::









::: {#cell-fig-kim2020-correlation-Cck-Eya2-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya2-P8
#| fig-cap: "Correlation between Cck and Eya2 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Cck", "Eya2"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya2 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya2-P8-1.png){#fig-kim2020-correlation-Cck-Eya2-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Eya2-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya2-P45
#| fig-cap: "Correlation between Cck and Eya2 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Cck", "Eya2"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya2 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya2-P45-1.png){#fig-kim2020-correlation-Cck-Eya2-P45 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Cck-Eya3-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya3-E11
#| fig-cap: "Correlation between Cck and Eya3 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Cck", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya3 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya3-E11-1.png){#fig-kim2020-correlation-Cck-Eya3-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Eya3-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya3-E12
#| fig-cap: "Correlation between Cck and Eya3 expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Cck", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya3 expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya3-E12-1.png){#fig-kim2020-correlation-Cck-Eya3-E12 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Eya3-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya3-E13
#| fig-cap: "Correlation between Cck and Eya3 expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Cck", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya3 expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya3-E13-1.png){#fig-kim2020-correlation-Cck-Eya3-E13 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Eya3-E14 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya3-E14
#| fig-cap: "Correlation between Cck and Eya3 expression at E14 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E14", cells = cells_to_check),
    c("Cck", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "E14"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya3 expression at E14 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya3-E14-1.png){#fig-kim2020-correlation-Cck-Eya3-E14 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Eya3-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya3-E15
#| fig-cap: "Correlation between Cck and Eya3 expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Cck", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya3 expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya3-E15-1.png){#fig-kim2020-correlation-Cck-Eya3-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Eya3-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya3-E16
#| fig-cap: "Correlation between Cck and Eya3 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Cck", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya3 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya3-E16-1.png){#fig-kim2020-correlation-Cck-Eya3-E16 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Eya3-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya3-E18
#| fig-cap: "Correlation between Cck and Eya3 expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Cck", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya3 expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya3-E18-1.png){#fig-kim2020-correlation-Cck-Eya3-E18 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Eya3-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya3-P8
#| fig-cap: "Correlation between Cck and Eya3 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Cck", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya3 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya3-P8-1.png){#fig-kim2020-correlation-Cck-Eya3-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Eya3-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya3-P45
#| fig-cap: "Correlation between Cck and Eya3 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Cck", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya3 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya3-P45-1.png){#fig-kim2020-correlation-Cck-Eya3-P45 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Cck-Eya4-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya4-E11
#| fig-cap: "Correlation between Cck and Eya4 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Cck", "Eya4"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya4 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya4-E11-1.png){#fig-kim2020-correlation-Cck-Eya4-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Eya4-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya4-E12
#| fig-cap: "Correlation between Cck and Eya4 expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Cck", "Eya4"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya4 expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya4-E12-1.png){#fig-kim2020-correlation-Cck-Eya4-E12 width=1152}
:::
:::





::: {#cell-fig-kim2020-correlation-Cck-Eya4-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya4-E15
#| fig-cap: "Correlation between Cck and Eya4 expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Cck", "Eya4"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya4 expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya4-E15-1.png){#fig-kim2020-correlation-Cck-Eya4-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Eya4-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya4-E16
#| fig-cap: "Correlation between Cck and Eya4 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Cck", "Eya4"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya4 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya4-E16-1.png){#fig-kim2020-correlation-Cck-Eya4-E16 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Cck-Eya4-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya4-P8
#| fig-cap: "Correlation between Cck and Eya4 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Cck", "Eya4"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya4 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya4-P8-1.png){#fig-kim2020-correlation-Cck-Eya4-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Eya4-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Eya4-P45
#| fig-cap: "Correlation between Cck and Eya4 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Cck", "Eya4"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya4 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Eya4-P45-1.png){#fig-kim2020-correlation-Cck-Eya4-P45 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Cck-Sox2-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Sox2-E11
#| fig-cap: "Correlation between Cck and Sox2 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Cck", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Sox2 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Sox2-E11-1.png){#fig-kim2020-correlation-Cck-Sox2-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Sox2-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Sox2-E12
#| fig-cap: "Correlation between Cck and Sox2 expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Cck", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Sox2 expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Sox2-E12-1.png){#fig-kim2020-correlation-Cck-Sox2-E12 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Sox2-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Sox2-E13
#| fig-cap: "Correlation between Cck and Sox2 expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Cck", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Sox2 expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Sox2-E13-1.png){#fig-kim2020-correlation-Cck-Sox2-E13 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Sox2-E14 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Sox2-E14
#| fig-cap: "Correlation between Cck and Sox2 expression at E14 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E14", cells = cells_to_check),
    c("Cck", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E14"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Sox2 expression at E14 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Sox2-E14-1.png){#fig-kim2020-correlation-Cck-Sox2-E14 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Sox2-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Sox2-E15
#| fig-cap: "Correlation between Cck and Sox2 expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Cck", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Sox2 expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Sox2-E15-1.png){#fig-kim2020-correlation-Cck-Sox2-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Sox2-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Sox2-E16
#| fig-cap: "Correlation between Cck and Sox2 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Cck", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Sox2 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Sox2-E16-1.png){#fig-kim2020-correlation-Cck-Sox2-E16 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Sox2-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Sox2-E18
#| fig-cap: "Correlation between Cck and Sox2 expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Cck", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Sox2 expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Sox2-E18-1.png){#fig-kim2020-correlation-Cck-Sox2-E18 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Sox2-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Sox2-P8
#| fig-cap: "Correlation between Cck and Sox2 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Cck", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Sox2 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Sox2-P8-1.png){#fig-kim2020-correlation-Cck-Sox2-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Sox2-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Sox2-P45
#| fig-cap: "Correlation between Cck and Sox2 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Cck", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Sox2 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Sox2-P45-1.png){#fig-kim2020-correlation-Cck-Sox2-P45 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Cck-Hlf-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Hlf-E11
#| fig-cap: "Correlation between Cck and Hlf expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Cck", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Hlf expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Hlf-E11-1.png){#fig-kim2020-correlation-Cck-Hlf-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Hlf-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Hlf-E12
#| fig-cap: "Correlation between Cck and Hlf expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Cck", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Hlf expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Hlf-E12-1.png){#fig-kim2020-correlation-Cck-Hlf-E12 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Hlf-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Hlf-E13
#| fig-cap: "Correlation between Cck and Hlf expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Cck", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Hlf expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Hlf-E13-1.png){#fig-kim2020-correlation-Cck-Hlf-E13 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Hlf-E14 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Hlf-E14
#| fig-cap: "Correlation between Cck and Hlf expression at E14 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E14", cells = cells_to_check),
    c("Cck", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E14"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Hlf expression at E14 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Hlf-E14-1.png){#fig-kim2020-correlation-Cck-Hlf-E14 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Hlf-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Hlf-E15
#| fig-cap: "Correlation between Cck and Hlf expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Cck", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Hlf expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Hlf-E15-1.png){#fig-kim2020-correlation-Cck-Hlf-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Hlf-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Hlf-E16
#| fig-cap: "Correlation between Cck and Hlf expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Cck", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Hlf expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Hlf-E16-1.png){#fig-kim2020-correlation-Cck-Hlf-E16 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Hlf-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Hlf-E18
#| fig-cap: "Correlation between Cck and Hlf expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Cck", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Hlf expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Hlf-E18-1.png){#fig-kim2020-correlation-Cck-Hlf-E18 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Hlf-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Hlf-P8
#| fig-cap: "Correlation between Cck and Hlf expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Cck", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Hlf expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Hlf-P8-1.png){#fig-kim2020-correlation-Cck-Hlf-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Hlf-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Hlf-P45
#| fig-cap: "Correlation between Cck and Hlf expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Cck", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Hlf expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Hlf-P45-1.png){#fig-kim2020-correlation-Cck-Hlf-P45 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Cck-Tshr-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Tshr-E11
#| fig-cap: "Correlation between Cck and Tshr expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Cck", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Tshr expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Tshr-E11-1.png){#fig-kim2020-correlation-Cck-Tshr-E11 width=1152}
:::
:::







::: {#cell-fig-kim2020-correlation-Cck-Tshr-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Tshr-E15
#| fig-cap: "Correlation between Cck and Tshr expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Cck", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Tshr expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Tshr-E15-1.png){#fig-kim2020-correlation-Cck-Tshr-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Tshr-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Tshr-E16
#| fig-cap: "Correlation between Cck and Tshr expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Cck", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Tshr expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Tshr-E16-1.png){#fig-kim2020-correlation-Cck-Tshr-E16 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Tshr-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Tshr-E18
#| fig-cap: "Correlation between Cck and Tshr expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Cck", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Tshr expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Tshr-E18-1.png){#fig-kim2020-correlation-Cck-Tshr-E18 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Tshr-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Tshr-P8
#| fig-cap: "Correlation between Cck and Tshr expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Cck", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Tshr expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Tshr-P8-1.png){#fig-kim2020-correlation-Cck-Tshr-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Tshr-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Tshr-P45
#| fig-cap: "Correlation between Cck and Tshr expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Cck", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Tshr expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Tshr-P45-1.png){#fig-kim2020-correlation-Cck-Tshr-P45 width=1152}
:::
:::











::: {#cell-fig-kim2020-correlation-Cck-Cckar-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Cckar-E15
#| fig-cap: "Correlation between Cck and Cckar expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Cck", "Cckar"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Cckar expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Cckar-E15-1.png){#fig-kim2020-correlation-Cck-Cckar-E15 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Cck-Cckar-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Cckar-E18
#| fig-cap: "Correlation between Cck and Cckar expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Cck", "Cckar"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Cckar expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Cckar-E18-1.png){#fig-kim2020-correlation-Cck-Cckar-E18 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Cck-Cckar-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Cckar-P45
#| fig-cap: "Correlation between Cck and Cckar expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Cck", "Cckar"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Cckar expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Cckar-P45-1.png){#fig-kim2020-correlation-Cck-Cckar-P45 width=1152}
:::
:::

















::: {#cell-fig-kim2020-correlation-Cck-Cckbr-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Cckbr-P8
#| fig-cap: "Correlation between Cck and Cckbr expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Cck", "Cckbr"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Cckbr expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Cckbr-P8-1.png){#fig-kim2020-correlation-Cck-Cckbr-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Cckbr-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Cckbr-P45
#| fig-cap: "Correlation between Cck and Cckbr expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Cck", "Cckbr"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Cckbr expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Cckbr-P45-1.png){#fig-kim2020-correlation-Cck-Cckbr-P45 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Cck-Gpr173-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Gpr173-E11
#| fig-cap: "Correlation between Cck and Gpr173 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Cck", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Gpr173 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Gpr173-E11-1.png){#fig-kim2020-correlation-Cck-Gpr173-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Gpr173-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Gpr173-E12
#| fig-cap: "Correlation between Cck and Gpr173 expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Cck", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Gpr173 expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Gpr173-E12-1.png){#fig-kim2020-correlation-Cck-Gpr173-E12 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Gpr173-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Gpr173-E13
#| fig-cap: "Correlation between Cck and Gpr173 expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Cck", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Gpr173 expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Gpr173-E13-1.png){#fig-kim2020-correlation-Cck-Gpr173-E13 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Gpr173-E14 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Gpr173-E14
#| fig-cap: "Correlation between Cck and Gpr173 expression at E14 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E14", cells = cells_to_check),
    c("Cck", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E14"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Gpr173 expression at E14 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Gpr173-E14-1.png){#fig-kim2020-correlation-Cck-Gpr173-E14 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Gpr173-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Gpr173-E15
#| fig-cap: "Correlation between Cck and Gpr173 expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Cck", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Gpr173 expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Gpr173-E15-1.png){#fig-kim2020-correlation-Cck-Gpr173-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Gpr173-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Gpr173-E16
#| fig-cap: "Correlation between Cck and Gpr173 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Cck", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Gpr173 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Gpr173-E16-1.png){#fig-kim2020-correlation-Cck-Gpr173-E16 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Gpr173-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Gpr173-E18
#| fig-cap: "Correlation between Cck and Gpr173 expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Cck", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Gpr173 expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Gpr173-E18-1.png){#fig-kim2020-correlation-Cck-Gpr173-E18 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Gpr173-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Gpr173-P8
#| fig-cap: "Correlation between Cck and Gpr173 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Cck", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Gpr173 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Gpr173-P8-1.png){#fig-kim2020-correlation-Cck-Gpr173-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Gpr173-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Gpr173-P45
#| fig-cap: "Correlation between Cck and Gpr173 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Cck", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Gpr173 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Gpr173-P45-1.png){#fig-kim2020-correlation-Cck-Gpr173-P45 width=1152}
:::
:::

















::: {#cell-fig-kim2020-correlation-Cck-Lhx3-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Lhx3-P8
#| fig-cap: "Correlation between Cck and Lhx3 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Cck", "Lhx3"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Lhx3 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Lhx3-P8-1.png){#fig-kim2020-correlation-Cck-Lhx3-P8 width=1152}
:::
:::





::: {#cell-fig-kim2020-correlation-Cck-Lhx4-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Lhx4-E11
#| fig-cap: "Correlation between Cck and Lhx4 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Cck", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Lhx4 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Lhx4-E11-1.png){#fig-kim2020-correlation-Cck-Lhx4-E11 width=1152}
:::
:::









::: {#cell-fig-kim2020-correlation-Cck-Lhx4-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Lhx4-E16
#| fig-cap: "Correlation between Cck and Lhx4 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Cck", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Lhx4 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Lhx4-E16-1.png){#fig-kim2020-correlation-Cck-Lhx4-E16 width=1152}
:::
:::









::: {#cell-fig-kim2020-correlation-Cck-Gata2-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Gata2-E11
#| fig-cap: "Correlation between Cck and Gata2 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Cck", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Gata2 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Gata2-E11-1.png){#fig-kim2020-correlation-Cck-Gata2-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Gata2-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Gata2-E12
#| fig-cap: "Correlation between Cck and Gata2 expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Cck", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Gata2 expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Gata2-E12-1.png){#fig-kim2020-correlation-Cck-Gata2-E12 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Gata2-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Gata2-E13
#| fig-cap: "Correlation between Cck and Gata2 expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Cck", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Gata2 expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Gata2-E13-1.png){#fig-kim2020-correlation-Cck-Gata2-E13 width=1152}
:::
:::





::: {#cell-fig-kim2020-correlation-Cck-Gata2-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Gata2-E16
#| fig-cap: "Correlation between Cck and Gata2 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Cck", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Gata2 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Gata2-E16-1.png){#fig-kim2020-correlation-Cck-Gata2-E16 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Cck-Gata2-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Gata2-P8
#| fig-cap: "Correlation between Cck and Gata2 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Cck", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Gata2 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Gata2-P8-1.png){#fig-kim2020-correlation-Cck-Gata2-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cck-Gata2-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cck-Gata2-P45
#| fig-cap: "Correlation between Cck and Gata2 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Cck", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Gata2 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cck-Gata2-P45-1.png){#fig-kim2020-correlation-Cck-Gata2-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Pitx1-Eya1-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Eya1-E10
#| fig-cap: "Correlation between Pitx1 and Eya1 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Pitx1", "Eya1"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Eya1 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Eya1-E10-1.png){#fig-kim2020-correlation-Pitx1-Eya1-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Pitx1-Eya1-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Eya1-E11
#| fig-cap: "Correlation between Pitx1 and Eya1 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Pitx1", "Eya1"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Eya1 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Eya1-E11-1.png){#fig-kim2020-correlation-Pitx1-Eya1-E11 width=1152}
:::
:::













::: {#cell-fig-kim2020-correlation-Pitx1-Eya1-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Eya1-P8
#| fig-cap: "Correlation between Pitx1 and Eya1 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Pitx1", "Eya1"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Eya1 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Eya1-P8-1.png){#fig-kim2020-correlation-Pitx1-Eya1-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Pitx1-Eya1-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Eya1-P45
#| fig-cap: "Correlation between Pitx1 and Eya1 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Pitx1", "Eya1"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Eya1 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Eya1-P45-1.png){#fig-kim2020-correlation-Pitx1-Eya1-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Pitx1-Eya2-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Eya2-E10
#| fig-cap: "Correlation between Pitx1 and Eya2 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Pitx1", "Eya2"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Eya2 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Eya2-E10-1.png){#fig-kim2020-correlation-Pitx1-Eya2-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Pitx1-Eya2-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Eya2-E11
#| fig-cap: "Correlation between Pitx1 and Eya2 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Pitx1", "Eya2"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Eya2 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Eya2-E11-1.png){#fig-kim2020-correlation-Pitx1-Eya2-E11 width=1152}
:::
:::













::: {#cell-fig-kim2020-correlation-Pitx1-Eya2-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Eya2-P8
#| fig-cap: "Correlation between Pitx1 and Eya2 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Pitx1", "Eya2"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Eya2 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Eya2-P8-1.png){#fig-kim2020-correlation-Pitx1-Eya2-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Pitx1-Eya2-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Eya2-P45
#| fig-cap: "Correlation between Pitx1 and Eya2 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Pitx1", "Eya2"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Eya2 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Eya2-P45-1.png){#fig-kim2020-correlation-Pitx1-Eya2-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Pitx1-Eya3-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Eya3-E10
#| fig-cap: "Correlation between Pitx1 and Eya3 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Pitx1", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Eya3 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Eya3-E10-1.png){#fig-kim2020-correlation-Pitx1-Eya3-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Pitx1-Eya3-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Eya3-E11
#| fig-cap: "Correlation between Pitx1 and Eya3 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Pitx1", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Eya3 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Eya3-E11-1.png){#fig-kim2020-correlation-Pitx1-Eya3-E11 width=1152}
:::
:::













::: {#cell-fig-kim2020-correlation-Pitx1-Eya3-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Eya3-P8
#| fig-cap: "Correlation between Pitx1 and Eya3 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Pitx1", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Eya3 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Eya3-P8-1.png){#fig-kim2020-correlation-Pitx1-Eya3-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Pitx1-Eya3-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Eya3-P45
#| fig-cap: "Correlation between Pitx1 and Eya3 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Pitx1", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Eya3 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Eya3-P45-1.png){#fig-kim2020-correlation-Pitx1-Eya3-P45 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Pitx1-Eya4-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Eya4-E11
#| fig-cap: "Correlation between Pitx1 and Eya4 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Pitx1", "Eya4"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Eya4 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Eya4-E11-1.png){#fig-kim2020-correlation-Pitx1-Eya4-E11 width=1152}
:::
:::













::: {#cell-fig-kim2020-correlation-Pitx1-Eya4-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Eya4-P8
#| fig-cap: "Correlation between Pitx1 and Eya4 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Pitx1", "Eya4"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Eya4 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Eya4-P8-1.png){#fig-kim2020-correlation-Pitx1-Eya4-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Pitx1-Eya4-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Eya4-P45
#| fig-cap: "Correlation between Pitx1 and Eya4 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Pitx1", "Eya4"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Eya4 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Eya4-P45-1.png){#fig-kim2020-correlation-Pitx1-Eya4-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Pitx1-Sox2-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Sox2-E10
#| fig-cap: "Correlation between Pitx1 and Sox2 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Pitx1", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Sox2 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Sox2-E10-1.png){#fig-kim2020-correlation-Pitx1-Sox2-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Pitx1-Sox2-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Sox2-E11
#| fig-cap: "Correlation between Pitx1 and Sox2 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Pitx1", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Sox2 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Sox2-E11-1.png){#fig-kim2020-correlation-Pitx1-Sox2-E11 width=1152}
:::
:::













::: {#cell-fig-kim2020-correlation-Pitx1-Sox2-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Sox2-P8
#| fig-cap: "Correlation between Pitx1 and Sox2 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Pitx1", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Sox2 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Sox2-P8-1.png){#fig-kim2020-correlation-Pitx1-Sox2-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Pitx1-Sox2-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Sox2-P45
#| fig-cap: "Correlation between Pitx1 and Sox2 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Pitx1", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Sox2 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Sox2-P45-1.png){#fig-kim2020-correlation-Pitx1-Sox2-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Pitx1-Hlf-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Hlf-E10
#| fig-cap: "Correlation between Pitx1 and Hlf expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Pitx1", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Hlf expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Hlf-E10-1.png){#fig-kim2020-correlation-Pitx1-Hlf-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Pitx1-Hlf-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Hlf-E11
#| fig-cap: "Correlation between Pitx1 and Hlf expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Pitx1", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Hlf expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Hlf-E11-1.png){#fig-kim2020-correlation-Pitx1-Hlf-E11 width=1152}
:::
:::













::: {#cell-fig-kim2020-correlation-Pitx1-Hlf-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Hlf-P8
#| fig-cap: "Correlation between Pitx1 and Hlf expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Pitx1", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Hlf expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Hlf-P8-1.png){#fig-kim2020-correlation-Pitx1-Hlf-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Pitx1-Hlf-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Hlf-P45
#| fig-cap: "Correlation between Pitx1 and Hlf expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Pitx1", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Hlf expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Hlf-P45-1.png){#fig-kim2020-correlation-Pitx1-Hlf-P45 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Pitx1-Tshr-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Tshr-E11
#| fig-cap: "Correlation between Pitx1 and Tshr expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Pitx1", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Tshr expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Tshr-E11-1.png){#fig-kim2020-correlation-Pitx1-Tshr-E11 width=1152}
:::
:::













::: {#cell-fig-kim2020-correlation-Pitx1-Tshr-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Tshr-P8
#| fig-cap: "Correlation between Pitx1 and Tshr expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Pitx1", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Tshr expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Tshr-P8-1.png){#fig-kim2020-correlation-Pitx1-Tshr-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Pitx1-Tshr-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Tshr-P45
#| fig-cap: "Correlation between Pitx1 and Tshr expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Pitx1", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Tshr expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Tshr-P45-1.png){#fig-kim2020-correlation-Pitx1-Tshr-P45 width=1152}
:::
:::



















::: {#cell-fig-kim2020-correlation-Pitx1-Cckar-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Cckar-P45
#| fig-cap: "Correlation between Pitx1 and Cckar expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Pitx1", "Cckar"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Cckar expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Cckar-P45-1.png){#fig-kim2020-correlation-Pitx1-Cckar-P45 width=1152}
:::
:::

















::: {#cell-fig-kim2020-correlation-Pitx1-Cckbr-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Cckbr-P8
#| fig-cap: "Correlation between Pitx1 and Cckbr expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Pitx1", "Cckbr"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Cckbr expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Cckbr-P8-1.png){#fig-kim2020-correlation-Pitx1-Cckbr-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Pitx1-Cckbr-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Cckbr-P45
#| fig-cap: "Correlation between Pitx1 and Cckbr expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Pitx1", "Cckbr"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Cckbr expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Cckbr-P45-1.png){#fig-kim2020-correlation-Pitx1-Cckbr-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Pitx1-Gpr173-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Gpr173-E10
#| fig-cap: "Correlation between Pitx1 and Gpr173 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Pitx1", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Gpr173 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Gpr173-E10-1.png){#fig-kim2020-correlation-Pitx1-Gpr173-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Pitx1-Gpr173-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Gpr173-E11
#| fig-cap: "Correlation between Pitx1 and Gpr173 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Pitx1", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Gpr173 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Gpr173-E11-1.png){#fig-kim2020-correlation-Pitx1-Gpr173-E11 width=1152}
:::
:::













::: {#cell-fig-kim2020-correlation-Pitx1-Gpr173-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Gpr173-P8
#| fig-cap: "Correlation between Pitx1 and Gpr173 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Pitx1", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Gpr173 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Gpr173-P8-1.png){#fig-kim2020-correlation-Pitx1-Gpr173-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Pitx1-Gpr173-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Gpr173-P45
#| fig-cap: "Correlation between Pitx1 and Gpr173 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Pitx1", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Gpr173 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Gpr173-P45-1.png){#fig-kim2020-correlation-Pitx1-Gpr173-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Pitx1-Lhx3-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Lhx3-E10
#| fig-cap: "Correlation between Pitx1 and Lhx3 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Pitx1", "Lhx3"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Lhx3 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Lhx3-E10-1.png){#fig-kim2020-correlation-Pitx1-Lhx3-E10 width=1152}
:::
:::















::: {#cell-fig-kim2020-correlation-Pitx1-Lhx3-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Lhx3-P8
#| fig-cap: "Correlation between Pitx1 and Lhx3 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Pitx1", "Lhx3"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Lhx3 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Lhx3-P8-1.png){#fig-kim2020-correlation-Pitx1-Lhx3-P8 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Pitx1-Lhx4-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Lhx4-E10
#| fig-cap: "Correlation between Pitx1 and Lhx4 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Pitx1", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Lhx4 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Lhx4-E10-1.png){#fig-kim2020-correlation-Pitx1-Lhx4-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Pitx1-Lhx4-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Lhx4-E11
#| fig-cap: "Correlation between Pitx1 and Lhx4 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Pitx1", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Lhx4 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Lhx4-E11-1.png){#fig-kim2020-correlation-Pitx1-Lhx4-E11 width=1152}
:::
:::

















::: {#cell-fig-kim2020-correlation-Pitx1-Gata2-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Gata2-E10
#| fig-cap: "Correlation between Pitx1 and Gata2 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Pitx1", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Gata2 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Gata2-E10-1.png){#fig-kim2020-correlation-Pitx1-Gata2-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Pitx1-Gata2-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Gata2-E11
#| fig-cap: "Correlation between Pitx1 and Gata2 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Pitx1", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Gata2 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Gata2-E11-1.png){#fig-kim2020-correlation-Pitx1-Gata2-E11 width=1152}
:::
:::













::: {#cell-fig-kim2020-correlation-Pitx1-Gata2-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Gata2-P8
#| fig-cap: "Correlation between Pitx1 and Gata2 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Pitx1", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Gata2 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Gata2-P8-1.png){#fig-kim2020-correlation-Pitx1-Gata2-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Pitx1-Gata2-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Pitx1-Gata2-P45
#| fig-cap: "Correlation between Pitx1 and Gata2 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Pitx1", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Gata2 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Pitx1-Gata2-P45-1.png){#fig-kim2020-correlation-Pitx1-Gata2-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Eya2-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Eya2-E10
#| fig-cap: "Correlation between Eya1 and Eya2 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Eya1", "Eya2"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya2 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Eya2-E10-1.png){#fig-kim2020-correlation-Eya1-Eya2-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Eya2-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Eya2-E11
#| fig-cap: "Correlation between Eya1 and Eya2 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya1", "Eya2"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya2 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Eya2-E11-1.png){#fig-kim2020-correlation-Eya1-Eya2-E11 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya1-Eya2-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Eya2-E13
#| fig-cap: "Correlation between Eya1 and Eya2 expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Eya1", "Eya2"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya2 expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Eya2-E13-1.png){#fig-kim2020-correlation-Eya1-Eya2-E13 width=1152}
:::
:::









::: {#cell-fig-kim2020-correlation-Eya1-Eya2-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Eya2-P8
#| fig-cap: "Correlation between Eya1 and Eya2 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya1", "Eya2"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya2 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Eya2-P8-1.png){#fig-kim2020-correlation-Eya1-Eya2-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Eya2-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Eya2-P45
#| fig-cap: "Correlation between Eya1 and Eya2 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya1", "Eya2"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya2 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Eya2-P45-1.png){#fig-kim2020-correlation-Eya1-Eya2-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Eya3-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Eya3-E10
#| fig-cap: "Correlation between Eya1 and Eya3 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Eya1", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya3 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Eya3-E10-1.png){#fig-kim2020-correlation-Eya1-Eya3-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Eya3-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Eya3-E11
#| fig-cap: "Correlation between Eya1 and Eya3 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya1", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya3 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Eya3-E11-1.png){#fig-kim2020-correlation-Eya1-Eya3-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Eya3-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Eya3-E12
#| fig-cap: "Correlation between Eya1 and Eya3 expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Eya1", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya3 expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Eya3-E12-1.png){#fig-kim2020-correlation-Eya1-Eya3-E12 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Eya3-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Eya3-E13
#| fig-cap: "Correlation between Eya1 and Eya3 expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Eya1", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya3 expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Eya3-E13-1.png){#fig-kim2020-correlation-Eya1-Eya3-E13 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Eya3-E14 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Eya3-E14
#| fig-cap: "Correlation between Eya1 and Eya3 expression at E14 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E14", cells = cells_to_check),
    c("Eya1", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "E14"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya3 expression at E14 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Eya3-E14-1.png){#fig-kim2020-correlation-Eya1-Eya3-E14 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Eya3-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Eya3-E15
#| fig-cap: "Correlation between Eya1 and Eya3 expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Eya1", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya3 expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Eya3-E15-1.png){#fig-kim2020-correlation-Eya1-Eya3-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Eya3-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Eya3-E16
#| fig-cap: "Correlation between Eya1 and Eya3 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Eya1", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya3 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Eya3-E16-1.png){#fig-kim2020-correlation-Eya1-Eya3-E16 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Eya3-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Eya3-E18
#| fig-cap: "Correlation between Eya1 and Eya3 expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Eya1", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya3 expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Eya3-E18-1.png){#fig-kim2020-correlation-Eya1-Eya3-E18 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Eya3-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Eya3-P8
#| fig-cap: "Correlation between Eya1 and Eya3 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya1", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya3 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Eya3-P8-1.png){#fig-kim2020-correlation-Eya1-Eya3-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Eya3-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Eya3-P45
#| fig-cap: "Correlation between Eya1 and Eya3 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya1", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya3 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Eya3-P45-1.png){#fig-kim2020-correlation-Eya1-Eya3-P45 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya1-Eya4-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Eya4-E11
#| fig-cap: "Correlation between Eya1 and Eya4 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya1", "Eya4"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya4 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Eya4-E11-1.png){#fig-kim2020-correlation-Eya1-Eya4-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Eya4-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Eya4-E12
#| fig-cap: "Correlation between Eya1 and Eya4 expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Eya1", "Eya4"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya4 expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Eya4-E12-1.png){#fig-kim2020-correlation-Eya1-Eya4-E12 width=1152}
:::
:::





::: {#cell-fig-kim2020-correlation-Eya1-Eya4-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Eya4-E15
#| fig-cap: "Correlation between Eya1 and Eya4 expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Eya1", "Eya4"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya4 expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Eya4-E15-1.png){#fig-kim2020-correlation-Eya1-Eya4-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Eya4-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Eya4-E16
#| fig-cap: "Correlation between Eya1 and Eya4 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Eya1", "Eya4"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya4 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Eya4-E16-1.png){#fig-kim2020-correlation-Eya1-Eya4-E16 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya1-Eya4-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Eya4-P8
#| fig-cap: "Correlation between Eya1 and Eya4 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya1", "Eya4"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya4 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Eya4-P8-1.png){#fig-kim2020-correlation-Eya1-Eya4-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Eya4-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Eya4-P45
#| fig-cap: "Correlation between Eya1 and Eya4 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya1", "Eya4"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya4 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Eya4-P45-1.png){#fig-kim2020-correlation-Eya1-Eya4-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Sox2-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Sox2-E10
#| fig-cap: "Correlation between Eya1 and Sox2 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Eya1", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Sox2 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Sox2-E10-1.png){#fig-kim2020-correlation-Eya1-Sox2-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Sox2-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Sox2-E11
#| fig-cap: "Correlation between Eya1 and Sox2 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya1", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Sox2 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Sox2-E11-1.png){#fig-kim2020-correlation-Eya1-Sox2-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Sox2-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Sox2-E12
#| fig-cap: "Correlation between Eya1 and Sox2 expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Eya1", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Sox2 expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Sox2-E12-1.png){#fig-kim2020-correlation-Eya1-Sox2-E12 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Sox2-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Sox2-E13
#| fig-cap: "Correlation between Eya1 and Sox2 expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Eya1", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Sox2 expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Sox2-E13-1.png){#fig-kim2020-correlation-Eya1-Sox2-E13 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Sox2-E14 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Sox2-E14
#| fig-cap: "Correlation between Eya1 and Sox2 expression at E14 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E14", cells = cells_to_check),
    c("Eya1", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E14"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Sox2 expression at E14 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Sox2-E14-1.png){#fig-kim2020-correlation-Eya1-Sox2-E14 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Sox2-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Sox2-E15
#| fig-cap: "Correlation between Eya1 and Sox2 expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Eya1", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Sox2 expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Sox2-E15-1.png){#fig-kim2020-correlation-Eya1-Sox2-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Sox2-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Sox2-E16
#| fig-cap: "Correlation between Eya1 and Sox2 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Eya1", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Sox2 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Sox2-E16-1.png){#fig-kim2020-correlation-Eya1-Sox2-E16 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Sox2-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Sox2-E18
#| fig-cap: "Correlation between Eya1 and Sox2 expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Eya1", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Sox2 expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Sox2-E18-1.png){#fig-kim2020-correlation-Eya1-Sox2-E18 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Sox2-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Sox2-P8
#| fig-cap: "Correlation between Eya1 and Sox2 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya1", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Sox2 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Sox2-P8-1.png){#fig-kim2020-correlation-Eya1-Sox2-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Sox2-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Sox2-P45
#| fig-cap: "Correlation between Eya1 and Sox2 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya1", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Sox2 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Sox2-P45-1.png){#fig-kim2020-correlation-Eya1-Sox2-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Hlf-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Hlf-E10
#| fig-cap: "Correlation between Eya1 and Hlf expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Eya1", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Hlf expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Hlf-E10-1.png){#fig-kim2020-correlation-Eya1-Hlf-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Hlf-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Hlf-E11
#| fig-cap: "Correlation between Eya1 and Hlf expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya1", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Hlf expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Hlf-E11-1.png){#fig-kim2020-correlation-Eya1-Hlf-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Hlf-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Hlf-E12
#| fig-cap: "Correlation between Eya1 and Hlf expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Eya1", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Hlf expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Hlf-E12-1.png){#fig-kim2020-correlation-Eya1-Hlf-E12 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Hlf-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Hlf-E13
#| fig-cap: "Correlation between Eya1 and Hlf expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Eya1", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Hlf expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Hlf-E13-1.png){#fig-kim2020-correlation-Eya1-Hlf-E13 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Hlf-E14 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Hlf-E14
#| fig-cap: "Correlation between Eya1 and Hlf expression at E14 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E14", cells = cells_to_check),
    c("Eya1", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E14"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Hlf expression at E14 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Hlf-E14-1.png){#fig-kim2020-correlation-Eya1-Hlf-E14 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Hlf-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Hlf-E15
#| fig-cap: "Correlation between Eya1 and Hlf expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Eya1", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Hlf expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Hlf-E15-1.png){#fig-kim2020-correlation-Eya1-Hlf-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Hlf-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Hlf-E16
#| fig-cap: "Correlation between Eya1 and Hlf expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Eya1", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Hlf expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Hlf-E16-1.png){#fig-kim2020-correlation-Eya1-Hlf-E16 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Hlf-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Hlf-E18
#| fig-cap: "Correlation between Eya1 and Hlf expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Eya1", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Hlf expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Hlf-E18-1.png){#fig-kim2020-correlation-Eya1-Hlf-E18 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Hlf-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Hlf-P8
#| fig-cap: "Correlation between Eya1 and Hlf expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya1", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Hlf expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Hlf-P8-1.png){#fig-kim2020-correlation-Eya1-Hlf-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Hlf-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Hlf-P45
#| fig-cap: "Correlation between Eya1 and Hlf expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya1", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Hlf expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Hlf-P45-1.png){#fig-kim2020-correlation-Eya1-Hlf-P45 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya1-Tshr-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Tshr-E11
#| fig-cap: "Correlation between Eya1 and Tshr expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya1", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Tshr expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Tshr-E11-1.png){#fig-kim2020-correlation-Eya1-Tshr-E11 width=1152}
:::
:::







::: {#cell-fig-kim2020-correlation-Eya1-Tshr-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Tshr-E15
#| fig-cap: "Correlation between Eya1 and Tshr expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Eya1", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Tshr expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Tshr-E15-1.png){#fig-kim2020-correlation-Eya1-Tshr-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Tshr-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Tshr-E16
#| fig-cap: "Correlation between Eya1 and Tshr expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Eya1", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Tshr expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Tshr-E16-1.png){#fig-kim2020-correlation-Eya1-Tshr-E16 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Tshr-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Tshr-E18
#| fig-cap: "Correlation between Eya1 and Tshr expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Eya1", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Tshr expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Tshr-E18-1.png){#fig-kim2020-correlation-Eya1-Tshr-E18 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Tshr-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Tshr-P8
#| fig-cap: "Correlation between Eya1 and Tshr expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya1", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Tshr expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Tshr-P8-1.png){#fig-kim2020-correlation-Eya1-Tshr-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Tshr-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Tshr-P45
#| fig-cap: "Correlation between Eya1 and Tshr expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya1", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Tshr expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Tshr-P45-1.png){#fig-kim2020-correlation-Eya1-Tshr-P45 width=1152}
:::
:::











::: {#cell-fig-kim2020-correlation-Eya1-Cckar-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Cckar-E15
#| fig-cap: "Correlation between Eya1 and Cckar expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Eya1", "Cckar"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Cckar expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Cckar-E15-1.png){#fig-kim2020-correlation-Eya1-Cckar-E15 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya1-Cckar-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Cckar-E18
#| fig-cap: "Correlation between Eya1 and Cckar expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Eya1", "Cckar"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Cckar expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Cckar-E18-1.png){#fig-kim2020-correlation-Eya1-Cckar-E18 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya1-Cckar-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Cckar-P45
#| fig-cap: "Correlation between Eya1 and Cckar expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya1", "Cckar"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Cckar expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Cckar-P45-1.png){#fig-kim2020-correlation-Eya1-Cckar-P45 width=1152}
:::
:::

















::: {#cell-fig-kim2020-correlation-Eya1-Cckbr-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Cckbr-P8
#| fig-cap: "Correlation between Eya1 and Cckbr expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya1", "Cckbr"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Cckbr expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Cckbr-P8-1.png){#fig-kim2020-correlation-Eya1-Cckbr-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Cckbr-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Cckbr-P45
#| fig-cap: "Correlation between Eya1 and Cckbr expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya1", "Cckbr"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Cckbr expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Cckbr-P45-1.png){#fig-kim2020-correlation-Eya1-Cckbr-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Gpr173-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Gpr173-E10
#| fig-cap: "Correlation between Eya1 and Gpr173 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Eya1", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Gpr173 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Gpr173-E10-1.png){#fig-kim2020-correlation-Eya1-Gpr173-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Gpr173-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Gpr173-E11
#| fig-cap: "Correlation between Eya1 and Gpr173 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya1", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Gpr173 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Gpr173-E11-1.png){#fig-kim2020-correlation-Eya1-Gpr173-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Gpr173-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Gpr173-E12
#| fig-cap: "Correlation between Eya1 and Gpr173 expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Eya1", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Gpr173 expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Gpr173-E12-1.png){#fig-kim2020-correlation-Eya1-Gpr173-E12 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Gpr173-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Gpr173-E13
#| fig-cap: "Correlation between Eya1 and Gpr173 expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Eya1", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Gpr173 expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Gpr173-E13-1.png){#fig-kim2020-correlation-Eya1-Gpr173-E13 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Gpr173-E14 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Gpr173-E14
#| fig-cap: "Correlation between Eya1 and Gpr173 expression at E14 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E14", cells = cells_to_check),
    c("Eya1", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E14"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Gpr173 expression at E14 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Gpr173-E14-1.png){#fig-kim2020-correlation-Eya1-Gpr173-E14 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Gpr173-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Gpr173-E15
#| fig-cap: "Correlation between Eya1 and Gpr173 expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Eya1", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Gpr173 expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Gpr173-E15-1.png){#fig-kim2020-correlation-Eya1-Gpr173-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Gpr173-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Gpr173-E16
#| fig-cap: "Correlation between Eya1 and Gpr173 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Eya1", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Gpr173 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Gpr173-E16-1.png){#fig-kim2020-correlation-Eya1-Gpr173-E16 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Gpr173-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Gpr173-E18
#| fig-cap: "Correlation between Eya1 and Gpr173 expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Eya1", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Gpr173 expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Gpr173-E18-1.png){#fig-kim2020-correlation-Eya1-Gpr173-E18 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Gpr173-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Gpr173-P8
#| fig-cap: "Correlation between Eya1 and Gpr173 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya1", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Gpr173 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Gpr173-P8-1.png){#fig-kim2020-correlation-Eya1-Gpr173-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Gpr173-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Gpr173-P45
#| fig-cap: "Correlation between Eya1 and Gpr173 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya1", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Gpr173 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Gpr173-P45-1.png){#fig-kim2020-correlation-Eya1-Gpr173-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Lhx3-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Lhx3-E10
#| fig-cap: "Correlation between Eya1 and Lhx3 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Eya1", "Lhx3"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Lhx3 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Lhx3-E10-1.png){#fig-kim2020-correlation-Eya1-Lhx3-E10 width=1152}
:::
:::















::: {#cell-fig-kim2020-correlation-Eya1-Lhx3-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Lhx3-P8
#| fig-cap: "Correlation between Eya1 and Lhx3 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya1", "Lhx3"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Lhx3 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Lhx3-P8-1.png){#fig-kim2020-correlation-Eya1-Lhx3-P8 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya1-Lhx4-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Lhx4-E10
#| fig-cap: "Correlation between Eya1 and Lhx4 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Eya1", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Lhx4 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Lhx4-E10-1.png){#fig-kim2020-correlation-Eya1-Lhx4-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Lhx4-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Lhx4-E11
#| fig-cap: "Correlation between Eya1 and Lhx4 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya1", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Lhx4 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Lhx4-E11-1.png){#fig-kim2020-correlation-Eya1-Lhx4-E11 width=1152}
:::
:::









::: {#cell-fig-kim2020-correlation-Eya1-Lhx4-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Lhx4-E16
#| fig-cap: "Correlation between Eya1 and Lhx4 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Eya1", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Lhx4 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Lhx4-E16-1.png){#fig-kim2020-correlation-Eya1-Lhx4-E16 width=1152}
:::
:::







::: {#cell-fig-kim2020-correlation-Eya1-Gata2-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Gata2-E10
#| fig-cap: "Correlation between Eya1 and Gata2 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Eya1", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Gata2 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Gata2-E10-1.png){#fig-kim2020-correlation-Eya1-Gata2-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Gata2-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Gata2-E11
#| fig-cap: "Correlation between Eya1 and Gata2 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya1", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Gata2 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Gata2-E11-1.png){#fig-kim2020-correlation-Eya1-Gata2-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Gata2-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Gata2-E12
#| fig-cap: "Correlation between Eya1 and Gata2 expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Eya1", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Gata2 expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Gata2-E12-1.png){#fig-kim2020-correlation-Eya1-Gata2-E12 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Gata2-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Gata2-E13
#| fig-cap: "Correlation between Eya1 and Gata2 expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Eya1", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Gata2 expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Gata2-E13-1.png){#fig-kim2020-correlation-Eya1-Gata2-E13 width=1152}
:::
:::





::: {#cell-fig-kim2020-correlation-Eya1-Gata2-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Gata2-E16
#| fig-cap: "Correlation between Eya1 and Gata2 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Eya1", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Gata2 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Gata2-E16-1.png){#fig-kim2020-correlation-Eya1-Gata2-E16 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya1-Gata2-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Gata2-P8
#| fig-cap: "Correlation between Eya1 and Gata2 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya1", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Gata2 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Gata2-P8-1.png){#fig-kim2020-correlation-Eya1-Gata2-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya1-Gata2-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya1-Gata2-P45
#| fig-cap: "Correlation between Eya1 and Gata2 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya1", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Gata2 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya1-Gata2-P45-1.png){#fig-kim2020-correlation-Eya1-Gata2-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya2-Eya3-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Eya3-E10
#| fig-cap: "Correlation between Eya2 and Eya3 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Eya2", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Eya3 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Eya3-E10-1.png){#fig-kim2020-correlation-Eya2-Eya3-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya2-Eya3-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Eya3-E11
#| fig-cap: "Correlation between Eya2 and Eya3 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya2", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Eya3 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Eya3-E11-1.png){#fig-kim2020-correlation-Eya2-Eya3-E11 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya2-Eya3-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Eya3-E13
#| fig-cap: "Correlation between Eya2 and Eya3 expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Eya2", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Eya3 expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Eya3-E13-1.png){#fig-kim2020-correlation-Eya2-Eya3-E13 width=1152}
:::
:::









::: {#cell-fig-kim2020-correlation-Eya2-Eya3-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Eya3-P8
#| fig-cap: "Correlation between Eya2 and Eya3 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya2", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Eya3 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Eya3-P8-1.png){#fig-kim2020-correlation-Eya2-Eya3-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya2-Eya3-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Eya3-P45
#| fig-cap: "Correlation between Eya2 and Eya3 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya2", "Eya3"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Eya3 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Eya3-P45-1.png){#fig-kim2020-correlation-Eya2-Eya3-P45 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya2-Eya4-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Eya4-E11
#| fig-cap: "Correlation between Eya2 and Eya4 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya2", "Eya4"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Eya4 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Eya4-E11-1.png){#fig-kim2020-correlation-Eya2-Eya4-E11 width=1152}
:::
:::













::: {#cell-fig-kim2020-correlation-Eya2-Eya4-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Eya4-P8
#| fig-cap: "Correlation between Eya2 and Eya4 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya2", "Eya4"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Eya4 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Eya4-P8-1.png){#fig-kim2020-correlation-Eya2-Eya4-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya2-Eya4-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Eya4-P45
#| fig-cap: "Correlation between Eya2 and Eya4 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya2", "Eya4"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Eya4 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Eya4-P45-1.png){#fig-kim2020-correlation-Eya2-Eya4-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya2-Sox2-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Sox2-E10
#| fig-cap: "Correlation between Eya2 and Sox2 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Eya2", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Sox2 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Sox2-E10-1.png){#fig-kim2020-correlation-Eya2-Sox2-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya2-Sox2-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Sox2-E11
#| fig-cap: "Correlation between Eya2 and Sox2 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya2", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Sox2 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Sox2-E11-1.png){#fig-kim2020-correlation-Eya2-Sox2-E11 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya2-Sox2-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Sox2-E13
#| fig-cap: "Correlation between Eya2 and Sox2 expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Eya2", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Sox2 expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Sox2-E13-1.png){#fig-kim2020-correlation-Eya2-Sox2-E13 width=1152}
:::
:::









::: {#cell-fig-kim2020-correlation-Eya2-Sox2-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Sox2-P8
#| fig-cap: "Correlation between Eya2 and Sox2 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya2", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Sox2 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Sox2-P8-1.png){#fig-kim2020-correlation-Eya2-Sox2-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya2-Sox2-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Sox2-P45
#| fig-cap: "Correlation between Eya2 and Sox2 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya2", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Sox2 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Sox2-P45-1.png){#fig-kim2020-correlation-Eya2-Sox2-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya2-Hlf-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Hlf-E10
#| fig-cap: "Correlation between Eya2 and Hlf expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Eya2", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Hlf expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Hlf-E10-1.png){#fig-kim2020-correlation-Eya2-Hlf-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya2-Hlf-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Hlf-E11
#| fig-cap: "Correlation between Eya2 and Hlf expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya2", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Hlf expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Hlf-E11-1.png){#fig-kim2020-correlation-Eya2-Hlf-E11 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya2-Hlf-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Hlf-E13
#| fig-cap: "Correlation between Eya2 and Hlf expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Eya2", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Hlf expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Hlf-E13-1.png){#fig-kim2020-correlation-Eya2-Hlf-E13 width=1152}
:::
:::









::: {#cell-fig-kim2020-correlation-Eya2-Hlf-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Hlf-P8
#| fig-cap: "Correlation between Eya2 and Hlf expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya2", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Hlf expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Hlf-P8-1.png){#fig-kim2020-correlation-Eya2-Hlf-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya2-Hlf-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Hlf-P45
#| fig-cap: "Correlation between Eya2 and Hlf expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya2", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Hlf expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Hlf-P45-1.png){#fig-kim2020-correlation-Eya2-Hlf-P45 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya2-Tshr-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Tshr-E11
#| fig-cap: "Correlation between Eya2 and Tshr expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya2", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Tshr expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Tshr-E11-1.png){#fig-kim2020-correlation-Eya2-Tshr-E11 width=1152}
:::
:::













::: {#cell-fig-kim2020-correlation-Eya2-Tshr-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Tshr-P8
#| fig-cap: "Correlation between Eya2 and Tshr expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya2", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Tshr expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Tshr-P8-1.png){#fig-kim2020-correlation-Eya2-Tshr-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya2-Tshr-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Tshr-P45
#| fig-cap: "Correlation between Eya2 and Tshr expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya2", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Tshr expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Tshr-P45-1.png){#fig-kim2020-correlation-Eya2-Tshr-P45 width=1152}
:::
:::



















::: {#cell-fig-kim2020-correlation-Eya2-Cckar-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Cckar-P45
#| fig-cap: "Correlation between Eya2 and Cckar expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya2", "Cckar"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Cckar expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Cckar-P45-1.png){#fig-kim2020-correlation-Eya2-Cckar-P45 width=1152}
:::
:::

















::: {#cell-fig-kim2020-correlation-Eya2-Cckbr-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Cckbr-P8
#| fig-cap: "Correlation between Eya2 and Cckbr expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya2", "Cckbr"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Cckbr expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Cckbr-P8-1.png){#fig-kim2020-correlation-Eya2-Cckbr-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya2-Cckbr-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Cckbr-P45
#| fig-cap: "Correlation between Eya2 and Cckbr expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya2", "Cckbr"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Cckbr expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Cckbr-P45-1.png){#fig-kim2020-correlation-Eya2-Cckbr-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya2-Gpr173-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Gpr173-E10
#| fig-cap: "Correlation between Eya2 and Gpr173 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Eya2", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Gpr173 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Gpr173-E10-1.png){#fig-kim2020-correlation-Eya2-Gpr173-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya2-Gpr173-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Gpr173-E11
#| fig-cap: "Correlation between Eya2 and Gpr173 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya2", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Gpr173 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Gpr173-E11-1.png){#fig-kim2020-correlation-Eya2-Gpr173-E11 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya2-Gpr173-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Gpr173-E13
#| fig-cap: "Correlation between Eya2 and Gpr173 expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Eya2", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Gpr173 expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Gpr173-E13-1.png){#fig-kim2020-correlation-Eya2-Gpr173-E13 width=1152}
:::
:::









::: {#cell-fig-kim2020-correlation-Eya2-Gpr173-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Gpr173-P8
#| fig-cap: "Correlation between Eya2 and Gpr173 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya2", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Gpr173 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Gpr173-P8-1.png){#fig-kim2020-correlation-Eya2-Gpr173-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya2-Gpr173-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Gpr173-P45
#| fig-cap: "Correlation between Eya2 and Gpr173 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya2", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Gpr173 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Gpr173-P45-1.png){#fig-kim2020-correlation-Eya2-Gpr173-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya2-Lhx3-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Lhx3-E10
#| fig-cap: "Correlation between Eya2 and Lhx3 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Eya2", "Lhx3"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Lhx3 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Lhx3-E10-1.png){#fig-kim2020-correlation-Eya2-Lhx3-E10 width=1152}
:::
:::















::: {#cell-fig-kim2020-correlation-Eya2-Lhx3-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Lhx3-P8
#| fig-cap: "Correlation between Eya2 and Lhx3 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya2", "Lhx3"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Lhx3 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Lhx3-P8-1.png){#fig-kim2020-correlation-Eya2-Lhx3-P8 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya2-Lhx4-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Lhx4-E10
#| fig-cap: "Correlation between Eya2 and Lhx4 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Eya2", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Lhx4 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Lhx4-E10-1.png){#fig-kim2020-correlation-Eya2-Lhx4-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya2-Lhx4-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Lhx4-E11
#| fig-cap: "Correlation between Eya2 and Lhx4 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya2", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Lhx4 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Lhx4-E11-1.png){#fig-kim2020-correlation-Eya2-Lhx4-E11 width=1152}
:::
:::

















::: {#cell-fig-kim2020-correlation-Eya2-Gata2-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Gata2-E10
#| fig-cap: "Correlation between Eya2 and Gata2 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Eya2", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Gata2 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Gata2-E10-1.png){#fig-kim2020-correlation-Eya2-Gata2-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya2-Gata2-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Gata2-E11
#| fig-cap: "Correlation between Eya2 and Gata2 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya2", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Gata2 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Gata2-E11-1.png){#fig-kim2020-correlation-Eya2-Gata2-E11 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya2-Gata2-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Gata2-E13
#| fig-cap: "Correlation between Eya2 and Gata2 expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Eya2", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Gata2 expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Gata2-E13-1.png){#fig-kim2020-correlation-Eya2-Gata2-E13 width=1152}
:::
:::









::: {#cell-fig-kim2020-correlation-Eya2-Gata2-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Gata2-P8
#| fig-cap: "Correlation between Eya2 and Gata2 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya2", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Gata2 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Gata2-P8-1.png){#fig-kim2020-correlation-Eya2-Gata2-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya2-Gata2-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya2-Gata2-P45
#| fig-cap: "Correlation between Eya2 and Gata2 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya2", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Gata2 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya2-Gata2-P45-1.png){#fig-kim2020-correlation-Eya2-Gata2-P45 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya3-Eya4-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Eya4-E11
#| fig-cap: "Correlation between Eya3 and Eya4 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya3", "Eya4"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Eya4 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Eya4-E11-1.png){#fig-kim2020-correlation-Eya3-Eya4-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Eya4-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Eya4-E12
#| fig-cap: "Correlation between Eya3 and Eya4 expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Eya3", "Eya4"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Eya4 expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Eya4-E12-1.png){#fig-kim2020-correlation-Eya3-Eya4-E12 width=1152}
:::
:::





::: {#cell-fig-kim2020-correlation-Eya3-Eya4-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Eya4-E15
#| fig-cap: "Correlation between Eya3 and Eya4 expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Eya3", "Eya4"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Eya4 expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Eya4-E15-1.png){#fig-kim2020-correlation-Eya3-Eya4-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Eya4-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Eya4-E16
#| fig-cap: "Correlation between Eya3 and Eya4 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Eya3", "Eya4"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Eya4 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Eya4-E16-1.png){#fig-kim2020-correlation-Eya3-Eya4-E16 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya3-Eya4-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Eya4-P8
#| fig-cap: "Correlation between Eya3 and Eya4 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya3", "Eya4"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Eya4 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Eya4-P8-1.png){#fig-kim2020-correlation-Eya3-Eya4-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Eya4-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Eya4-P45
#| fig-cap: "Correlation between Eya3 and Eya4 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya3", "Eya4"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Eya4 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Eya4-P45-1.png){#fig-kim2020-correlation-Eya3-Eya4-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Sox2-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Sox2-E10
#| fig-cap: "Correlation between Eya3 and Sox2 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Eya3", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Sox2 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Sox2-E10-1.png){#fig-kim2020-correlation-Eya3-Sox2-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Sox2-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Sox2-E11
#| fig-cap: "Correlation between Eya3 and Sox2 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya3", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Sox2 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Sox2-E11-1.png){#fig-kim2020-correlation-Eya3-Sox2-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Sox2-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Sox2-E12
#| fig-cap: "Correlation between Eya3 and Sox2 expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Eya3", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Sox2 expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Sox2-E12-1.png){#fig-kim2020-correlation-Eya3-Sox2-E12 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Sox2-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Sox2-E13
#| fig-cap: "Correlation between Eya3 and Sox2 expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Eya3", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Sox2 expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Sox2-E13-1.png){#fig-kim2020-correlation-Eya3-Sox2-E13 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Sox2-E14 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Sox2-E14
#| fig-cap: "Correlation between Eya3 and Sox2 expression at E14 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E14", cells = cells_to_check),
    c("Eya3", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E14"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Sox2 expression at E14 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Sox2-E14-1.png){#fig-kim2020-correlation-Eya3-Sox2-E14 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Sox2-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Sox2-E15
#| fig-cap: "Correlation between Eya3 and Sox2 expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Eya3", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Sox2 expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Sox2-E15-1.png){#fig-kim2020-correlation-Eya3-Sox2-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Sox2-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Sox2-E16
#| fig-cap: "Correlation between Eya3 and Sox2 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Eya3", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Sox2 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Sox2-E16-1.png){#fig-kim2020-correlation-Eya3-Sox2-E16 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Sox2-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Sox2-E18
#| fig-cap: "Correlation between Eya3 and Sox2 expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Eya3", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Sox2 expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Sox2-E18-1.png){#fig-kim2020-correlation-Eya3-Sox2-E18 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Sox2-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Sox2-P8
#| fig-cap: "Correlation between Eya3 and Sox2 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya3", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Sox2 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Sox2-P8-1.png){#fig-kim2020-correlation-Eya3-Sox2-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Sox2-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Sox2-P45
#| fig-cap: "Correlation between Eya3 and Sox2 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya3", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Sox2 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Sox2-P45-1.png){#fig-kim2020-correlation-Eya3-Sox2-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Hlf-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Hlf-E10
#| fig-cap: "Correlation between Eya3 and Hlf expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Eya3", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Hlf expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Hlf-E10-1.png){#fig-kim2020-correlation-Eya3-Hlf-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Hlf-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Hlf-E11
#| fig-cap: "Correlation between Eya3 and Hlf expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya3", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Hlf expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Hlf-E11-1.png){#fig-kim2020-correlation-Eya3-Hlf-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Hlf-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Hlf-E12
#| fig-cap: "Correlation between Eya3 and Hlf expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Eya3", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Hlf expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Hlf-E12-1.png){#fig-kim2020-correlation-Eya3-Hlf-E12 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Hlf-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Hlf-E13
#| fig-cap: "Correlation between Eya3 and Hlf expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Eya3", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Hlf expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Hlf-E13-1.png){#fig-kim2020-correlation-Eya3-Hlf-E13 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Hlf-E14 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Hlf-E14
#| fig-cap: "Correlation between Eya3 and Hlf expression at E14 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E14", cells = cells_to_check),
    c("Eya3", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E14"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Hlf expression at E14 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Hlf-E14-1.png){#fig-kim2020-correlation-Eya3-Hlf-E14 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Hlf-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Hlf-E15
#| fig-cap: "Correlation between Eya3 and Hlf expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Eya3", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Hlf expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Hlf-E15-1.png){#fig-kim2020-correlation-Eya3-Hlf-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Hlf-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Hlf-E16
#| fig-cap: "Correlation between Eya3 and Hlf expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Eya3", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Hlf expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Hlf-E16-1.png){#fig-kim2020-correlation-Eya3-Hlf-E16 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Hlf-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Hlf-E18
#| fig-cap: "Correlation between Eya3 and Hlf expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Eya3", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Hlf expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Hlf-E18-1.png){#fig-kim2020-correlation-Eya3-Hlf-E18 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Hlf-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Hlf-P8
#| fig-cap: "Correlation between Eya3 and Hlf expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya3", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Hlf expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Hlf-P8-1.png){#fig-kim2020-correlation-Eya3-Hlf-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Hlf-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Hlf-P45
#| fig-cap: "Correlation between Eya3 and Hlf expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya3", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Hlf expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Hlf-P45-1.png){#fig-kim2020-correlation-Eya3-Hlf-P45 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya3-Tshr-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Tshr-E11
#| fig-cap: "Correlation between Eya3 and Tshr expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya3", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Tshr expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Tshr-E11-1.png){#fig-kim2020-correlation-Eya3-Tshr-E11 width=1152}
:::
:::







::: {#cell-fig-kim2020-correlation-Eya3-Tshr-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Tshr-E15
#| fig-cap: "Correlation between Eya3 and Tshr expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Eya3", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Tshr expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Tshr-E15-1.png){#fig-kim2020-correlation-Eya3-Tshr-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Tshr-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Tshr-E16
#| fig-cap: "Correlation between Eya3 and Tshr expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Eya3", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Tshr expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Tshr-E16-1.png){#fig-kim2020-correlation-Eya3-Tshr-E16 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Tshr-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Tshr-E18
#| fig-cap: "Correlation between Eya3 and Tshr expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Eya3", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Tshr expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Tshr-E18-1.png){#fig-kim2020-correlation-Eya3-Tshr-E18 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Tshr-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Tshr-P8
#| fig-cap: "Correlation between Eya3 and Tshr expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya3", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Tshr expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Tshr-P8-1.png){#fig-kim2020-correlation-Eya3-Tshr-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Tshr-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Tshr-P45
#| fig-cap: "Correlation between Eya3 and Tshr expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya3", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Tshr expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Tshr-P45-1.png){#fig-kim2020-correlation-Eya3-Tshr-P45 width=1152}
:::
:::











::: {#cell-fig-kim2020-correlation-Eya3-Cckar-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Cckar-E15
#| fig-cap: "Correlation between Eya3 and Cckar expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Eya3", "Cckar"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Cckar expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Cckar-E15-1.png){#fig-kim2020-correlation-Eya3-Cckar-E15 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya3-Cckar-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Cckar-E18
#| fig-cap: "Correlation between Eya3 and Cckar expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Eya3", "Cckar"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Cckar expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Cckar-E18-1.png){#fig-kim2020-correlation-Eya3-Cckar-E18 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya3-Cckar-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Cckar-P45
#| fig-cap: "Correlation between Eya3 and Cckar expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya3", "Cckar"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Cckar expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Cckar-P45-1.png){#fig-kim2020-correlation-Eya3-Cckar-P45 width=1152}
:::
:::

















::: {#cell-fig-kim2020-correlation-Eya3-Cckbr-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Cckbr-P8
#| fig-cap: "Correlation between Eya3 and Cckbr expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya3", "Cckbr"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Cckbr expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Cckbr-P8-1.png){#fig-kim2020-correlation-Eya3-Cckbr-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Cckbr-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Cckbr-P45
#| fig-cap: "Correlation between Eya3 and Cckbr expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya3", "Cckbr"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Cckbr expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Cckbr-P45-1.png){#fig-kim2020-correlation-Eya3-Cckbr-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Gpr173-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Gpr173-E10
#| fig-cap: "Correlation between Eya3 and Gpr173 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Eya3", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Gpr173 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Gpr173-E10-1.png){#fig-kim2020-correlation-Eya3-Gpr173-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Gpr173-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Gpr173-E11
#| fig-cap: "Correlation between Eya3 and Gpr173 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya3", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Gpr173 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Gpr173-E11-1.png){#fig-kim2020-correlation-Eya3-Gpr173-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Gpr173-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Gpr173-E12
#| fig-cap: "Correlation between Eya3 and Gpr173 expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Eya3", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Gpr173 expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Gpr173-E12-1.png){#fig-kim2020-correlation-Eya3-Gpr173-E12 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Gpr173-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Gpr173-E13
#| fig-cap: "Correlation between Eya3 and Gpr173 expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Eya3", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Gpr173 expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Gpr173-E13-1.png){#fig-kim2020-correlation-Eya3-Gpr173-E13 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Gpr173-E14 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Gpr173-E14
#| fig-cap: "Correlation between Eya3 and Gpr173 expression at E14 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E14", cells = cells_to_check),
    c("Eya3", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E14"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Gpr173 expression at E14 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Gpr173-E14-1.png){#fig-kim2020-correlation-Eya3-Gpr173-E14 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Gpr173-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Gpr173-E15
#| fig-cap: "Correlation between Eya3 and Gpr173 expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Eya3", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Gpr173 expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Gpr173-E15-1.png){#fig-kim2020-correlation-Eya3-Gpr173-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Gpr173-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Gpr173-E16
#| fig-cap: "Correlation between Eya3 and Gpr173 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Eya3", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Gpr173 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Gpr173-E16-1.png){#fig-kim2020-correlation-Eya3-Gpr173-E16 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Gpr173-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Gpr173-E18
#| fig-cap: "Correlation between Eya3 and Gpr173 expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Eya3", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Gpr173 expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Gpr173-E18-1.png){#fig-kim2020-correlation-Eya3-Gpr173-E18 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Gpr173-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Gpr173-P8
#| fig-cap: "Correlation between Eya3 and Gpr173 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya3", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Gpr173 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Gpr173-P8-1.png){#fig-kim2020-correlation-Eya3-Gpr173-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Gpr173-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Gpr173-P45
#| fig-cap: "Correlation between Eya3 and Gpr173 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya3", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Gpr173 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Gpr173-P45-1.png){#fig-kim2020-correlation-Eya3-Gpr173-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Lhx3-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Lhx3-E10
#| fig-cap: "Correlation between Eya3 and Lhx3 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Eya3", "Lhx3"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Lhx3 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Lhx3-E10-1.png){#fig-kim2020-correlation-Eya3-Lhx3-E10 width=1152}
:::
:::















::: {#cell-fig-kim2020-correlation-Eya3-Lhx3-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Lhx3-P8
#| fig-cap: "Correlation between Eya3 and Lhx3 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya3", "Lhx3"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Lhx3 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Lhx3-P8-1.png){#fig-kim2020-correlation-Eya3-Lhx3-P8 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya3-Lhx4-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Lhx4-E10
#| fig-cap: "Correlation between Eya3 and Lhx4 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Eya3", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Lhx4 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Lhx4-E10-1.png){#fig-kim2020-correlation-Eya3-Lhx4-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Lhx4-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Lhx4-E11
#| fig-cap: "Correlation between Eya3 and Lhx4 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya3", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Lhx4 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Lhx4-E11-1.png){#fig-kim2020-correlation-Eya3-Lhx4-E11 width=1152}
:::
:::









::: {#cell-fig-kim2020-correlation-Eya3-Lhx4-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Lhx4-E16
#| fig-cap: "Correlation between Eya3 and Lhx4 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Eya3", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Lhx4 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Lhx4-E16-1.png){#fig-kim2020-correlation-Eya3-Lhx4-E16 width=1152}
:::
:::







::: {#cell-fig-kim2020-correlation-Eya3-Gata2-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Gata2-E10
#| fig-cap: "Correlation between Eya3 and Gata2 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Eya3", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Gata2 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Gata2-E10-1.png){#fig-kim2020-correlation-Eya3-Gata2-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Gata2-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Gata2-E11
#| fig-cap: "Correlation between Eya3 and Gata2 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya3", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Gata2 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Gata2-E11-1.png){#fig-kim2020-correlation-Eya3-Gata2-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Gata2-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Gata2-E12
#| fig-cap: "Correlation between Eya3 and Gata2 expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Eya3", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Gata2 expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Gata2-E12-1.png){#fig-kim2020-correlation-Eya3-Gata2-E12 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Gata2-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Gata2-E13
#| fig-cap: "Correlation between Eya3 and Gata2 expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Eya3", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Gata2 expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Gata2-E13-1.png){#fig-kim2020-correlation-Eya3-Gata2-E13 width=1152}
:::
:::





::: {#cell-fig-kim2020-correlation-Eya3-Gata2-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Gata2-E16
#| fig-cap: "Correlation between Eya3 and Gata2 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Eya3", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Gata2 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Gata2-E16-1.png){#fig-kim2020-correlation-Eya3-Gata2-E16 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya3-Gata2-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Gata2-P8
#| fig-cap: "Correlation between Eya3 and Gata2 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya3", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Gata2 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Gata2-P8-1.png){#fig-kim2020-correlation-Eya3-Gata2-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya3-Gata2-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya3-Gata2-P45
#| fig-cap: "Correlation between Eya3 and Gata2 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya3", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Gata2 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya3-Gata2-P45-1.png){#fig-kim2020-correlation-Eya3-Gata2-P45 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya4-Sox2-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Sox2-E11
#| fig-cap: "Correlation between Eya4 and Sox2 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya4", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Sox2 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Sox2-E11-1.png){#fig-kim2020-correlation-Eya4-Sox2-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya4-Sox2-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Sox2-E12
#| fig-cap: "Correlation between Eya4 and Sox2 expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Eya4", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Sox2 expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Sox2-E12-1.png){#fig-kim2020-correlation-Eya4-Sox2-E12 width=1152}
:::
:::





::: {#cell-fig-kim2020-correlation-Eya4-Sox2-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Sox2-E15
#| fig-cap: "Correlation between Eya4 and Sox2 expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Eya4", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Sox2 expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Sox2-E15-1.png){#fig-kim2020-correlation-Eya4-Sox2-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya4-Sox2-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Sox2-E16
#| fig-cap: "Correlation between Eya4 and Sox2 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Eya4", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Sox2 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Sox2-E16-1.png){#fig-kim2020-correlation-Eya4-Sox2-E16 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya4-Sox2-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Sox2-P8
#| fig-cap: "Correlation between Eya4 and Sox2 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya4", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Sox2 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Sox2-P8-1.png){#fig-kim2020-correlation-Eya4-Sox2-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya4-Sox2-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Sox2-P45
#| fig-cap: "Correlation between Eya4 and Sox2 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya4", "Sox2"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Sox2 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Sox2-P45-1.png){#fig-kim2020-correlation-Eya4-Sox2-P45 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya4-Hlf-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Hlf-E11
#| fig-cap: "Correlation between Eya4 and Hlf expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya4", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Hlf expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Hlf-E11-1.png){#fig-kim2020-correlation-Eya4-Hlf-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya4-Hlf-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Hlf-E12
#| fig-cap: "Correlation between Eya4 and Hlf expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Eya4", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Hlf expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Hlf-E12-1.png){#fig-kim2020-correlation-Eya4-Hlf-E12 width=1152}
:::
:::





::: {#cell-fig-kim2020-correlation-Eya4-Hlf-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Hlf-E15
#| fig-cap: "Correlation between Eya4 and Hlf expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Eya4", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Hlf expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Hlf-E15-1.png){#fig-kim2020-correlation-Eya4-Hlf-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya4-Hlf-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Hlf-E16
#| fig-cap: "Correlation between Eya4 and Hlf expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Eya4", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Hlf expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Hlf-E16-1.png){#fig-kim2020-correlation-Eya4-Hlf-E16 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya4-Hlf-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Hlf-P8
#| fig-cap: "Correlation between Eya4 and Hlf expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya4", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Hlf expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Hlf-P8-1.png){#fig-kim2020-correlation-Eya4-Hlf-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya4-Hlf-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Hlf-P45
#| fig-cap: "Correlation between Eya4 and Hlf expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya4", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Hlf expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Hlf-P45-1.png){#fig-kim2020-correlation-Eya4-Hlf-P45 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya4-Tshr-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Tshr-E11
#| fig-cap: "Correlation between Eya4 and Tshr expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya4", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Tshr expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Tshr-E11-1.png){#fig-kim2020-correlation-Eya4-Tshr-E11 width=1152}
:::
:::







::: {#cell-fig-kim2020-correlation-Eya4-Tshr-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Tshr-E15
#| fig-cap: "Correlation between Eya4 and Tshr expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Eya4", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Tshr expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Tshr-E15-1.png){#fig-kim2020-correlation-Eya4-Tshr-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya4-Tshr-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Tshr-E16
#| fig-cap: "Correlation between Eya4 and Tshr expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Eya4", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Tshr expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Tshr-E16-1.png){#fig-kim2020-correlation-Eya4-Tshr-E16 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya4-Tshr-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Tshr-P8
#| fig-cap: "Correlation between Eya4 and Tshr expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya4", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Tshr expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Tshr-P8-1.png){#fig-kim2020-correlation-Eya4-Tshr-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya4-Tshr-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Tshr-P45
#| fig-cap: "Correlation between Eya4 and Tshr expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya4", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Tshr expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Tshr-P45-1.png){#fig-kim2020-correlation-Eya4-Tshr-P45 width=1152}
:::
:::











::: {#cell-fig-kim2020-correlation-Eya4-Cckar-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Cckar-E15
#| fig-cap: "Correlation between Eya4 and Cckar expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Eya4", "Cckar"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Cckar expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Cckar-E15-1.png){#fig-kim2020-correlation-Eya4-Cckar-E15 width=1152}
:::
:::







::: {#cell-fig-kim2020-correlation-Eya4-Cckar-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Cckar-P45
#| fig-cap: "Correlation between Eya4 and Cckar expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya4", "Cckar"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Cckar expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Cckar-P45-1.png){#fig-kim2020-correlation-Eya4-Cckar-P45 width=1152}
:::
:::

















::: {#cell-fig-kim2020-correlation-Eya4-Cckbr-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Cckbr-P8
#| fig-cap: "Correlation between Eya4 and Cckbr expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya4", "Cckbr"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Cckbr expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Cckbr-P8-1.png){#fig-kim2020-correlation-Eya4-Cckbr-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya4-Cckbr-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Cckbr-P45
#| fig-cap: "Correlation between Eya4 and Cckbr expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya4", "Cckbr"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Cckbr expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Cckbr-P45-1.png){#fig-kim2020-correlation-Eya4-Cckbr-P45 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya4-Gpr173-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Gpr173-E11
#| fig-cap: "Correlation between Eya4 and Gpr173 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya4", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Gpr173 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Gpr173-E11-1.png){#fig-kim2020-correlation-Eya4-Gpr173-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya4-Gpr173-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Gpr173-E12
#| fig-cap: "Correlation between Eya4 and Gpr173 expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Eya4", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Gpr173 expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Gpr173-E12-1.png){#fig-kim2020-correlation-Eya4-Gpr173-E12 width=1152}
:::
:::





::: {#cell-fig-kim2020-correlation-Eya4-Gpr173-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Gpr173-E15
#| fig-cap: "Correlation between Eya4 and Gpr173 expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Eya4", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Gpr173 expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Gpr173-E15-1.png){#fig-kim2020-correlation-Eya4-Gpr173-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya4-Gpr173-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Gpr173-E16
#| fig-cap: "Correlation between Eya4 and Gpr173 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Eya4", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Gpr173 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Gpr173-E16-1.png){#fig-kim2020-correlation-Eya4-Gpr173-E16 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya4-Gpr173-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Gpr173-P8
#| fig-cap: "Correlation between Eya4 and Gpr173 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya4", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Gpr173 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Gpr173-P8-1.png){#fig-kim2020-correlation-Eya4-Gpr173-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya4-Gpr173-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Gpr173-P45
#| fig-cap: "Correlation between Eya4 and Gpr173 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya4", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Gpr173 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Gpr173-P45-1.png){#fig-kim2020-correlation-Eya4-Gpr173-P45 width=1152}
:::
:::

















::: {#cell-fig-kim2020-correlation-Eya4-Lhx3-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Lhx3-P8
#| fig-cap: "Correlation between Eya4 and Lhx3 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya4", "Lhx3"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Lhx3 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Lhx3-P8-1.png){#fig-kim2020-correlation-Eya4-Lhx3-P8 width=1152}
:::
:::





::: {#cell-fig-kim2020-correlation-Eya4-Lhx4-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Lhx4-E11
#| fig-cap: "Correlation between Eya4 and Lhx4 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya4", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Lhx4 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Lhx4-E11-1.png){#fig-kim2020-correlation-Eya4-Lhx4-E11 width=1152}
:::
:::









::: {#cell-fig-kim2020-correlation-Eya4-Lhx4-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Lhx4-E16
#| fig-cap: "Correlation between Eya4 and Lhx4 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Eya4", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Lhx4 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Lhx4-E16-1.png){#fig-kim2020-correlation-Eya4-Lhx4-E16 width=1152}
:::
:::









::: {#cell-fig-kim2020-correlation-Eya4-Gata2-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Gata2-E11
#| fig-cap: "Correlation between Eya4 and Gata2 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Eya4", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Gata2 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Gata2-E11-1.png){#fig-kim2020-correlation-Eya4-Gata2-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya4-Gata2-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Gata2-E12
#| fig-cap: "Correlation between Eya4 and Gata2 expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Eya4", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Gata2 expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Gata2-E12-1.png){#fig-kim2020-correlation-Eya4-Gata2-E12 width=1152}
:::
:::







::: {#cell-fig-kim2020-correlation-Eya4-Gata2-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Gata2-E16
#| fig-cap: "Correlation between Eya4 and Gata2 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Eya4", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Gata2 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Gata2-E16-1.png){#fig-kim2020-correlation-Eya4-Gata2-E16 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Eya4-Gata2-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Gata2-P8
#| fig-cap: "Correlation between Eya4 and Gata2 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Eya4", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Gata2 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Gata2-P8-1.png){#fig-kim2020-correlation-Eya4-Gata2-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Eya4-Gata2-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Eya4-Gata2-P45
#| fig-cap: "Correlation between Eya4 and Gata2 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Eya4", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Gata2 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Eya4-Gata2-P45-1.png){#fig-kim2020-correlation-Eya4-Gata2-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Hlf-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Hlf-E10
#| fig-cap: "Correlation between Sox2 and Hlf expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Sox2", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Hlf expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Hlf-E10-1.png){#fig-kim2020-correlation-Sox2-Hlf-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Hlf-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Hlf-E11
#| fig-cap: "Correlation between Sox2 and Hlf expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Sox2", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Hlf expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Hlf-E11-1.png){#fig-kim2020-correlation-Sox2-Hlf-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Hlf-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Hlf-E12
#| fig-cap: "Correlation between Sox2 and Hlf expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Sox2", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Hlf expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Hlf-E12-1.png){#fig-kim2020-correlation-Sox2-Hlf-E12 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Hlf-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Hlf-E13
#| fig-cap: "Correlation between Sox2 and Hlf expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Sox2", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Hlf expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Hlf-E13-1.png){#fig-kim2020-correlation-Sox2-Hlf-E13 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Hlf-E14 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Hlf-E14
#| fig-cap: "Correlation between Sox2 and Hlf expression at E14 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E14", cells = cells_to_check),
    c("Sox2", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E14"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Hlf expression at E14 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Hlf-E14-1.png){#fig-kim2020-correlation-Sox2-Hlf-E14 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Hlf-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Hlf-E15
#| fig-cap: "Correlation between Sox2 and Hlf expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Sox2", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Hlf expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Hlf-E15-1.png){#fig-kim2020-correlation-Sox2-Hlf-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Hlf-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Hlf-E16
#| fig-cap: "Correlation between Sox2 and Hlf expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Sox2", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Hlf expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Hlf-E16-1.png){#fig-kim2020-correlation-Sox2-Hlf-E16 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Hlf-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Hlf-E18
#| fig-cap: "Correlation between Sox2 and Hlf expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Sox2", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Hlf expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Hlf-E18-1.png){#fig-kim2020-correlation-Sox2-Hlf-E18 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Hlf-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Hlf-P8
#| fig-cap: "Correlation between Sox2 and Hlf expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Sox2", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Hlf expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Hlf-P8-1.png){#fig-kim2020-correlation-Sox2-Hlf-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Hlf-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Hlf-P45
#| fig-cap: "Correlation between Sox2 and Hlf expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Sox2", "Hlf"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Hlf expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Hlf-P45-1.png){#fig-kim2020-correlation-Sox2-Hlf-P45 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Sox2-Tshr-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Tshr-E11
#| fig-cap: "Correlation between Sox2 and Tshr expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Sox2", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Tshr expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Tshr-E11-1.png){#fig-kim2020-correlation-Sox2-Tshr-E11 width=1152}
:::
:::







::: {#cell-fig-kim2020-correlation-Sox2-Tshr-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Tshr-E15
#| fig-cap: "Correlation between Sox2 and Tshr expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Sox2", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Tshr expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Tshr-E15-1.png){#fig-kim2020-correlation-Sox2-Tshr-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Tshr-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Tshr-E16
#| fig-cap: "Correlation between Sox2 and Tshr expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Sox2", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Tshr expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Tshr-E16-1.png){#fig-kim2020-correlation-Sox2-Tshr-E16 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Tshr-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Tshr-E18
#| fig-cap: "Correlation between Sox2 and Tshr expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Sox2", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Tshr expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Tshr-E18-1.png){#fig-kim2020-correlation-Sox2-Tshr-E18 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Tshr-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Tshr-P8
#| fig-cap: "Correlation between Sox2 and Tshr expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Sox2", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Tshr expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Tshr-P8-1.png){#fig-kim2020-correlation-Sox2-Tshr-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Tshr-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Tshr-P45
#| fig-cap: "Correlation between Sox2 and Tshr expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Sox2", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Tshr expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Tshr-P45-1.png){#fig-kim2020-correlation-Sox2-Tshr-P45 width=1152}
:::
:::











::: {#cell-fig-kim2020-correlation-Sox2-Cckar-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Cckar-E15
#| fig-cap: "Correlation between Sox2 and Cckar expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Sox2", "Cckar"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Cckar expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Cckar-E15-1.png){#fig-kim2020-correlation-Sox2-Cckar-E15 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Sox2-Cckar-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Cckar-E18
#| fig-cap: "Correlation between Sox2 and Cckar expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Sox2", "Cckar"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Cckar expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Cckar-E18-1.png){#fig-kim2020-correlation-Sox2-Cckar-E18 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Sox2-Cckar-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Cckar-P45
#| fig-cap: "Correlation between Sox2 and Cckar expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Sox2", "Cckar"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Cckar expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Cckar-P45-1.png){#fig-kim2020-correlation-Sox2-Cckar-P45 width=1152}
:::
:::

















::: {#cell-fig-kim2020-correlation-Sox2-Cckbr-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Cckbr-P8
#| fig-cap: "Correlation between Sox2 and Cckbr expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Sox2", "Cckbr"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Cckbr expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Cckbr-P8-1.png){#fig-kim2020-correlation-Sox2-Cckbr-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Cckbr-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Cckbr-P45
#| fig-cap: "Correlation between Sox2 and Cckbr expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Sox2", "Cckbr"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Cckbr expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Cckbr-P45-1.png){#fig-kim2020-correlation-Sox2-Cckbr-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Gpr173-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Gpr173-E10
#| fig-cap: "Correlation between Sox2 and Gpr173 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Sox2", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Gpr173 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Gpr173-E10-1.png){#fig-kim2020-correlation-Sox2-Gpr173-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Gpr173-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Gpr173-E11
#| fig-cap: "Correlation between Sox2 and Gpr173 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Sox2", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Gpr173 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Gpr173-E11-1.png){#fig-kim2020-correlation-Sox2-Gpr173-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Gpr173-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Gpr173-E12
#| fig-cap: "Correlation between Sox2 and Gpr173 expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Sox2", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Gpr173 expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Gpr173-E12-1.png){#fig-kim2020-correlation-Sox2-Gpr173-E12 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Gpr173-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Gpr173-E13
#| fig-cap: "Correlation between Sox2 and Gpr173 expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Sox2", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Gpr173 expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Gpr173-E13-1.png){#fig-kim2020-correlation-Sox2-Gpr173-E13 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Gpr173-E14 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Gpr173-E14
#| fig-cap: "Correlation between Sox2 and Gpr173 expression at E14 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E14", cells = cells_to_check),
    c("Sox2", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E14"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Gpr173 expression at E14 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Gpr173-E14-1.png){#fig-kim2020-correlation-Sox2-Gpr173-E14 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Gpr173-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Gpr173-E15
#| fig-cap: "Correlation between Sox2 and Gpr173 expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Sox2", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Gpr173 expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Gpr173-E15-1.png){#fig-kim2020-correlation-Sox2-Gpr173-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Gpr173-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Gpr173-E16
#| fig-cap: "Correlation between Sox2 and Gpr173 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Sox2", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Gpr173 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Gpr173-E16-1.png){#fig-kim2020-correlation-Sox2-Gpr173-E16 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Gpr173-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Gpr173-E18
#| fig-cap: "Correlation between Sox2 and Gpr173 expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Sox2", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Gpr173 expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Gpr173-E18-1.png){#fig-kim2020-correlation-Sox2-Gpr173-E18 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Gpr173-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Gpr173-P8
#| fig-cap: "Correlation between Sox2 and Gpr173 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Sox2", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Gpr173 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Gpr173-P8-1.png){#fig-kim2020-correlation-Sox2-Gpr173-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Gpr173-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Gpr173-P45
#| fig-cap: "Correlation between Sox2 and Gpr173 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Sox2", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Gpr173 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Gpr173-P45-1.png){#fig-kim2020-correlation-Sox2-Gpr173-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Lhx3-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Lhx3-E10
#| fig-cap: "Correlation between Sox2 and Lhx3 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Sox2", "Lhx3"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Lhx3 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Lhx3-E10-1.png){#fig-kim2020-correlation-Sox2-Lhx3-E10 width=1152}
:::
:::















::: {#cell-fig-kim2020-correlation-Sox2-Lhx3-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Lhx3-P8
#| fig-cap: "Correlation between Sox2 and Lhx3 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Sox2", "Lhx3"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Lhx3 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Lhx3-P8-1.png){#fig-kim2020-correlation-Sox2-Lhx3-P8 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Sox2-Lhx4-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Lhx4-E10
#| fig-cap: "Correlation between Sox2 and Lhx4 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Sox2", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Lhx4 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Lhx4-E10-1.png){#fig-kim2020-correlation-Sox2-Lhx4-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Lhx4-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Lhx4-E11
#| fig-cap: "Correlation between Sox2 and Lhx4 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Sox2", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Lhx4 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Lhx4-E11-1.png){#fig-kim2020-correlation-Sox2-Lhx4-E11 width=1152}
:::
:::









::: {#cell-fig-kim2020-correlation-Sox2-Lhx4-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Lhx4-E16
#| fig-cap: "Correlation between Sox2 and Lhx4 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Sox2", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Lhx4 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Lhx4-E16-1.png){#fig-kim2020-correlation-Sox2-Lhx4-E16 width=1152}
:::
:::







::: {#cell-fig-kim2020-correlation-Sox2-Gata2-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Gata2-E10
#| fig-cap: "Correlation between Sox2 and Gata2 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Sox2", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Gata2 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Gata2-E10-1.png){#fig-kim2020-correlation-Sox2-Gata2-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Gata2-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Gata2-E11
#| fig-cap: "Correlation between Sox2 and Gata2 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Sox2", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Gata2 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Gata2-E11-1.png){#fig-kim2020-correlation-Sox2-Gata2-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Gata2-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Gata2-E12
#| fig-cap: "Correlation between Sox2 and Gata2 expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Sox2", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Gata2 expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Gata2-E12-1.png){#fig-kim2020-correlation-Sox2-Gata2-E12 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Gata2-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Gata2-E13
#| fig-cap: "Correlation between Sox2 and Gata2 expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Sox2", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Gata2 expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Gata2-E13-1.png){#fig-kim2020-correlation-Sox2-Gata2-E13 width=1152}
:::
:::





::: {#cell-fig-kim2020-correlation-Sox2-Gata2-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Gata2-E16
#| fig-cap: "Correlation between Sox2 and Gata2 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Sox2", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Gata2 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Gata2-E16-1.png){#fig-kim2020-correlation-Sox2-Gata2-E16 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Sox2-Gata2-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Gata2-P8
#| fig-cap: "Correlation between Sox2 and Gata2 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Sox2", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Gata2 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Gata2-P8-1.png){#fig-kim2020-correlation-Sox2-Gata2-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Sox2-Gata2-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Sox2-Gata2-P45
#| fig-cap: "Correlation between Sox2 and Gata2 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Sox2", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Gata2 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Sox2-Gata2-P45-1.png){#fig-kim2020-correlation-Sox2-Gata2-P45 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Hlf-Tshr-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Tshr-E11
#| fig-cap: "Correlation between Hlf and Tshr expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Hlf", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Tshr expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Tshr-E11-1.png){#fig-kim2020-correlation-Hlf-Tshr-E11 width=1152}
:::
:::







::: {#cell-fig-kim2020-correlation-Hlf-Tshr-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Tshr-E15
#| fig-cap: "Correlation between Hlf and Tshr expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Hlf", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Tshr expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Tshr-E15-1.png){#fig-kim2020-correlation-Hlf-Tshr-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Hlf-Tshr-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Tshr-E16
#| fig-cap: "Correlation between Hlf and Tshr expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Hlf", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Tshr expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Tshr-E16-1.png){#fig-kim2020-correlation-Hlf-Tshr-E16 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Hlf-Tshr-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Tshr-E18
#| fig-cap: "Correlation between Hlf and Tshr expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Hlf", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Tshr expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Tshr-E18-1.png){#fig-kim2020-correlation-Hlf-Tshr-E18 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Hlf-Tshr-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Tshr-P8
#| fig-cap: "Correlation between Hlf and Tshr expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Hlf", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Tshr expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Tshr-P8-1.png){#fig-kim2020-correlation-Hlf-Tshr-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Hlf-Tshr-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Tshr-P45
#| fig-cap: "Correlation between Hlf and Tshr expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Hlf", "Tshr"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Tshr expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Tshr-P45-1.png){#fig-kim2020-correlation-Hlf-Tshr-P45 width=1152}
:::
:::











::: {#cell-fig-kim2020-correlation-Hlf-Cckar-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Cckar-E15
#| fig-cap: "Correlation between Hlf and Cckar expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Hlf", "Cckar"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Cckar expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Cckar-E15-1.png){#fig-kim2020-correlation-Hlf-Cckar-E15 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Hlf-Cckar-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Cckar-E18
#| fig-cap: "Correlation between Hlf and Cckar expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Hlf", "Cckar"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Cckar expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Cckar-E18-1.png){#fig-kim2020-correlation-Hlf-Cckar-E18 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Hlf-Cckar-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Cckar-P45
#| fig-cap: "Correlation between Hlf and Cckar expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Hlf", "Cckar"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Cckar expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Cckar-P45-1.png){#fig-kim2020-correlation-Hlf-Cckar-P45 width=1152}
:::
:::

















::: {#cell-fig-kim2020-correlation-Hlf-Cckbr-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Cckbr-P8
#| fig-cap: "Correlation between Hlf and Cckbr expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Hlf", "Cckbr"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Cckbr expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Cckbr-P8-1.png){#fig-kim2020-correlation-Hlf-Cckbr-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Hlf-Cckbr-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Cckbr-P45
#| fig-cap: "Correlation between Hlf and Cckbr expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Hlf", "Cckbr"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Cckbr expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Cckbr-P45-1.png){#fig-kim2020-correlation-Hlf-Cckbr-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Hlf-Gpr173-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Gpr173-E10
#| fig-cap: "Correlation between Hlf and Gpr173 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Hlf", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Gpr173 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Gpr173-E10-1.png){#fig-kim2020-correlation-Hlf-Gpr173-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Hlf-Gpr173-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Gpr173-E11
#| fig-cap: "Correlation between Hlf and Gpr173 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Hlf", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Gpr173 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Gpr173-E11-1.png){#fig-kim2020-correlation-Hlf-Gpr173-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Hlf-Gpr173-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Gpr173-E12
#| fig-cap: "Correlation between Hlf and Gpr173 expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Hlf", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Gpr173 expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Gpr173-E12-1.png){#fig-kim2020-correlation-Hlf-Gpr173-E12 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Hlf-Gpr173-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Gpr173-E13
#| fig-cap: "Correlation between Hlf and Gpr173 expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Hlf", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Gpr173 expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Gpr173-E13-1.png){#fig-kim2020-correlation-Hlf-Gpr173-E13 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Hlf-Gpr173-E14 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Gpr173-E14
#| fig-cap: "Correlation between Hlf and Gpr173 expression at E14 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E14", cells = cells_to_check),
    c("Hlf", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E14"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Gpr173 expression at E14 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Gpr173-E14-1.png){#fig-kim2020-correlation-Hlf-Gpr173-E14 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Hlf-Gpr173-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Gpr173-E15
#| fig-cap: "Correlation between Hlf and Gpr173 expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Hlf", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Gpr173 expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Gpr173-E15-1.png){#fig-kim2020-correlation-Hlf-Gpr173-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Hlf-Gpr173-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Gpr173-E16
#| fig-cap: "Correlation between Hlf and Gpr173 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Hlf", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Gpr173 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Gpr173-E16-1.png){#fig-kim2020-correlation-Hlf-Gpr173-E16 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Hlf-Gpr173-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Gpr173-E18
#| fig-cap: "Correlation between Hlf and Gpr173 expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Hlf", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Gpr173 expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Gpr173-E18-1.png){#fig-kim2020-correlation-Hlf-Gpr173-E18 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Hlf-Gpr173-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Gpr173-P8
#| fig-cap: "Correlation between Hlf and Gpr173 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Hlf", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Gpr173 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Gpr173-P8-1.png){#fig-kim2020-correlation-Hlf-Gpr173-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Hlf-Gpr173-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Gpr173-P45
#| fig-cap: "Correlation between Hlf and Gpr173 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Hlf", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Gpr173 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Gpr173-P45-1.png){#fig-kim2020-correlation-Hlf-Gpr173-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Hlf-Lhx3-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Lhx3-E10
#| fig-cap: "Correlation between Hlf and Lhx3 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Hlf", "Lhx3"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Lhx3 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Lhx3-E10-1.png){#fig-kim2020-correlation-Hlf-Lhx3-E10 width=1152}
:::
:::















::: {#cell-fig-kim2020-correlation-Hlf-Lhx3-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Lhx3-P8
#| fig-cap: "Correlation between Hlf and Lhx3 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Hlf", "Lhx3"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Lhx3 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Lhx3-P8-1.png){#fig-kim2020-correlation-Hlf-Lhx3-P8 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Hlf-Lhx4-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Lhx4-E10
#| fig-cap: "Correlation between Hlf and Lhx4 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Hlf", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Lhx4 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Lhx4-E10-1.png){#fig-kim2020-correlation-Hlf-Lhx4-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Hlf-Lhx4-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Lhx4-E11
#| fig-cap: "Correlation between Hlf and Lhx4 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Hlf", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Lhx4 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Lhx4-E11-1.png){#fig-kim2020-correlation-Hlf-Lhx4-E11 width=1152}
:::
:::









::: {#cell-fig-kim2020-correlation-Hlf-Lhx4-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Lhx4-E16
#| fig-cap: "Correlation between Hlf and Lhx4 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Hlf", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Lhx4 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Lhx4-E16-1.png){#fig-kim2020-correlation-Hlf-Lhx4-E16 width=1152}
:::
:::







::: {#cell-fig-kim2020-correlation-Hlf-Gata2-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Gata2-E10
#| fig-cap: "Correlation between Hlf and Gata2 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Hlf", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Gata2 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Gata2-E10-1.png){#fig-kim2020-correlation-Hlf-Gata2-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Hlf-Gata2-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Gata2-E11
#| fig-cap: "Correlation between Hlf and Gata2 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Hlf", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Gata2 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Gata2-E11-1.png){#fig-kim2020-correlation-Hlf-Gata2-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Hlf-Gata2-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Gata2-E12
#| fig-cap: "Correlation between Hlf and Gata2 expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Hlf", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Gata2 expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Gata2-E12-1.png){#fig-kim2020-correlation-Hlf-Gata2-E12 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Hlf-Gata2-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Gata2-E13
#| fig-cap: "Correlation between Hlf and Gata2 expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Hlf", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Gata2 expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Gata2-E13-1.png){#fig-kim2020-correlation-Hlf-Gata2-E13 width=1152}
:::
:::





::: {#cell-fig-kim2020-correlation-Hlf-Gata2-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Gata2-E16
#| fig-cap: "Correlation between Hlf and Gata2 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Hlf", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Gata2 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Gata2-E16-1.png){#fig-kim2020-correlation-Hlf-Gata2-E16 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Hlf-Gata2-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Gata2-P8
#| fig-cap: "Correlation between Hlf and Gata2 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Hlf", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Gata2 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Gata2-P8-1.png){#fig-kim2020-correlation-Hlf-Gata2-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Hlf-Gata2-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Hlf-Gata2-P45
#| fig-cap: "Correlation between Hlf and Gata2 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Hlf", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Gata2 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Hlf-Gata2-P45-1.png){#fig-kim2020-correlation-Hlf-Gata2-P45 width=1152}
:::
:::











::: {#cell-fig-kim2020-correlation-Tshr-Cckar-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshr-Cckar-E15
#| fig-cap: "Correlation between Tshr and Cckar expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Tshr", "Cckar"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Cckar expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshr-Cckar-E15-1.png){#fig-kim2020-correlation-Tshr-Cckar-E15 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Tshr-Cckar-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshr-Cckar-E18
#| fig-cap: "Correlation between Tshr and Cckar expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Tshr", "Cckar"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Cckar expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshr-Cckar-E18-1.png){#fig-kim2020-correlation-Tshr-Cckar-E18 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Tshr-Cckar-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshr-Cckar-P45
#| fig-cap: "Correlation between Tshr and Cckar expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Tshr", "Cckar"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Cckar expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshr-Cckar-P45-1.png){#fig-kim2020-correlation-Tshr-Cckar-P45 width=1152}
:::
:::

















::: {#cell-fig-kim2020-correlation-Tshr-Cckbr-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshr-Cckbr-P8
#| fig-cap: "Correlation between Tshr and Cckbr expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Tshr", "Cckbr"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Cckbr expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshr-Cckbr-P8-1.png){#fig-kim2020-correlation-Tshr-Cckbr-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Tshr-Cckbr-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshr-Cckbr-P45
#| fig-cap: "Correlation between Tshr and Cckbr expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Tshr", "Cckbr"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Cckbr expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshr-Cckbr-P45-1.png){#fig-kim2020-correlation-Tshr-Cckbr-P45 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Tshr-Gpr173-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshr-Gpr173-E11
#| fig-cap: "Correlation between Tshr and Gpr173 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Tshr", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Gpr173 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshr-Gpr173-E11-1.png){#fig-kim2020-correlation-Tshr-Gpr173-E11 width=1152}
:::
:::







::: {#cell-fig-kim2020-correlation-Tshr-Gpr173-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshr-Gpr173-E15
#| fig-cap: "Correlation between Tshr and Gpr173 expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Tshr", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Gpr173 expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshr-Gpr173-E15-1.png){#fig-kim2020-correlation-Tshr-Gpr173-E15 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Tshr-Gpr173-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshr-Gpr173-E16
#| fig-cap: "Correlation between Tshr and Gpr173 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Tshr", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Gpr173 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshr-Gpr173-E16-1.png){#fig-kim2020-correlation-Tshr-Gpr173-E16 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Tshr-Gpr173-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshr-Gpr173-E18
#| fig-cap: "Correlation between Tshr and Gpr173 expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Tshr", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Gpr173 expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshr-Gpr173-E18-1.png){#fig-kim2020-correlation-Tshr-Gpr173-E18 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Tshr-Gpr173-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshr-Gpr173-P8
#| fig-cap: "Correlation between Tshr and Gpr173 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Tshr", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Gpr173 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshr-Gpr173-P8-1.png){#fig-kim2020-correlation-Tshr-Gpr173-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Tshr-Gpr173-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshr-Gpr173-P45
#| fig-cap: "Correlation between Tshr and Gpr173 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Tshr", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Gpr173 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshr-Gpr173-P45-1.png){#fig-kim2020-correlation-Tshr-Gpr173-P45 width=1152}
:::
:::

















::: {#cell-fig-kim2020-correlation-Tshr-Lhx3-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshr-Lhx3-P8
#| fig-cap: "Correlation between Tshr and Lhx3 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Tshr", "Lhx3"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Lhx3 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshr-Lhx3-P8-1.png){#fig-kim2020-correlation-Tshr-Lhx3-P8 width=1152}
:::
:::





::: {#cell-fig-kim2020-correlation-Tshr-Lhx4-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshr-Lhx4-E11
#| fig-cap: "Correlation between Tshr and Lhx4 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Tshr", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Lhx4 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshr-Lhx4-E11-1.png){#fig-kim2020-correlation-Tshr-Lhx4-E11 width=1152}
:::
:::









::: {#cell-fig-kim2020-correlation-Tshr-Lhx4-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshr-Lhx4-E16
#| fig-cap: "Correlation between Tshr and Lhx4 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Tshr", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Lhx4 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshr-Lhx4-E16-1.png){#fig-kim2020-correlation-Tshr-Lhx4-E16 width=1152}
:::
:::









::: {#cell-fig-kim2020-correlation-Tshr-Gata2-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshr-Gata2-E11
#| fig-cap: "Correlation between Tshr and Gata2 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Tshr", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Gata2 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshr-Gata2-E11-1.png){#fig-kim2020-correlation-Tshr-Gata2-E11 width=1152}
:::
:::









::: {#cell-fig-kim2020-correlation-Tshr-Gata2-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshr-Gata2-E16
#| fig-cap: "Correlation between Tshr and Gata2 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Tshr", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Gata2 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshr-Gata2-E16-1.png){#fig-kim2020-correlation-Tshr-Gata2-E16 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Tshr-Gata2-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshr-Gata2-P8
#| fig-cap: "Correlation between Tshr and Gata2 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Tshr", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Gata2 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshr-Gata2-P8-1.png){#fig-kim2020-correlation-Tshr-Gata2-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Tshr-Gata2-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Tshr-Gata2-P45
#| fig-cap: "Correlation between Tshr and Gata2 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Tshr", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Gata2 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Tshr-Gata2-P45-1.png){#fig-kim2020-correlation-Tshr-Gata2-P45 width=1152}
:::
:::



















::: {#cell-fig-kim2020-correlation-Cckar-Cckbr-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cckar-Cckbr-P45
#| fig-cap: "Correlation between Cckar and Cckbr expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Cckar", "Cckbr"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cckar and Cckbr expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cckar-Cckbr-P45-1.png){#fig-kim2020-correlation-Cckar-Cckbr-P45 width=1152}
:::
:::











::: {#cell-fig-kim2020-correlation-Cckar-Gpr173-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cckar-Gpr173-E15
#| fig-cap: "Correlation between Cckar and Gpr173 expression at E15 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E15", cells = cells_to_check),
    c("Cckar", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cckar and Gpr173 expression at E15 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cckar-Gpr173-E15-1.png){#fig-kim2020-correlation-Cckar-Gpr173-E15 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Cckar-Gpr173-E18 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cckar-Gpr173-E18
#| fig-cap: "Correlation between Cckar and Gpr173 expression at E18 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E18", cells = cells_to_check),
    c("Cckar", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "E18"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cckar and Gpr173 expression at E18 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cckar-Gpr173-E18-1.png){#fig-kim2020-correlation-Cckar-Gpr173-E18 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Cckar-Gpr173-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cckar-Gpr173-P45
#| fig-cap: "Correlation between Cckar and Gpr173 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Cckar", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cckar and Gpr173 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cckar-Gpr173-P45-1.png){#fig-kim2020-correlation-Cckar-Gpr173-P45 width=1152}
:::
:::



























































::: {#cell-fig-kim2020-correlation-Cckar-Gata2-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cckar-Gata2-P45
#| fig-cap: "Correlation between Cckar and Gata2 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Cckar", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cckar and Gata2 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cckar-Gata2-P45-1.png){#fig-kim2020-correlation-Cckar-Gata2-P45 width=1152}
:::
:::

















::: {#cell-fig-kim2020-correlation-Cckbr-Gpr173-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cckbr-Gpr173-P8
#| fig-cap: "Correlation between Cckbr and Gpr173 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Cckbr", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cckbr and Gpr173 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cckbr-Gpr173-P8-1.png){#fig-kim2020-correlation-Cckbr-Gpr173-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cckbr-Gpr173-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cckbr-Gpr173-P45
#| fig-cap: "Correlation between Cckbr and Gpr173 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Cckbr", "Gpr173"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cckbr and Gpr173 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cckbr-Gpr173-P45-1.png){#fig-kim2020-correlation-Cckbr-Gpr173-P45 width=1152}
:::
:::

















::: {#cell-fig-kim2020-correlation-Cckbr-Lhx3-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cckbr-Lhx3-P8
#| fig-cap: "Correlation between Cckbr and Lhx3 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Cckbr", "Lhx3"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cckbr and Lhx3 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cckbr-Lhx3-P8-1.png){#fig-kim2020-correlation-Cckbr-Lhx3-P8 width=1152}
:::
:::







































::: {#cell-fig-kim2020-correlation-Cckbr-Gata2-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cckbr-Gata2-P8
#| fig-cap: "Correlation between Cckbr and Gata2 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Cckbr", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cckbr and Gata2 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cckbr-Gata2-P8-1.png){#fig-kim2020-correlation-Cckbr-Gata2-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Cckbr-Gata2-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Cckbr-Gata2-P45
#| fig-cap: "Correlation between Cckbr and Gata2 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Cckbr", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cckbr and Gata2 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Cckbr-Gata2-P45-1.png){#fig-kim2020-correlation-Cckbr-Gata2-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Gpr173-Lhx3-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Gpr173-Lhx3-E10
#| fig-cap: "Correlation between Gpr173 and Lhx3 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Gpr173", "Lhx3"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Gpr173 and Lhx3 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Gpr173-Lhx3-E10-1.png){#fig-kim2020-correlation-Gpr173-Lhx3-E10 width=1152}
:::
:::















::: {#cell-fig-kim2020-correlation-Gpr173-Lhx3-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Gpr173-Lhx3-P8
#| fig-cap: "Correlation between Gpr173 and Lhx3 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Gpr173", "Lhx3"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Gpr173 and Lhx3 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Gpr173-Lhx3-P8-1.png){#fig-kim2020-correlation-Gpr173-Lhx3-P8 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Gpr173-Lhx4-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Gpr173-Lhx4-E10
#| fig-cap: "Correlation between Gpr173 and Lhx4 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Gpr173", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Gpr173 and Lhx4 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Gpr173-Lhx4-E10-1.png){#fig-kim2020-correlation-Gpr173-Lhx4-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Gpr173-Lhx4-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Gpr173-Lhx4-E11
#| fig-cap: "Correlation between Gpr173 and Lhx4 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Gpr173", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Gpr173 and Lhx4 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Gpr173-Lhx4-E11-1.png){#fig-kim2020-correlation-Gpr173-Lhx4-E11 width=1152}
:::
:::









::: {#cell-fig-kim2020-correlation-Gpr173-Lhx4-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Gpr173-Lhx4-E16
#| fig-cap: "Correlation between Gpr173 and Lhx4 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Gpr173", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Gpr173 and Lhx4 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Gpr173-Lhx4-E16-1.png){#fig-kim2020-correlation-Gpr173-Lhx4-E16 width=1152}
:::
:::







::: {#cell-fig-kim2020-correlation-Gpr173-Gata2-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Gpr173-Gata2-E10
#| fig-cap: "Correlation between Gpr173 and Gata2 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Gpr173", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Gpr173 and Gata2 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Gpr173-Gata2-E10-1.png){#fig-kim2020-correlation-Gpr173-Gata2-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Gpr173-Gata2-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Gpr173-Gata2-E11
#| fig-cap: "Correlation between Gpr173 and Gata2 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Gpr173", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Gpr173 and Gata2 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Gpr173-Gata2-E11-1.png){#fig-kim2020-correlation-Gpr173-Gata2-E11 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Gpr173-Gata2-E12 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Gpr173-Gata2-E12
#| fig-cap: "Correlation between Gpr173 and Gata2 expression at E12 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E12", cells = cells_to_check),
    c("Gpr173", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E12"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Gpr173 and Gata2 expression at E12 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Gpr173-Gata2-E12-1.png){#fig-kim2020-correlation-Gpr173-Gata2-E12 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Gpr173-Gata2-E13 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Gpr173-Gata2-E13
#| fig-cap: "Correlation between Gpr173 and Gata2 expression at E13 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E13", cells = cells_to_check),
    c("Gpr173", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E13"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Gpr173 and Gata2 expression at E13 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Gpr173-Gata2-E13-1.png){#fig-kim2020-correlation-Gpr173-Gata2-E13 width=1152}
:::
:::





::: {#cell-fig-kim2020-correlation-Gpr173-Gata2-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Gpr173-Gata2-E16
#| fig-cap: "Correlation between Gpr173 and Gata2 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Gpr173", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Gpr173 and Gata2 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Gpr173-Gata2-E16-1.png){#fig-kim2020-correlation-Gpr173-Gata2-E16 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Gpr173-Gata2-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Gpr173-Gata2-P8
#| fig-cap: "Correlation between Gpr173 and Gata2 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Gpr173", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Gpr173 and Gata2 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Gpr173-Gata2-P8-1.png){#fig-kim2020-correlation-Gpr173-Gata2-P8 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Gpr173-Gata2-P45 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Gpr173-Gata2-P45
#| fig-cap: "Correlation between Gpr173 and Gata2 expression at P45 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P45", cells = cells_to_check),
    c("Gpr173", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "P45"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Gpr173 and Gata2 expression at P45 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Gpr173-Gata2-P45-1.png){#fig-kim2020-correlation-Gpr173-Gata2-P45 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Lhx3-Lhx4-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Lhx3-Lhx4-E10
#| fig-cap: "Correlation between Lhx3 and Lhx4 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Lhx3", "Lhx4"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Lhx3 and Lhx4 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Lhx3-Lhx4-E10-1.png){#fig-kim2020-correlation-Lhx3-Lhx4-E10 width=1152}
:::
:::



















::: {#cell-fig-kim2020-correlation-Lhx3-Gata2-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Lhx3-Gata2-E10
#| fig-cap: "Correlation between Lhx3 and Gata2 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Lhx3", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Lhx3 and Gata2 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Lhx3-Gata2-E10-1.png){#fig-kim2020-correlation-Lhx3-Gata2-E10 width=1152}
:::
:::















::: {#cell-fig-kim2020-correlation-Lhx3-Gata2-P8 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Lhx3-Gata2-P8
#| fig-cap: "Correlation between Lhx3 and Gata2 expression at P8 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "P8", cells = cells_to_check),
    c("Lhx3", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "P8"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Lhx3 and Gata2 expression at P8 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Lhx3-Gata2-P8-1.png){#fig-kim2020-correlation-Lhx3-Gata2-P8 width=1152}
:::
:::



::: {#cell-fig-kim2020-correlation-Lhx4-Gata2-E10 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Lhx4-Gata2-E10
#| fig-cap: "Correlation between Lhx4 and Gata2 expression at E10 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E10", cells = cells_to_check),
    c("Lhx4", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Lhx4 and Gata2 expression at E10 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Lhx4-Gata2-E10-1.png){#fig-kim2020-correlation-Lhx4-Gata2-E10 width=1152}
:::
:::

::: {#cell-fig-kim2020-correlation-Lhx4-Gata2-E11 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Lhx4-Gata2-E11
#| fig-cap: "Correlation between Lhx4 and Gata2 expression at E11 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E11", cells = cells_to_check),
    c("Lhx4", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E11"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Lhx4 and Gata2 expression at E11 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Lhx4-Gata2-E11-1.png){#fig-kim2020-correlation-Lhx4-Gata2-E11 width=1152}
:::
:::









::: {#cell-fig-kim2020-correlation-Lhx4-Gata2-E16 .cell}

```{.r .cell-code .hidden}
#| label: fig-kim2020-correlation-Lhx4-Gata2-E16
#| fig-cap: "Correlation between Lhx4 and Gata2 expression at E16 in Kim et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.kim |> subset(Age == "E16", cells = cells_to_check),
    c("Lhx4", "Gata2"),
    "predicted.id",
    cells_to_check,
    age = "E16"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Lhx4 and Gata2 expression at E16 in Kim et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-kim2020-correlation-Lhx4-Gata2-E16-1.png){#fig-kim2020-correlation-Lhx4-Gata2-E16 width=1152}
:::
:::












#### Dotplots





::: {#cell-fig-dotplot-dendrogram-genes-kim2020 .cell}

```{.r .cell-code .hidden}
#| label: fig-dotplot-dendrogram-genes-kim2020
#| fig-cap: "Dotplot of selected genes in hypothalamus across different developmental stages in the Kim DW et al. (2020) dataset. Cells are colored by expression level.  Note the distinct localization patterns of each gene."
#| fig-width: 6
#| fig-height: 6
Idents(srt.kim) <- "Age"

DotPlot_scCustom(seurat_object = subset(srt.kim, cells = cells_to_check), colors_use = viridis(n = 30, alpha = .75, direction = -1, option = "E"), features = genes_to_check, flip_axes = T, x_lab_rotate = TRUE, dot.scale = 15)
```

::: {.cell-output-display}
![Dotplot of selected genes in hypothalamus across different developmental stages in the Kim DW et al. (2020) dataset. Cells are colored by expression level.  Note the distinct localization patterns of each gene.](01-de_test-focus_pars_tub_files/figure-html/fig-dotplot-dendrogram-genes-kim2020-1.png){#fig-dotplot-dendrogram-genes-kim2020 width=576}
:::
:::





### Romanov et al. 2020 dataset





::: {#cell-fig-violin-gene-interactions-romanov2020 .cell}

```{.r .cell-code .hidden}
#| label: fig-violin-gene-interactions-romanov2020
#| fig-cap: Gene expression of projected Pars Tuberalis clusters in hypothalamus across different developmental stages in the Romanov et al. (2020) dataset.
#| fig-width: 12
#| fig-height: 12
# Violin plots for srt.romanov.pub (split by its original clusters, 'wtree')
Idents(srt.romanov.pub) <- "wtree"
cells_to_check <- WhichCells(srt.romanov.pub, idents = c("38", "42", "45") %>% .[. %in% levels(srt.romanov.pub)])
dat_romanov <- srt.romanov.pub@assays$RNA@counts[, cells_to_check]
genes_to_check <- c(
  "Tshb",
  "Cck",
  "Pitx1",
  "Eya1",
  "Eya2",
  "Eya3",
  "Eya4",
  "Sox2",
  "Hlf",
  "Tshr",
  "Cckar",
  "Cckbr",
  "Gpr173",
  "Foxl2",
  "Lhx3",
  "Lhx4",
  "Pit1",
  "Gata2"
) %>%
  .[. %in% (srt.romanov.pub@assays$RNA@counts[, cells_to_check] |> rowSums() %>% .[. > 3] %>% names())]


dat_romanov <- dat_romanov[genes_to_check, cells_to_check] %>% as.data.frame() %>% t()

dat_romanov <- dat_romanov |> as.data.frame() |> mutate(Age = srt.romanov.pub[["Age"]][rownames(dat_romanov), ])

# Complete gene pairs
gene_pairs <- combn(
  genes_to_check, 2
) |>
  as_tibble() |>
  array_tree(margin = 2)

ages <- srt.romanov.pub$Age |>
  table() %>%
  .[. > 0] %>%
  names()

# Violin plots of individual gene expression split by cluster
Idents(srt.romanov.pub) <- "Age" # Reset identity back to Age

Stacked_VlnPlot(subset(srt.romanov.pub, cells = cells_to_check), features = genes_to_check, split.by = "wtree", colors_use = srt.romanov.pub@misc$wtree_Colour_Pal[c("38", "42", "45")])
```

::: {.cell-output .cell-output-stderr .hidden}

```
Warning: Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
Groups with fewer than two datapoints have been dropped.
ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.
```


:::

::: {.cell-output-display}
![Gene expression of projected Pars Tuberalis clusters in hypothalamus across different developmental stages in the Romanov et al. (2020) dataset.](01-de_test-focus_pars_tub_files/figure-html/fig-violin-gene-interactions-romanov2020-1.png){#fig-violin-gene-interactions-romanov2020 width=1152}
:::
:::

::: {.cell}

````{.python .cell-code .hidden}
# | label: generate-romanov2020-correlation-plots
# | output: false
# | echo: false

import pandas as pd

# Access R variables in Python using py$
gene_pairs = r.gene_pairs
ages = r.ages
dat_romanov = r.dat_romanov  # Assuming you have dat_romanov in R

# Convert gene_pairs and data to a Pandas DataFrame
gene_pairs = pd.DataFrame(gene_pairs).T
dat_romanov = pd.DataFrame(dat_romanov)


def generate_plot_chunk(genes, age):
    # Filter the data for the current age
    age_data = dat_romanov[dat_romanov["Age"] == age]

    # Check if both genes have at least two unique values in the age subset
    gene1_values = age_data[genes[0]].unique()
    gene2_values = age_data[genes[1]].unique()

    if len(gene1_values) >= 2 and len(gene2_values) >= 2:
        # Generate the plot chunk
        chunk_text = f"""```{{r}}
#| label: fig-romanov2020-correlation-{genes[0]}-{genes[1]}-{age}
#| fig-cap: "Correlation between {genes[0]} and {genes[1]} expression at {age} in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "{age}", cells = cells_to_check),
    c("{genes[0]}", "{genes[1]}"),
    "wtree",
    cells_to_check,
    age = "{age}"
)
```"""
        return chunk_text
    else:
        # Skip this gene pair for this age if one or both genes don't have enough values
        return ""


# Generate the Quarto markdown content
qmd_content = "#### Gene Correlation Plots\n\n"

for _, row in gene_pairs.iterrows():
    for age in ages:
        chunk = generate_plot_chunk([row[0], row[1]], age)
        if chunk:
            qmd_content += chunk + "\n\n"

# Write the content to a .qmd file
with open("gene_correlation_plots_romanov2020.qmd", "w") as f:
    f.write(qmd_content)
````

::: {.cell-output .cell-output-stdout}

```
78459
```


:::
:::


#### Gene Correlation Plots

::: {#cell-fig-romanov2020-correlation-Tshb-Cck-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Cck-P2
#| fig-cap: "Correlation between Tshb and Cck expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Tshb", "Cck"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Cck expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Cck-P2-1.png){#fig-romanov2020-correlation-Tshb-Cck-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Cck-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Cck-P10
#| fig-cap: "Correlation between Tshb and Cck expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Tshb", "Cck"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Cck expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Cck-P10-1.png){#fig-romanov2020-correlation-Tshb-Cck-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Cck-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Cck-P23
#| fig-cap: "Correlation between Tshb and Cck expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Tshb", "Cck"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Cck expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Cck-P23-1.png){#fig-romanov2020-correlation-Tshb-Cck-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Pitx1-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Pitx1-P2
#| fig-cap: "Correlation between Tshb and Pitx1 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Tshb", "Pitx1"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Pitx1 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Pitx1-P2-1.png){#fig-romanov2020-correlation-Tshb-Pitx1-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Pitx1-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Pitx1-P10
#| fig-cap: "Correlation between Tshb and Pitx1 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Tshb", "Pitx1"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Pitx1 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Pitx1-P10-1.png){#fig-romanov2020-correlation-Tshb-Pitx1-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Pitx1-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Pitx1-P23
#| fig-cap: "Correlation between Tshb and Pitx1 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Tshb", "Pitx1"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Pitx1 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Pitx1-P23-1.png){#fig-romanov2020-correlation-Tshb-Pitx1-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Eya1-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Eya1-P2
#| fig-cap: "Correlation between Tshb and Eya1 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Tshb", "Eya1"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Eya1 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Eya1-P2-1.png){#fig-romanov2020-correlation-Tshb-Eya1-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Eya1-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Eya1-P10
#| fig-cap: "Correlation between Tshb and Eya1 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Tshb", "Eya1"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Eya1 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Eya1-P10-1.png){#fig-romanov2020-correlation-Tshb-Eya1-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Eya1-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Eya1-P23
#| fig-cap: "Correlation between Tshb and Eya1 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Tshb", "Eya1"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Eya1 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Eya1-P23-1.png){#fig-romanov2020-correlation-Tshb-Eya1-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Eya2-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Eya2-P2
#| fig-cap: "Correlation between Tshb and Eya2 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Tshb", "Eya2"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Eya2 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Eya2-P2-1.png){#fig-romanov2020-correlation-Tshb-Eya2-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Eya2-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Eya2-P23
#| fig-cap: "Correlation between Tshb and Eya2 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Tshb", "Eya2"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Eya2 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Eya2-P23-1.png){#fig-romanov2020-correlation-Tshb-Eya2-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Eya3-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Eya3-P2
#| fig-cap: "Correlation between Tshb and Eya3 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Tshb", "Eya3"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Eya3 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Eya3-P2-1.png){#fig-romanov2020-correlation-Tshb-Eya3-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Eya3-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Eya3-P23
#| fig-cap: "Correlation between Tshb and Eya3 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Tshb", "Eya3"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Eya3 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Eya3-P23-1.png){#fig-romanov2020-correlation-Tshb-Eya3-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Eya4-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Eya4-P2
#| fig-cap: "Correlation between Tshb and Eya4 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Tshb", "Eya4"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Eya4 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Eya4-P2-1.png){#fig-romanov2020-correlation-Tshb-Eya4-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Eya4-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Eya4-P10
#| fig-cap: "Correlation between Tshb and Eya4 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Tshb", "Eya4"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Eya4 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Eya4-P10-1.png){#fig-romanov2020-correlation-Tshb-Eya4-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Eya4-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Eya4-P23
#| fig-cap: "Correlation between Tshb and Eya4 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Tshb", "Eya4"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Eya4 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Eya4-P23-1.png){#fig-romanov2020-correlation-Tshb-Eya4-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Sox2-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Sox2-P2
#| fig-cap: "Correlation between Tshb and Sox2 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Tshb", "Sox2"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Sox2 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Sox2-P2-1.png){#fig-romanov2020-correlation-Tshb-Sox2-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Sox2-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Sox2-P10
#| fig-cap: "Correlation between Tshb and Sox2 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Tshb", "Sox2"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Sox2 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Sox2-P10-1.png){#fig-romanov2020-correlation-Tshb-Sox2-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Sox2-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Sox2-P23
#| fig-cap: "Correlation between Tshb and Sox2 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Tshb", "Sox2"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Sox2 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Sox2-P23-1.png){#fig-romanov2020-correlation-Tshb-Sox2-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Hlf-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Hlf-P2
#| fig-cap: "Correlation between Tshb and Hlf expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Tshb", "Hlf"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Hlf expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Hlf-P2-1.png){#fig-romanov2020-correlation-Tshb-Hlf-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Hlf-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Hlf-P10
#| fig-cap: "Correlation between Tshb and Hlf expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Tshb", "Hlf"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Hlf expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Hlf-P10-1.png){#fig-romanov2020-correlation-Tshb-Hlf-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Hlf-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Hlf-P23
#| fig-cap: "Correlation between Tshb and Hlf expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Tshb", "Hlf"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Hlf expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Hlf-P23-1.png){#fig-romanov2020-correlation-Tshb-Hlf-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Tshr-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Tshr-P2
#| fig-cap: "Correlation between Tshb and Tshr expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Tshb", "Tshr"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Tshr expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Tshr-P2-1.png){#fig-romanov2020-correlation-Tshb-Tshr-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Tshr-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Tshr-P10
#| fig-cap: "Correlation between Tshb and Tshr expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Tshb", "Tshr"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Tshr expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Tshr-P10-1.png){#fig-romanov2020-correlation-Tshb-Tshr-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Tshr-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Tshr-P23
#| fig-cap: "Correlation between Tshb and Tshr expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Tshb", "Tshr"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Tshr expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Tshr-P23-1.png){#fig-romanov2020-correlation-Tshb-Tshr-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Gpr173-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Gpr173-P2
#| fig-cap: "Correlation between Tshb and Gpr173 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Tshb", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Gpr173 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Gpr173-P2-1.png){#fig-romanov2020-correlation-Tshb-Gpr173-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Gpr173-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Gpr173-P10
#| fig-cap: "Correlation between Tshb and Gpr173 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Tshb", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Gpr173 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Gpr173-P10-1.png){#fig-romanov2020-correlation-Tshb-Gpr173-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Gpr173-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Gpr173-P23
#| fig-cap: "Correlation between Tshb and Gpr173 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Tshb", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Gpr173 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Gpr173-P23-1.png){#fig-romanov2020-correlation-Tshb-Gpr173-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Lhx3-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Lhx3-P2
#| fig-cap: "Correlation between Tshb and Lhx3 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Tshb", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Lhx3 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Lhx3-P2-1.png){#fig-romanov2020-correlation-Tshb-Lhx3-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Lhx3-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Lhx3-P10
#| fig-cap: "Correlation between Tshb and Lhx3 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Tshb", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Lhx3 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Lhx3-P10-1.png){#fig-romanov2020-correlation-Tshb-Lhx3-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Lhx3-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Lhx3-P23
#| fig-cap: "Correlation between Tshb and Lhx3 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Tshb", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Lhx3 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Lhx3-P23-1.png){#fig-romanov2020-correlation-Tshb-Lhx3-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Gata2-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Gata2-P2
#| fig-cap: "Correlation between Tshb and Gata2 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Tshb", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Gata2 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Gata2-P2-1.png){#fig-romanov2020-correlation-Tshb-Gata2-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Gata2-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Gata2-P10
#| fig-cap: "Correlation between Tshb and Gata2 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Tshb", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Gata2 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Gata2-P10-1.png){#fig-romanov2020-correlation-Tshb-Gata2-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshb-Gata2-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshb-Gata2-P23
#| fig-cap: "Correlation between Tshb and Gata2 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Tshb", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshb and Gata2 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshb-Gata2-P23-1.png){#fig-romanov2020-correlation-Tshb-Gata2-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Pitx1-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Pitx1-P2
#| fig-cap: "Correlation between Cck and Pitx1 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Cck", "Pitx1"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Pitx1 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Pitx1-P2-1.png){#fig-romanov2020-correlation-Cck-Pitx1-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Pitx1-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Pitx1-P10
#| fig-cap: "Correlation between Cck and Pitx1 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Cck", "Pitx1"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Pitx1 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Pitx1-P10-1.png){#fig-romanov2020-correlation-Cck-Pitx1-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Pitx1-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Pitx1-P23
#| fig-cap: "Correlation between Cck and Pitx1 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Cck", "Pitx1"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Pitx1 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Pitx1-P23-1.png){#fig-romanov2020-correlation-Cck-Pitx1-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Eya1-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Eya1-P2
#| fig-cap: "Correlation between Cck and Eya1 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Cck", "Eya1"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya1 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Eya1-P2-1.png){#fig-romanov2020-correlation-Cck-Eya1-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Eya1-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Eya1-P10
#| fig-cap: "Correlation between Cck and Eya1 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Cck", "Eya1"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya1 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Eya1-P10-1.png){#fig-romanov2020-correlation-Cck-Eya1-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Eya1-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Eya1-P23
#| fig-cap: "Correlation between Cck and Eya1 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Cck", "Eya1"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya1 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Eya1-P23-1.png){#fig-romanov2020-correlation-Cck-Eya1-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Eya2-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Eya2-P2
#| fig-cap: "Correlation between Cck and Eya2 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Cck", "Eya2"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya2 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Eya2-P2-1.png){#fig-romanov2020-correlation-Cck-Eya2-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Eya2-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Eya2-P23
#| fig-cap: "Correlation between Cck and Eya2 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Cck", "Eya2"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya2 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Eya2-P23-1.png){#fig-romanov2020-correlation-Cck-Eya2-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Eya3-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Eya3-P2
#| fig-cap: "Correlation between Cck and Eya3 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Cck", "Eya3"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya3 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Eya3-P2-1.png){#fig-romanov2020-correlation-Cck-Eya3-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Eya3-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Eya3-P23
#| fig-cap: "Correlation between Cck and Eya3 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Cck", "Eya3"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya3 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Eya3-P23-1.png){#fig-romanov2020-correlation-Cck-Eya3-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Eya4-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Eya4-P2
#| fig-cap: "Correlation between Cck and Eya4 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Cck", "Eya4"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya4 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Eya4-P2-1.png){#fig-romanov2020-correlation-Cck-Eya4-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Eya4-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Eya4-P10
#| fig-cap: "Correlation between Cck and Eya4 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Cck", "Eya4"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya4 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Eya4-P10-1.png){#fig-romanov2020-correlation-Cck-Eya4-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Eya4-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Eya4-P23
#| fig-cap: "Correlation between Cck and Eya4 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Cck", "Eya4"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Eya4 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Eya4-P23-1.png){#fig-romanov2020-correlation-Cck-Eya4-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Sox2-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Sox2-P2
#| fig-cap: "Correlation between Cck and Sox2 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Cck", "Sox2"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Sox2 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Sox2-P2-1.png){#fig-romanov2020-correlation-Cck-Sox2-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Sox2-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Sox2-P10
#| fig-cap: "Correlation between Cck and Sox2 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Cck", "Sox2"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Sox2 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Sox2-P10-1.png){#fig-romanov2020-correlation-Cck-Sox2-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Sox2-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Sox2-P23
#| fig-cap: "Correlation between Cck and Sox2 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Cck", "Sox2"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Sox2 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Sox2-P23-1.png){#fig-romanov2020-correlation-Cck-Sox2-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Hlf-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Hlf-P2
#| fig-cap: "Correlation between Cck and Hlf expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Cck", "Hlf"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Hlf expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Hlf-P2-1.png){#fig-romanov2020-correlation-Cck-Hlf-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Hlf-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Hlf-P10
#| fig-cap: "Correlation between Cck and Hlf expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Cck", "Hlf"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Hlf expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Hlf-P10-1.png){#fig-romanov2020-correlation-Cck-Hlf-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Hlf-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Hlf-P23
#| fig-cap: "Correlation between Cck and Hlf expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Cck", "Hlf"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Hlf expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Hlf-P23-1.png){#fig-romanov2020-correlation-Cck-Hlf-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Tshr-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Tshr-P2
#| fig-cap: "Correlation between Cck and Tshr expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Cck", "Tshr"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Tshr expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Tshr-P2-1.png){#fig-romanov2020-correlation-Cck-Tshr-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Tshr-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Tshr-P10
#| fig-cap: "Correlation between Cck and Tshr expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Cck", "Tshr"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Tshr expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Tshr-P10-1.png){#fig-romanov2020-correlation-Cck-Tshr-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Tshr-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Tshr-P23
#| fig-cap: "Correlation between Cck and Tshr expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Cck", "Tshr"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Tshr expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Tshr-P23-1.png){#fig-romanov2020-correlation-Cck-Tshr-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Gpr173-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Gpr173-P2
#| fig-cap: "Correlation between Cck and Gpr173 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Cck", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Gpr173 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Gpr173-P2-1.png){#fig-romanov2020-correlation-Cck-Gpr173-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Gpr173-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Gpr173-P10
#| fig-cap: "Correlation between Cck and Gpr173 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Cck", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Gpr173 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Gpr173-P10-1.png){#fig-romanov2020-correlation-Cck-Gpr173-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Gpr173-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Gpr173-P23
#| fig-cap: "Correlation between Cck and Gpr173 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Cck", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Gpr173 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Gpr173-P23-1.png){#fig-romanov2020-correlation-Cck-Gpr173-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Lhx3-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Lhx3-P2
#| fig-cap: "Correlation between Cck and Lhx3 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Cck", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Lhx3 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Lhx3-P2-1.png){#fig-romanov2020-correlation-Cck-Lhx3-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Lhx3-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Lhx3-P10
#| fig-cap: "Correlation between Cck and Lhx3 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Cck", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Lhx3 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Lhx3-P10-1.png){#fig-romanov2020-correlation-Cck-Lhx3-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Lhx3-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Lhx3-P23
#| fig-cap: "Correlation between Cck and Lhx3 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Cck", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Lhx3 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Lhx3-P23-1.png){#fig-romanov2020-correlation-Cck-Lhx3-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Gata2-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Gata2-P2
#| fig-cap: "Correlation between Cck and Gata2 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Cck", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Gata2 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Gata2-P2-1.png){#fig-romanov2020-correlation-Cck-Gata2-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Gata2-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Gata2-P10
#| fig-cap: "Correlation between Cck and Gata2 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Cck", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Gata2 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Gata2-P10-1.png){#fig-romanov2020-correlation-Cck-Gata2-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Cck-Gata2-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Cck-Gata2-P23
#| fig-cap: "Correlation between Cck and Gata2 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Cck", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Cck and Gata2 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Cck-Gata2-P23-1.png){#fig-romanov2020-correlation-Cck-Gata2-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Eya1-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Eya1-P2
#| fig-cap: "Correlation between Pitx1 and Eya1 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Pitx1", "Eya1"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Eya1 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Eya1-P2-1.png){#fig-romanov2020-correlation-Pitx1-Eya1-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Eya1-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Eya1-P10
#| fig-cap: "Correlation between Pitx1 and Eya1 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Pitx1", "Eya1"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Eya1 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Eya1-P10-1.png){#fig-romanov2020-correlation-Pitx1-Eya1-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Eya1-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Eya1-P23
#| fig-cap: "Correlation between Pitx1 and Eya1 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Pitx1", "Eya1"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Eya1 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Eya1-P23-1.png){#fig-romanov2020-correlation-Pitx1-Eya1-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Eya2-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Eya2-P2
#| fig-cap: "Correlation between Pitx1 and Eya2 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Pitx1", "Eya2"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Eya2 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Eya2-P2-1.png){#fig-romanov2020-correlation-Pitx1-Eya2-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Eya2-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Eya2-P23
#| fig-cap: "Correlation between Pitx1 and Eya2 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Pitx1", "Eya2"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Eya2 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Eya2-P23-1.png){#fig-romanov2020-correlation-Pitx1-Eya2-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Eya3-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Eya3-P2
#| fig-cap: "Correlation between Pitx1 and Eya3 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Pitx1", "Eya3"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Eya3 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Eya3-P2-1.png){#fig-romanov2020-correlation-Pitx1-Eya3-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Eya3-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Eya3-P23
#| fig-cap: "Correlation between Pitx1 and Eya3 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Pitx1", "Eya3"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Eya3 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Eya3-P23-1.png){#fig-romanov2020-correlation-Pitx1-Eya3-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Eya4-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Eya4-P2
#| fig-cap: "Correlation between Pitx1 and Eya4 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Pitx1", "Eya4"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Eya4 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Eya4-P2-1.png){#fig-romanov2020-correlation-Pitx1-Eya4-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Eya4-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Eya4-P10
#| fig-cap: "Correlation between Pitx1 and Eya4 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Pitx1", "Eya4"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Eya4 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Eya4-P10-1.png){#fig-romanov2020-correlation-Pitx1-Eya4-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Eya4-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Eya4-P23
#| fig-cap: "Correlation between Pitx1 and Eya4 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Pitx1", "Eya4"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Eya4 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Eya4-P23-1.png){#fig-romanov2020-correlation-Pitx1-Eya4-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Sox2-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Sox2-P2
#| fig-cap: "Correlation between Pitx1 and Sox2 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Pitx1", "Sox2"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Sox2 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Sox2-P2-1.png){#fig-romanov2020-correlation-Pitx1-Sox2-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Sox2-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Sox2-P10
#| fig-cap: "Correlation between Pitx1 and Sox2 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Pitx1", "Sox2"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Sox2 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Sox2-P10-1.png){#fig-romanov2020-correlation-Pitx1-Sox2-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Sox2-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Sox2-P23
#| fig-cap: "Correlation between Pitx1 and Sox2 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Pitx1", "Sox2"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Sox2 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Sox2-P23-1.png){#fig-romanov2020-correlation-Pitx1-Sox2-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Hlf-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Hlf-P2
#| fig-cap: "Correlation between Pitx1 and Hlf expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Pitx1", "Hlf"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Hlf expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Hlf-P2-1.png){#fig-romanov2020-correlation-Pitx1-Hlf-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Hlf-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Hlf-P10
#| fig-cap: "Correlation between Pitx1 and Hlf expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Pitx1", "Hlf"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Hlf expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Hlf-P10-1.png){#fig-romanov2020-correlation-Pitx1-Hlf-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Hlf-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Hlf-P23
#| fig-cap: "Correlation between Pitx1 and Hlf expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Pitx1", "Hlf"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Hlf expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Hlf-P23-1.png){#fig-romanov2020-correlation-Pitx1-Hlf-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Tshr-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Tshr-P2
#| fig-cap: "Correlation between Pitx1 and Tshr expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Pitx1", "Tshr"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Tshr expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Tshr-P2-1.png){#fig-romanov2020-correlation-Pitx1-Tshr-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Tshr-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Tshr-P10
#| fig-cap: "Correlation between Pitx1 and Tshr expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Pitx1", "Tshr"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Tshr expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Tshr-P10-1.png){#fig-romanov2020-correlation-Pitx1-Tshr-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Tshr-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Tshr-P23
#| fig-cap: "Correlation between Pitx1 and Tshr expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Pitx1", "Tshr"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Tshr expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Tshr-P23-1.png){#fig-romanov2020-correlation-Pitx1-Tshr-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Gpr173-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Gpr173-P2
#| fig-cap: "Correlation between Pitx1 and Gpr173 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Pitx1", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Gpr173 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Gpr173-P2-1.png){#fig-romanov2020-correlation-Pitx1-Gpr173-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Gpr173-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Gpr173-P10
#| fig-cap: "Correlation between Pitx1 and Gpr173 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Pitx1", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Gpr173 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Gpr173-P10-1.png){#fig-romanov2020-correlation-Pitx1-Gpr173-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Gpr173-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Gpr173-P23
#| fig-cap: "Correlation between Pitx1 and Gpr173 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Pitx1", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Gpr173 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Gpr173-P23-1.png){#fig-romanov2020-correlation-Pitx1-Gpr173-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Lhx3-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Lhx3-P2
#| fig-cap: "Correlation between Pitx1 and Lhx3 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Pitx1", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Lhx3 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Lhx3-P2-1.png){#fig-romanov2020-correlation-Pitx1-Lhx3-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Lhx3-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Lhx3-P10
#| fig-cap: "Correlation between Pitx1 and Lhx3 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Pitx1", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Lhx3 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Lhx3-P10-1.png){#fig-romanov2020-correlation-Pitx1-Lhx3-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Lhx3-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Lhx3-P23
#| fig-cap: "Correlation between Pitx1 and Lhx3 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Pitx1", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Lhx3 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Lhx3-P23-1.png){#fig-romanov2020-correlation-Pitx1-Lhx3-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Gata2-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Gata2-P2
#| fig-cap: "Correlation between Pitx1 and Gata2 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Pitx1", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Gata2 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Gata2-P2-1.png){#fig-romanov2020-correlation-Pitx1-Gata2-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Gata2-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Gata2-P10
#| fig-cap: "Correlation between Pitx1 and Gata2 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Pitx1", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Gata2 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Gata2-P10-1.png){#fig-romanov2020-correlation-Pitx1-Gata2-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Pitx1-Gata2-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Pitx1-Gata2-P23
#| fig-cap: "Correlation between Pitx1 and Gata2 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Pitx1", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Pitx1 and Gata2 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Pitx1-Gata2-P23-1.png){#fig-romanov2020-correlation-Pitx1-Gata2-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya1-Eya2-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya1-Eya2-P2
#| fig-cap: "Correlation between Eya1 and Eya2 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya1", "Eya2"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya2 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya1-Eya2-P2-1.png){#fig-romanov2020-correlation-Eya1-Eya2-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya1-Eya2-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya1-Eya2-P23
#| fig-cap: "Correlation between Eya1 and Eya2 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya1", "Eya2"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya2 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya1-Eya2-P23-1.png){#fig-romanov2020-correlation-Eya1-Eya2-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya1-Eya3-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya1-Eya3-P2
#| fig-cap: "Correlation between Eya1 and Eya3 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya1", "Eya3"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya3 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya1-Eya3-P2-1.png){#fig-romanov2020-correlation-Eya1-Eya3-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya1-Eya3-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya1-Eya3-P23
#| fig-cap: "Correlation between Eya1 and Eya3 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya1", "Eya3"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya3 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya1-Eya3-P23-1.png){#fig-romanov2020-correlation-Eya1-Eya3-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya1-Eya4-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya1-Eya4-P2
#| fig-cap: "Correlation between Eya1 and Eya4 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya1", "Eya4"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya4 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya1-Eya4-P2-1.png){#fig-romanov2020-correlation-Eya1-Eya4-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya1-Eya4-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya1-Eya4-P10
#| fig-cap: "Correlation between Eya1 and Eya4 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Eya1", "Eya4"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya4 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya1-Eya4-P10-1.png){#fig-romanov2020-correlation-Eya1-Eya4-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya1-Eya4-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya1-Eya4-P23
#| fig-cap: "Correlation between Eya1 and Eya4 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya1", "Eya4"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Eya4 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya1-Eya4-P23-1.png){#fig-romanov2020-correlation-Eya1-Eya4-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya1-Sox2-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya1-Sox2-P2
#| fig-cap: "Correlation between Eya1 and Sox2 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya1", "Sox2"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Sox2 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya1-Sox2-P2-1.png){#fig-romanov2020-correlation-Eya1-Sox2-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya1-Sox2-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya1-Sox2-P10
#| fig-cap: "Correlation between Eya1 and Sox2 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Eya1", "Sox2"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Sox2 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya1-Sox2-P10-1.png){#fig-romanov2020-correlation-Eya1-Sox2-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya1-Sox2-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya1-Sox2-P23
#| fig-cap: "Correlation between Eya1 and Sox2 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya1", "Sox2"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Sox2 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya1-Sox2-P23-1.png){#fig-romanov2020-correlation-Eya1-Sox2-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya1-Hlf-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya1-Hlf-P2
#| fig-cap: "Correlation between Eya1 and Hlf expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya1", "Hlf"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Hlf expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya1-Hlf-P2-1.png){#fig-romanov2020-correlation-Eya1-Hlf-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya1-Hlf-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya1-Hlf-P10
#| fig-cap: "Correlation between Eya1 and Hlf expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Eya1", "Hlf"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Hlf expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya1-Hlf-P10-1.png){#fig-romanov2020-correlation-Eya1-Hlf-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya1-Hlf-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya1-Hlf-P23
#| fig-cap: "Correlation between Eya1 and Hlf expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya1", "Hlf"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Hlf expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya1-Hlf-P23-1.png){#fig-romanov2020-correlation-Eya1-Hlf-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya1-Tshr-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya1-Tshr-P2
#| fig-cap: "Correlation between Eya1 and Tshr expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya1", "Tshr"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Tshr expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya1-Tshr-P2-1.png){#fig-romanov2020-correlation-Eya1-Tshr-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya1-Tshr-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya1-Tshr-P10
#| fig-cap: "Correlation between Eya1 and Tshr expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Eya1", "Tshr"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Tshr expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya1-Tshr-P10-1.png){#fig-romanov2020-correlation-Eya1-Tshr-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya1-Tshr-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya1-Tshr-P23
#| fig-cap: "Correlation between Eya1 and Tshr expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya1", "Tshr"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Tshr expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya1-Tshr-P23-1.png){#fig-romanov2020-correlation-Eya1-Tshr-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya1-Gpr173-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya1-Gpr173-P2
#| fig-cap: "Correlation between Eya1 and Gpr173 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya1", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Gpr173 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya1-Gpr173-P2-1.png){#fig-romanov2020-correlation-Eya1-Gpr173-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya1-Gpr173-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya1-Gpr173-P10
#| fig-cap: "Correlation between Eya1 and Gpr173 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Eya1", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Gpr173 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya1-Gpr173-P10-1.png){#fig-romanov2020-correlation-Eya1-Gpr173-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya1-Gpr173-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya1-Gpr173-P23
#| fig-cap: "Correlation between Eya1 and Gpr173 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya1", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Gpr173 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya1-Gpr173-P23-1.png){#fig-romanov2020-correlation-Eya1-Gpr173-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya1-Lhx3-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya1-Lhx3-P2
#| fig-cap: "Correlation between Eya1 and Lhx3 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya1", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Lhx3 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya1-Lhx3-P2-1.png){#fig-romanov2020-correlation-Eya1-Lhx3-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya1-Lhx3-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya1-Lhx3-P10
#| fig-cap: "Correlation between Eya1 and Lhx3 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Eya1", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Lhx3 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya1-Lhx3-P10-1.png){#fig-romanov2020-correlation-Eya1-Lhx3-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya1-Lhx3-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya1-Lhx3-P23
#| fig-cap: "Correlation between Eya1 and Lhx3 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya1", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Lhx3 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya1-Lhx3-P23-1.png){#fig-romanov2020-correlation-Eya1-Lhx3-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya1-Gata2-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya1-Gata2-P2
#| fig-cap: "Correlation between Eya1 and Gata2 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya1", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Gata2 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya1-Gata2-P2-1.png){#fig-romanov2020-correlation-Eya1-Gata2-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya1-Gata2-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya1-Gata2-P10
#| fig-cap: "Correlation between Eya1 and Gata2 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Eya1", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Gata2 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya1-Gata2-P10-1.png){#fig-romanov2020-correlation-Eya1-Gata2-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya1-Gata2-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya1-Gata2-P23
#| fig-cap: "Correlation between Eya1 and Gata2 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya1", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya1 and Gata2 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya1-Gata2-P23-1.png){#fig-romanov2020-correlation-Eya1-Gata2-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya2-Eya3-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya2-Eya3-P2
#| fig-cap: "Correlation between Eya2 and Eya3 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya2", "Eya3"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Eya3 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya2-Eya3-P2-1.png){#fig-romanov2020-correlation-Eya2-Eya3-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya2-Eya3-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya2-Eya3-P23
#| fig-cap: "Correlation between Eya2 and Eya3 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya2", "Eya3"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Eya3 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya2-Eya3-P23-1.png){#fig-romanov2020-correlation-Eya2-Eya3-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya2-Eya4-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya2-Eya4-E15
#| fig-cap: "Correlation between Eya2 and Eya4 expression at E15 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "E15", cells = cells_to_check),
    c("Eya2", "Eya4"),
    "wtree",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Eya4 expression at E15 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya2-Eya4-E15-1.png){#fig-romanov2020-correlation-Eya2-Eya4-E15 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya2-Eya4-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya2-Eya4-P2
#| fig-cap: "Correlation between Eya2 and Eya4 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya2", "Eya4"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Eya4 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya2-Eya4-P2-1.png){#fig-romanov2020-correlation-Eya2-Eya4-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya2-Eya4-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya2-Eya4-P23
#| fig-cap: "Correlation between Eya2 and Eya4 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya2", "Eya4"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Eya4 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya2-Eya4-P23-1.png){#fig-romanov2020-correlation-Eya2-Eya4-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya2-Sox2-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya2-Sox2-E15
#| fig-cap: "Correlation between Eya2 and Sox2 expression at E15 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "E15", cells = cells_to_check),
    c("Eya2", "Sox2"),
    "wtree",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Sox2 expression at E15 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya2-Sox2-E15-1.png){#fig-romanov2020-correlation-Eya2-Sox2-E15 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya2-Sox2-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya2-Sox2-P2
#| fig-cap: "Correlation between Eya2 and Sox2 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya2", "Sox2"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Sox2 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya2-Sox2-P2-1.png){#fig-romanov2020-correlation-Eya2-Sox2-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya2-Sox2-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya2-Sox2-P23
#| fig-cap: "Correlation between Eya2 and Sox2 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya2", "Sox2"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Sox2 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya2-Sox2-P23-1.png){#fig-romanov2020-correlation-Eya2-Sox2-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya2-Hlf-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya2-Hlf-P2
#| fig-cap: "Correlation between Eya2 and Hlf expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya2", "Hlf"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Hlf expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya2-Hlf-P2-1.png){#fig-romanov2020-correlation-Eya2-Hlf-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya2-Hlf-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya2-Hlf-P23
#| fig-cap: "Correlation between Eya2 and Hlf expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya2", "Hlf"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Hlf expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya2-Hlf-P23-1.png){#fig-romanov2020-correlation-Eya2-Hlf-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya2-Tshr-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya2-Tshr-P2
#| fig-cap: "Correlation between Eya2 and Tshr expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya2", "Tshr"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Tshr expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya2-Tshr-P2-1.png){#fig-romanov2020-correlation-Eya2-Tshr-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya2-Tshr-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya2-Tshr-P23
#| fig-cap: "Correlation between Eya2 and Tshr expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya2", "Tshr"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Tshr expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya2-Tshr-P23-1.png){#fig-romanov2020-correlation-Eya2-Tshr-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya2-Gpr173-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya2-Gpr173-P2
#| fig-cap: "Correlation between Eya2 and Gpr173 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya2", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Gpr173 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya2-Gpr173-P2-1.png){#fig-romanov2020-correlation-Eya2-Gpr173-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya2-Gpr173-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya2-Gpr173-P23
#| fig-cap: "Correlation between Eya2 and Gpr173 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya2", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Gpr173 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya2-Gpr173-P23-1.png){#fig-romanov2020-correlation-Eya2-Gpr173-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya2-Lhx3-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya2-Lhx3-P2
#| fig-cap: "Correlation between Eya2 and Lhx3 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya2", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Lhx3 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya2-Lhx3-P2-1.png){#fig-romanov2020-correlation-Eya2-Lhx3-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya2-Lhx3-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya2-Lhx3-P23
#| fig-cap: "Correlation between Eya2 and Lhx3 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya2", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Lhx3 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya2-Lhx3-P23-1.png){#fig-romanov2020-correlation-Eya2-Lhx3-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya2-Gata2-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya2-Gata2-P2
#| fig-cap: "Correlation between Eya2 and Gata2 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya2", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Gata2 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya2-Gata2-P2-1.png){#fig-romanov2020-correlation-Eya2-Gata2-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya2-Gata2-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya2-Gata2-P23
#| fig-cap: "Correlation between Eya2 and Gata2 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya2", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya2 and Gata2 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya2-Gata2-P23-1.png){#fig-romanov2020-correlation-Eya2-Gata2-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya3-Eya4-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya3-Eya4-P2
#| fig-cap: "Correlation between Eya3 and Eya4 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya3", "Eya4"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Eya4 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya3-Eya4-P2-1.png){#fig-romanov2020-correlation-Eya3-Eya4-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya3-Eya4-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya3-Eya4-P23
#| fig-cap: "Correlation between Eya3 and Eya4 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya3", "Eya4"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Eya4 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya3-Eya4-P23-1.png){#fig-romanov2020-correlation-Eya3-Eya4-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya3-Sox2-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya3-Sox2-P2
#| fig-cap: "Correlation between Eya3 and Sox2 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya3", "Sox2"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Sox2 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya3-Sox2-P2-1.png){#fig-romanov2020-correlation-Eya3-Sox2-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya3-Sox2-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya3-Sox2-P23
#| fig-cap: "Correlation between Eya3 and Sox2 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya3", "Sox2"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Sox2 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya3-Sox2-P23-1.png){#fig-romanov2020-correlation-Eya3-Sox2-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya3-Hlf-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya3-Hlf-P2
#| fig-cap: "Correlation between Eya3 and Hlf expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya3", "Hlf"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Hlf expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya3-Hlf-P2-1.png){#fig-romanov2020-correlation-Eya3-Hlf-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya3-Hlf-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya3-Hlf-P23
#| fig-cap: "Correlation between Eya3 and Hlf expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya3", "Hlf"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Hlf expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya3-Hlf-P23-1.png){#fig-romanov2020-correlation-Eya3-Hlf-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya3-Tshr-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya3-Tshr-P2
#| fig-cap: "Correlation between Eya3 and Tshr expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya3", "Tshr"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Tshr expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya3-Tshr-P2-1.png){#fig-romanov2020-correlation-Eya3-Tshr-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya3-Tshr-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya3-Tshr-P23
#| fig-cap: "Correlation between Eya3 and Tshr expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya3", "Tshr"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Tshr expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya3-Tshr-P23-1.png){#fig-romanov2020-correlation-Eya3-Tshr-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya3-Gpr173-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya3-Gpr173-P2
#| fig-cap: "Correlation between Eya3 and Gpr173 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya3", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Gpr173 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya3-Gpr173-P2-1.png){#fig-romanov2020-correlation-Eya3-Gpr173-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya3-Gpr173-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya3-Gpr173-P23
#| fig-cap: "Correlation between Eya3 and Gpr173 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya3", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Gpr173 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya3-Gpr173-P23-1.png){#fig-romanov2020-correlation-Eya3-Gpr173-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya3-Lhx3-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya3-Lhx3-P2
#| fig-cap: "Correlation between Eya3 and Lhx3 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya3", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Lhx3 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya3-Lhx3-P2-1.png){#fig-romanov2020-correlation-Eya3-Lhx3-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya3-Lhx3-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya3-Lhx3-P23
#| fig-cap: "Correlation between Eya3 and Lhx3 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya3", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Lhx3 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya3-Lhx3-P23-1.png){#fig-romanov2020-correlation-Eya3-Lhx3-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya3-Gata2-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya3-Gata2-P2
#| fig-cap: "Correlation between Eya3 and Gata2 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya3", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Gata2 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya3-Gata2-P2-1.png){#fig-romanov2020-correlation-Eya3-Gata2-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya3-Gata2-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya3-Gata2-P23
#| fig-cap: "Correlation between Eya3 and Gata2 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya3", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya3 and Gata2 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya3-Gata2-P23-1.png){#fig-romanov2020-correlation-Eya3-Gata2-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya4-Sox2-E15 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya4-Sox2-E15
#| fig-cap: "Correlation between Eya4 and Sox2 expression at E15 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "E15", cells = cells_to_check),
    c("Eya4", "Sox2"),
    "wtree",
    cells_to_check,
    age = "E15"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Sox2 expression at E15 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya4-Sox2-E15-1.png){#fig-romanov2020-correlation-Eya4-Sox2-E15 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya4-Sox2-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya4-Sox2-P2
#| fig-cap: "Correlation between Eya4 and Sox2 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya4", "Sox2"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Sox2 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya4-Sox2-P2-1.png){#fig-romanov2020-correlation-Eya4-Sox2-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya4-Sox2-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya4-Sox2-P10
#| fig-cap: "Correlation between Eya4 and Sox2 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Eya4", "Sox2"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Sox2 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya4-Sox2-P10-1.png){#fig-romanov2020-correlation-Eya4-Sox2-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya4-Sox2-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya4-Sox2-P23
#| fig-cap: "Correlation between Eya4 and Sox2 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya4", "Sox2"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Sox2 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya4-Sox2-P23-1.png){#fig-romanov2020-correlation-Eya4-Sox2-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya4-Hlf-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya4-Hlf-P2
#| fig-cap: "Correlation between Eya4 and Hlf expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya4", "Hlf"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Hlf expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya4-Hlf-P2-1.png){#fig-romanov2020-correlation-Eya4-Hlf-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya4-Hlf-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya4-Hlf-P10
#| fig-cap: "Correlation between Eya4 and Hlf expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Eya4", "Hlf"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Hlf expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya4-Hlf-P10-1.png){#fig-romanov2020-correlation-Eya4-Hlf-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya4-Hlf-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya4-Hlf-P23
#| fig-cap: "Correlation between Eya4 and Hlf expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya4", "Hlf"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Hlf expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya4-Hlf-P23-1.png){#fig-romanov2020-correlation-Eya4-Hlf-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya4-Tshr-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya4-Tshr-P2
#| fig-cap: "Correlation between Eya4 and Tshr expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya4", "Tshr"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Tshr expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya4-Tshr-P2-1.png){#fig-romanov2020-correlation-Eya4-Tshr-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya4-Tshr-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya4-Tshr-P10
#| fig-cap: "Correlation between Eya4 and Tshr expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Eya4", "Tshr"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Tshr expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya4-Tshr-P10-1.png){#fig-romanov2020-correlation-Eya4-Tshr-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya4-Tshr-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya4-Tshr-P23
#| fig-cap: "Correlation between Eya4 and Tshr expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya4", "Tshr"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Tshr expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya4-Tshr-P23-1.png){#fig-romanov2020-correlation-Eya4-Tshr-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya4-Gpr173-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya4-Gpr173-P2
#| fig-cap: "Correlation between Eya4 and Gpr173 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya4", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Gpr173 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya4-Gpr173-P2-1.png){#fig-romanov2020-correlation-Eya4-Gpr173-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya4-Gpr173-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya4-Gpr173-P10
#| fig-cap: "Correlation between Eya4 and Gpr173 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Eya4", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Gpr173 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya4-Gpr173-P10-1.png){#fig-romanov2020-correlation-Eya4-Gpr173-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya4-Gpr173-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya4-Gpr173-P23
#| fig-cap: "Correlation between Eya4 and Gpr173 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya4", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Gpr173 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya4-Gpr173-P23-1.png){#fig-romanov2020-correlation-Eya4-Gpr173-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya4-Lhx3-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya4-Lhx3-P2
#| fig-cap: "Correlation between Eya4 and Lhx3 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya4", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Lhx3 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya4-Lhx3-P2-1.png){#fig-romanov2020-correlation-Eya4-Lhx3-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya4-Lhx3-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya4-Lhx3-P10
#| fig-cap: "Correlation between Eya4 and Lhx3 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Eya4", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Lhx3 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya4-Lhx3-P10-1.png){#fig-romanov2020-correlation-Eya4-Lhx3-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya4-Lhx3-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya4-Lhx3-P23
#| fig-cap: "Correlation between Eya4 and Lhx3 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya4", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Lhx3 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya4-Lhx3-P23-1.png){#fig-romanov2020-correlation-Eya4-Lhx3-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya4-Gata2-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya4-Gata2-P2
#| fig-cap: "Correlation between Eya4 and Gata2 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Eya4", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Gata2 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya4-Gata2-P2-1.png){#fig-romanov2020-correlation-Eya4-Gata2-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya4-Gata2-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya4-Gata2-P10
#| fig-cap: "Correlation between Eya4 and Gata2 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Eya4", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Gata2 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya4-Gata2-P10-1.png){#fig-romanov2020-correlation-Eya4-Gata2-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Eya4-Gata2-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Eya4-Gata2-P23
#| fig-cap: "Correlation between Eya4 and Gata2 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Eya4", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Eya4 and Gata2 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Eya4-Gata2-P23-1.png){#fig-romanov2020-correlation-Eya4-Gata2-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Sox2-Hlf-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Sox2-Hlf-P2
#| fig-cap: "Correlation between Sox2 and Hlf expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Sox2", "Hlf"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Hlf expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Sox2-Hlf-P2-1.png){#fig-romanov2020-correlation-Sox2-Hlf-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Sox2-Hlf-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Sox2-Hlf-P10
#| fig-cap: "Correlation between Sox2 and Hlf expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Sox2", "Hlf"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Hlf expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Sox2-Hlf-P10-1.png){#fig-romanov2020-correlation-Sox2-Hlf-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Sox2-Hlf-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Sox2-Hlf-P23
#| fig-cap: "Correlation between Sox2 and Hlf expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Sox2", "Hlf"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Hlf expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Sox2-Hlf-P23-1.png){#fig-romanov2020-correlation-Sox2-Hlf-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Sox2-Tshr-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Sox2-Tshr-P2
#| fig-cap: "Correlation between Sox2 and Tshr expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Sox2", "Tshr"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Tshr expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Sox2-Tshr-P2-1.png){#fig-romanov2020-correlation-Sox2-Tshr-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Sox2-Tshr-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Sox2-Tshr-P10
#| fig-cap: "Correlation between Sox2 and Tshr expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Sox2", "Tshr"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Tshr expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Sox2-Tshr-P10-1.png){#fig-romanov2020-correlation-Sox2-Tshr-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Sox2-Tshr-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Sox2-Tshr-P23
#| fig-cap: "Correlation between Sox2 and Tshr expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Sox2", "Tshr"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Tshr expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Sox2-Tshr-P23-1.png){#fig-romanov2020-correlation-Sox2-Tshr-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Sox2-Gpr173-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Sox2-Gpr173-P2
#| fig-cap: "Correlation between Sox2 and Gpr173 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Sox2", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Gpr173 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Sox2-Gpr173-P2-1.png){#fig-romanov2020-correlation-Sox2-Gpr173-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Sox2-Gpr173-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Sox2-Gpr173-P10
#| fig-cap: "Correlation between Sox2 and Gpr173 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Sox2", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Gpr173 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Sox2-Gpr173-P10-1.png){#fig-romanov2020-correlation-Sox2-Gpr173-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Sox2-Gpr173-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Sox2-Gpr173-P23
#| fig-cap: "Correlation between Sox2 and Gpr173 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Sox2", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Gpr173 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Sox2-Gpr173-P23-1.png){#fig-romanov2020-correlation-Sox2-Gpr173-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Sox2-Lhx3-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Sox2-Lhx3-P2
#| fig-cap: "Correlation between Sox2 and Lhx3 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Sox2", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Lhx3 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Sox2-Lhx3-P2-1.png){#fig-romanov2020-correlation-Sox2-Lhx3-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Sox2-Lhx3-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Sox2-Lhx3-P10
#| fig-cap: "Correlation between Sox2 and Lhx3 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Sox2", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Lhx3 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Sox2-Lhx3-P10-1.png){#fig-romanov2020-correlation-Sox2-Lhx3-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Sox2-Lhx3-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Sox2-Lhx3-P23
#| fig-cap: "Correlation between Sox2 and Lhx3 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Sox2", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Lhx3 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Sox2-Lhx3-P23-1.png){#fig-romanov2020-correlation-Sox2-Lhx3-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Sox2-Gata2-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Sox2-Gata2-P2
#| fig-cap: "Correlation between Sox2 and Gata2 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Sox2", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Gata2 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Sox2-Gata2-P2-1.png){#fig-romanov2020-correlation-Sox2-Gata2-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Sox2-Gata2-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Sox2-Gata2-P10
#| fig-cap: "Correlation between Sox2 and Gata2 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Sox2", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Gata2 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Sox2-Gata2-P10-1.png){#fig-romanov2020-correlation-Sox2-Gata2-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Sox2-Gata2-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Sox2-Gata2-P23
#| fig-cap: "Correlation between Sox2 and Gata2 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Sox2", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Sox2 and Gata2 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Sox2-Gata2-P23-1.png){#fig-romanov2020-correlation-Sox2-Gata2-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Hlf-Tshr-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Hlf-Tshr-P2
#| fig-cap: "Correlation between Hlf and Tshr expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Hlf", "Tshr"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Tshr expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Hlf-Tshr-P2-1.png){#fig-romanov2020-correlation-Hlf-Tshr-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Hlf-Tshr-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Hlf-Tshr-P10
#| fig-cap: "Correlation between Hlf and Tshr expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Hlf", "Tshr"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Tshr expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Hlf-Tshr-P10-1.png){#fig-romanov2020-correlation-Hlf-Tshr-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Hlf-Tshr-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Hlf-Tshr-P23
#| fig-cap: "Correlation between Hlf and Tshr expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Hlf", "Tshr"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Tshr expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Hlf-Tshr-P23-1.png){#fig-romanov2020-correlation-Hlf-Tshr-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Hlf-Gpr173-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Hlf-Gpr173-P2
#| fig-cap: "Correlation between Hlf and Gpr173 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Hlf", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Gpr173 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Hlf-Gpr173-P2-1.png){#fig-romanov2020-correlation-Hlf-Gpr173-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Hlf-Gpr173-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Hlf-Gpr173-P10
#| fig-cap: "Correlation between Hlf and Gpr173 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Hlf", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Gpr173 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Hlf-Gpr173-P10-1.png){#fig-romanov2020-correlation-Hlf-Gpr173-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Hlf-Gpr173-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Hlf-Gpr173-P23
#| fig-cap: "Correlation between Hlf and Gpr173 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Hlf", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Gpr173 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Hlf-Gpr173-P23-1.png){#fig-romanov2020-correlation-Hlf-Gpr173-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Hlf-Lhx3-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Hlf-Lhx3-P2
#| fig-cap: "Correlation between Hlf and Lhx3 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Hlf", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Lhx3 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Hlf-Lhx3-P2-1.png){#fig-romanov2020-correlation-Hlf-Lhx3-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Hlf-Lhx3-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Hlf-Lhx3-P10
#| fig-cap: "Correlation between Hlf and Lhx3 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Hlf", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Lhx3 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Hlf-Lhx3-P10-1.png){#fig-romanov2020-correlation-Hlf-Lhx3-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Hlf-Lhx3-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Hlf-Lhx3-P23
#| fig-cap: "Correlation between Hlf and Lhx3 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Hlf", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Lhx3 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Hlf-Lhx3-P23-1.png){#fig-romanov2020-correlation-Hlf-Lhx3-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Hlf-Gata2-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Hlf-Gata2-P2
#| fig-cap: "Correlation between Hlf and Gata2 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Hlf", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Gata2 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Hlf-Gata2-P2-1.png){#fig-romanov2020-correlation-Hlf-Gata2-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Hlf-Gata2-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Hlf-Gata2-P10
#| fig-cap: "Correlation between Hlf and Gata2 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Hlf", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Gata2 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Hlf-Gata2-P10-1.png){#fig-romanov2020-correlation-Hlf-Gata2-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Hlf-Gata2-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Hlf-Gata2-P23
#| fig-cap: "Correlation between Hlf and Gata2 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Hlf", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Hlf and Gata2 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Hlf-Gata2-P23-1.png){#fig-romanov2020-correlation-Hlf-Gata2-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshr-Gpr173-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshr-Gpr173-P2
#| fig-cap: "Correlation between Tshr and Gpr173 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Tshr", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Gpr173 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshr-Gpr173-P2-1.png){#fig-romanov2020-correlation-Tshr-Gpr173-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshr-Gpr173-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshr-Gpr173-P10
#| fig-cap: "Correlation between Tshr and Gpr173 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Tshr", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Gpr173 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshr-Gpr173-P10-1.png){#fig-romanov2020-correlation-Tshr-Gpr173-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshr-Gpr173-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshr-Gpr173-P23
#| fig-cap: "Correlation between Tshr and Gpr173 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Tshr", "Gpr173"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Gpr173 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshr-Gpr173-P23-1.png){#fig-romanov2020-correlation-Tshr-Gpr173-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshr-Lhx3-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshr-Lhx3-P2
#| fig-cap: "Correlation between Tshr and Lhx3 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Tshr", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Lhx3 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshr-Lhx3-P2-1.png){#fig-romanov2020-correlation-Tshr-Lhx3-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshr-Lhx3-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshr-Lhx3-P10
#| fig-cap: "Correlation between Tshr and Lhx3 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Tshr", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Lhx3 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshr-Lhx3-P10-1.png){#fig-romanov2020-correlation-Tshr-Lhx3-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshr-Lhx3-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshr-Lhx3-P23
#| fig-cap: "Correlation between Tshr and Lhx3 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Tshr", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Lhx3 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshr-Lhx3-P23-1.png){#fig-romanov2020-correlation-Tshr-Lhx3-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshr-Gata2-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshr-Gata2-P2
#| fig-cap: "Correlation between Tshr and Gata2 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Tshr", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Gata2 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshr-Gata2-P2-1.png){#fig-romanov2020-correlation-Tshr-Gata2-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshr-Gata2-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshr-Gata2-P10
#| fig-cap: "Correlation between Tshr and Gata2 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Tshr", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Gata2 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshr-Gata2-P10-1.png){#fig-romanov2020-correlation-Tshr-Gata2-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Tshr-Gata2-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Tshr-Gata2-P23
#| fig-cap: "Correlation between Tshr and Gata2 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Tshr", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Tshr and Gata2 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Tshr-Gata2-P23-1.png){#fig-romanov2020-correlation-Tshr-Gata2-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Gpr173-Lhx3-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Gpr173-Lhx3-P2
#| fig-cap: "Correlation between Gpr173 and Lhx3 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Gpr173", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Gpr173 and Lhx3 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Gpr173-Lhx3-P2-1.png){#fig-romanov2020-correlation-Gpr173-Lhx3-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Gpr173-Lhx3-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Gpr173-Lhx3-P10
#| fig-cap: "Correlation between Gpr173 and Lhx3 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Gpr173", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Gpr173 and Lhx3 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Gpr173-Lhx3-P10-1.png){#fig-romanov2020-correlation-Gpr173-Lhx3-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Gpr173-Lhx3-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Gpr173-Lhx3-P23
#| fig-cap: "Correlation between Gpr173 and Lhx3 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Gpr173", "Lhx3"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Gpr173 and Lhx3 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Gpr173-Lhx3-P23-1.png){#fig-romanov2020-correlation-Gpr173-Lhx3-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Gpr173-Gata2-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Gpr173-Gata2-P2
#| fig-cap: "Correlation between Gpr173 and Gata2 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Gpr173", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Gpr173 and Gata2 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Gpr173-Gata2-P2-1.png){#fig-romanov2020-correlation-Gpr173-Gata2-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Gpr173-Gata2-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Gpr173-Gata2-P10
#| fig-cap: "Correlation between Gpr173 and Gata2 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Gpr173", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Gpr173 and Gata2 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Gpr173-Gata2-P10-1.png){#fig-romanov2020-correlation-Gpr173-Gata2-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Gpr173-Gata2-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Gpr173-Gata2-P23
#| fig-cap: "Correlation between Gpr173 and Gata2 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Gpr173", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Gpr173 and Gata2 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Gpr173-Gata2-P23-1.png){#fig-romanov2020-correlation-Gpr173-Gata2-P23 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Lhx3-Gata2-P2 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Lhx3-Gata2-P2
#| fig-cap: "Correlation between Lhx3 and Gata2 expression at P2 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P2", cells = cells_to_check),
    c("Lhx3", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P2"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Lhx3 and Gata2 expression at P2 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Lhx3-Gata2-P2-1.png){#fig-romanov2020-correlation-Lhx3-Gata2-P2 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Lhx3-Gata2-P10 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Lhx3-Gata2-P10
#| fig-cap: "Correlation between Lhx3 and Gata2 expression at P10 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P10", cells = cells_to_check),
    c("Lhx3", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P10"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Lhx3 and Gata2 expression at P10 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Lhx3-Gata2-P10-1.png){#fig-romanov2020-correlation-Lhx3-Gata2-P10 width=1152}
:::
:::

::: {#cell-fig-romanov2020-correlation-Lhx3-Gata2-P23 .cell}

```{.r .cell-code .hidden}
#| label: fig-romanov2020-correlation-Lhx3-Gata2-P23
#| fig-cap: "Correlation between Lhx3 and Gata2 expression at P23 in Romanov et al. 2020 dataset"
#| fig-width: 12
#| fig-height: 8

generate_correlation_plots(
    srt.romanov.pub |> subset(Age == "P23", cells = cells_to_check),
    c("Lhx3", "Gata2"),
    "wtree",
    cells_to_check,
    age = "P23"
)
```

::: {.cell-output .cell-output-stderr .hidden}

```
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_xsidebin()` using `bins = 30`. Pick better value with `binwidth`.
`stat_ysidebin()` using `bins = 30`. Pick better value with `binwidth`.
```


:::

::: {.cell-output-display}
![Correlation between Lhx3 and Gata2 expression at P23 in Romanov et al. 2020 dataset](01-de_test-focus_pars_tub_files/figure-html/fig-romanov2020-correlation-Lhx3-Gata2-P23-1.png){#fig-romanov2020-correlation-Lhx3-Gata2-P23 width=1152}
:::
:::






#### Dotplots





::: {#cell-fig-dotplot-dendrogram-genes-romanov2020 .cell}

```{.r .cell-code .hidden}
#| label: fig-dotplot-dendrogram-genes-romanov2020
#| fig-cap: "Dotplot of selected genes in hypothalamus across different developmental stages in the Romanov et al. (2020) dataset. Cells are colored by expression level.  Note the distinct localization patterns of each gene."
#| fig-width: 6
#| fig-height: 6
Idents(srt.kim) <- "Age"

DotPlot_scCustom(seurat_object = subset(srt.romanov.pub, cells = cells_to_check), colors_use = viridis(n = 30, alpha = .75, direction = -1, option = "E"), features = genes_to_check, flip_axes = T, x_lab_rotate = TRUE, dot.scale = 15)
```

::: {.cell-output-display}
![Dotplot of selected genes in hypothalamus across different developmental stages in the Romanov et al. (2020) dataset. Cells are colored by expression level.  Note the distinct localization patterns of each gene.](01-de_test-focus_pars_tub_files/figure-html/fig-dotplot-dendrogram-genes-romanov2020-1.png){#fig-dotplot-dendrogram-genes-romanov2020 width=576}
:::
:::





## Barplot of Sox2 and Tshr expression in hypothalamus across different developmental stages





::: {#cell-fig-sox2-tshr-bargraph .cell}

```{.r .cell-code .hidden}
#| label: fig-sox2-tshr-bargraph
#| fig-width: 7
#| #| fig-height: 6
# Create a new variable that combines Sox2_pos and Tshr_pos
df <- sbs_mtx %>%
  mutate(Age = factor(Age, levels = c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45"), ordered = TRUE), VarComb = case_when(
    Sox2_pos & Tshr_pos ~ "Sox2+/Tshr+",
    Sox2_pos & !Tshr_pos ~ "Sox2+/Tshr-",
    !Sox2_pos & Tshr_pos ~ "Sox2-/Tshr+",
    !Sox2_pos & !Tshr_pos ~ "Sox2-/Tshr-"
  ))

# Calculate counts and proportions for each category
df_counts <- df %>%
  group_by(Age, VarComb) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n))
```

::: {.cell-output .cell-output-stderr .hidden}

```
`summarise()` has grouped output by 'Age'. You can override using the `.groups`
argument.
```


:::

```{.r .cell-code .hidden}
#| label: fig-sox2-tshr-bargraph
#| fig-width: 7
#| #| fig-height: 6
# Calculate the total counts for each category
df_total_counts <- df %>%
  group_by(Age) %>%
  summarise(total_n = n())

# Create a vector of Age names ordered by prop for Sox2+/Tshr+ cases for each Age
ordered_Ages <- df_counts %>%
  filter(VarComb == "Sox2+/Tshr-") %>%
  arrange(desc(prop)) %>%
  pull(Age)

# Reorder the factor levels of Age
df_counts$Age <- factor(df_counts$Age, levels = ordered_Ages)
df_total_counts$Age <- factor(df_total_counts$Age, levels = ordered_Ages)

# Reorder the factor levels of combinations
df_counts$VarComb <- factor(df_counts$VarComb, levels = rev(c("Sox2+/Tshr-", "Sox2+/Tshr+", "Sox2-/Tshr+", "Sox2-/Tshr-")))

# Create a stacked bar plot
ggplot(df_counts, aes(x = Age, y = n, fill = VarComb)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("Sox2+/Tshr+" = "purple", "Sox2+/Tshr-" = "red3", "Sox2-/Tshr+" = "royalblue", "Sox2-/Tshr-" = "grey50")) +
  labs(x = "Hypothalamic subAge", y = "Number of cells", fill = "Expression") +
  scale_x_discrete(labels = c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45"))

ggplot(df_counts, aes(x = Age, y = prop, fill = VarComb)) +
  geom_bar(stat = "identity", color = "black", position = "fill") +
  scale_fill_manual(values = c("Sox2+/Tshr-" = "red3", "Sox2+/Tshr+" = "purple", "Sox2-/Tshr+" = "royalblue", "Sox2-/Tshr-" = "grey50")) +
  labs(x = "Hypothalamic subAge", y = "Proportion of cells", fill = "Expression") +
  scale_x_discrete(labels = c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45"))

ggplot(df_counts |> filter(!VarComb %in% c("Sox2-/Tshr-", "Sox2+/Tshr-")), aes(x = Age, y = n, fill = VarComb)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("Sox2+/Tshr+" = "purple", "Sox2+/Tshr-" = "red3", "Sox2-/Tshr+" = "royalblue", "Sox2-/Tshr-" = "grey50")) +
  labs(x = "Hypothalamic subAge", y = "Number of cells", fill = "Expression") +
  scale_x_discrete(labels = c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45"))

ggplot(df_counts |> filter(!VarComb %in% c("Sox2-/Tshr-", "Sox2+/Tshr-")), aes(x = Age, y = prop, fill = VarComb)) +
  geom_bar(stat = "identity", color = "black", position = "fill") +
  scale_fill_manual(values = c("Sox2+/Tshr-" = "red3", "Sox2+/Tshr+" = "purple", "Sox2-/Tshr+" = "royalblue", "Sox2-/Tshr-" = "grey50")) +
  labs(x = "Hypothalamic subAge", y = "Proportion of cells", fill = "Expression") +
  scale_x_discrete(labels = c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45"))
```

::: {.cell-output-display}
![](01-de_test-focus_pars_tub_files/figure-html/fig-sox2-tshr-bargraph-1.png){#fig-sox2-tshr-bargraph-1 width=672}
:::

::: {.cell-output-display}
![](01-de_test-focus_pars_tub_files/figure-html/fig-sox2-tshr-bargraph-2.png){#fig-sox2-tshr-bargraph-2 width=672}
:::

::: {.cell-output-display}
![](01-de_test-focus_pars_tub_files/figure-html/fig-sox2-tshr-bargraph-3.png){#fig-sox2-tshr-bargraph-3 width=672}
:::

::: {.cell-output-display}
![](01-de_test-focus_pars_tub_files/figure-html/fig-sox2-tshr-bargraph-4.png){#fig-sox2-tshr-bargraph-4 width=672}
:::
:::





# Calculate and plot chi2 test of independence between Sox2 and Cckbr expression in hypothalamus across different developmental stages





::: {.cell}

```{.r .cell-code .hidden}
#| label: get-goi-sox2-Cckbr
sbs_mtx <- GetAssayData(object = srt.kim, layer = "counts", assay = "RNA")[c("Sox2", "Cckbr"), ] %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  select(Sox2, Cckbr) %>%
  dplyr::bind_cols(srt.kim@meta.data) %>%
  select(Age, Sox2, Cckbr) %>%
  mutate(
    Sox2_pos = Sox2 > 0,
    Cckbr_pos = Cckbr > 0
  )

sbs_mtx %>% skimr::skim()
```

::: {.cell-output-display}

Table: Data summary

|                         |           |
|:------------------------|:----------|
|Name                     |Piped data |
|Number of rows           |128006     |
|Number of columns        |5          |
|_______________________  |           |
|Column type frequency:   |           |
|factor                   |1          |
|logical                  |2          |
|numeric                  |2          |
|________________________ |           |
|Group variables          |None       |


**Variable type: factor**

|skim_variable | n_missing| complete_rate|ordered | n_unique|top_counts                                    |
|:-------------|---------:|-------------:|:-------|--------:|:---------------------------------------------|
|Age           |         0|             1|FALSE   |       12|P45: 17025, E16: 16781, E14: 16543, P8: 13677 |


**Variable type: logical**

|skim_variable | n_missing| complete_rate| mean|count                   |
|:-------------|---------:|-------------:|----:|:-----------------------|
|Sox2_pos      |         0|             1| 0.14|FAL: 110343, TRU: 17663 |
|Cckbr_pos     |         0|             1| 0.00|FAL: 127679, TRU: 327   |


**Variable type: numeric**

|skim_variable | n_missing| complete_rate| mean|   sd| p0| p25| p50| p75| p100|hist  |
|:-------------|---------:|-------------:|----:|----:|--:|---:|---:|---:|----:|:-----|
|Sox2          |         0|             1| 0.27| 1.21|  0|   0|   0|   0|   82|▇▁▁▁▁ |
|Cckbr         |         0|             1| 0.00| 0.06|  0|   0|   0|   0|    7|▇▁▁▁▁ |


:::
:::

::: {#cell-fig-sox2-Cckbr-stats .cell}

```{.r .cell-code .hidden}
#| label: fig-sox2-Cckbr-stats
#| fig-width: 8
#| fig-height: 24
write_csv(sbs_mtx, here(tables_dir, "Sox2-Cckbr-expression-status-between-Ages-on-evaluation-datasets.csv"))


# plot
grouped_ggpiestats(
  data = sbs_mtx,
  x = Cckbr_pos,
  y = Sox2_pos,
  grouping.var = Age,
  perc.k = 1,
  package = "ggsci",
  palette = "category10_d3",
  title.text = "Sox2 specification of Cckbr-positive hypothalamic development",
  caption.text = "Asterisks denote results from proportion tests; \n***: p < 0.001, ns: non-significant",
  plotgrid.args = list(nrow = 8)
)
```

::: {.cell-output-display}
![](01-de_test-focus_pars_tub_files/figure-html/fig-sox2-Cckbr-stats-1.png){#fig-sox2-Cckbr-stats width=768}
:::
:::





## Barplot of Sox2 and Cckbr expression in hypothalamus across different developmental stages





::: {#cell-fig-sox2-Cckbr-bargraph .cell}

```{.r .cell-code .hidden}
#| label: fig-sox2-Cckbr-bargraph
#| fig-width: 7
#| fig-height: 6
# Create a new variable that combines Sox2_pos and Cckbr_pos
df <- sbs_mtx %>%
  mutate(Age = factor(Age, levels = c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45"), ordered = TRUE), VarComb = case_when(
    Sox2_pos & Cckbr_pos ~ "Sox2+/Cckbr+",
    Sox2_pos & !Cckbr_pos ~ "Sox2+/Cckbr-",
    !Sox2_pos & Cckbr_pos ~ "Sox2-/Cckbr+",
    !Sox2_pos & !Cckbr_pos ~ "Sox2-/Cckbr-"
  ))

# Calculate counts and proportions for each category
df_counts <- df %>%
  group_by(Age, VarComb) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n))
```

::: {.cell-output .cell-output-stderr .hidden}

```
`summarise()` has grouped output by 'Age'. You can override using the `.groups`
argument.
```


:::

```{.r .cell-code .hidden}
#| label: fig-sox2-Cckbr-bargraph
#| fig-width: 7
#| fig-height: 6
# Calculate the total counts for each category
df_total_counts <- df %>%
  group_by(Age) %>%
  summarise(total_n = n())

# Create a vector of Age names ordered by prop for Sox2+/Cckbr+ cases for each Age
ordered_Ages <- df_counts %>%
  filter(VarComb == "Sox2+/Cckbr-") %>%
  arrange(desc(prop)) %>%
  pull(Age)

# Reorder the factor levels of Age
df_counts$Age <- factor(df_counts$Age, levels = ordered_Ages)
df_total_counts$Age <- factor(df_total_counts$Age, levels = ordered_Ages)

# Reorder the factor levels of combinations
df_counts$VarComb <- factor(df_counts$VarComb, levels = rev(c("Sox2+/Cckbr-", "Sox2+/Cckbr+", "Sox2-/Cckbr+", "Sox2-/Cckbr-")))

# Create a stacked bar plot
ggplot(df_counts, aes(x = Age, y = n, fill = VarComb)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("Sox2+/Cckbr+" = "yellow", "Sox2+/Cckbr-" = "red3", "Sox2-/Cckbr+" = "green", "Sox2-/Cckbr-" = "grey50")) +
  labs(x = "Hypothalamic subAge", y = "Number of cells", fill = "Expression") +
  scale_x_discrete(labels = c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45"))

ggplot(df_counts, aes(x = Age, y = prop, fill = VarComb)) +
  geom_bar(stat = "identity", color = "black", position = "fill") +
  scale_fill_manual(values = c("Sox2+/Cckbr-" = "red3", "Sox2+/Cckbr+" = "yellow", "Sox2-/Cckbr+" = "green", "Sox2-/Cckbr-" = "grey50")) +
  labs(x = "Hypothalamic subAge", y = "Proportion of cells", fill = "Expression") +
  scale_x_discrete(labels = c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45"))

ggplot(df_counts |> filter(!VarComb %in% c("Sox2-/Cckbr-", "Sox2+/Cckbr-")), aes(x = Age, y = n, fill = VarComb)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("Sox2+/Cckbr+" = "yellow", "Sox2+/Cckbr-" = "red3", "Sox2-/Cckbr+" = "green", "Sox2-/Cckbr-" = "grey50")) +
  labs(x = "Hypothalamic subAge", y = "Number of cells", fill = "Expression") +
  scale_x_discrete(labels = c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45"))

ggplot(df_counts |> filter(!VarComb %in% c("Sox2-/Cckbr-", "Sox2+/Cckbr-")), aes(x = Age, y = prop, fill = VarComb)) +
  geom_bar(stat = "identity", color = "black", position = "fill") +
  scale_fill_manual(values = c("Sox2+/Cckbr-" = "red3", "Sox2+/Cckbr+" = "yellow", "Sox2-/Cckbr+" = "green", "Sox2-/Cckbr-" = "grey50")) +
  labs(x = "Hypothalamic subAge", y = "Proportion of cells", fill = "Expression") +
  scale_x_discrete(labels = c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45"))
```

::: {.cell-output-display}
![](01-de_test-focus_pars_tub_files/figure-html/fig-sox2-Cckbr-bargraph-1.png){#fig-sox2-Cckbr-bargraph-1 width=672}
:::

::: {.cell-output-display}
![](01-de_test-focus_pars_tub_files/figure-html/fig-sox2-Cckbr-bargraph-2.png){#fig-sox2-Cckbr-bargraph-2 width=672}
:::

::: {.cell-output-display}
![](01-de_test-focus_pars_tub_files/figure-html/fig-sox2-Cckbr-bargraph-3.png){#fig-sox2-Cckbr-bargraph-3 width=672}
:::

::: {.cell-output-display}
![](01-de_test-focus_pars_tub_files/figure-html/fig-sox2-Cckbr-bargraph-4.png){#fig-sox2-Cckbr-bargraph-4 width=672}
:::
:::





# Calculate and plot chi2 test of independence between Cckbr and Tshr expression in hypothalamus across different developmental stages





::: {.cell}

```{.r .cell-code .hidden}
#| label: get-goi-Cckbr-tshr
sbs_mtx <- GetAssayData(object = srt.kim, layer = "counts", assay = "RNA")[c("Cckbr", "Tshr"), ] %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  select(Cckbr, Tshr) %>%
  dplyr::bind_cols(srt.kim@meta.data) %>%
  select(Age, Cckbr, Tshr) %>%
  mutate(
    Cckbr_pos = Cckbr > 0,
    Tshr_pos = Tshr > 0
  )

sbs_mtx %>% skimr::skim()
```

::: {.cell-output-display}

Table: Data summary

|                         |           |
|:------------------------|:----------|
|Name                     |Piped data |
|Number of rows           |128006     |
|Number of columns        |5          |
|_______________________  |           |
|Column type frequency:   |           |
|factor                   |1          |
|logical                  |2          |
|numeric                  |2          |
|________________________ |           |
|Group variables          |None       |


**Variable type: factor**

|skim_variable | n_missing| complete_rate|ordered | n_unique|top_counts                                    |
|:-------------|---------:|-------------:|:-------|--------:|:---------------------------------------------|
|Age           |         0|             1|FALSE   |       12|P45: 17025, E16: 16781, E14: 16543, P8: 13677 |


**Variable type: logical**

|skim_variable | n_missing| complete_rate| mean|count                 |
|:-------------|---------:|-------------:|----:|:---------------------|
|Cckbr_pos     |         0|             1|    0|FAL: 127679, TRU: 327 |
|Tshr_pos      |         0|             1|    0|FAL: 127714, TRU: 292 |


**Variable type: numeric**

|skim_variable | n_missing| complete_rate| mean|   sd| p0| p25| p50| p75| p100|hist  |
|:-------------|---------:|-------------:|----:|----:|--:|---:|---:|---:|----:|:-----|
|Cckbr         |         0|             1|    0| 0.06|  0|   0|   0|   0|    7|▇▁▁▁▁ |
|Tshr          |         0|             1|    0| 0.06|  0|   0|   0|   0|    5|▇▁▁▁▁ |


:::
:::

::: {#cell-fig-Cckbr-tshr-stats .cell}

```{.r .cell-code .hidden}
#| label: fig-Cckbr-tshr-stats
#| fig-width: 8
#| fig-height: 24
write_csv(sbs_mtx, here(tables_dir, "Cckbr-Tshr-expression-status-between-Ages-on-evaluation-datasets.csv"))


# plot
grouped_ggpiestats(
  data = sbs_mtx,
  x = Cckbr_pos,
  y = Tshr_pos,
  grouping.var = Age,
  perc.k = 1,
  package = "ggsci",
  palette = "category10_d3",
  title.text = "Cckbr specification of Tshr-positive hypothalamic development",
  caption.text = "Asterisks denote results from proportion tests; \n***: p < 0.001, ns: non-significant",
  plotgrid.args = list(nrow = 8)
)
```

::: {.cell-output-display}
![](01-de_test-focus_pars_tub_files/figure-html/fig-Cckbr-tshr-stats-1.png){#fig-Cckbr-tshr-stats width=768}
:::
:::





## Barplot of Cckbr and Tshr expression in hypothalamus across different developmental stages





::: {#cell-fig-Cckbr-tshr-bargraph .cell}

```{.r .cell-code .hidden}
#| label: fig-Cckbr-tshr-bargraph
#| fig-width: 7
#| fig-height: 6
# Create a new variable that combines Cckbr_pos and Tshr_pos
df <- sbs_mtx %>%
  mutate(Age = factor(Age, levels = c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45"), ordered = TRUE), VarComb = case_when(
    Cckbr_pos & Tshr_pos ~ "Cckbr+/Tshr+",
    Cckbr_pos & !Tshr_pos ~ "Cckbr+/Tshr-",
    !Cckbr_pos & Tshr_pos ~ "Cckbr-/Tshr+",
    !Cckbr_pos & !Tshr_pos ~ "Cckbr-/Tshr-"
  ))

# Calculate counts and proportions for each category
df_counts <- df %>%
  group_by(Age, VarComb) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n))
```

::: {.cell-output .cell-output-stderr .hidden}

```
`summarise()` has grouped output by 'Age'. You can override using the `.groups`
argument.
```


:::

```{.r .cell-code .hidden}
#| label: fig-Cckbr-tshr-bargraph
#| fig-width: 7
#| fig-height: 6
# Calculate the total counts for each category
df_total_counts <- df %>%
  group_by(Age) %>%
  summarise(total_n = n())

# Create a vector of Age names ordered by prop for Cckbr+/Tshr+ cases for each Age
ordered_Ages <- df_counts %>%
  filter(VarComb == "Cckbr+/Tshr-") %>%
  arrange(desc(prop)) %>%
  pull(Age)

# Reorder the factor levels of Age
df_counts$Age <- factor(df_counts$Age, levels = ordered_Ages)
df_total_counts$Age <- factor(df_total_counts$Age, levels = ordered_Ages)

# Reorder the factor levels of combinations
df_counts$VarComb <- factor(df_counts$VarComb, levels = rev(c("Cckbr+/Tshr-", "Cckbr+/Tshr+", "Cckbr-/Tshr+", "Cckbr-/Tshr-")))

# Create a stacked bar plot
ggplot(df_counts, aes(x = Age, y = n, fill = VarComb)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("Cckbr+/Tshr+" = "cyan", "Cckbr+/Tshr-" = "green", "Cckbr-/Tshr+" = "royalblue", "Cckbr-/Tshr-" = "grey50")) +
  labs(x = "Hypothalamic subAge", y = "Number of cells", fill = "Expression") +
  scale_x_discrete(labels = c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45"))

ggplot(df_counts, aes(x = Age, y = prop, fill = VarComb)) +
  geom_bar(stat = "identity", color = "black", position = "fill") +
  scale_fill_manual(values = c("Cckbr+/Tshr-" = "green", "Cckbr+/Tshr+" = "cyan", "Cckbr-/Tshr+" = "royalblue", "Cckbr-/Tshr-" = "grey50")) +
  labs(x = "Hypothalamic subAge", y = "Proportion of cells", fill = "Expression") +
  scale_x_discrete(labels = c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45"))

ggplot(df_counts |> filter(!VarComb %in% c("Cckbr-/Tshr-", "Cckbr+/Tshr-")), aes(x = Age, y = n, fill = VarComb)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("Cckbr+/Tshr+" = "cyan", "Cckbr+/Tshr-" = "green", "Cckbr-/Tshr+" = "royalblue", "Cckbr-/Tshr-" = "grey50")) +
  labs(x = "Hypothalamic subAge", y = "Number of cells", fill = "Expression") +
  scale_x_discrete(labels = c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45"))

ggplot(df_counts |> filter(!VarComb %in% c("Cckbr-/Tshr-", "Cckbr+/Tshr-")), aes(x = Age, y = prop, fill = VarComb)) +
  geom_bar(stat = "identity", color = "black", position = "fill") +
  scale_fill_manual(values = c("Cckbr+/Tshr-" = "green", "Cckbr+/Tshr+" = "cyan", "Cckbr-/Tshr+" = "royalblue", "Cckbr-/Tshr-" = "grey50")) +
  labs(x = "Hypothalamic subAge", y = "Proportion of cells", fill = "Expression") +
  scale_x_discrete(labels = c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45"))
```

::: {.cell-output-display}
![](01-de_test-focus_pars_tub_files/figure-html/fig-Cckbr-tshr-bargraph-1.png){#fig-Cckbr-tshr-bargraph-1 width=672}
:::

::: {.cell-output-display}
![](01-de_test-focus_pars_tub_files/figure-html/fig-Cckbr-tshr-bargraph-2.png){#fig-Cckbr-tshr-bargraph-2 width=672}
:::

::: {.cell-output-display}
![](01-de_test-focus_pars_tub_files/figure-html/fig-Cckbr-tshr-bargraph-3.png){#fig-Cckbr-tshr-bargraph-3 width=672}
:::

::: {.cell-output-display}
![](01-de_test-focus_pars_tub_files/figure-html/fig-Cckbr-tshr-bargraph-4.png){#fig-Cckbr-tshr-bargraph-4 width=672}
:::
:::





# Calculate and plot chi2 test of independence between Sox2 and Gpr173 expression in hypothalamus across different developmental stages





::: {.cell}

```{.r .cell-code .hidden}
#| label: get-goi-sox2-Gpr173
sbs_mtx <- GetAssayData(object = srt.kim, layer = "counts", assay = "RNA")[c("Sox2", "Gpr173"), ] %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  select(Sox2, Gpr173) %>%
  dplyr::bind_cols(srt.kim@meta.data) %>%
  select(Age, Sox2, Gpr173) %>%
  mutate(
    Sox2_pos = Sox2 > 0,
    Gpr173_pos = Gpr173 > 0
  )

sbs_mtx %>% skimr::skim()
```

::: {.cell-output-display}

Table: Data summary

|                         |           |
|:------------------------|:----------|
|Name                     |Piped data |
|Number of rows           |128006     |
|Number of columns        |5          |
|_______________________  |           |
|Column type frequency:   |           |
|factor                   |1          |
|logical                  |2          |
|numeric                  |2          |
|________________________ |           |
|Group variables          |None       |


**Variable type: factor**

|skim_variable | n_missing| complete_rate|ordered | n_unique|top_counts                                    |
|:-------------|---------:|-------------:|:-------|--------:|:---------------------------------------------|
|Age           |         0|             1|FALSE   |       12|P45: 17025, E16: 16781, E14: 16543, P8: 13677 |


**Variable type: logical**

|skim_variable | n_missing| complete_rate| mean|count                   |
|:-------------|---------:|-------------:|----:|:-----------------------|
|Sox2_pos      |         0|             1| 0.14|FAL: 110343, TRU: 17663 |
|Gpr173_pos    |         0|             1| 0.02|FAL: 125963, TRU: 2043  |


**Variable type: numeric**

|skim_variable | n_missing| complete_rate| mean|   sd| p0| p25| p50| p75| p100|hist  |
|:-------------|---------:|-------------:|----:|----:|--:|---:|---:|---:|----:|:-----|
|Sox2          |         0|             1| 0.27| 1.21|  0|   0|   0|   0|   82|▇▁▁▁▁ |
|Gpr173        |         0|             1| 0.02| 0.13|  0|   0|   0|   0|    4|▇▁▁▁▁ |


:::
:::

::: {#cell-fig-sox2-Gpr173-stats .cell}

```{.r .cell-code .hidden}
#| label: fig-sox2-Gpr173-stats
#| fig-width: 8
#| fig-height: 24
write_csv(sbs_mtx, here(tables_dir, "Sox2-Gpr173-expression-status-between-Ages-on-evaluation-datasets.csv"))


# plot
grouped_ggpiestats(
  data = sbs_mtx,
  x = Gpr173_pos,
  y = Sox2_pos,
  grouping.var = Age,
  perc.k = 1,
  package = "ggsci",
  palette = "category10_d3",
  title.text = "Sox2 specification of Gpr173-positive hypothalamic development",
  caption.text = "Asterisks denote results from proportion tests; \n***: p < 0.001, ns: non-significant",
  plotgrid.args = list(nrow = 8)
)
```

::: {.cell-output-display}
![](01-de_test-focus_pars_tub_files/figure-html/fig-sox2-Gpr173-stats-1.png){#fig-sox2-Gpr173-stats width=768}
:::
:::





## Barplot of Sox2 and Gpr173 expression in hypothalamus across different developmental stages





::: {#cell-fig-sox2-Gpr173-bargraph .cell}

```{.r .cell-code .hidden}
#| label: fig-sox2-Gpr173-bargraph
#| fig-width: 7
#| fig-height: 6
# Create a new variable that combines Sox2_pos and Gpr173_pos
df <- sbs_mtx %>%
  mutate(Age = factor(Age, levels = c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45"), ordered = TRUE), VarComb = case_when(
    Sox2_pos & Gpr173_pos ~ "Sox2+/Gpr173+",
    Sox2_pos & !Gpr173_pos ~ "Sox2+/Gpr173-",
    !Sox2_pos & Gpr173_pos ~ "Sox2-/Gpr173+",
    !Sox2_pos & !Gpr173_pos ~ "Sox2-/Gpr173-"
  ))

# Calculate counts and proportions for each category
df_counts <- df %>%
  group_by(Age, VarComb) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n))
```

::: {.cell-output .cell-output-stderr .hidden}

```
`summarise()` has grouped output by 'Age'. You can override using the `.groups`
argument.
```


:::

```{.r .cell-code .hidden}
#| label: fig-sox2-Gpr173-bargraph
#| fig-width: 7
#| fig-height: 6
# Calculate the total counts for each category
df_total_counts <- df %>%
  group_by(Age) %>%
  summarise(total_n = n())

# Create a vector of Age names ordered by prop for Sox2+/Gpr173+ cases for each Age
ordered_Ages <- df_counts %>%
  filter(VarComb == "Sox2+/Gpr173-") %>%
  arrange(desc(prop)) %>%
  pull(Age)

# Reorder the factor levels of Age
df_counts$Age <- factor(df_counts$Age, levels = ordered_Ages)
df_total_counts$Age <- factor(df_total_counts$Age, levels = ordered_Ages)

# Reorder the factor levels of combinations
df_counts$VarComb <- factor(df_counts$VarComb, levels = rev(c("Sox2+/Gpr173-", "Sox2+/Gpr173+", "Sox2-/Gpr173+", "Sox2-/Gpr173-")))

# Create a stacked bar plot
ggplot(df_counts, aes(x = Age, y = n, fill = VarComb)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("Sox2+/Gpr173+" = "deeppink", "Sox2+/Gpr173-" = "red3", "Sox2-/Gpr173+" = "orchid", "Sox2-/Gpr173-" = "grey50")) +
  labs(x = "Hypothalamic subAge", y = "Number of cells", fill = "Expression") +
  scale_x_discrete(labels = c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45"))

ggplot(df_counts, aes(x = Age, y = prop, fill = VarComb)) +
  geom_bar(stat = "identity", color = "black", position = "fill") +
  scale_fill_manual(values = c("Sox2+/Gpr173-" = "red3", "Sox2+/Gpr173+" = "deeppink", "Sox2-/Gpr173+" = "orchid", "Sox2-/Gpr173-" = "grey50")) +
  labs(x = "Hypothalamic subAge", y = "Proportion of cells", fill = "Expression") +
  scale_x_discrete(labels = c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45"))

ggplot(df_counts |> filter(!VarComb %in% c("Sox2-/Gpr173-", "Sox2+/Gpr173-")), aes(x = Age, y = n, fill = VarComb)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("Sox2+/Gpr173+" = "deeppink", "Sox2+/Gpr173-" = "red3", "Sox2-/Gpr173+" = "orchid", "Sox2-/Gpr173-" = "grey50")) +
  labs(x = "Hypothalamic subAge", y = "Number of cells", fill = "Expression") +
  scale_x_discrete(labels = c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45"))

ggplot(df_counts |> filter(!VarComb %in% c("Sox2-/Gpr173-", "Sox2+/Gpr173-")), aes(x = Age, y = prop, fill = VarComb)) +
  geom_bar(stat = "identity", color = "black", position = "fill") +
  scale_fill_manual(values = c("Sox2+/Gpr173-" = "red3", "Sox2+/Gpr173+" = "deeppink", "Sox2-/Gpr173+" = "orchid", "Sox2-/Gpr173-" = "grey50")) +
  labs(x = "Hypothalamic subAge", y = "Proportion of cells", fill = "Expression") +
  scale_x_discrete(labels = c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45"))
```

::: {.cell-output-display}
![](01-de_test-focus_pars_tub_files/figure-html/fig-sox2-Gpr173-bargraph-1.png){#fig-sox2-Gpr173-bargraph-1 width=672}
:::

::: {.cell-output-display}
![](01-de_test-focus_pars_tub_files/figure-html/fig-sox2-Gpr173-bargraph-2.png){#fig-sox2-Gpr173-bargraph-2 width=672}
:::

::: {.cell-output-display}
![](01-de_test-focus_pars_tub_files/figure-html/fig-sox2-Gpr173-bargraph-3.png){#fig-sox2-Gpr173-bargraph-3 width=672}
:::

::: {.cell-output-display}
![](01-de_test-focus_pars_tub_files/figure-html/fig-sox2-Gpr173-bargraph-4.png){#fig-sox2-Gpr173-bargraph-4 width=672}
:::
:::





# Calculate and plot chi2 test of independence between Gpr173 and Tshr expression in hypothalamus across different developmental stages





::: {.cell}

```{.r .cell-code .hidden}
#| label: get-goi-Gpr173-tshr
sbs_mtx <- GetAssayData(object = srt.kim, layer = "counts", assay = "RNA")[c("Gpr173", "Tshr"), ] %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  select(Gpr173, Tshr) %>%
  dplyr::bind_cols(srt.kim@meta.data) %>%
  select(Age, Gpr173, Tshr) %>%
  mutate(
    Gpr173_pos = Gpr173 > 0,
    Tshr_pos = Tshr > 0
  )

sbs_mtx %>% skimr::skim()
```

::: {.cell-output-display}

Table: Data summary

|                         |           |
|:------------------------|:----------|
|Name                     |Piped data |
|Number of rows           |128006     |
|Number of columns        |5          |
|_______________________  |           |
|Column type frequency:   |           |
|factor                   |1          |
|logical                  |2          |
|numeric                  |2          |
|________________________ |           |
|Group variables          |None       |


**Variable type: factor**

|skim_variable | n_missing| complete_rate|ordered | n_unique|top_counts                                    |
|:-------------|---------:|-------------:|:-------|--------:|:---------------------------------------------|
|Age           |         0|             1|FALSE   |       12|P45: 17025, E16: 16781, E14: 16543, P8: 13677 |


**Variable type: logical**

|skim_variable | n_missing| complete_rate| mean|count                  |
|:-------------|---------:|-------------:|----:|:----------------------|
|Gpr173_pos    |         0|             1| 0.02|FAL: 125963, TRU: 2043 |
|Tshr_pos      |         0|             1| 0.00|FAL: 127714, TRU: 292  |


**Variable type: numeric**

|skim_variable | n_missing| complete_rate| mean|   sd| p0| p25| p50| p75| p100|hist  |
|:-------------|---------:|-------------:|----:|----:|--:|---:|---:|---:|----:|:-----|
|Gpr173        |         0|             1| 0.02| 0.13|  0|   0|   0|   0|    4|▇▁▁▁▁ |
|Tshr          |         0|             1| 0.00| 0.06|  0|   0|   0|   0|    5|▇▁▁▁▁ |


:::
:::

::: {#cell-fig-Gpr173-tshr-stats .cell}

```{.r .cell-code .hidden}
#| label: fig-Gpr173-tshr-stats
#| fig-width: 8
#| fig-height: 24
write_csv(sbs_mtx, here(tables_dir, "Gpr173-Tshr-expression-status-between-Ages-on-evaluation-datasets.csv"))


# plot
grouped_ggpiestats(
  data = sbs_mtx,
  x = Gpr173_pos,
  y = Tshr_pos,
  grouping.var = Age,
  perc.k = 1,
  package = "ggsci",
  palette = "category10_d3",
  title.text = "Gpr173 specification of Tshr-positive hypothalamic development",
  caption.text = "Asterisks denote results from proportion tests; \n***: p < 0.001, ns: non-significant",
  plotgrid.args = list(nrow = 8)
)
```

::: {.cell-output-display}
![](01-de_test-focus_pars_tub_files/figure-html/fig-Gpr173-tshr-stats-1.png){#fig-Gpr173-tshr-stats width=768}
:::
:::





## Barplot of Gpr173 and Tshr expression in hypothalamus across different developmental stages





::: {#cell-fig-Gpr173-tshr-bargraph .cell}

```{.r .cell-code .hidden}
#| label: fig-Gpr173-tshr-bargraph
#| fig-width: 7
#| fig-height: 6
# Create a new variable that combines Gpr173_pos and Tshr_pos
df <- sbs_mtx %>%
  mutate(Age = factor(Age, levels = c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45"), ordered = TRUE), VarComb = case_when(
    Gpr173_pos & Tshr_pos ~ "Gpr173+/Tshr+",
    Gpr173_pos & !Tshr_pos ~ "Gpr173+/Tshr-",
    !Gpr173_pos & Tshr_pos ~ "Gpr173-/Tshr+",
    !Gpr173_pos & !Tshr_pos ~ "Gpr173-/Tshr-"
  ))

# Calculate counts and proportions for each category
df_counts <- df %>%
  group_by(Age, VarComb) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n))
```

::: {.cell-output .cell-output-stderr .hidden}

```
`summarise()` has grouped output by 'Age'. You can override using the `.groups`
argument.
```


:::

```{.r .cell-code .hidden}
#| label: fig-Gpr173-tshr-bargraph
#| fig-width: 7
#| fig-height: 6
# Calculate the total counts for each category
df_total_counts <- df %>%
  group_by(Age) %>%
  summarise(total_n = n())

# Create a vector of Age names ordered by prop for Gpr173+/Tshr+ cases for each Age
ordered_Ages <- df_counts %>%
  filter(VarComb == "Gpr173+/Tshr-") %>%
  arrange(desc(prop)) %>%
  pull(Age)

# Reorder the factor levels of Age
df_counts$Age <- factor(df_counts$Age, levels = ordered_Ages)
df_total_counts$Age <- factor(df_total_counts$Age, levels = ordered_Ages)

# Reorder the factor levels of combinations
df_counts$VarComb <- factor(df_counts$VarComb, levels = rev(c("Gpr173+/Tshr-", "Gpr173+/Tshr+", "Gpr173-/Tshr+", "Gpr173-/Tshr-")))

# Create a stacked bar plot
ggplot(df_counts, aes(x = Age, y = n, fill = VarComb)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("Gpr173+/Tshr+" = "magenta4", "Gpr173+/Tshr-" = "orchid", "Gpr173-/Tshr+" = "royalblue", "Gpr173-/Tshr-" = "grey50")) +
  labs(x = "Hypothalamic subAge", y = "Number of cells", fill = "Expression") +
  scale_x_discrete(labels = c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45"))

ggplot(df_counts, aes(x = Age, y = prop, fill = VarComb)) +
  geom_bar(stat = "identity", color = "black", position = "fill") +
  scale_fill_manual(values = c("Gpr173+/Tshr-" = "orchid", "Gpr173+/Tshr+" = "magenta4", "Gpr173-/Tshr+" = "royalblue", "Gpr173-/Tshr-" = "grey50")) +
  labs(x = "Hypothalamic subAge", y = "Proportion of cells", fill = "Expression") +
  scale_x_discrete(labels = c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45"))

ggplot(df_counts |> filter(!VarComb %in% c("Gpr173-/Tshr-", "Gpr173+/Tshr-")), aes(x = Age, y = n, fill = VarComb)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("Gpr173+/Tshr+" = "magenta4", "Gpr173+/Tshr-" = "orchid", "Gpr173-/Tshr+" = "royalblue", "Gpr173-/Tshr-" = "grey50")) +
  labs(x = "Hypothalamic subAge", y = "Number of cells", fill = "Expression") +
  scale_x_discrete(labels = c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45"))

ggplot(df_counts |> filter(!VarComb %in% c("Gpr173-/Tshr-", "Gpr173+/Tshr-")), aes(x = Age, y = prop, fill = VarComb)) +
  geom_bar(stat = "identity", color = "black", position = "fill") +
  scale_fill_manual(values = c("Gpr173+/Tshr-" = "orchid", "Gpr173+/Tshr+" = "magenta4", "Gpr173-/Tshr+" = "royalblue", "Gpr173-/Tshr-" = "grey50")) +
  labs(x = "Hypothalamic subAge", y = "Proportion of cells", fill = "Expression") +
  scale_x_discrete(labels = c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45"))
```

::: {.cell-output-display}
![](01-de_test-focus_pars_tub_files/figure-html/fig-Gpr173-tshr-bargraph-1.png){#fig-Gpr173-tshr-bargraph-1 width=672}
:::

::: {.cell-output-display}
![](01-de_test-focus_pars_tub_files/figure-html/fig-Gpr173-tshr-bargraph-2.png){#fig-Gpr173-tshr-bargraph-2 width=672}
:::

::: {.cell-output-display}
![](01-de_test-focus_pars_tub_files/figure-html/fig-Gpr173-tshr-bargraph-3.png){#fig-Gpr173-tshr-bargraph-3 width=672}
:::

::: {.cell-output-display}
![](01-de_test-focus_pars_tub_files/figure-html/fig-Gpr173-tshr-bargraph-4.png){#fig-Gpr173-tshr-bargraph-4 width=672}
:::
:::

::: {.cell}

```{.r .cell-code .hidden}
sessioninfo::session_info()
```

::: {.cell-output .cell-output-stdout}

```
─ Session info ───────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.4.2 (2024-10-31)
 os       Ubuntu 22.04.5 LTS
 system   x86_64, linux-gnu
 ui       X11
 language en_US:en
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       Etc/UTC
 date     2025-01-30
 pandoc   3.1.11.1 @ /home/etretiakov/micromamba/bin/ (via rmarkdown)

─ Packages ───────────────────────────────────────────────────────────────────
 package          * version     date (UTC) lib source
 abind              1.4-8       2024-09-12 [2] RSPM
 base64enc          0.1-3       2015-07-28 [2] RSPM (R 4.4.0)
 BayesFactor        0.9.12-4.7  2024-01-24 [2] RSPM
 bayestestR         0.15.0      2024-10-17 [2] RSPM (R 4.4.0)
 beeswarm           0.4.0       2021-06-01 [2] RSPM (R 4.4.0)
 BiocManager        1.30.25     2024-08-28 [2] RSPM (R 4.4.0)
 bit                4.5.0       2024-09-20 [2] RSPM
 bit64              4.5.2       2024-09-22 [2] RSPM
 circlize           0.4.16      2024-12-04 [2] Github (jokergoo/circlize@9b21578)
 cli                3.6.3       2024-06-21 [2] RSPM (R 4.4.0)
 cluster            2.1.6       2023-12-01 [2] CRAN (R 4.4.2)
 coda               0.19-4.1    2024-01-31 [2] RSPM
 codetools          0.2-20      2024-03-31 [2] CRAN (R 4.4.2)
 colorspace         2.1-1       2024-07-26 [2] RSPM (R 4.4.0)
 correlation        0.8.6       2024-10-26 [2] RSPM (R 4.4.0)
 cowplot          * 1.1.3       2024-01-22 [2] RSPM
 crayon             1.5.3       2024-06-20 [2] RSPM (R 4.4.0)
 data.table         1.16.2      2024-10-10 [2] RSPM
 datawizard         0.13.0.17   2024-12-04 [2] Github (easystats/datawizard@25f8ec4)
 deldir             2.0-4       2024-02-28 [2] RSPM (R 4.4.0)
 digest             0.6.37      2024-08-19 [2] RSPM (R 4.4.0)
 dotCall64          1.2         2024-10-04 [2] RSPM
 dplyr            * 1.1.4       2023-11-17 [2] RSPM (R 4.4.0)
 effectsize         0.8.9       2024-07-03 [2] RSPM (R 4.4.0)
 emmeans            1.10.5      2024-10-14 [2] RSPM
 estimability       1.5.1       2024-05-12 [2] RSPM (R 4.4.0)
 evaluate           1.0.1       2024-10-10 [2] RSPM (R 4.4.0)
 fansi              1.0.6       2023-12-08 [2] RSPM (R 4.4.0)
 farver             2.1.2       2024-05-13 [2] RSPM (R 4.4.0)
 fastDummies        1.7.4       2024-08-16 [2] RSPM
 fastmap            1.2.0       2024-05-15 [2] RSPM (R 4.4.0)
 fitdistrplus       1.2-1       2024-07-12 [2] RSPM (R 4.4.0)
 forcats          * 1.0.0       2023-01-29 [2] RSPM
 future           * 1.34.0      2024-07-29 [2] RSPM
 future.apply       1.11.3      2024-10-27 [2] RSPM
 generics           0.1.3       2022-07-05 [2] RSPM (R 4.4.0)
 ggbeeswarm         0.7.2       2024-12-04 [2] Github (eclarke/ggbeeswarm@14ef76c)
 ggplot2          * 3.5.1       2024-04-23 [2] RSPM (R 4.4.0)
 ggprism            1.0.5       2024-12-04 [2] Github (csdaw/ggprism@b6e6c0e)
 ggrastr            1.0.2       2024-12-04 [2] Github (VPetukhov/ggrastr@50ca3e0)
 ggrepel            0.9.6.9999  2024-12-04 [2] Github (slowkow/ggrepel@e72a66d)
 ggridges           0.5.6       2024-01-23 [2] RSPM
 ggsci              3.2.0       2024-12-04 [2] Github (nanxstats/ggsci@b5bf1fd)
 ggside             0.3.1.9999  2024-12-04 [2] Github (jtlandis/ggside@47c24a4)
 ggstatsplot      * 0.12.5.9000 2024-12-04 [2] Github (IndrajeetPatil/ggstatsplot@d312b9f)
 GlobalOptions      0.1.2       2020-06-10 [2] RSPM (R 4.4.0)
 globals            0.16.3      2024-03-08 [2] RSPM
 glue               1.8.0       2024-09-30 [2] RSPM (R 4.4.0)
 goftest            1.2-3       2021-10-07 [2] RSPM
 gridExtra          2.3         2017-09-09 [2] RSPM
 gtable             0.3.6       2024-10-25 [2] RSPM (R 4.4.0)
 here             * 1.0.1       2020-12-13 [2] RSPM
 hexbin           * 1.28.5      2024-11-13 [2] RSPM (R 4.4.0)
 hms                1.1.3       2023-03-21 [2] RSPM
 htmltools          0.5.8.1     2024-04-04 [2] RSPM (R 4.4.0)
 htmlwidgets        1.6.4       2023-12-06 [2] RSPM (R 4.4.0)
 httpuv             1.6.15      2024-03-26 [2] RSPM (R 4.4.0)
 httr               1.4.7       2023-08-15 [2] RSPM (R 4.4.0)
 ica                1.0-3       2022-07-08 [2] RSPM
 igraph             2.1.1       2024-10-19 [2] RSPM (R 4.4.0)
 insight            1.0.0.2     2024-12-04 [2] Github (easystats/insight@8e78b12)
 irlba              2.3.5.1     2022-10-03 [2] RSPM
 janitor            2.2.0.9000  2024-12-04 [2] Github (sfirke/janitor@6ee7919)
 jsonlite           1.8.9       2024-09-20 [2] RSPM (R 4.4.0)
 KernSmooth         2.23-24     2024-05-17 [2] CRAN (R 4.4.2)
 knitr              1.49        2024-11-08 [2] RSPM
 labeling           0.4.3       2023-08-29 [2] RSPM (R 4.4.0)
 later              1.4.1       2024-11-27 [2] RSPM (R 4.4.0)
 lattice            0.22-6      2024-03-20 [2] CRAN (R 4.4.2)
 lazyeval           0.2.2       2019-03-15 [2] RSPM (R 4.4.0)
 leiden             0.4.3.1     2023-11-17 [2] RSPM
 lifecycle          1.0.4       2023-11-07 [2] RSPM (R 4.4.0)
 listenv            0.9.1       2024-01-29 [2] RSPM
 lmtest             0.9-40      2022-03-21 [2] RSPM (R 4.4.0)
 lubridate        * 1.9.3       2023-09-27 [2] RSPM
 magrittr         * 2.0.3       2022-03-30 [2] RSPM (R 4.4.0)
 MASS               7.3-61      2024-06-13 [2] CRAN (R 4.4.2)
 Matrix             1.7-1       2024-10-18 [2] CRAN (R 4.4.2)
 MatrixModels       0.5-3       2023-11-06 [2] RSPM
 matrixStats        1.4.1       2024-09-08 [2] RSPM (R 4.4.0)
 mgcv               1.9-1       2023-12-21 [2] CRAN (R 4.4.2)
 mime               0.12        2021-09-28 [2] RSPM (R 4.4.0)
 miniUI             0.1.1.1     2018-05-18 [2] RSPM
 multcomp           1.4-26      2024-07-18 [2] RSPM
 munsell            0.5.1       2024-04-01 [2] RSPM (R 4.4.0)
 mvtnorm            1.3-2       2024-11-04 [2] RSPM
 nlme               3.1-166     2024-08-14 [2] CRAN (R 4.4.2)
 paletteer          1.6.0       2024-01-21 [2] RSPM
 parallelly         1.39.0      2024-11-07 [2] RSPM
 parameters         0.24.0.3    2024-12-04 [2] Github (easystats/parameters@eff54e5)
 patchwork        * 1.3.0.9000  2024-12-04 [2] Github (thomasp85/patchwork@2695a9f)
 pbapply            1.7-2       2023-06-27 [2] RSPM
 pillar             1.9.0       2023-03-22 [2] RSPM (R 4.4.0)
 pkgconfig          2.0.3       2019-09-22 [2] RSPM (R 4.4.0)
 plotly             4.10.4      2024-01-13 [2] RSPM
 plyr               1.8.9       2023-10-02 [2] RSPM
 png                0.1-8       2022-11-29 [2] RSPM
 polyclip           1.10-7      2024-07-23 [2] RSPM (R 4.4.0)
 prismatic          1.1.2       2024-04-10 [2] RSPM
 progressr          0.15.1      2024-11-22 [2] RSPM
 promises           1.3.2       2024-11-28 [2] RSPM (R 4.4.0)
 purrr            * 1.0.2       2023-08-10 [2] RSPM (R 4.4.0)
 R.methodsS3        1.8.2       2022-06-13 [2] RSPM (R 4.4.0)
 R.oo               1.27.0      2024-11-01 [2] RSPM (R 4.4.0)
 R.utils            2.12.3      2023-11-18 [2] RSPM (R 4.4.0)
 R6                 2.5.1       2021-08-19 [2] RSPM (R 4.4.0)
 RANN               2.6.2       2024-08-25 [2] RSPM (R 4.4.0)
 RColorBrewer     * 1.1-3       2022-04-03 [2] RSPM
 Rcpp               1.0.13-1    2024-11-02 [2] RSPM (R 4.4.0)
 RcppAnnoy          0.0.22      2024-01-23 [2] RSPM
 RcppHNSW           0.6.0       2024-02-04 [2] RSPM
 readr            * 2.1.5       2024-01-10 [2] RSPM
 rematch2           2.1.2       2020-05-01 [2] RSPM
 remotes            2.5.0       2024-03-17 [2] RSPM
 repr               1.1.7       2024-03-22 [2] RSPM
 reshape2           1.4.4       2020-04-09 [2] RSPM
 reticulate         1.40.0.9000 2024-12-04 [2] Github (rstudio/reticulate@61f0fa4)
 rhdf5              2.50.0      2024-10-29 [2] RSPM (R 4.4.2)
 rhdf5filters       1.18.0      2024-10-29 [2] RSPM (R 4.4.2)
 Rhdf5lib           1.28.0      2024-10-29 [2] RSPM (R 4.4.2)
 rlang              1.1.4       2024-06-04 [2] RSPM (R 4.4.0)
 rmarkdown          2.29        2024-11-04 [2] RSPM
 ROCR               1.0-11      2020-05-02 [2] RSPM (R 4.4.0)
 rprojroot          2.0.4       2023-11-05 [2] RSPM (R 4.4.0)
 RSpectra           0.16-2      2024-07-18 [2] RSPM
 rstudioapi         0.17.1      2024-10-22 [2] RSPM
 rsvd               1.0.5       2021-04-16 [2] RSPM (R 4.4.0)
 Rtsne              0.17        2023-12-07 [2] RSPM (R 4.4.0)
 sandwich           3.1-1       2024-09-15 [2] RSPM
 scales             1.3.0       2023-11-28 [2] RSPM (R 4.4.0)
 scattermore        1.2         2023-06-12 [2] RSPM
 scCustomize      * 3.0.1       2025-01-30 [2] Github (samuel-marsh/scCustomize@3299b95)
 schard             0.0.1       2024-12-04 [2] Github (cellgeni/schard@c22b46d)
 sctransform        0.4.1       2023-10-19 [2] RSPM
 sessioninfo        1.2.2       2021-12-06 [2] RSPM
 Seurat           * 5.1.0       2024-12-04 [2] Github (satijalab/seurat@1549dcb)
 SeuratObject     * 5.0.99.9001 2024-12-04 [2] Github (satijalab/seurat-object@42e53ba)
 SeuratWrappers   * 0.4.0       2024-12-04 [2] Github (satijalab/seurat-wrappers@a1eb0d8)
 shape              1.4.6.1     2024-02-23 [2] RSPM
 shiny              1.9.1       2024-08-01 [2] RSPM (R 4.4.0)
 skimr            * 2.1.5       2024-12-04 [2] Github (ropensci/skimr@d5126aa)
 snakecase          0.11.1      2023-08-27 [2] RSPM (R 4.4.0)
 sp               * 2.1-4       2024-04-30 [2] RSPM
 spam               2.11-0      2024-10-03 [2] RSPM
 spatstat.data      3.1-4       2024-11-15 [2] RSPM (R 4.4.0)
 spatstat.explore   3.3-3       2024-10-22 [2] RSPM
 spatstat.geom      3.3-4       2024-11-18 [2] RSPM (R 4.4.0)
 spatstat.random    3.3-2       2024-09-18 [2] RSPM (R 4.4.0)
 spatstat.sparse    3.1-0       2024-06-21 [2] RSPM
 spatstat.univar    3.1-1       2024-11-05 [2] RSPM (R 4.4.0)
 spatstat.utils     3.1-1       2024-11-03 [2] RSPM (R 4.4.0)
 statsExpressions   1.6.1       2024-10-31 [2] RSPM
 stringi            1.8.4       2024-05-06 [2] RSPM (R 4.4.0)
 stringr          * 1.5.1       2023-11-14 [2] RSPM (R 4.4.0)
 survival           3.7-0       2024-06-05 [2] CRAN (R 4.4.2)
 tensor             1.5         2012-05-05 [2] RSPM
 TH.data            1.1-2       2023-04-17 [2] RSPM (R 4.4.0)
 tibble           * 3.2.1       2023-03-20 [2] RSPM (R 4.4.0)
 tidyr            * 1.3.1       2024-01-24 [2] RSPM (R 4.4.0)
 tidyselect         1.2.1       2024-03-11 [2] RSPM (R 4.4.0)
 tidyverse        * 2.0.0.9000  2024-12-04 [2] Github (tidyverse/tidyverse@c06a3c9)
 timechange         0.3.0       2024-01-18 [2] RSPM
 tzdb               0.4.0       2023-05-12 [2] RSPM
 utf8               1.2.4       2023-10-22 [2] RSPM (R 4.4.0)
 uwot               0.2.2       2024-04-21 [2] RSPM (R 4.4.0)
 vctrs              0.6.5       2023-12-01 [2] RSPM (R 4.4.0)
 vipor              0.4.7       2023-12-18 [2] RSPM (R 4.4.0)
 viridis          * 0.6.5       2024-01-29 [2] RSPM
 viridisLite      * 0.4.2       2023-05-02 [2] RSPM (R 4.4.0)
 vroom              1.6.5       2023-12-05 [2] RSPM
 withr              3.0.2       2024-10-28 [2] RSPM (R 4.4.0)
 xfun               0.49        2024-10-31 [2] RSPM (R 4.4.0)
 xtable             1.8-4       2019-04-21 [2] RSPM (R 4.4.0)
 yaml               2.3.10      2024-07-26 [2] RSPM
 zeallot            0.1.0       2018-01-28 [2] RSPM
 zoo                1.8-12      2023-04-13 [2] RSPM (R 4.4.0)

 [1] /home/etretiakov/R/x86_64-pc-linux-gnu-library/4.4
 [2] /opt/R/4.4.2/lib/R/library

─ Python configuration ───────────────────────────────────────────────────────
 python:         /opt/conda/bin/python
 libpython:      /opt/conda/lib/libpython3.12.so
 pythonhome:     /opt/conda:/opt/conda
 version:        3.12.8 | packaged by conda-forge | (main, Dec  5 2024, 14:24:40) [GCC 13.3.0]
 numpy:          /opt/conda/lib/python3.12/site-packages/numpy
 numpy_version:  1.26.4
 
 NOTE: Python version was forced by RETICULATE_PYTHON

──────────────────────────────────────────────────────────────────────────────
```


:::
:::
