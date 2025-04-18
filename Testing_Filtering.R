---
title: "Testing_Filtering"
output: html_document
date: "2025-04-17"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
```{r}
#including sanjana's script here, and working on plotting things in the subsequent cells
library(ArchR)
library(Giotto)

addArchRGenome("hg38") 
addArchRThreads(threads = 16)

arrow_file <-"/projectnb/paxlab/carlos/ArchR_Gastric/ArrowFiles/D01942_NG05827.arrow"
position_files <- "/projectnb/paxlab/carlos/ArchR_Gastric/ArrowFiles/tissue_positions_list.csv"

# Load ArchR project
proj <- ArchRProject(ArrowFiles = arrow_file, outputDirectory = "ArchR_Output")

# Load spatial data
spatial_data <- fread(position_files)
setnames(spatial_data, old = c("V1", "V2", "V3", "V4", "V5", "V6"), 
         new = c("barcode", "in_tissue", "array_row", "array_col", "x_spatial", "y_spatial"))
run_id <- 'D01942_NG05827'
spatial_data[, barcode := paste0(run_id, "#", barcode, "-1")]
# Get on-tissue barcodes (16,314)
filtered_barcodes <- spatial_data[in_tissue == 1]$barcode

proj_filtered <- proj[proj$cellNames %in% filtered_barcodes]

# Get the gene score matrix
geneScoreMatrix <- getMatrixFromProject(proj_filtered, useMatrix = "GeneScoreMatrix")
genescore_matrix <- geneScoreMatrix@assays@data$GeneScoreMatrix

# Get gene names from annotations
gene_ano <- getGeneAnnotation(proj_filtered)
gene_ids <- gene_ano$genes$gene_id
gene_symbols <- gene_ano$genes$symbol

# Assign gene names to sparse matrix
rownames(genescore_matrix) <- rowData(geneScoreMatrix)$name 

proj_filtered <- addIterativeLSI( proj_filtered,
                                  name = "IterativeLSI",
                                  force = TRUE
)

proj_filtered <- addImputeWeights(
  ArchRProj = proj_filtered,
  reducedDims = "IterativeLSI"
)
## added this line
spatial_coords <- spatial_data[barcode %in% proj_filtered$cellNames, .(barcode, x_spatial, y_spatial)]

archrImputedMatrix <- imputeMatrix(
  mat = genescore_matrix,
  imputeWeights = getImputeWeights(proj_filtered)
)

giotto_obj <- createGiottoObject(
  expression = list(
    spatial_atac = list(
      raw = genescore_matrix,
      archr_imputed = archrImputedMatrix
    )
  ),
  spatial_locs = spatial_coords,
  instructions = createGiottoInstructions(python_path = NULL)
)
```
```{r}
proj <- addDoubletScores(
  input = proj,
  k = 10,                          # Number of neighbors
  knnMethod = "UMAP",             # or "LSI" if no UMAP yet
  LSIMethod = 1
)

```
```{r}
# Basic QC thresholds (adjust to your dataset as needed)
min_tss <- 8
min_frags <- 1000
max_doublet_score <- 1.5  # Optional, only if you've run addDoubletScores()

# Make sure doublet scores are computed if you're going to use them
# proj <- addDoubletScores(proj)

# Apply filtering using QC metrics
proj_filtered_qc <- subsetArchRProject(
  ArchRProj = proj,
  cells = proj$cellNames[
    proj$TSSEnrichment >= min_tss &
    proj$nFrags >= min_frags &
    proj$DoubletEnrichment <= max_doublet_score
  ],
  outputDirectory = "ArchR_Output_Filtered_QC"
)

```
```{r}
proj_filtered_qc
```

```{r}
d <- plotTSSEnrichment(proj_filtered_qc, returnDF = TRUE)
head(d)
```
```{r}
d_ <- getCellColData(proj_filtered_qc, select = c("TSSEnrichment", "nFrags")) %>% 
  as.data.frame() %>% 
  mutate(log10_nFrags = log10(nFrags))

```

```{r}
ggplot(d_, aes(x = log10_nFrags, y = TSSEnrichment)) +
  geom_point(alpha = 0.3, size = 0.5, color = "dodgerblue") +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", color = NA, alpha = 0.4) +
  geom_vline(xintercept = 3, linetype = "dashed") +
  geom_hline(yintercept = 10, linetype = "dashed") +
  scale_fill_viridis_c() +
  labs(
    title = "TSS Enrichment vs Log10 Unique Fragments",
    x = "Log10(Unique Fragments)",
    y = "TSS Enrichment"
  ) +
  theme_minimal()
```
```{r}
proj_filtered <- addDoubletScores(
  input = proj_filtered,
  k = 10,                          # Number of neighbors
  knnMethod = "UMAP",             # or "LSI" if no UMAP yet
  LSIMethod = 1
)

```
```{r}

# Basic QC thresholds (adjust to your dataset as needed)
min_tss <- 8
min_frags <- 1000
max_doublet_score <- 1.5  # Optional, only if you've run addDoubletScores()

# Make sure doublet scores are computed if you're going to use them
# proj <- addDoubletScores(proj)

# Apply filtering using QC metrics
proj_filtered_wqc <- subsetArchRProject(
  ArchRProj = proj_filtered,
  cells = proj_filtered$cellNames[
    proj_filtered$TSSEnrichment >= min_tss &
    proj_filtered$nFrags >= min_frags &
    proj_filtered$DoubletEnrichment <= max_doublet_score
  ],
  outputDirectory ="ArchR_Output_",
  force = TRUE
)


```
```{r}
d__ <- getCellColData(proj_filtered_wqc, select = c("TSSEnrichment", "nFrags")) %>% 
  as.data.frame() %>% 
  mutate(log10_nFrags = log10(nFrags))

```
```{r}
d__
```
```{r}
ggplot(d__, aes(x = log10_nFrags, y = TSSEnrichment)) +
  geom_point(alpha = 0.3, size = 0.5, color = "dodgerblue") +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", color = NA, alpha = 0.4) +
  geom_vline(xintercept = 3, linetype = "dashed") +
  geom_hline(yintercept = 10, linetype = "dashed") +
  scale_fill_viridis_c() +
  labs(
    title = "TSS Enrichment vs Log10 Unique Fragments",
    x = "Log10(Unique Fragments)",
    y = "TSS Enrichment"
  ) +
  theme_minimal()
```
```{r}
# Basic QC thresholds (adjust to your dataset as needed)
min_tss <- 4
min_frags <- 1000
max_doublet_score <- 1.5  # Optional, only if you've run addDoubletScores()

# Make sure doublet scores are computed if you're going to use them
# proj <- addDoubletScores(proj)

# Apply filtering using QC metrics
proj_filtered_wqc <- subsetArchRProject(
  ArchRProj = proj_filtered,
  cells = proj_filtered$cellNames[
    proj_filtered$TSSEnrichment >= min_tss &
    proj_filtered$nFrags >= min_frags &
    proj_filtered$DoubletEnrichment <= max_doublet_score
  ],
  outputDirectory ="ArchR_Output_",
  force = TRUE
)
d__ <- getCellColData(proj_filtered_wqc, select = c("TSSEnrichment", "nFrags")) %>% 
  as.data.frame() %>% 
  mutate(log10_nFrags = log10(nFrags))

ggplot(d__, aes(x = log10_nFrags, y = TSSEnrichment)) +
  geom_point(alpha = 0.3, size = 0.5, color = "dodgerblue") +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", color = NA, alpha = 0.4) +
  geom_vline(xintercept = 3, linetype = "dashed") +
  geom_hline(yintercept = 10, linetype = "dashed") +
  scale_fill_viridis_c() +
  labs(
    title = "TSS Enrichment vs Log10 Unique Fragments",
    x = "Log10(Unique Fragments)",
    y = "TSS Enrichment"
  ) +
  theme_minimal()
```


```{r}
# Basic QC thresholds (adjust to your dataset as needed)
min_tss <- 8
min_frags <- 0
max_doublet_score <- 1.5  # Optional, only if you've run addDoubletScores()

# Make sure doublet scores are computed if you're going to use them
# proj <- addDoubletScores(proj)

# Apply filtering using QC metrics
proj_filtered_wqc <- subsetArchRProject(
  ArchRProj = proj_filtered,
  cells = proj_filtered$cellNames[
    proj_filtered$TSSEnrichment >= min_tss &
    proj_filtered$nFrags >= min_frags &
    proj_filtered$DoubletEnrichment <= max_doublet_score
  ],
  outputDirectory ="ArchR_Output_",
  force = TRUE
)
d__ <- getCellColData(proj_filtered_wqc, select = c("TSSEnrichment", "nFrags")) %>% 
  as.data.frame() %>% 
  mutate(log10_nFrags = log10(nFrags))

ggplot(d__, aes(x = log10_nFrags, y = TSSEnrichment)) +
  geom_point(alpha = 0.3, size = 0.5, color = "dodgerblue") +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", color = NA, alpha = 0.4) +
  geom_vline(xintercept = 3, linetype = "dashed") +
  geom_hline(yintercept = 10, linetype = "dashed") +
  scale_fill_viridis_c() +
  labs(
    title = "TSS Enrichment vs Log10 Unique Fragments",
    x = "Log10(Unique Fragments)",
    y = "TSS Enrichment"
  ) +
  theme_minimal()
```
```{r}
# Basic QC thresholds (adjust to your dataset as needed)
min_tss <- 8
min_frags <- 1000
max_doublet_score <- 1000  # Optional, only if you've run addDoubletScores()

# Make sure doublet scores are computed if you're going to use them
# proj <- addDoubletScores(proj)

# Apply filtering using QC metrics
proj_filtered_wqc <- subsetArchRProject(
  ArchRProj = proj_filtered,
  cells = proj_filtered$cellNames[
    proj_filtered$TSSEnrichment >= min_tss &
    proj_filtered$nFrags >= min_frags &
    proj_filtered$DoubletEnrichment <= max_doublet_score
  ],
  outputDirectory ="ArchR_Output_",
  force = TRUE
)
d__ <- getCellColData(proj_filtered_wqc, select = c("TSSEnrichment", "nFrags")) %>% 
  as.data.frame() %>% 
  mutate(log10_nFrags = log10(nFrags))

ggplot(d__, aes(x = log10_nFrags, y = TSSEnrichment)) +
  geom_point(alpha = 0.3, size = 0.5, color = "dodgerblue") +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", color = NA, alpha = 0.4) +
  geom_vline(xintercept = 3, linetype = "dashed") +
  geom_hline(yintercept = 10, linetype = "dashed") +
  scale_fill_viridis_c() +
  labs(
    title = "TSS Enrichment vs Log10 Unique Fragments",
    x = "Log10(Unique Fragments)",
    y = "TSS Enrichment"
  ) +
  theme_minimal()
```

```{r}
# Get the gene score matrix from the QC-filtered project
geneScoreMatrix_wqc <- getMatrixFromProject(proj_filtered_wqc, useMatrix = "GeneScoreMatrix")
genescore_matrix_wqc <- geneScoreMatrix_wqc@assays@data$GeneScoreMatrix

# Set gene names
rownames(genescore_matrix_wqc) <- rowData(geneScoreMatrix_wqc)$name

# Impute matrix using weights from the new filtered project
proj_filtered_wqc <- addIterativeLSI(
  ArchRProj = proj_filtered_wqc,
  name = "IterativeLSI",
  force = TRUE
)

proj_filtered_wqc <- addImputeWeights(
  ArchRProj = proj_filtered_wqc,
  reducedDims = "IterativeLSI"
)

archrImputedMatrix_wqc <- imputeMatrix(
  mat = genescore_matrix_wqc,
  imputeWeights = getImputeWeights(proj_filtered_wqc)
)

# Get new spatial coords (match barcodes to the filtered project)
spatial_coords_wqc <- spatial_data[barcode %in% proj_filtered_wqc$cellNames, .(barcode, x_spatial, y_spatial)]

# Create Giotto object
giotto_obj_wqc <- createGiottoObject(
  expression = list(
    spatial_atac = list(
      raw = genescore_matrix_wqc,
      archr_imputed = archrImputedMatrix_wqc
    )
  ),
  spatial_locs = spatial_coords_wqc,
  instructions = createGiottoInstructions(python_path = NULL)
)

```
```{r}
d__
#~15k when treshold is set to 0, similar when set to 2, closer to 14k with 4

```
