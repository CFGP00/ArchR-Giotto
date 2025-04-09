library(ArchR)
library(Giotto)

addArchRGenome("hg38") 
addArchRThreads(threads = 16)

arrow_file <-".../Dries_259_ArchRProject/ArrowFiles/D01942_NG05827.arrow"
position_files <- ".../tissue_positions_list.csv"

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

# Create spatial co-ords with only the barcode,x and y values
spatial_coords <- spatial_data[barcode %in% proj_filtered$cellNames, .(barcode, x_spatial, y_spatial)]

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
