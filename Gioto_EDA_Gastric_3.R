# Load Required Libraries
library(Giotto)
library(data.table)
library(ggplot2)

# Define Directories
output_dir <- "/projectnb/paxlab/carlos/Giotto_Results/"

# Ensure Output Directory Exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define File Path (GeneScoreMatrix from ArchR)
gene_score_file <- "/projectnb/paxlab/carlos/ArchR_Gastric/TFMatrix_cleaned_.csv"

# Load Data
gene_scores <- fread(gene_score_file, data.table = FALSE)

# Convert rownames back
rownames(gene_scores) <- gene_scores[,1]
gene_scores <- gene_scores[,-1]  

# Create Giotto Object
giotto_obj <- createGiottoObject(expression = gene_scores)

# Save the Giotto object before any modifications
#saveGiotto(giotto_obj, save_dir = output_dir)

print("Giotto object loaded and saved! Ready for preprocessing.")

# Aggregate, Normalize, and Filter Giotto Datta:

#giotto_obj <- calculateOverlapRaster(giotto_obj
                                     #spatial_info = "cell",
                                     #feat_info = "rna"
                         #            )

#giotto_obj <- overlapToMatrix(giotto_obj)
go <- giotto_obj
go <- addSpatialCentroidLocations(go)

#filter/normalize

filterDistributions(go, detection = "feats")
filterDistributions(go, detection = "cells")

go <- filterGiotto(go, feat_det_in_min_cells = 100,
                   min_det_feats_per_cell = 20,
                   expression_threshold = 1)

go <- normalizeGiotto(go)

#getting stats
go <- addStatistics(go)

# Step 6: Dimension Reduction

# Calculate highly variable features
go <- calculateHVF(gobject = go)
cat(fDataDT(go)[, sum(hvf == "yes")], "hvf found")

# Only 18 hvf found -> better to use ALL genes -> feats_to_use = NULL
go <- runPCA(
  gobject = go,
  spat_unit = "cell",
  expression_values = "scaled",
  feats_to_use = NULL,
  scale_unit = FALSE,
  center = FALSE
)

#PCA part

# Visualize Screeplot and PCA
screePlot(go,
          ncp = 20,
          save_param = list(
            save_name = "sg_screePlot"
          )
)
plotPCA(gobject = go,
        spat_unit = "cell",
        dim_reduction_name = "pca",
        dim1_to_use = 1,
        dim2_to_use = 2
)

#Plotting TSNE and UMAP:

# Run and Plot tSNE and UMAP
go <- runtSNE(go,
              dimensions_to_use = 1:10,
              spat_unit = "cell",
              check_duplicates = FALSE
)
go <- runUMAP(go,
              dimensions_to_use = 1:10,
              spat_unit = "cell"
)
plotTSNE(go,
         point_size = 0.01,
         save_param = list(
           save_name = "sg_tSNE"
         )
)
plotUMAP(go,
         point_size = 0.01,
         save_param = list(
           save_name = "sg_UMAP"
         )
)

# Part 7: Clustering
# Clustering and UMAP cluster visualization
go <- createNearestNetwork(go,
                           dimensions_to_use = 1:10,
                           k = 10,
                           spat_unit = "cell"
)
go <- doLeidenCluster(go,
                      resolution = 0.25,
                      n_iterations = 100,
                      spat_unit = "cell"
)
# Plot Leiden clusters onto UMAP
plotUMAP(
  gobject = go,
  spat_unit = "cell",
  cell_color = "leiden_clus",
  show_legend = FALSE,
  point_size = 0.01,
  point_shape = "no_border",
  save_param = list(save_name = "sg_umap_leiden")
)

#spatial 7.2 spatial leiden clustering:

# Plot Leiden clusters onto spatial image plot
my_spatPlot <- spatPlot2D(
  gobject = go,
  spat_unit = "cell",
  cell_color = "leiden_clus",
  point_size = 0.4,
  point_shape = "no_border",
  show_legend = TRUE,
  #image_name = gimg,
  save_param = list(
    save_name = "sg_spat_leiden",
    base_width = 15,
    base_height = 15
  )
)

# 8 Cell Type Marker Gene Detection

# Identify gene markers_gini per cluster- changed method to existing method
markers_gini <- findMarkers_one_vs_all(
  gobject = go,
  method = "gini",
  expression_values = "normalized",
  cluster_column = "leiden_clus",
  min_feats = 1, rank_score = 2
)

# Display details about the marker genes
#markers_gini[, head(.SD, 2), by = "cluster", .SDcols = colnames(markers_gini)]

markers_gini[, head(.SD, 2), by = "cluster", .SDcols = colnames(markers_gini)]

# Violinplots to show marker expression
topgenes_gini <- unique(markers_gini[, head(.SD, 2), by = "cluster"])

violinPlot(go, 
           feats = topgenes_gini$feats[1:10], 
           cluster_column = "leiden_clus")

violinPlot(go, 
           feats = topgenes_gini$feats[11:20], 
           cluster_column = "leiden_clus")

# Known markers_gini to Annotate Giotto - these might change ? ?
selected_genes <- c(
  "My12", "G6pc", "Ppp1r1a", "Grik5", "Hsd11b2", "Rhbg", "Mapk11",
  "Egf17", "Gpr55", "Acsm2", "Tpm2", "D1c1", "Shisa3",
  "Tspan2", "Sox17", "Eef2", "Cd79b", "Ctss", "Serpina1f", "Cyp51"
)

cell_metadata <- pDataDT(go)
cluster_order <- unique(cell_metadata$leiden_clus)

# Plot markers_gini to clusters heatmap
plotMetaDataHeatmap(go,
                    expression_values = "scaled",
                    metadata_cols = c("leiden_clus"),
                    selected_feats = selected_genes,
                    custom_feat_order = rev(selected_genes),
                    custom_cluster_order = cluster_order
)
# spatial Gene Expression Patterns
# 9.1
plotStatDelaunayNetwork(gobject = sg, 
                        maximum_distance = 250)

sg <- createSpatialNetwork(
  gobject = go, 
  minimum_k = 2,
  maximum_distance_delaunay = 250
)

sg <- createSpatialNetwork(
  gobject = go, 
  minimum_k = 2,
  method = "kNN", 
  k = 10
)
#9.2

km_spatialfeats <- binSpect(go)

spatFeatPlot2D(sg,
               expression_values = "scaled",
               feats = km_spatialfeats[1:4]$feats,
               point_shape = "no_border",
               show_network = FALSE, 
               network_color = "lightgrey",
               point_size = 0.5,
               cow_n_col = 2
)
#9.3
rank_spatialfeats <- binSpect(go, 
                              bin_method = "rank")

spatFeatPlot2D(go,
               expression_values = "scaled",
               feats = rank_spatialfeats[1:4]$feats,
               point_shape = "no_border",
               show_network = FALSE, 
               network_color = "lightgrey", 
               point_size = 0.5,
               cow_n_col = 2
)