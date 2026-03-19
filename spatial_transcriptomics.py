# Language: Python

# Install necessary packages if not already installed:
# pip install scanpy squidpy matplotlib seaborn

import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt

# ------------------------------
# Load example spatial transcriptomics dataset (Visium)
# ------------------------------
adata = sq.datasets.visium_hne_adata()  # includes spatial coordinates

# ------------------------------
# Preprocessing
# ------------------------------
sc.pp.filter_cells(adata, min_genes=200)       # Filter cells with too few genes
sc.pp.filter_genes(adata, min_cells=3)        # Filter genes detected in very few cells
sc.pp.normalize_total(adata, target_sum=1e4)   # Normalize counts per cell
sc.pp.log1p(adata)                              # Log-transform data
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)  # Select top variable genes
sc.pp.scale(adata, max_value=10)               # Scale data

# ------------------------------
# Dimensionality reduction and clustering
# ------------------------------
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_neighbors=10)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)  # clustering

# ------------------------------
# Visualization: UMAP with clusters
# ------------------------------
sc.pl.umap(adata, color='leiden', title='Spatial Transcriptomics Clusters')

# ------------------------------
# Spatial plotting
# ------------------------------
sc.pl.spatial(adata, img_key="hires", color="leiden", title='Spatial Clusters on Tissue', size=1.5)

# ------------------------------
# Spatial Neighborhood Enrichment Analysis (Squidpy)
# ------------------------------
sq.gr.spatial_neighbors(adata)            # Build spatial neighbors graph
sq.gr.nhood_enrichment(adata, cluster_key="leiden")  # Test enrichment between clusters
sq.pl.nhood_enrichment(adata, cluster_key="leiden", cmap="viridis")

# ------------------------------
# Gene-specific spatial plotting
# ------------------------------
genes_of_interest = ["PLAC8", "APOE", "GFAP"]  # replace with relevant genes in your dataset
for gene in genes_of_interest:
    sc.pl.spatial(adata, color=gene, cmap='magma', size=1.5, title=f"Spatial expression of {gene}")

# ------------------------------
# Example: Moran's I for spatial autocorrelation of gene expression
# ------------------------------
sq.gr.spatial_autocorr(adata, mode='moran', genes=genes_of_interest)
sc.pl.spatial(adata, color=genes_of_interest[0], cmap='coolwarm', size=1.5, title=f"Moran's I {genes_of_interest[0]}")

# ------------------------------
# Summary:
# - This workflow loads Visium spatial transcriptomics data.
# - Normalizes, logs, and scales counts.
# - Performs PCA, clustering, and UMAP visualization.
# - Plots clusters spatially on tissue slides.
# - Performs spatial neighborhood enrichment and autocorrelation.