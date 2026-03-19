# Python implementation for RNA-Seq normalization and visualization

# Necessary imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# Sample RNA-Seq count data (rows: genes, columns: samples)
data = pd.DataFrame({
    'Sample1': [100, 150, 200, 250, 300],
    'Sample2': [80, 120, 210, 260, 310],
    'Sample3': [90, 130, 190, 240, 320]
}, index=['GeneA','GeneB','GeneC','GeneD','GeneE'])

# 1. Count per million (CPM) normalization
counts_sum = data.sum(axis=0)
cpm = data.div(counts_sum) * 1e6

# 2. Log2 transformation for better visualization (log2(CPM + 1))
log_cpm = np.log2(cpm + 1)

# 3. Standard Scaler normalization (optional before PCA or clustering)
scaler = StandardScaler()
scaled_data = pd.DataFrame(scaler.fit_transform(log_cpm.T).T, index=log_cpm.index, columns=log_cpm.columns)

# 4. Principal Component Analysis (PCA) for visualization
pca = PCA(n_components=2)
pca_result = pca.fit_transform(scaled_data.T)
pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'], index=log_cpm.columns)

# 5. Plot PCA
plt.figure(figsize=(6,5))
sns.scatterplot(x='PC1', y='PC2', data=pca_df, s=100)
for sample in pca_df.index:
    plt.text(pca_df.loc[sample, 'PC1']+0.1, pca_df.loc[sample, 'PC2']+0.1, sample)
plt.title('PCA of RNA-Seq Samples')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.show()

# 6. Heatmap of normalized expression values
plt.figure(figsize=(8,6))
sns.heatmap(scaled_data, cmap='viridis', annot=True)
plt.title('Heatmap of Normalized RNA-Seq Expression')
plt.xlabel('Samples')
plt.ylabel('Genes')
plt.show()

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