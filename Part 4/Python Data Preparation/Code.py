# pip install numpy h5py anndata scanpy pandas
import numpy, h5py, anndata, scanpy, pandas
print("NumPy:", numpy.__version__)
print("h5py:", h5py.__version__)
print("anndata:", anndata.__version__)
print("scanpy:", scanpy.__version__)
print("pandas:", pandas.__version__)


import anndata as ad
import pandas as pd

adata = ad.read_h5ad("C:/Esmaeil/irAEsProject/Backup/Part 4/0_Data/GSE253720_Biopsy_RNA_Final.h5ad")

# Extract raw counts matrix
expr = pd.DataFrame(adata.layers["counts"].toarray(), index=adata.obs_names, columns=adata.var_names)
# Metadata
obs = adata.obs
var = adata.var

# Save to CSVs (or optionally RDS later)
expr.to_csv("C:/Esmaeil/irAEsProject/Backup/Part 4/0_Data/Python Data Preparation/expression.csv")
obs.to_csv("C:/Esmaeil/irAEsProject/Backup/Part 4/0_Data/Python Data Preparation/metadata_obs.csv")
var.to_csv("C:/Esmaeil/irAEsProject/Backup/Part 4/0_Data/Python Data Preparation/metadata_var.csv")


print("✅ File loaded successfully!")
print("Shape:", adata.shape)
print("Observations (cells):", len(adata.obs))
print("Variables (genes):", len(adata.var))
print("Metadata columns:", list(adata.obs.columns[:10]))