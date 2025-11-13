# pip install numpy h5py anndata scanpy pandas
import numpy, h5py, anndata, scanpy, pandas
print("NumPy:", numpy.__version__)
print("h5py:", h5py.__version__)
print("anndata:", anndata.__version__)
print("scanpy:", scanpy.__version__)
print("pandas:", pandas.__version__)



# Load the .h5ad file
adata = ad.read_h5ad("C:/Esmaeil/irAEsProject/Backup/Part 4/0_Data/GSE253720_Biopsy_RNA_Final.h5ad")



# Extract raw counts matrix
expr = pd.DataFrame(adata.layers["counts"].toarray(), index=adata.obs_names, columns=adata.var_names)
# Check if all values are integers to confirm that this is the raw counts matrix, not normalized or log-transformed data
has_fraction = ((expr.values % 1) != 0).any()
print("Are there any real decimal numbers in expr (e.g., 2.6, 7.9)?", has_fraction)



# For additional confirmation, we extract the default matrix (adata.X) to show that it contains decimals (non-integer values), indicating that it is normalized/log-transformed and not the raw counts matrix
expr2 = adata.to_df()
has_fraction = ((expr2.values % 1) != 0).any()
print("Are there any real decimal numbers in expr2 (e.g., 2.6, 7.9)?", has_fraction)



# So, based on this comparison, expr contains the raw counts that we want, while expr2 contains the normalized data. Therefore, we should save expr, the raw counts matrix.
expr.to_csv("C:/Esmaeil/irAEsProject/Backup/Part 4/0_Data/Python Data Preparation/expression.csv")



# And save Metadata
obs = adata.obs
obs.to_csv("C:/Esmaeil/irAEsProject/Backup/Part 4/0_Data/Python Data Preparation/metadata_obs.csv")

var = adata.var
var.to_csv("C:/Esmaeil/irAEsProject/Backup/Part 4/0_Data/Python Data Preparation/metadata_var.csv")



# print summary
print("✅ File loaded successfully!")
print("Shape:", adata.shape)
print("Observations (cells):", len(adata.obs))
print("Variables (genes):", len(adata.var))
print("Metadata columns:", list(adata.obs.columns[:10]))