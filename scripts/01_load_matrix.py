import pandas as pd
import scipy.io
import anndata as ad

RAW = "data/raw/GSE129785_scATAC-TME-TCells"

# cell metadata from authors (tab-separated with header)
obs = pd.read_csv(f"{RAW}.cell_barcodes.txt.gz", sep="\t")
obs.index = obs["Group_Barcode"].values

# peak names, split into genomic coordinates
var = pd.read_csv(f"{RAW}.peaks.txt.gz", sep="\t")
var.index = var["Feature"].values
var[["chrom", "start", "end"]] = var["Feature"].str.split("_", expand=True)
var["start"] = var["start"].astype(int)
var["end"] = var["end"].astype(int)

# load sparse count matrix (uncompressed despite GEO's .gz label)
X = scipy.io.mmread(f"{RAW}.mtx").tocsr()

# orient as cells x peaks
if X.shape == (len(var), len(obs)):
    X = X.T.tocsr()
assert X.shape == (len(obs), len(var)), f"Unexpected shape {X.shape}"

# assemble and add sample annotations
adata = ad.AnnData(X=X.astype("float32"), obs=obs, var=var)
adata.obs["patient"] = adata.obs["Group"].str.split("_").str[0]
adata.obs["timepoint"] = adata.obs["Group"].str.split("_").str[-1]

# save and summarize
adata.write_h5ad("data/processed/tcells_raw.h5ad", compression="gzip")
print(adata)
print(adata.obs["timepoint"].value_counts())
