import pandas as pd

qname = pd.read_csv(snakemake.input["qname"], sep="\t", names=["QNAME"])
cb = pd.read_csv(snakemake.input["cb"], sep="\t", names=["CB"])
ub = pd.read_csv(snakemake.input["ub"], sep="\t", names=["UB"])
reads_df = pd.concat([qname, cb, ub], axis=1).drop_duplicates()
reads_df.to_csv(snakemake.output[0], sep="\t", index=False)

qname = pd.read_csv(snakemake.input["path_qname"], sep="\t", names=["QNAME"])
yp = pd.read_csv(snakemake.input["yp"], sep="\t", names=["YP"])
yp_df = pd.concat([qname, yp], axis=1).drop_duplicates()
yp_df.to_csv(snakemake.output[1], sep="\t", index=False)

df = yp_df.merge(reads_df, on="QNAME")
df["CB"] = df.CB.str.strip("CB:Z:") # remove leading CB:Z: so we can merge
df["UB"] = df.UB.str.strip("UB:Z:") # remove leading UB:Z: so we can merge
df.to_csv(snakemake.output[2], sep="\t", index=False)
