import pandas as pd
import json


cols = ["Length", "reads", "pct_reads", "cum_reads", "cum_pct_reads", "bases", "pct_bases", "cum_bases", "cum_pct_bases"]
df = pd.read_csv(snakemake.input[0], sep="\t", comment="#", names=cols)
metadata = {}
metadata["max_readlength"] = int(df["Length"].max())
metadata["sjdbOverhang"] = int(df["Length"].max()-1)
with open(snakemake.output[0], "w") as write_file:
    json.dump(metadata, write_file)
