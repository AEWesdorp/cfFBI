import glob
import pandas as pd

#configfile: "config/config.yaml"


# unit name = unique name?
units = (
    pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str}, skipinitialspace=True)
    .set_index(["sample_name"], drop=False)
    .sort_index()
)

def get_fq(wildcards):
    #no trimming, use raw reads
    u = units.loc[wildcards.sample_name]['path_to_R1_R2'] # .str.strip()
    return u
