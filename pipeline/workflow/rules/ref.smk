# unit name = unique name?
units = (
    pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "ref_genome": str})
    .set_index(["sample_name"], drop=False)
    .sort_index()
)

def get_reference(sample_name):
	ref_genome = units.loc[sample_name]['ref_genome'].values[0]
	return ref_genome
