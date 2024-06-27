import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# Load data
foal_nonMT_IS = pd.read_csv("../../output/04_fragmentomics/foal_cohort_EquCabAll_host_nonMT_motif.csv")
metadata = pd.read_csv("../../output/02_tables/02_data_merged/summary_basic_stats.csv")

# Rename and select relevant columns
metadata = metadata.rename(columns={'patientID': 'sample_id'})[['sample_id', 'FID', 'nSIRS_class']]

# Merge metadata with foal_nonMT_IS
foal_nonMT_IS_meta = pd.merge(metadata, foal_nonMT_IS, on='sample_id')

# Data processing
foal_nonMT_IS_meta = (
    foal_nonMT_IS_meta[~foal_nonMT_IS_meta['sample_id'].str.contains('^[P/N]')]
    .assign(
        EndMotif=lambda df: df['EndMotif'].str[:1],
        read=lambda df: df['side'].str[:2]
    )
    .query("read == 'R1'")
    .assign(
        nSIRS_class=lambda df: pd.Categorical(df['nSIRS_class'], categories=['H', 'nS-', 'sS-', 'S+'])
    )
)

# Group and summarize data
grouped_df = (
    foal_nonMT_IS_meta
    .groupby(['nSIRS_class', 'FID', 'read', 'EndMotif'])
    .agg(sCount=('Count', 'sum'))
    .reset_index()
)

# Calculate nsCount
grouped_df['nsCount'] = grouped_df.groupby(['nSIRS_class', 'FID', 'read'])['sCount'].transform(lambda x: np.log10(x / x.sum() * 4))

# Remove NA and zero values from nsCount
grouped_df = grouped_df[(~grouped_df['nsCount'].isna()) & (grouped_df['nsCount'] != 0)]

# Create boxplot
boxprops = dict(linewidth=2.5)
medianprops = dict(linewidth=2.5, color='black')
plt.boxplot([grouped_df[grouped_df['nSIRS_class'] == cls]['nsCount'] for cls in grouped_df['nSIRS_class'].unique()],
            positions=np.arange(1, len(grouped_df['nSIRS_class'].unique()) + 1),
            labels=grouped_df['nSIRS_class'].unique(),
            patch_artist=True, boxprops=boxprops, medianprops=medianprops)

# Customize plot aesthetics
plt.title('Boxplot of nsCount ("normalized, log10") by nSIRS_class and EndMotif') #this is the end-goal
plt.xlabel('nSIRS_class')
plt.ylabel('nsCount')
plt.xticks(rotation=45)
plt.tight_layout()

# Show plot
plt.show()

### @Liting: these counts are correct!
grouped_df

### @Liting
# There are some errors here, namely that the x-axis should be the EndMotif
# And the color should be the nSIRS_class

#As you suggest, we could add:
# - significant testing
# - split y-axis for each end <-- for now I just plot R1, but we have the same for R2!
# - add dots within the box as well <-- yess I agree 100%
