import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['mathtext.rm'] = 'Arial'
plt.rcParams['font.family'] = ['Arial']
SMALL_SIZE = 5
MEDIUM_SIZE = 6
BIGGER_SIZE = 7
plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE, frameon=False)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
sns.set(rc={'font.family': 'Arial'}, font_scale=1.0)
sns.set_style("whitegrid", {'axes.grid': False})

if __name__ == '__main__':
    df = pd.read_csv("../../output/02_tables/02_data_merged/EquAllRS_08_merged_domains_kindom.csv")
    TOTAL_NUMBER_OF_FOALS_plus1 = 32
    # Read CSV
    foal_nonMT_IS_meta = pd.read_csv("../../output/04_fragmentomics/library_prep_comparison_EquCab3_nonMT_IS_meta.csv")
    # Convert some columns into factors (categorical variables)
    foal_nonMT_IS_meta = foal_nonMT_IS_meta.apply(
        lambda x: pd.factorize(x)[0] if x.name not in ['TLEN', 'Count', 'side', 'EndMotif_tmp', 'sample_id'] else x)
    # Take absolute value of TLEN
    foal_nonMT_IS_meta['TLEN'] = foal_nonMT_IS_meta['TLEN'].abs()
    # Filter rows where TLEN is less than 300 and read is 'R1'
    foal_nonMT_IS_meta = foal_nonMT_IS_meta[
        (foal_nonMT_IS_meta['TLEN'] < 300) & (foal_nonMT_IS_meta['side'].str.startswith('R1'))]
    # Group by sample_id, TLEN, read and calculate sum Count
    summary_df = foal_nonMT_IS_meta.groupby(['sample_id', 'TLEN'])['Count'].sum().reset_index()
    # Group by read and calculate nCount
    summary_df['nCount'] = summary_df['Count'].transform(lambda x: x / x.sum() * 100)


    # Load the data
    foal_nonMT_IS_meta = pd.read_csv("../../output/04_fragmentomics/library_prep_comparison_EquCab3_nonMT_IS_meta.csv")
    foal_nonMT_IS_meta = foal_nonMT_IS_meta.apply(
        lambda x: pd.factorize(x)[0] if x.name not in ['TLEN', 'Count', 'side', 'EndMotif_tmp', 'sample_id'] else x)
    foal_nonMT_IS_meta['TLEN'] = foal_nonMT_IS_meta['TLEN'].abs()
    foal_nonMT_IS_meta = foal_nonMT_IS_meta[
        (foal_nonMT_IS_meta['TLEN'] < 300) & (foal_nonMT_IS_meta['side'].str.startswith('R1'))]
    summary_df = foal_nonMT_IS_meta.groupby(['sample_id', 'TLEN'])['Count'].sum().reset_index()
    summary_df['nCount'] = summary_df['Count'].transform(lambda x: x / x.sum() * 100)
    m = {'2105643Mod': 'Moderate', '2105643Small': 'Small', '2105643SmallBBSS': 'Small + extra bead'}
    summary_df['sample_id'].replace(m, inplace=True)

    zymo = ['Pseudomonas', 'Escherichia', 'Salmonella', 'Limosilactobacillus', 'Enterococcus', 'Staphylococcus',
            'Listeria', 'Bacillus', 'Saccharomyces', 'Cryptococcus'][::-1]
    zymo_species = ['Pseudomonas aeruginosa', 'Escherichia coli', 'Salmonella enterica',
                    'Limosilactobacillus fermentum',
                    'Enterococcus faecalis', 'Staphylococcus aureus', 'Listeria monocytogenes', 'Bacillus spizizenii',
                    'Saccharomyces cerevisiae', 'Cryptococcus neoformans']
    per_sample_species = []
    per_sample = []
    for PC_index in [1, 2, 3]:
        PC1 = pd.read_csv(
            f"../../output/01_pipeline/foal_cohort_EquCabAll/results/bracken_output/after_host_mapping/PC{PC_index}_RS_conf0.8.output",
            sep='\t', index_col=0)
        PC1['name'] = PC1.index
        PC1['Genus'] = PC1['name'].apply(lambda x: x.split(" ")[0])
        genus = PC1.groupby('Genus').sum(numeric_only=True)['fraction_total_reads']
        genus_high_PC1 = genus[genus > 0.00]
        per_sample.append(genus_high_PC1[zymo].values)
        per_sample_species.append(PC1.loc[zymo_species]['fraction_total_reads'].values)

    output = pd.DataFrame(np.array(per_sample), columns=zymo, index=['PC1', 'PC2', 'PC3'])
    outputT = output.T
    outputT.insert(loc=0, column='Expected', value=[0.02] * 2 + [0.12] * 8)
    output2 = outputT.T
    output2.T[::-1].to_csv('../../output/03_microbial/source_data/Fig1F.csv')

    output_species = pd.DataFrame(np.array(per_sample_species), columns=zymo_species, index=['PC1', 'PC2', 'PC3'])
    output_speciesT = output_species.T
    output_speciesT.insert(loc=0, column='Expected', value=[0.02] * 2 + [0.12] * 8)
    output_speciesT2 = output_speciesT.T
    output_speciesT2.T.to_csv('../../output/03_microbial/source_data/SupplTables_PC_species.csv')

    # Create a single figure with GridSpec for different width subplots
    fig = plt.figure(figsize=(15, 5))
    gs = GridSpec(1, 3, width_ratios=[1, 0.6, 0.6])  # Adjust width_ratios to make the middle plot narrower

    # Plot 1
    ax1 = fig.add_subplot(gs[0])
    sns.lineplot(ax=ax1, x='TLEN', y='nCount', hue='sample_id',
                 data=summary_df[summary_df['sample_id'] != '2105643ModBBSS'], palette='Greys')
    ax1.legend(title='')  # Remove the column name from the legend title
    ax1.set_ylabel('Fraction')
    ax1.set_xlabel('cfDNA length (bp)')
    ax1.spines[['right', 'top']].set_visible(False)

    # Plot 2
    ax2 = fig.add_subplot(gs[1])
    spikein_fraction = df[['R1_02_50mer', 'R1_02_100mer', 'R1_02_150mer', 'R1_01_raw_fastq']].apply(
        lambda x: x / x['R1_01_raw_fastq'], axis=1).drop(columns=['R1_01_raw_fastq'])
    spikein_fraction_enrichment = spikein_fraction.div(spikein_fraction.sum(axis=1), axis=0) / (1 / 3) * 100
    sns.boxplot(ax=ax2, data=spikein_fraction_enrichment, color='white', fliersize=0,
                boxprops=dict(facecolor='white', edgecolor='black'),
                medianprops=dict(color='black'),
                whiskerprops=dict(color='black'),
                capprops=dict(color='black'),
                flierprops=dict(markerfacecolor='black', markeredgecolor='black'))
    ax2.spines[['right', 'top']].set_visible(False)
    ax2.set_ylabel('Read length enrichment (%)')
    ax2.set_xticklabels(['50 bp\nspike-in', '100 bp\nspike-in', '150 bp\nspike-in'])

    # Plot 3
    ax3 = fig.add_subplot(gs[2])
    output2.plot(kind='bar', stacked=True, colormap='twilight_shifted', linewidth=0.5, ax=ax3)
    ax3.set_ylabel('Fraction')
    handles, labels = ax3.get_legend_handles_labels()
    ax3.spines[['right', 'top']].set_visible(False)
    ax3.legend(reversed(handles), reversed(labels), bbox_to_anchor=(1, 1))

    plt.tight_layout()
    output_path = '../../output-figures/Fig1_combined_size_modified.pdf'
    plt.savefig(output_path)
    plt.savefig(output_path + '.png')
    plt.show()

    ## Extra figure not showed
    output_path = '../../output/03_microbial/figures/Fig1E_extended.pdf'

    spikein_fraction = df[['R1_02_50mer', 'R1_02_100mer', 'R1_02_150mer', 'R1_01_raw_fastq']].apply(
        lambda x: x / x['R1_01_raw_fastq'], axis=1).drop(columns=['R1_01_raw_fastq'])
    # Enrichment compared to 1/3
    spikein_fraction_enrichment = spikein_fraction.div(spikein_fraction.sum(axis=1), axis=0) / (1 / 3) * 100
    ## 4 Negative Controls does not have spike-in and count as NaN in the plotting (not included in average)
    # sns.boxplot( data=spikein_fraction_enrichment, color = 'black', fill=False, fliersize=0,)
    sns.boxplot(data=spikein_fraction_enrichment, color='white', fliersize=0, )
    spikein_fraction_enrichment['prep_n'] = df['prep_n'].copy()

    # For dots per batch
    for n, color in zip([1, 2, 3], ['orange', 'green', 'lightblue', ]):
        sub_df1 = spikein_fraction_enrichment[spikein_fraction_enrichment['prep_n'] == f'Srsly prep {n}']
        k1234 = sub_df1.drop(columns='prep_n')
        ax = sns.stripplot(data=k1234, color=color, linewidth=0, label=f'Srsly prep {n}')
    # Black
    sns.stripplot(data=spikein_fraction_enrichment[-3:], color='white', edgecolor='black', linewidth=1, label='PCs')

    ax.spines[['right', 'top']].set_visible(False)
    ax.set_ylabel('Read length enrichment (%)')
    ax.set_xticklabels(['50 bp\nspike-in', '100 bp\nspike-in', '150 bp\nspike-in', ])
    ax.legend(bbox_to_anchor=(1, 1))
    plt.tight_layout()
    # plt.savefig(output_path)
    # plt.savefig(output_path + '.png')
    # plt.show()
