import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
import statsmodels.api as sm
from scipy.stats import f
import os
import statannotations
from statannotations.Annotator import Annotator
import matplotlib.colors as mcolors

## Print versions
print("numpy", np.__version__)
print("pandas", pd.__version__)
print("seaborn", sns.__version__)
print("scipy", scipy.__version__)
print("statsmodels", sm.__version__)
print("statannotations", statannotations.__version__)

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
    stats = pd.read_csv("../../output/02_tables/02_data_merged/all_basic_stats.csv")
    df = pd.read_csv("../../output/02_tables/02_data_merged/EquAllRS_08_merged_domains_kindom.csv")
#    df_genus = pd.read_csv("../../output/02_tables/02_data_merged/EquAllRS_08_merged_genera.csv")
    df_raw_seq = df.loc[:,
                 ['FID', 'R1_01_raw_fastq', 'R1_02_spike_in_unmapped', 'R1_03_uniq_fastq', 'R1_04_fastp_fastq',
                  'R1_05_adapt_remov_fastq', ]]
    df_raw_seq.set_index('FID', drop=True, inplace=True)

    df_raw_seq_diff = df_raw_seq.diff(axis=1).apply(lambda x: -x)
    df_raw_seq_diff['Duplicated reads'] = df_raw_seq_diff['R1_03_uniq_fastq']
    df_raw_seq_diff['Low quality reads'] = df_raw_seq_diff[['R1_04_fastp_fastq', 'R1_05_adapt_remov_fastq']].sum(axis=1)
    df_raw_seq_diff['Spike in reads'] = df_raw_seq_diff['R1_02_spike_in_unmapped']
    df_raw_seq_diff = df_raw_seq_diff.drop(columns=df_raw_seq_diff.loc[:, df_raw_seq_diff.columns.str.contains('R1')])
    df_raw_seq_diff['Good quality reads'] = df_raw_seq['R1_01_raw_fastq'] - df_raw_seq_diff.sum(axis=1)
    df_raw_seq_diff = df_raw_seq_diff[df_raw_seq_diff.columns[::-1]].sort_index()
    df_raw_seq_diff.to_csv('../../output/03_microbial/source_data/SupplFig202AB_test.csv')
    outpath = '../../output/03_microbial/figures/SupplFig202AB_test.pdf'

    fig, axes = plt.subplots(1, 2, figsize=(12, 6), width_ratios=[1, 2], sharey=True)
    # ['#B3E2CD','#FDCDAC']
    df_raw_seq_diff.plot(kind='barh', stacked=True, ax=axes[0], width=0.9, colormap='Pastel2', linewidth=0,
                         legend=False)
    # Don't show the series name
    axes[0].set_ylabel('')
    axes[0].set_xlabel('Read count')
    axes[0].set_xlim(0, 70000000)

    # Normalize
    df_raw_seq_diff_normalized = df_raw_seq_diff.div(df_raw_seq_diff.sum(axis=1), axis=0) * 100
    df_raw_seq_diff_normalized.plot(kind='barh', stacked=True, ax=axes[1], width=0.9, colormap='Pastel2', linewidth=0)
    axes[1].set_xlabel('Read ratio (%)')
    plt.gca().invert_yaxis()
    axes[1].legend(bbox_to_anchor=(1, 1))
    plt.tight_layout()
    plt.savefig(outpath)
    plt.savefig(outpath + '.png', dpi=300)

    plt.show()


    df2 = stats[['FID','Mapped to host genomes','unclassified_read_count','classified_read_count']]
    df2.set_index('FID', drop = True, inplace = True)
    df2.to_csv('../../output/03_microbial/source_data/SupplFig202CD_test.csv')

    outpath = '../../output/03_microbial/figures/SupplFig202CD-not_test.pdf'
    colors = ['#93c2ad', '#597569', '#c0fce1']  # Example colors
    # Create a colormap from the list of colors
    cmap_name = 'custom_cmap'
    n_bins = len(colors)  # Number of colors
    cus_cmap = mcolors.LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)
    cus_cmap2 = mcolors.LinearSegmentedColormap.from_list(cmap_name, colors[1:], N=n_bins)

    hatches = ['\\', '\\\\\\', '']  # Define hatching patterns
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), width_ratios=[1, 2], sharey=True)
    df2 = df2.rename(columns={'unclassified_read_count': 'Unclassified', 'classified_read_count': 'Classified'})
    bars1 = df2.plot(kind='barh',
                     stacked=True,
                     ax=axes[0],
                     width=0.9,
                     colormap=cus_cmap,
                     linewidth=1,
                     legend=True, )
    for bar, hatch in zip(bars1.containers, hatches):
        for patch in bar:
            patch.set_hatch(hatch)
            patch.set_edgecolor(colors[0])

    # Don't show the series name
    axes[1].set_ylabel('')
    axes[0].set_ylabel('')

    axes[1].set_xlabel('Read count')
    df4 = df2[['Unclassified', 'Classified']].div(df2[['Unclassified', 'Classified']].sum(axis=1), axis=0) * 100
    axes[0].set_xlim(0, 70000000)

    bars4 = df4.plot(kind='barh',
                     stacked=True,
                     ax=axes[1],
                     width=0.9,
                     colormap=cus_cmap2,
                     legend=False)
    for bar, hatch in zip(bars4.containers, hatches[1:]):
        for patch in bar:
            patch.set_hatch(hatch)
            patch.set_edgecolor(colors[0])
    plt.gca().invert_yaxis()
    axes[1].legend(bbox_to_anchor=(1, 1))
    plt.tight_layout()
    plt.savefig(outpath)
    plt.savefig(outpath + '.png', dpi=300)
    plt.show()