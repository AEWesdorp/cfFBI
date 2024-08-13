import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
from scipy.stats import f
from scipy.stats import mannwhitneyu
import statsmodels.api as sm
import os
import statannotations
from statannotations.Annotator import Annotator
import matplotlib.colors as mcolors
import skbio
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa
import matplotlib.pyplot as plt
from scipy.spatial.distance import braycurtis, pdist, squareform
from scipy.stats import entropy
import numpy as np
import pandas as pd
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt


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
sns.set_style("whitegrid")

import pandas as pd
import numpy as np


def rarefy(df, target_depth):
    rarefied_df = pd.DataFrame(index=df.index, columns=df.columns)

    for sample in df.index:
        counts = df.loc[sample]
        if counts.sum() >= target_depth:
            # Corrected rarefaction process
            rarefied_counts = np.random.choice(
                counts.index, size=target_depth, replace=True, p=counts / counts.sum()
            )
            rarefied_counts = pd.Series(rarefied_counts).value_counts()
            rarefied_df.loc[sample] = rarefied_counts.reindex(df.columns, fill_value=0)
        else:
            rarefied_df.loc[sample] = counts

    return rarefied_df


# Sample data for illustration


# Function to calculate alpha diversity
def calculate_alpha_diversity(group, abundance_label):
    # All: 'Exact abundance'
    # Rarefy: 'Exact Rarefy abundance'
    abundance_per_species = group[abundance_label]
    species_counts = group['Species'].unique()
    # print(species_counts)
    # print(species_counts == abundance_per_species)

    # Species richness (S)
    species_richness = len(abundance_per_species)

    # Shannon diversity index (H)

    shannon_index = entropy(abundance_per_species, base=np.e)

    # Simpson's diversity index (D)
    # total_count = group[abundance_label].sum()
    # simpson_index = 1 - sum((count / total_count) ** 2 for count in abundance_per_species)

    return pd.Series({
        'Species Richness': species_richness,
        'Shannon Index': shannon_index,
        # 'Simpson Index': simpson_index
    })



def import_data():
    df = pd.read_csv("../../output/02_tables/02_data_merged/EquAllRS_08_merged_domains_kindom.csv")
    TOTAL_NUMBER_OF_FOALS_plus1 = 32
    TOTAL_NUMBER_OF_FOAL_SAMPLES_plus1 = 33

    return df

def preprocessing(df):
    df[['Bacterial fraction (with contam)', 'Fungal fraction (with contam)', 'Viral fraction (with contam)',
        'to_drop']] = df[['Bacteria', 'Fungi', 'Viruses', 'total_filtered_effective_reads']].apply(
        lambda row: row / row['total_filtered_effective_reads'], axis=1)
    df['Foal MT cfDNA fraction'] = df['MT_reads'].astype(int) / df['Mapped to host genomes']
    df['Foal MT cfDNA conc. in plasma (ng/ml)'] = df['Foal MT cfDNA fraction'] * df['cfDNA yield (ng/ml)']
    df.rename(columns={'cfDNA yield (ng/ml)': 'Total cfDNA conc. in plasma (ng/ml)'}, inplace=True)
    return df

def exclude_contam(df, basic):
    THRESHOLD = 10
    labels = basic[['FID', 'nSIRS_class', 'nSIRS_score']].copy()
    labels.index = labels['FID']
    taxid_dict = pd.read_pickle("../../output/02_tables/03_intermediate/taxid_dict.pickle.gz")
    bac = pd.read_csv('../../output/02_tables/02_data_merged/EquAllRS_08_bacteria_species_for_decontam_taxid.csv')
    bac = bac.rename(columns=taxid_dict, )
    total_count = bac.sum(axis=1)
    bac_normalized = bac.apply(lambda x: x / total_count, axis=0)
    bac_index = pd.read_csv('../../output/02_tables/02_data_merged/EquAllRS_08_FID_for_decontam.csv')
    bac_normalized['FID'] = bac_index['FID']
    bac_normalized.index = bac_index['FID']
    bac['Total QC reads'] = total_count
    bac_index = pd.read_csv('../../output/02_tables/02_data_merged/EquAllRS_08_FID_for_decontam.csv')
    bac['FID'] = bac_index['FID']
    bac.index = bac_index['FID']
    bac['nSIRS_class'] = labels['nSIRS_class']
    bac['nSIRS_score'] = labels['nSIRS_score']

    # Read contaminants
    contaminants_table = pd.read_csv('../../output/02_tables/03_intermediate/decontam_species_joint.csv')
    contaminants = contaminants_table[contaminants_table['Contaminants'] == True]['name'].values
    contaminants_iso = contaminants_table[contaminants_table['Isolation related contaminant (p<0.25)'] == True][
        'name'].values
    contaminants_prep = contaminants_table[contaminants_table['Library prep related contaminant (p<0.25)'] == True][
        'name'].values

    contaminants = list(contaminants_iso) + list(contaminants_prep)
    contaminants_bacterial = set(list(contaminants_iso) + list(contaminants_prep)).intersection(bac.columns)

    bac_is_contam = bac[list(contaminants_bacterial)]
    is_contam = pd.Series(bac_is_contam.sum(axis=1), name='Bacterial contaminants read counts')
    is_contam = is_contam.reset_index()
    df = df.merge(is_contam, on='FID', how='left')

    bac_not_contam = bac.drop(columns=list(contaminants_bacterial))
    bac_not_contam_bac = bac_not_contam.iloc[:, :-5]
    not_contam = pd.Series(bac_not_contam_bac.sum(axis=1), name='Bacterial reads')
    not_contam = not_contam.reset_index()
    df = df.merge(not_contam, on='FID', how='left')

    df['Bacterial fraction'] = df['Bacterial reads'] / df['total_filtered_effective_reads']

    return df, bac


def get_bray_curtis_matrix_from_long_table(df, variable='Exact Rarefied abundance'):
    # Pivot the table to have species as columns and IDs as rows
    pivot_df = df.pivot(index='FID', columns='Species', values=variable).fillna(0)

    # Calculate Bray-Curtis dissimilarity matrix
    bray_curtis_matrix = pd.DataFrame(squareform(pdist(pivot_df, metric='braycurtis')),
                                      index=pivot_df.index,
                                      columns=pivot_df.index)

    # Print the results
    # print(bray_curtis_matrix)
    return bray_curtis_matrix, pivot_df


def get_jaccard_index_from_long_table(df, variable='Exact Rarefied abundance'):
    # Pivot the table to have species as columns and IDs as rows
    pivot_df = df.pivot(index='FID', columns='Species', values=variable).fillna(0)

    # Convert abundance data to binary presence/absence data
    binary_df = pivot_df.applymap(lambda x: 1 if x > 0 else 0)

    # Calculate Jaccard index dissimilarity matrix
    jaccard_index_matrix = pd.DataFrame(squareform(pdist(binary_df, metric='jaccard')),
                                        index=binary_df.index,
                                        columns=binary_df.index)

    return jaccard_index_matrix, binary_df


def plot_bray_curtis_clustermap(bray_curtis_matrix_log, row_colors, out_path):
    g = sns.clustermap(bray_curtis_matrix_log,
                       figsize=(14, 14),
                       row_colors=row_colors,
                       col_colors=row_colors,
                       yticklabels=False,
                       dendrogram_ratio=(0.05, 0),
                       cbar_pos=(0.00, 0.8, 0.05, 0.18),
                       cmap='Greys_r'
                       )
    # Get the reordered row indices
    row_order = g.dendrogram_row.reordered_ind
    col_order = g.dendrogram_col.reordered_ind

    # Map the nSIRS_class to colors according to the clustered order
    ordered_row_colors = df_rarefied_no_sS['nSIRS_class'].iloc[row_order].map(color_map).values

    # Apply the correct row colors
    g.ax_heatmap.set_yticks(np.arange(len(row_order)) + 0.5)
    g.ax_heatmap.set_yticklabels(bray_curtis_matrix_log.index[row_order])

    g.ax_heatmap.set_xticks(np.arange(len(col_order)) + 0.5)
    g.ax_heatmap.set_xticklabels(bray_curtis_matrix_log.columns[col_order])

    # Rotate and adjust the column labels
    plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)

    # Add the row colors
    for label in categories:
        g.ax_row_dendrogram.bar(0, 0, color=color_map[label], label=label, linewidth=0)
    g.ax_row_dendrogram.legend(title='nSIRS_class', loc='center', ncol=len(categories), bbox_to_anchor=(0.5, -0.1))
#    plt.tight_layout()
    plt.savefig(out_path)
    plt.savefig(out_path + '.png',dpi=300)

    plt.show()

if __name__ == '__main__':
    df = import_data()
    df = preprocessing(df)
    basic = pd.read_csv('../../output/02_tables/02_data_merged/all_basic_stats.csv')
    df, bac = exclude_contam(df, basic)

    custom_order = ['H', 'nS-', 'S+', ]
    df['nSIRS_class'] = pd.Categorical(df['nSIRS_class'], categories=custom_order, ordered=True)
    custom_order2 = ['Y', 'N']

    df['Lived'] = pd.Categorical(df['Lived'], categories=custom_order2, ordered=True, )
    replacement = {'Y': 'S+ Lived', 'N': 'S+ Died'}
    df['Lived'].replace(replacement, inplace=True)
    # df['S+ Survival'] = df['Lived'].apply(get_survival_S)
    septic_group = df[df['nSIRS_class'] == 'S+'].copy()

    long_table = pd.read_csv('../../output/02_tables/04_source_data/contaminantFree_bacteria_long_inc_batched.csv')
    df_original = bac.iloc[:, :bac.columns.get_loc('all_filtered_reads_exc_listed_here') + 1]
    # Get minimum filtered read count in all samples
    min_read_count = int(df_original.sum(axis=1).min())
    max_read_count = int(df_original.sum(axis=1).max())
    # Subsample to same amount of total reads as the minimal read count
    df_rarefied = rarefy(df_original, min_read_count)
    print(df_rarefied.isnull().sum().sum())
    df_rarefied = df_rarefied.astype(float)
    df_rarefied[['FID', 'nSIRS_class', 'nSIRS_score']] = bac[['FID', 'nSIRS_class', 'nSIRS_score']].copy()
    df_rarefied_no_sS = df_rarefied[df_rarefied['nSIRS_class'] != 'sS-']
    df_rarefied_melted = df_rarefied.melt(id_vars=['FID', 'nSIRS_class','nSIRS_score'], var_name='Species', value_name='Exact Rarefied abundance')
    df_rarefied_long_table = long_table.merge( right = df_rarefied_melted, how  = 'left')

    # Take Non-contaminants, with exact count more than 10
    n = 0
    sub1_rarefied = df_rarefied_long_table[(df_rarefied_long_table['Contaminant'] == False) & \
                                           (df_rarefied_long_table['Exact_count_morethaninc_10'] == True) & \
                                           (df_rarefied_long_table['Species'] != 'all_filtered_reads_exc_listed_here') & \
#                                           (df_rarefied_long_table['nSIRS_class'] != 'sS-') &\
                                           (df_rarefied_long_table['Exact Rarefied abundance'] >= n)
                                           ]

    alpha_diversity_results = sub1_rarefied.groupby(['FID', 'nSIRS_class']).apply(
        lambda x: calculate_alpha_diversity(x, 'Exact Rarefied abundance')).reset_index()
    alpha_diversity_results_with_sS =alpha_diversity_results.copy()

    alpha_diversity_results['nSIRS_class'] = pd.Categorical(alpha_diversity_results['nSIRS_class'],
                                                            categories=['H', 'nS-', 'S+',])  # 'sS-'
    alpha_diversity_results_with_sS['nSIRS_class'] = pd.Categorical(alpha_diversity_results['nSIRS_class'],
                                                            categories=['H', 'nS-', 'S+','sS-'])

    alpha_melted = alpha_diversity_results.melt(id_vars=['FID', 'nSIRS_class'],
                                                value_vars=['Species Richness', 'Shannon Index', ],
                                                var_name='measurement',
                                                value_name='value')  # 'Simpson Index'

    alpha_melted['nSIRS_class'] = pd.Categorical(alpha_melted['nSIRS_class'], categories=['H', 'nS-', 'S+', ]) # 'sS-'
    # Plot
    custom_palette = ['#6066B6', '#7CA2C2', '#A84750', ] #'grey']

    fig, axes = plt.subplots(2, 2, figsize=(5, 5))

    sns.boxplot(x='nSIRS_class', y='Bacterial fraction', data=df[df['FID'] != 'F06'], palette=custom_palette,
                showfliers=False, ax=axes[0][0])
    sns.stripplot(x='nSIRS_class', y='Bacterial fraction', data=df[df['FID'] != 'F06'], palette=custom_palette,
                  ax=axes[0][0], linewidth=1)
    axes[0][0].set_ylim([0, 0.00050])
    ## Add outlier F06 manually at the level it is (0.05)

    sns.boxplot(x='Lived', y='Bacterial fraction', data=septic_group[septic_group['FID'] != 'F06'],
                palette=['#e35f6c', '#85363e'], showfliers=False, ax=axes[1][0])
    sns.stripplot(x='Lived', y='Bacterial fraction', data=septic_group[septic_group['FID'] != 'F06'],
                  palette=['#e35f6c', '#85363e'], ax=axes[1][0], linewidth=1)
    axes[1][0].set_ylim([0, 0.00058])
    axes[1][0].set_ylim([0, 0.00058])

    ## Add outlier F06 manually at the level it is (0.05)

    # plt.suptitle(f'At least {n} observation, Rarefied')
    sns.boxplot(x='nSIRS_class', y='Species Richness', data=alpha_diversity_results, ax=axes[0][1],
                palette=custom_palette, showfliers=False)
    sns.stripplot(x='nSIRS_class',
                  y='Species Richness',
                  data=alpha_diversity_results,
                  ax=axes[0][1],
                  palette=custom_palette,
                  linewidth=1)
    sns.boxplot(x='nSIRS_class', y='Shannon Index', data=alpha_diversity_results, ax=axes[1][1], palette=custom_palette,
                showfliers=False)
    sns.stripplot(x='nSIRS_class', y='Shannon Index', data=alpha_diversity_results, ax=axes[1][1],
                  palette=custom_palette, linewidth=1)

    axes[0][0].spines[['right', 'top']].set_visible(False)
    axes[0][1].spines[['right', 'top']].set_visible(False)
    axes[1][0].spines[['right', 'top']].set_visible(False)
    axes[1][1].spines[['right', 'top']].set_visible(False)
    axes[0][0].set_xlabel('')
    axes[0][1].set_xlabel('')
    axes[1][0].set_xlabel('')
    axes[1][1].set_xlabel('')

    ## Add annotation
    for df1, feature_y, category_x, ax, pairs in zip(
            [df, septic_group, alpha_diversity_results, alpha_diversity_results],
            ['Bacterial fraction', 'Bacterial fraction', 'Species Richness', 'Shannon Index'],
            ['nSIRS_class', 'Lived', 'nSIRS_class', 'nSIRS_class'],
            [axes[0][0], axes[1][0], axes[0][1], axes[1][1]],
            [[('H', 'S+'), ('nS-', 'S+')], [('S+ Lived', 'S+ Died')], [('H', 'S+'), ('nS-', 'S+')],
             [('H', 'S+'), ('nS-', 'S+')]]
            ):
        annotator = Annotator(ax, pairs, data=df1, x=category_x, y=feature_y,
                              order=df1[category_x].cat.categories)  # hide_non_significant=True
        annotator.configure(test='Mann-Whitney', text_format='star', loc='outside', verbose=0)
        annotator.apply_and_annotate()

    plt.tight_layout()
    plt.savefig(f'../../output-figures/Fig3_alpha_diversity_atleast{n}_rarefied_no_sS_2by2_sig_test.pdf')
#    plt.savefig(f'../../output-figures/Fig3_alpha_diversity_atleast{n}_rarefied_no_sS_2by2_sig_test.png')
    plt.show()


    # Beta diversity:
    sub1_excl_sS_rarefied = sub1_rarefied[sub1_rarefied['nSIRS_class'] != 'sS-']
    sub1_excl_sS_rarefied['Exact Rarefied abundance (log)'] = sub1_excl_sS_rarefied['Exact Rarefied abundance'].apply(
        lambda x: np.log(x) if x != 0 else np.log(x + 0.0000001))
    bray_curtis_matrix_log, pivot_df3 = get_bray_curtis_matrix_from_long_table(sub1_excl_sS_rarefied,
                                                                               'Exact Rarefied abundance (log)')
    categories = ['H', 'nS-', 'S+', ] #'sS-'
    palette = custom_palette
    color_map = dict(zip(categories, palette))
    # Map the nSIRS_class to colors
    row_colors = df_rarefied_no_sS['nSIRS_class'].map(color_map)
    out_beta = '../../output-figures/Fig3_clustering_bray_curtis_matrix_no_sS_log_c0.0000001_tight_layout_test.pdf'
    plot_bray_curtis_clustermap(bray_curtis_matrix_log, row_colors, out_beta)




    # Suppl Figs Age as a co-variant
    df2 = df.copy()
    custom_palette_with_sS = ['#6066B6', '#7CA2C2', '#A84750','grey' ]
    custom_order = ['H', 'nS-', 'S+', 'sS-']
    df2['nSIRS_class'] = pd.Categorical(df2['nSIRS_class'], categories=custom_order, ordered=True)
    df2['nSIRS_class'] = df2['nSIRS_class'].fillna('sS-')
    df2['Foal age at presentation (hours)'] = basic['Foal age at presentation (hours)'].copy()
    df3 = df2[df2['FID'] != 'F06']

    # Suppl Fig 6A: Age and bacterial fraction scatterplot
    fig, ax = plt.subplots(figsize=(5, 4))
    sns.scatterplot(x='Foal age at presentation (hours)',
                    y='Bacterial fraction',
                    hue='nSIRS_class',
                    data=df3,
                    palette=custom_palette_with_sS,
                    ax=ax)
    ax.set_ylim(0, 0.0005)
    ax.set_xlim(0, 160)
    plt.tight_layout()
    plt.savefig('../../output-figures/SupplFig_foal_age_bacterial_load.pdf')
    plt.show()
    # Suppl Fig 6B
    fig, ax = plt.subplots(figsize=(6, 4))
    print('Outlier is not plotted: \n', df2[(df2['FID'] == 'F06')][['Foal age at presentation (hours)', 'Bacterial fraction']])
    df3['Age'] = pd.cut(df3['Foal age at presentation (hours)'].copy(),
                        bins=[-float('inf'), 24, float('inf')],
                        labels=['<1 day old', '>1 day old'])
    sns.boxplot(x='Age', y='Bacterial fraction', data=df3, palette=custom_palette_with_sS, showfliers=False, hue='nSIRS_class',
                ax=ax)
    sns.stripplot(x='Age', y='Bacterial fraction', data=df3, palette=custom_palette_with_sS, linewidth=1,
                  ax=ax, hue='nSIRS_class', dodge=True)

    ax.spines[['right', 'top']].set_visible(False)
    ax.set_xlabel('')
    ax.legend(bbox_to_anchor=(1.1, 1))
    plt.tight_layout()
    plt.savefig('../../output-figures/SupplFig6B_age_group.pdf')
    plt.show()
    # Suppl Figure Species richness 6C
    alpha_diversity_results_with_sS['Foal age at presentation (hours)'] = basic['Foal age at presentation (hours)'][:-7].copy()
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.scatterplot('Foal age at presentation (hours)', 'Species Richness',
                    data=alpha_diversity_results_with_sS, hue='nSIRS_class',
                    palette=custom_palette_with_sS, ax=ax)
    ax.set_xlim(0, 150)
    ax.legend(bbox_to_anchor=(1.3, 1))
    plt.tight_layout()
    plt.savefig('../../output-figures/SupplFig6C_Age_Species_richness_test.pdf')
    plt.show()

    # Suppl Figure 6D Shannon index
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.scatterplot('Foal age at presentation (hours)', 'Shannon Index',
                    data=alpha_diversity_results_with_sS, hue='nSIRS_class',
                    palette=custom_palette_with_sS, ax=ax)
    ax.set_ylim(0, 3)
    ax.set_xlim(0, 150)
    ax.legend(bbox_to_anchor=(1.3, 1))
    plt.tight_layout()
    plt.savefig('../../output-figures/SupplFig6D_Age_shannon_index_test.pdf')
    plt.show()
    # Export dataframes
    df[['FID','Bacterial fraction','nSIRS_class','Lived']].to_csv('../../output/02_tables/04_source_data/Fig3_bacterial_load_test.csv')
    alpha_diversity_results_with_sS[['FID', 'nSIRS_class','Shannon Index','Species Richness']].to_csv('../../output/02_tables/04_source_data/Fig3_diversity_measurements_test.csv', index = False)
