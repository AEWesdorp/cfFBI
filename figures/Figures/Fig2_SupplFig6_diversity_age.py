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
from scipy.stats import pearsonr
from matplotlib.colors import ListedColormap



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

def get_top_n_relative_abundance(sub1):
    species_df = sub1.pivot(index='FID', columns='Species', values=['Relative abundance'], )
    species_df.columns = species_df.columns.get_level_values(1)
    species_df = species_df.loc[:, species_df.sum(axis=0) > 0]
    top_n = 45
    # Sort the DataFrame by the column in descending order
    sorted_df = sub1.sort_values(by='Relative abundance', ascending=False)
    # Get the top 20 values
    top_values = sorted_df.head(top_n)
    # Extract the indices of the top 20 values
    top_indices = top_values.index
    # Extract top species
    species_top_n = []
    for x in top_indices:
        species_top_n.append(sub1.loc[x, 'Species'])
        species_top_n = list(set(species_top_n))
#        print(len(species_top_n))
        # Non repeat 20 to relative abundance species

    species_top_n.sort()
    return species_top_n

def get_top_species_names(sub1, species_df, species_top_n):
    top20_df = species_df[species_top_n]
    top20_df = top20_df[top20_df.sum().sort_values(ascending=False).index]
    not_top_20_species = []
    for x in species_df.columns:
        if x not in species_top_n:
            not_top_20_species.append(x)
    not_top20_df = species_df[not_top_20_species]
    mask_morethan10 = sub1.pivot(index='FID', columns='Species', values=['Exact_count_morethaninc_10'], )
    mask_morethan10.columns = mask_morethan10.columns.get_level_values(1)
    not_top20_df_mask_morethan10 = mask_morethan10[not_top_20_species]
    more_than_10_sum_per_foal = not_top20_df[not_top20_df_mask_morethan10].fillna(0).sum(axis=1)
    less_than_10_sum_per_foal = not_top20_df[~not_top20_df_mask_morethan10].fillna(0).sum(axis=1)
    top20_df.sort_index(axis=1, inplace=True)
    top20_df['Others 1'] = more_than_10_sum_per_foal
    top20_df['Others 2'] = less_than_10_sum_per_foal
    color_per_species_dict = {}
    for i, species_x in enumerate(species_df.columns):
        color_per_species_dict[species_x] = (0.50196, 0.50196, 0.50196, 0.3)
    for i, species_x in enumerate(top20_df.columns):
        if species_x.startswith('Others 1'):
            color_per_species_dict[species_x] = (0.50196, 0.50196, 0.50196, 1.0)
        elif species_x.startswith('Others 2'):
            color_per_species_dict[species_x] = (0.50196, 0.50196, 0.50196, 0.3)
        else:
            color_per_species_dict[species_x] = plt.cm.Spectral(i / len(top20_df.columns))
    return color_per_species_dict, not_top_20_species, more_than_10_sum_per_foal, less_than_10_sum_per_foal, top20_df, not_top20_df_mask_morethan10


def plot_species_overview(sub1,
                          species_top_n,
                          outpath,
                          not_top20_df_mask_morethan10,
                          figsize=(10, 5)):

    grouped = sub1.groupby('nSIRS_class')
#    grouped_by_FID = sub1.groupby('FID').sum(numeric_only=True)
    # Calculate the number of FIDs in each group. Use for the subplot y axis shape.

    # Style
    # sns.set_style("whitegrid",{'axes.grid': True})
    group_sizes = grouped['FID'].nunique()
    height_ratios = group_sizes / group_sizes.min()

    fig, axes = plt.subplots(grouped.sum(numeric_only=True).shape[0], 3, figsize=figsize,
                             gridspec_kw={'width_ratios': [2, 3.5, 4.5],
                                          'height_ratios': height_ratios,
                                          'wspace': 0.1,
                                          'hspace': 0.1}, )
    for i, (group_name, group_data) in enumerate(grouped):
        ax0 = axes[i, 0]
        ax1 = axes[i, 1]
        ax2 = axes[i, 2]
        # Filter the data for this group
        group_categories_FID_sequence = sorted(group_data['FID'].unique(), key=lambda x: int(x[-2:]))
        # Based on the index of grouped categories
        species_df = group_data.pivot(index='FID', columns='Species', values=['Relative abundance']).copy()
        species_df.columns = species_df.columns.get_level_values(1)
        # species_df = species_df.loc[:, species_df.sum(axis = 0 ) > 0  ]
        species_df = species_df.reindex(group_categories_FID_sequence)
        species_df = species_df.sort_index(axis=1)

        ## Subplot (a)
        ax0.barh(y=group_categories_FID_sequence, width=species_df.sum(axis=1), color='black')
        ax0.set_xscale('log')
        ax0.set_xlim(0.0000005, 0.1)

        ax0.spines[['right', 'top']].set_visible(False)
        if i != 2:
            ax0.set_xticklabels([])
        if i == 2:
            ax0.set_xlabel('Bacterial fraction')

        # Species_top_n is sorted by alphabetical order
        top20_df = species_df[species_top_n].copy()
        not_top20_df = species_df[not_top_20_species].copy()

        more_than_10_sum_per_foal = not_top20_df[not_top20_df_mask_morethan10].fillna(0).sum(axis=1)
        less_than_10_sum_per_foal = not_top20_df[~not_top20_df_mask_morethan10].fillna(0).sum(axis=1)
        # print(more_than_10_sum_per_foal)
        # print(less_than_10_sum_per_foal)
        top20_df['Others 1'] = more_than_10_sum_per_foal
        top20_df['Others 2'] = less_than_10_sum_per_foal
        top20_df_normalized = top20_df.fillna(0).apply(lambda x: x / top20_df.fillna(0).sum(axis=1))

        # Subplot (b)
        ax1.clear()
        # THis order should be the same as the order for the next plot...?!

        # Cmap is sorted by the top20_df order (of the first batch)
        # Relative between species --- show the top 20 and then the rest in 2 colors
        alist = [color_per_species_dict[x] for x in top20_df.columns]
        custom_cmap = ListedColormap(alist, name='custom_cmap')
        stackedbar = top20_df_normalized.plot(kind='barh', stacked=True, colormap=custom_cmap, linewidth=0.8, ax=ax1,
                                              label=top20_df_normalized.columns)  #

        ax1.set_xlim(0, 1.0)
        #    ax1.set_title(f"Group: {group_name}")
        ax1.set_ylabel('')
        ax1.spines[['right', 'top']].set_visible(False)
        ax1.legend().set_visible(False)
        if i != 2:
            ax1.set_xticklabels([])
        if i == 2:
            ax1.set_xlabel('Relative bacterial composition')
        # Suplot (c) - Boxplot and scatter plot
        group_data = group_data[group_data['Relative abundance'] > 0]
        group_data = group_data[group_data['Exact abundance'] > 0]
        box = sns.boxplot(x='Relative abundance',
                          y='FID',
                          data=group_data.sort_values('FID'),
                          boxprops=dict(facecolor='none', edgecolor='black', linewidth=0.8),  # was 10
                          medianprops=dict(color='black', linewidth=0.8),
                          whiskerprops=dict(color='black', linewidth=0.8),
                          capprops=dict(color='black', linewidth=0.8),
                          flierprops=dict(markerfacecolor='black', marker='o'),
                          notch=False,
                          showfliers=False,
                          ax=ax2,

                          )
        jitter_strength_y = 0.07  # Adjust this value for y-axis jitter
        for j, species_1 in enumerate(species_df.columns):
            if species_1 in top20_df.columns:
                pass
            elif species_df[species_1].sum() > 0:
                scatter = ax2.scatter(species_df[species_1],
                                      np.arange(len(group_categories_FID_sequence)) + np.random.normal(0,
                                                                                                       jitter_strength_y,
                                                                                                       len(group_categories_FID_sequence)),
                                      color=color_per_species_dict[species_1],
                                      s=10,
                                      )

        for j, species_top in enumerate(top20_df.columns[:-2]):
            scatter = ax2.scatter(species_df[species_top],
                                  range(len(group_categories_FID_sequence)) + np.random.normal(0, jitter_strength_y,
                                                                                               len(group_categories_FID_sequence)),
                                  color=color_per_species_dict[species_top],
                                  s=20,
                                  linewidths=0.8,
                                  edgecolor='black',
                                  )
        ax2.spines[['right', 'top']].set_visible(False)

        ax2.set_yticks(ticks=range(len(group_categories_FID_sequence)), labels=group_categories_FID_sequence)
        ax2.set_ylabel('')

        if i != 2:
            ax2.set_xticklabels([])

        ax2.set_xscale('log')
        ax2.invert_yaxis()
        ax2.set_xlim(0.00000001, 0.2)
        #    ax2.set_title(f"Group: {group_name}")
        print(f'Row {i} done')

    # plt.tight_layout()
    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, bbox_to_anchor=(1.01, 1.0), loc='center left')  # ,
    # Save and show the figure
    plt.savefig(outpath, bbox_inches='tight')
    plt.savefig(outpath + '.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_alpha_diversity(df, septic_group, outpath, custom_palette=['#6066B6', '#7CA2C2', '#A84750', ]):
    #    custom_palette = ['#6066B6', '#7CA2C2', '#A84750', ] #'grey']

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
    sns.boxplot(x='nSIRS_class', y='Shannon Index', data=alpha_diversity_results, ax=axes[1][1],
                palette=custom_palette,
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
    plt.savefig(outpath)
    plt.show()


def plot_age_covariant(df, alpha_diversity_results_with_sS):
    df2 = df.copy()
    custom_palette_with_sS = ['#6066B6', '#7CA2C2', '#A84750','grey' ]
    custom_order = ['H', 'nS-', 'S+', 'sS-']
    df2['nSIRS_class'] = pd.Categorical(df2['nSIRS_class'], categories=custom_order, ordered=True)
    df2['nSIRS_class'] = df2['nSIRS_class'].fillna('sS-')
    df2['Foal age at presentation (hours)'] = basic['Foal age at presentation (hours)'].copy()
    df3 = df2[df2['FID'] != 'F06']

    # Suppl Fig 6A: Age and bacterial fraction scatterplot
    df3= df3.dropna(subset=['Foal age at presentation (hours)', 'Bacterial fraction'])
    x = df3['Foal age at presentation (hours)'].values
    y = df3['Bacterial fraction'].values
    corr_coeff, _ = pearsonr(x, y)
    slope, intercept = np.polyfit(x, y, 1)

    fig, ax = plt.subplots(figsize=(5, 4))
    plt.text(0.05, 0.95, f'Pearson r = {corr_coeff:.2f}', ha='left', va='center', transform=plt.gca().transAxes,
             fontsize=7)
    plt.plot(x, slope * x + intercept, color='black', linewidth=1)  # Line of best fit
    sns.scatterplot(x='Foal age at presentation (hours)',
                    y='Bacterial fraction',
                    hue='nSIRS_class',
                    data=df3,
                    palette=custom_palette_with_sS,
                    ax=ax)
    ax.set_ylim(0, 0.0005)
    ax.set_xlim(0, 160)
    plt.tight_layout()
    plt.savefig('../../output-figures/SupplFig5A_foal_age_bacterial_load_without_F06.pdf')
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
    plt.savefig('../../output-figures/SupplFig5B_age_group.pdf')
    plt.show()
    # Suppl Figure Species richness 6C
    alpha_diversity_results_with_sS['Foal age at presentation (hours)'] = basic['Foal age at presentation (hours)'][:-7].copy()
    alpha_diversity_results_with_sS= alpha_diversity_results_with_sS.dropna(subset=['Foal age at presentation (hours)', 'Species Richness','Shannon Index'])

    fig, ax = plt.subplots(figsize=(6, 4))
    x = alpha_diversity_results_with_sS['Foal age at presentation (hours)'].values
    y = alpha_diversity_results_with_sS['Species Richness'].values
    corr_coeff, _ = pearsonr(x, y)
    slope, intercept = np.polyfit(x, y, 1)
    plt.plot(x, slope * x + intercept, color='black', linewidth=1)  # Line of best fit
    plt.text(0.05, 0.95, f'Pearson r = {corr_coeff:.2f}', ha='left', va='center', transform=plt.gca().transAxes,
             fontsize=7)

    sns.scatterplot('Foal age at presentation (hours)', 'Species Richness',
                    data=alpha_diversity_results_with_sS, hue='nSIRS_class',
                    palette=custom_palette_with_sS, ax=ax)
    ax.set_xlim(0, 160)
    ax.legend(bbox_to_anchor=(1.3, 1))
    plt.tight_layout()
    plt.savefig('../../output-figures/SupplFig5C_Age_Species_richness_test.pdf')
    plt.show()

    # Suppl Figure 6D Shannon index
    fig, ax = plt.subplots(figsize=(6, 4))
    x = alpha_diversity_results_with_sS['Foal age at presentation (hours)'].values
    y = alpha_diversity_results_with_sS['Shannon Index'].values
    corr_coeff, _ = pearsonr(x, y)
    slope, intercept = np.polyfit(x, y, 1)
    plt.plot(x, slope * x + intercept, color='black', linewidth=1)  # Line of best fit
    plt.text(0.05, 0.95, f'Pearson r = {corr_coeff:.2f}', ha='left', va='center', transform=plt.gca().transAxes,
             fontsize=7)
    sns.scatterplot('Foal age at presentation (hours)', 'Shannon Index',
                    data=alpha_diversity_results_with_sS, hue='nSIRS_class',
                    palette=custom_palette_with_sS, ax=ax)
    ax.set_ylim(0, 3)
    ax.set_xlim(0, 160)
    ax.legend(bbox_to_anchor=(1.3, 1))
    plt.tight_layout()
    plt.savefig('../../output-figures/SupplFig5D_Age_shannon_index_test.pdf')
    plt.show()
    # Export dataframes
    df[['FID','Bacterial fraction','nSIRS_class','Lived']].to_csv('../../output/02_tables/04_source_data/Fig3_bacterial_load_test.csv')
    alpha_diversity_results_with_sS[['FID', 'nSIRS_class','Shannon Index','Species Richness']].to_csv('../../output/02_tables/04_source_data/Fig3_diversity_measurements_test.csv', index = False)


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

    long_table = pd.read_csv('../../output/02_tables/04_source_data/contaminantFree_bacteria_long.csv')
    df_original = bac.iloc[:, :bac.columns.get_loc('all_filtered_reads_exc_listed_here') + 1]

    # Fig 2A-C
    df_new = df_original.copy()
    df_new[['FID', 'nSIRS_class', 'nSIRS_score']] = bac[['FID', 'nSIRS_class', 'nSIRS_score']].copy()
    df_new = df_new[df_new['nSIRS_class'] != 'sS-']
    df_original_melted = df_new.melt(id_vars=['FID', 'nSIRS_class','nSIRS_score'], var_name='Species', value_name='Exact Rarefied abundance')
    df_original_long_table = long_table.merge( right = df_original_melted, how  = 'left')
    sub1 = df_original_long_table[(long_table['Contaminant'] == False) &  (long_table['Species'] != 'all_filtered_reads_exc_listed_here') ]
    species_top_n = get_top_n_relative_abundance(sub1)
    species_df = sub1.pivot(index='FID', columns='Species', values=['Relative abundance'], )
    species_df.columns = species_df.columns.get_level_values(1)
    species_df = species_df.loc[:, species_df.sum(axis=0) > 0]
    (color_per_species_dict,
     not_top_20_species,
     more_than_10_sum_per_foal,
     less_than_10_sum_per_foal,
     top20_df,
     not_top20_df_mask_morethan10) = get_top_species_names(
        sub1, species_df, species_top_n)

    outpath = f'../../output-figures/Fig2ABC.pdf'
    plotting = plot_species_overview(sub1,
                          species_top_n,
                          outpath,
                          not_top20_df_mask_morethan10,
                          figsize=(10, 5))


    # Get minimum filtered read count in all samples
    min_read_count = int(df_original.sum(numeric_only = True, axis=1).min())
    max_read_count = int(df_original.sum(numeric_only = True, axis=1).max())
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
    alpha_diversity_results_with_sS['nSIRS_class'] = pd.Categorical(alpha_diversity_results['nSIRS_class'],
                                                            categories=['H', 'nS-', 'S+','sS-'])
    alpha_diversity_results_with_sS.to_csv('../../output/02_tables/04_source_data/alpha_diversity_results_with_sS.csv')

    alpha_diversity_results['nSIRS_class'] = pd.Categorical(alpha_diversity_results['nSIRS_class'],
                                                            categories=['H', 'nS-', 'S+',])  # 'sS-'
    alpha_melted = alpha_diversity_results.melt(id_vars=['FID', 'nSIRS_class'],
                                                value_vars=['Species Richness', 'Shannon Index', ],
                                                var_name='measurement',
                                                value_name='value')  # 'Simpson Index'

    alpha_melted['nSIRS_class'] = pd.Categorical(alpha_melted['nSIRS_class'], categories=['H', 'nS-', 'S+', ]) # 'sS-'
    # Plot
    outpath = f'../../output-figures/Fig2DEFG_alpha_diversity_atleast{n}_rarefied.pdf'
    plot_alpha_diversity(df, septic_group, outpath, custom_palette=['#6066B6', '#7CA2C2', '#A84750', ])


    # Beta diversity:
    sub1_excl_sS_rarefied = sub1_rarefied[sub1_rarefied['nSIRS_class'] != 'sS-']
    sub1_excl_sS_rarefied['Exact Rarefied abundance (log)'] = sub1_excl_sS_rarefied['Exact Rarefied abundance'].apply(
        lambda x: np.log(x) if x != 0 else np.log(x + 0.0000001))
    bray_curtis_matrix_log, pivot_df3 = get_bray_curtis_matrix_from_long_table(sub1_excl_sS_rarefied,
                                                                               'Exact Rarefied abundance (log)')
    categories = ['H', 'nS-', 'S+', ] #'sS-'
    palette = ['#6066B6', '#7CA2C2', '#A84750', ]
    color_map = dict(zip(categories, palette))
    # Map the nSIRS_class to colors
    row_colors = df_rarefied_no_sS['nSIRS_class'].map(color_map)
    out_beta = '../../output-figures/Fig6_clustering_bray_curtis_matrix.pdf'
    # Plotting
    plot_bray_curtis_clustermap(bray_curtis_matrix_log, row_colors, out_beta)

    # Suppl Figs Age as a co-variant
    plot_age_covariant(df, alpha_diversity_results_with_sS)