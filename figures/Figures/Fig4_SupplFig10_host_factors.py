import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
# import scipy
# from scipy.stats import f
# from scipy.stats import mannwhitneyu
# import statsmodels.api as sm
import os
import statannotations
from statannotations.Annotator import Annotator
import matplotlib.colors as mcolors


## Local conda env: py310
## Print versions
print("numpy", np.__version__) # 1.26.0
print("pandas", pd.__version__) # 1.5.3
print("seaborn", sns.__version__) # 0.11.2
print("statannotations", statannotations.__version__) # 0.6.0

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
sns.set_style("whitegrid", )# {'axes.grid': False}


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

def exclude_contam(df):
    THRESHOLD = 10
    basic = pd.read_csv('../../output/02_tables/02_data_merged/all_basic_stats.csv')
    labels = basic[['FID', 'nSIRS_class', 'nSIRS_score']].copy()
    # basic['nSIRS_score']
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

    return df

def plot_grouped_bar(categories, analysis, base_on_df, custom_palettes, out_path, dont_plot_outlier=False, figsize=(9, 6), test_loc = 'outside'):
    # Create subplots
    fig, axes = plt.subplots(len(categories), len(analysis), figsize=figsize)
    # Loop over categories and features
    for i, (category_x, df_x, custom_palette) in enumerate(zip(categories, base_on_df, custom_palettes)):
        for j, feature_y in enumerate(analysis):

            if len(categories) > 1 or len(analysis) > 1:
                ax = axes[i][j]
            else:
                ax = axes
            sns.boxplot(y=feature_y, x=category_x, data=df_x, ax=ax, palette=custom_palette, showfliers=False)
            if dont_plot_outlier:
                ylims = ax.get_ylim()
                sns.stripplot(y=feature_y, x=category_x, data=df_x, palette=custom_palette, ax=ax, jitter=0.2,
                              linewidth=1, )
                ax.set(ylim=ylims)
            else:
                sns.stripplot(y=feature_y, x=category_x, data=df_x, palette=custom_palette, ax=ax, jitter=0.2,
                              linewidth=1, )
            ax.spines[['right', 'top']].set_visible(False)
            ax.set_xlabel('')

            # Using a package...
            pairs = []
            for k in range(len(df[category_x].cat.categories)):
                for m in range(k + 1, len(df[category_x].cat.categories)):
                    pairs.append((df[category_x].cat.categories[k], df[category_x].cat.categories[m]))
            if ('H', 'nS-') in pairs:
                pairs.remove(('H', 'nS-'))

            annotator = Annotator(ax, pairs, data=df, x=category_x, y=feature_y,
                                  order=df[category_x].cat.categories)  # hide_non_significant=True
            annotator.configure(test='Mann-Whitney', text_format='star', loc=test_loc, verbose=0)
            annotator.apply_and_annotate()

    plt.tight_layout()
    plt.savefig(out_path)
    plt.savefig(out_path + '.png', dpi=300)
    plt.show()


def plot_grouped_bar_with_hue(categories, hues, grouped_df, custom_palettes, out_path):
    fig, axes = plt.subplots(len(categories), len(hues), figsize=(8, 6))
    # Loop over categories and features
    for j, category_x in enumerate(categories):
        for i, (feature_y, custom_palette) in enumerate(zip(hues, custom_palettes)):
            if feature_y == 'Lived':
                sub_r1 = grouped_df[(grouped_df['read'] == category_x) & (grouped_df['nSIRS_class'] == 'S+')].copy()
            elif feature_y == 'nSIRS_class':
                sub_r1 = grouped_df[(grouped_df['read'] == category_x) & (grouped_df['nSIRS_class'] != 'sS-')].copy()
                sub_r1['nSIRS_class'] = pd.Categorical(sub_r1['nSIRS_class'], categories=['H', 'nS-', 'S+'])
            else:
                print('Double check the sub-selection of sub_r1')
                exit()
            sns.boxplot(x='EndMotif', y='nsCount', hue=feature_y, data=sub_r1, linewidth=1, palette=custom_palette,
                        fliersize=0, ax=axes[i][j])
            swarm = sns.stripplot(x='EndMotif', y='nsCount', hue=feature_y, data=sub_r1, jitter=0.3, dodge=True,
                                  linewidth=1, palette=custom_palette, size=4, ax=axes[i][j])
            # boxplot.set_ylabel
            axes[i][j].set_ylabel("Relative fraction (log10)")
            axes[i][j].set_xlabel(f"{sub_r1['End'].unique()[0]} motif")

            axes[i][j].spines[['right', 'top']].set_visible(False)

            if j == 1:
                axes[i][j].legend(bbox_to_anchor=(1.2, 1))
            else:
                axes[i][j].legend().set_visible(False)
    plt.tight_layout()
    plt.savefig(out_path)
    plt.savefig(out_path + '.png', dpi=300)
    plt.show()


if __name__ == '__main__':
    df = import_data()
    df = preprocessing(df)
    df = exclude_contam(df)

    # 1. Plot cfDNA conc. and MT DNA fraction

    custom_order = ['H', 'nS-', 'S+']
    df['nSIRS_class'] = pd.Categorical(df['nSIRS_class'], categories=custom_order, ordered=True)
    custom_order2 = ['Y', 'N']

    df['Lived'] = pd.Categorical(df['Lived'], categories=custom_order2, ordered=True, )
    replacement = {'Y': 'S+ Lived', 'N': 'S+ Died'}
    df['Lived'].replace(replacement, inplace=True)

    septic_group = df[df['nSIRS_class'] == 'S+'].copy()

    fig2_data = df[['FID', 'nSIRS_class', 'patientID', 'Total cfDNA conc. in plasma (ng/ml)', 'Foal MT cfDNA conc. in plasma (ng/ml)',
                    'Foal MT cfDNA fraction', 'Bacterial fraction', ]][:32]
    fig2_data.to_csv('../../output/03_microbial/source_data/Fig2ABCD.csv')

    # Define your data
    analysis = ['Total cfDNA conc. in plasma (ng/ml)', 'Foal MT cfDNA fraction', 'Bacterial fraction', ]
    categories = ['nSIRS_class', 'Lived']
    # Ensure `df` and `septic_group` are defined earlier in your script.
    base_on_df = [df, septic_group]
    custom_palettes = [['#6066B6', '#7CA2C2', '#A84750'], ['#e35f6c', '#85363e']]
    out_path1 = "../../output/03_microbial/figures/Fig102_exclude_outlier_GRID_test.pdf"
    plot_grouped_bar(categories, analysis, base_on_df, custom_palettes, out_path1, dont_plot_outlier=True)
    out_path2 = "../../output/03_microbial/figures/Fig102_include_outlier_GRID_test.pdf"
    plot_grouped_bar(categories, analysis, base_on_df, custom_palettes, out_path2, dont_plot_outlier=False)


    # 2. Plot end motifs
    # Load data
    foal_nonMT_IS = pd.read_csv("../../output/04_fragmentomics/foal_cohort_EquCabAll_host_nonMT_motif.csv")
    metadata = pd.read_csv("../../output/02_tables/02_data_merged/summary_basic_stats.csv")
    complete = pd.read_csv("../../output/02_tables/02_data_merged/all_basic_stats.csv")

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
        # .query("read == 'R1'")
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
    grouped_df['nsCount'] = grouped_df.groupby(['nSIRS_class', 'FID', 'read'])['sCount'].transform(
        lambda x: np.log10(x / x.sum() / 0.25))
    # Remove NA and zero values from nsCount
    grouped_df = grouped_df[(~grouped_df['nsCount'].isna()) & (grouped_df['nsCount'] != 0)]
    def map_x_value(x):
        d = {"R1": "5' end", "R2": "3' end"}
        return d[x]
    grouped_df['End'] = grouped_df['read'].apply(lambda x: map_x_value(x))
    grouped_df = grouped_df.merge(df[['FID', 'Lived']], how='left', on='FID')
    complete = complete[['FID', 'Lived']]
    complete.index = complete.FID
    grouped_df['Lived'] = grouped_df['FID'].apply(lambda x: complete.to_dict()['Lived'][x])
    replacement = {'Y': 'S+ Lived', 'N': 'S+ Died'}
    grouped_df['Lived'].replace(replacement, inplace=True)
    custom_palettes = [['#6066B6', '#7CA2C2', '#A84750'], ['#e35f6c', '#85363e']]
    categories = ['R1', 'R2']
    hues = ['nSIRS_class', 'Lived']
    out_path3 = "../../output/03_microbial/figures/Fig2CDGH_EndMotif_GRID_test.pdf"
    plot_grouped_bar_with_hue(categories, hues, grouped_df, custom_palettes, out_path3)
    export_end_motif_pivot = grouped_df.pivot(index=['FID'], columns=['EndMotif', 'End'], values=['sCount', 'nsCount'])
    out_path4 = "../../output/02_tables/04_source_data/EndMotif_table_pivot.csv"
    export_end_motif_pivot.to_csv(out_path4)
    out_path5 = "../../output/02_tables/04_source_data/EndMotif_table_long.csv"
    grouped_df.to_csv(out_path5)
    # 3. Plot covariance

    order_prep_n = ['Srsly prep 1', 'Srsly prep 2', 'Srsly prep 3']
    order_strip_n = ['Strip 1', 'Strip 2', 'Strip 3', 'Strip 4', 'Strip 5', ]
    df['prep_n'] = pd.Categorical(df['prep_n'], categories=order_prep_n, ordered=True)
    df['strip_n'] = pd.Categorical(df['strip_n'], categories=order_strip_n, ordered=True, )
    df['isolation_batch'] = df['isolation_batch'].fillna(0).astype('int32')
    df['isolation_batch'] = pd.Categorical(df['isolation_batch'], categories=[1, 2, 3, 4, 5, 6], ordered=True, )
    df['isolation_by'] = pd.Categorical(df['isolation_by'], categories=['E', 'N'], ordered=True, )

    df.rename(mapper={'DNA concentration (elution total volume = 28ul-N, 35ul-E)': 'DNA conc. (ng/ul)'}, axis=1,
              inplace=True, )
    categories = ['prep_n',
                  'strip_n',
                  'isolation_batch', ]
    analysis = ['Total cfDNA conc. in plasma (ng/ml)',
                'Foal MT cfDNA fraction',
                'Bacterial fraction', ]
    base_on_df2 = [df,df,df]
    custom_palettes2 = ['Greys', 'Greys', 'Greys']
    out_path4 = "../../output/03_microbial/figures/SupplFig203_batch_effect_include_outlier_UPDATE_test.pdf"
    plot_grouped_bar(categories, analysis, base_on_df2, custom_palettes2, out_path4, dont_plot_outlier=False, figsize=(12, 12), test_loc= 'inside')

    # 4. extra analysis -- comparing MT cfDNA conc. in plasma and MT cfDNA fraction
    analysis = ['Total cfDNA conc. in plasma (ng/ml)', 'Foal MT cfDNA conc. in plasma (ng/ml)',
                'Foal MT cfDNA fraction', ]
    categories = ['nSIRS_class', 'Lived']
    base_on_df = [df, septic_group]
    custom_palettes = [['#6066B6', '#7CA2C2', '#A84750'], ['#e35f6c', '#85363e']]
    out_path5 = "../../output/03_microbial/figures/Extra Figure two analysis of MT.pdf"
    plot_grouped_bar(categories, analysis, base_on_df, custom_palettes, out_path5, dont_plot_outlier=False)

