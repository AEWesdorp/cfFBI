import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, LogFormatter

# Set font size, style
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

import re


def replace_dots(input_string):
    # Replace single dots with spaces, ensuring not to touch consecutive dots
    processed_string = re.sub(r'(?<!\.)\.(?!\.)', ' ', input_string)

    # Handle consecutive dots: surround them with spaces, maintaining one dot
    processed_string = re.sub(r'\.{2,}', lambda x: x.group(0)[0] + ' ', processed_string)
    # processed_string = re.sub(r'-', ' ', processed_string)

    return processed_string


def get_contaminants(contypes, path_inputs, path_outputs, species_with_at_least_one_count):
    contaminations = []
    for contype, path_x, path_output in zip(contypes, path_inputs, path_outputs):
        # Read CSV
        X = pd.read_csv(path_x)
        X.columns = [f'contamination type {contype}']
        con = pd.read_csv(path_output, index_col=0)
        con_names = con.index.to_series().apply(lambda x: replace_dots(x))
        con.index = con_names.values
        con = con.add_suffix(f'_con{contype}')
        contaminations.append(con)

    for i in range(len(contaminations)):
        if i == 0:
            contam_merged = contaminations[i]
        else:
            contam_merged = contam_merged.join(contaminations[i])

    # Change taxid to names
    taxid_dict = pd.read_pickle("../../output/02_tables/03_intermediate/taxid_dict.pickle.gz")

    # Convert list to names
    names_of_species_with_at_least_one_count = [taxid_dict[taxid] for taxid in species_with_at_least_one_count.values]

    # Subset contam table to only include species with at least 1 count
    contam_merged.drop(index=['all_filtered_reads_exc_listed_here'], inplace=True)
    contam_merged = contam_merged.reset_index()
    contam_merged['taxid'] = contam_merged['index'].apply(lambda x: x.split('X')[-1])
    contam_merged.drop(index=contam_merged[contam_merged['taxid'] == ''].index, inplace=True)
    contam_merged['name'] = contam_merged['taxid'].apply(lambda x: taxid_dict[x])
    contam_merged.index = contam_merged['name'].copy()
    contam_merged.drop(columns=['index'], inplace=True)
    # Subset to species with at least 1 count
    contam_merged = contam_merged.loc[names_of_species_with_at_least_one_count]
    return contam_merged


## Relationship of prevalence/frequency and p-value
from matplotlib.lines import Line2D


##### Leave out black list because they are genus level. --> /hpc/compgen/projects/cf-spi/cf-spi/analysis/ewesdorp/GenusBlacklist_PMID30497919.txt
def scatter_pvalue_freq(df, species_of_interest_list, colors, out_path):
    temp = contam_merged_prep_5.loc[species_of_interest_list]

    ## Also color Actinobacillus equlli

    # TODO: Color human (other skin microbe) , K. phaffi  ?
    # TODO: plot

    figsize = (15, 8)
    fig = plt.figure(figsize=figsize)
    plots_per_row = 2
    colors = colors
    legend_handles = [Line2D([0], [0], marker='^', color='w', markerfacecolor=color, markersize=10, label=species)
                      for color, species in zip(colors, species_of_interest_list)]

    for i, contype in zip(range(0, len(contypes)), contypes):
        ax = fig.add_subplot((len(contypes) // plots_per_row) + 1, plots_per_row, i + 1)
        ax.scatter(df[f'p_con{contype}'], df[f'freq_con{contype}'], s=2, c='k', )
        # Plotting known contaminants
        ax.scatter(temp[f'p_con{contype}'], temp[f'freq_con{contype}'], s=60, marker='^', c=colors, )
        ax.set_yscale('log')
        ax.set_title(f'Testing {contype} contaminants')
        ax.set_ylabel('Sum Frequency')
        ax.set_xlabel('p-value')
    #    ax.legend(species_of_interest_list,bbox_to_anchor= (1,1))
    fig.legend(handles=legend_handles, loc='upper right', bbox_to_anchor=(0.95, 0.4
                                                                          ))

    plt.tight_layout()
    plt.savefig(out_path)
    plt.savefig(out_path + '.png')

    plt.show()



if __name__ == '__main__':
    raw = pd.read_csv('../../output/02_tables/03_intermediate/contamA_iso_all_taxid_BATCHED.csv',index_col = 0)
    raw1 = raw.drop(columns = ['all_filtered_reads_exc_listed_here'])
    species_with_at_least_one_count = raw1.loc[:, raw1.sum(axis = 0) > 0].columns
    print("Species with at least 1 count:", species_with_at_least_one_count.shape)

    # Contaminant type
    contypes = ['Isolation related', 'Library prep related']
    path_inputs = [
        '../../output/02_tables/03_intermediate/ALT_contamA_input_volume_NO_batch.csv',
        '../../output/02_tables/03_intermediate/contamB_srsly_input_volume_NO_batch.csv',
    ]
    path_outputs = [
        '../../output/02_tables/03_intermediate/decontam_output_ALT_conA_no_batch_X_taxid.csv',
        '../../output/02_tables/03_intermediate/decontam_output_concB_no_batch_X_taxid.csv',
    ]
    contam_merged = get_contaminants(contypes, path_inputs, path_outputs, species_with_at_least_one_count)

    # PARAMETERS
    p_N = [f'p_con{x}' for x in contypes]
    prev_N = [f'prev_con{x}' for x in contypes]
    ## COLLECTING INDEXES
    indexes_high_prev = []
    for prev in prev_N:
        indexes_high_prev += list(contam_merged[contam_merged[prev] > 5].index)

    contam_merged_prep_5 = contam_merged.loc[contam_merged.index.isin(indexes_high_prev)]
    contam_merged_prep_5_plot = contam_merged.loc[contam_merged.index.isin(indexes_high_prev), p_N]

    # 1. Take p values of 2 tests
    contam_final_p = contam_merged_prep_5[
        ['p_conIsolation related', 'p_conLibrary prep related', 'freq_conIsolation related']].copy()

    # 2. sort the two columns
    contam_final_p_sorted = contam_final_p.sort_values(
        by=['freq_conIsolation related', 'p_conLibrary prep related', 'p_conIsolation related'], ascending=True)
    contam_final_p_sorted = contam_final_p_sorted.sort_values(by=['freq_conIsolation related'], ascending=False)

    # 3. test if either one is smaller than cutoff
    CUTOFF = 0.25
    # contam_true = contam_final_p.copy()
    # contam_true
    contam_final_p_sorted['Isolation related contaminant (p<0.25)'] = contam_final_p_sorted[
        'p_conIsolation related'].apply(lambda x: x < CUTOFF)
    contam_final_p_sorted['Library prep related contaminant (p<0.25)'] = contam_final_p_sorted[
        'p_conLibrary prep related'].apply(lambda x: x < CUTOFF)
    contam_final_p_sorted['Contaminants'] = contam_final_p_sorted['Isolation related contaminant (p<0.25)'] + \
                                            contam_final_p_sorted['Library prep related contaminant (p<0.25)']

    #### This table is used in pathogen_table_new.py. Should be integrated there.
    contam_final_p_sorted.to_csv('../../output/02_tables/03_intermediate/decontam_species_joint.csv')

    # Plot correlation of samples Suppl Fig D-E
    concA = pd.read_csv('../../output/02_tables/03_intermediate/ALT_contamA_input_volume_NO_batch.csv')
    concB = pd.read_csv('../../output/02_tables/03_intermediate/contamB_srsly_input_volume_NO_batch.csv')
    all = pd.read_csv('../../output/02_tables/02_data_merged/EquAllRS_08_all_species_for_decontam_taxid.csv')
    all.head()
    all_index = pd.read_csv('../../output/02_tables/02_data_merged/EquAllRS_08_FID_for_decontam.csv')
    all.index = all_index['FID']
    taxid_dict = pd.read_pickle("../../output/02_tables/03_intermediate/taxid_dict.pickle.gz")
    all = all.rename(columns=taxid_dict, )
    total_read_count =  all.sum(axis = 1)
    all_normalized = all.apply(lambda x: x / total_read_count)


    fig, axes = plt.subplots(1, 2, figsize=(6, 3))
    axes[0].scatter(x=concA.values, y=all_normalized['Homo sapiens'].values, s=8, c='#C68C74')
    axes[0].set_xscale('log')
    axes[0].set_yscale('log')
    axes[0].set_xlabel('1 / Input volume Srsly (µl)')
    axes[0].set_ylabel('Frequency')

    axes[1].scatter(x=concB.values, y=all_normalized['Komagataella phaffii'].values, s=8, c='red')
    axes[1].set_xscale('log')
    axes[1].set_yscale('log')

    axes[1].set_xlabel('Final Yield')
    axes[1].set_ylabel('Frequency')
    axes[0].set_title('Homo sapiens')
    axes[1].set_title('Komagataella phaffii')

    plt.savefig('../../output-figures/Suppl203D-E_examples_correlation_test.pdf')
    plt.show()

    # Plot p-value distribution Suppl Fig F-G
    figsize = (15, 8)
    fig = plt.figure(figsize=figsize)
    bins = 50
    plots_per_row = 2
    # Arbitary p-values cutoff
    CUTOFF = 0.25
    for i, contype in zip(range(0, len(contypes)), contypes):
        ax = fig.add_subplot((len(contypes) // plots_per_row) + 1, plots_per_row, i + 1)
        ########## Plot only species that occur at least in 5 samples. #####
        subdf = contam_merged_prep_5[f'p_con{contype}']
        ax.hist(subdf, range=(0, 1), bins=bins, color='k')
        ax.axvline(CUTOFF, color='#A84750', linewidth=1)
        ax.set_title(f'Testing {contype} contaminations')
        ax.set_xlabel('p-value')
        ax.set_ylabel('amount of species (n)')
    plt.tight_layout()
    plt.savefig(f"../../output-figures/Suppl203F_G_p-value-{CUTOFF}-{bins}bins-YIELD_test.pdf")
    plt.savefig(f"../../output-figures/Suppl203F_G_p-value-{CUTOFF}-{bins}bins-YIELD_test.pdf.png")
    plt.show()

    ## Suppl Fig H - I scatter plot p-value and
    # Blue: human
    # Green: Cutibacterium acnes
    # Red: K. phaffi
    # Yellow: Staphy
    # purple,orange
    out_path_scatter = '../../output-figures/Suppl203_H-I_test.pdf'
    colors = ['red','#c68c74','#a84950','lightgreen']
    #colors = ['red', 'blue','lightgreen','darkorange']
    # blue: '#6066B6'
    #### A. equuli: IT is NOT contaminant --> Change shape.
    species_of_interest = ['Komagataella phaffii', 'Homo sapiens','Cutibacterium acnes', 'Actinobacillus equuli']
    scatter_pvalue_freq(contam_merged_prep_5, species_of_interest, colors, out_path_scatter)


