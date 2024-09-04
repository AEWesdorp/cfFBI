# Instruction to self: Use py310 env

import pandas as pd


def is_contam(value, contaminants_iso):
    if value in contaminants_iso:
        return True
    else:
        return False


def is_higher_than_double_threshold(row, max_dict):
    species = row['Species']
    background = max_dict[species]
    if row['Relative abundance'] > background:
        return True
    else:
        return False


def is_higher_than_single_threshold(row, max_dict, nSIRS_class):
    species = row['Species']
    background = max_dict[(nSIRS_class, species)]
    if row['Relative abundance'] > background:
        return True
    else:
        return False
if __name__ == '__main__':
    # P-value cutoff
    CUTOFF = 0.25
    # Prevalence of contaminants
    EXIST_IN_AT_LEAST_X_SAMPLES = 6
    # Minimal counts of reads in each genera to be called as elevated
    THRESHOLD = 10
    basic = pd.read_csv('../../output/02_tables/02_data_merged/all_basic_stats.csv')
    labels = basic[['FID','nSIRS_class','nSIRS_score']].copy()
    labels.index = labels['FID']
    taxid_dict = pd.read_pickle( "../../output/02_tables/03_intermediate/taxid_dict.pickle.gz")
    bac = pd.read_csv('../../output/02_tables/02_data_merged/EquAllRS_08_bacteria_species_for_decontam_taxid.csv')
    bac = bac.rename(columns = taxid_dict, )
    total_count = bac.sum(axis = 1)
    bac_normalized = bac.apply(lambda x: x / total_count, axis = 0)
    bac_index = pd.read_csv('../../output/02_tables/02_data_merged/EquAllRS_08_FID_for_decontam.csv')
    bac_normalized['FID'] = bac_index['FID']
    bac_normalized.index = bac_index['FID']
    bac['Total QC reads'] = total_count
    bac['FID'] = bac_index['FID']
    bac.index = bac_index['FID']
    bac['nSIRS_class'] = labels['nSIRS_class']
    bac['nSIRS_score'] = labels['nSIRS_score']


    # Read contaminants
    # Gather results from decontam.
    # Name Contaminant type
    contypes = ['Isolation related', 'Isolation related batched', 'Library prep related',
                'Library prep related batched']
    # Gather volume paths
    path_inputs = [
        '../../output/02_tables/03_intermediate/ALT_contamA_input_volume_NO_batch.csv',
        '../../output/02_tables/03_intermediate/ALT_contamA_input_volume_BATCHED.csv',
        '../../output/02_tables/03_intermediate/contamB_srsly_input_volume_NO_batch.csv',
        '../../output/02_tables/03_intermediate/contamB_srsly_input_volume_BATCHED.csv',
    ]
    # Gather decontam output
    path_outputs = [
        '../../output/02_tables/03_intermediate/decontam_output_ALT_conA_no_batch_X_taxid.csv',
        '../../output/02_tables/03_intermediate/decontam_output_ALT_contamA_BATCHED_X_taxid.csv',
        '../../output/02_tables/03_intermediate/decontam_output_concB_no_batch_X_taxid.csv',
        '../../output/02_tables/03_intermediate/decontam_output_contamB_BATCHED_X_taxid.csv',

    ]
    contaminations = []
    for contype, path_x, path_output in zip(contypes, path_inputs, path_outputs):
        # Read CSV
        X = pd.read_csv(path_x)
        X.columns = [f'contamination type {contype}']
        con = pd.read_csv(path_output, index_col=0)
        con = con.add_suffix(f'_con{contype}')
        contaminations.append(con)
    for i in range(len(contaminations)):
        if i == 0:
            contam_merged = contaminations[i]
        else:
            contam_merged = contam_merged.join(contaminations[i])

    # Convert list to names
    raw = pd.read_csv('../../output/02_tables/03_intermediate/contamA_iso_all_taxid_BATCHED.csv', index_col=0)
    raw1 = raw.drop(columns=['all_filtered_reads_exc_listed_here'])
    species_with_at_least_one_count = raw1.loc[:, raw1.sum(axis=0) > 0].columns
    print("Species with at least 1 count:", species_with_at_least_one_count.shape[0])
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

    # Subset to species with at least 6 appearance in all samples
    prev_N = [f'prev_con{x}' for x in contypes]

    ## COLLECTING INDEXES of species with prevalence at least X.
    indexes_high_prev = []
    for prev in prev_N:
        indexes_high_prev += list(contam_merged[contam_merged[prev] >= EXIST_IN_AT_LEAST_X_SAMPLES].index)

    contam_merged_filter_prevalence = contam_merged.loc[contam_merged.index.isin(indexes_high_prev)]
    # Subset columns
    # contaminants_table_1 = contam_merged_filter_prevalence[
    #     ['p_conIsolation related', 'p_conIsolation related batched', 'p_conLibrary prep related',
    #      'p_conLibrary prep related batched', 'freq_conIsolation related']].copy()
    contaminants_table_1 = contam_merged_filter_prevalence[
        ['p_conIsolation related', 'p_conLibrary prep related','freq_conIsolation related']].copy()

    # 2. sort the two columns
    contaminants_table = contaminants_table_1.sort_values(
        by=['freq_conIsolation related', 'p_conLibrary prep related', 'p_conIsolation related'], ascending=True)
    contaminants_table = contaminants_table.sort_values(by=['freq_conIsolation related'], ascending=False)
    # 3. test if either one is smaller than cutoff
    contaminants_table['Isolation related contaminant (p<0.25)'] = contaminants_table[
        'p_conIsolation related'].apply(lambda x: x < CUTOFF)
    # contaminants_table['Isolation batch related contaminant (p<0.25)'] = contaminants_table[
    #     'p_conIsolation related batched'].apply(lambda x: x < CUTOFF)
    contaminants_table['Library prep related contaminant (p<0.25)'] = contaminants_table[
        'p_conLibrary prep related'].apply(lambda x: x < CUTOFF)
    # contaminants_table['Library prep batch related contaminant (p<0.25)'] = contaminants_table[
    #     'p_conLibrary prep related batched'].apply(lambda x: x < CUTOFF)
    contaminants_table['Contaminants'] = contaminants_table[
                                                            'Isolation related contaminant (p<0.25)'] + \
                                                        contaminants_table[
                                                            'Library prep related contaminant (p<0.25)']
    # contaminants_table['Contaminants batched'] = contaminants_table[
    #                                                     'Isolation batch related contaminant (p<0.25)'] + \
    #                                                 contaminants_table[
    #                                                     'Library prep batch related contaminant (p<0.25)']
#    contaminants_table['Contaminants'] = contaminants_table['Contaminants batched'] + contaminants_table[
#        'Contaminants not batched']

    contaminants_table.to_csv('../../output/02_tables/03_intermediate/decontam_species_joint_not_batched.csv')
    contaminants_table = contaminants_table.reset_index()
    contaminants = contaminants_table[contaminants_table['Contaminants'] == True]['name'].values
    contaminants_iso = contaminants_table[contaminants_table['Isolation related contaminant (p<0.25)'] == True]['name'].values
    contaminants_prep = contaminants_table[contaminants_table['Library prep related contaminant (p<0.25)'] == True]['name'].values
#    contaminants_iso_batched = contaminants_table[contaminants_table['Isolation batch related contaminant (p<0.25)'] == True]['name'].values
#    contaminants_prep_batched = contaminants_table[contaminants_table['Library prep batch related contaminant (p<0.25)'] == True]['name'].values


    # Melt bacteria count per sample per species
    bac_melted = bac.melt(id_vars=['FID','nSIRS_class','nSIRS_score','Total QC reads'], var_name = 'Species', value_name = 'Exact abundance' )
    bac_normalized_melted = bac_normalized.melt(id_vars=['FID'], var_name = 'Species', value_name = 'Relative abundance' )

    # Merge normalized and not normalized (Could be use to validate if normalized is correct with total_count)
    merged_df = pd.merge(bac_melted, bac_normalized_melted, on=['Species', 'FID'])
    merged_df['Genus'] = merged_df['Species'].apply(lambda x: x.split(' ')[0])
    merged_df['Isolation related contaminant (p<0.25)'] = merged_df['Species'].apply(lambda x : is_contam(x, contaminants_iso))
#    merged_df['Isolation batch related contaminant (p<0.25)'] = merged_df['Species'].apply(lambda x : is_contam(x, contaminants_iso_batched))
    merged_df['Library prep related contaminant (p<0.25)'] = merged_df['Species'].apply(lambda x : is_contam(x, contaminants_prep))
#    merged_df['Library prep batch related contaminant (p<0.25)'] = merged_df['Species'].apply(lambda x : is_contam(x, contaminants_prep_batched))

    merged_df['Contaminant'] = merged_df['Species'].apply(lambda x : is_contam(x, contaminants))
    # Add column more than 10

    merged_df[f'Exact_count_morethaninc_{THRESHOLD}'] =  merged_df['Exact abundance'].apply(lambda x: x >= THRESHOLD)
    merged_df[f'Exact_count_morethaninc_{2}'] =  merged_df['Exact abundance'].apply(lambda x: x >= 2)

    # Getting max per nSIRS class
    max_per_nSIRS_class = merged_df.groupby(['nSIRS_class','Species'])['Relative abundance'].max()
    max_per_nSIRS_class_dict = max_per_nSIRS_class.to_dict()
    max_per_nSIRS_class.head()
    # Get max per combined group (different data structure compared to above)
    H_nS_df  = merged_df[merged_df['nSIRS_class'].isin(['nS-','H'])]
    max_H_nS = H_nS_df.groupby('Species')['Relative abundance'].max()
    max_H_nS_dict = max_H_nS.to_dict()
    # Get max per combined group (different data structure compared to above) --> For testing the reverse
    S_nS_df  = merged_df[merged_df['nSIRS_class'].isin(['S+','nS-'])]
    max_nS_S = S_nS_df.groupby('Species')['Relative abundance'].max()
    max_nS_S_dict = max_nS_S.to_dict()


    merged_df['higher_than_H'] = merged_df.apply(lambda row: is_higher_than_single_threshold(row, max_per_nSIRS_class_dict, 'H'),
                                                 axis=1)
    merged_df['higher_than_nS-'] = merged_df.apply(
        lambda row: is_higher_than_single_threshold(row, max_per_nSIRS_class_dict, 'nS-'), axis=1)
    merged_df['higher_than_S+'] = merged_df.apply(lambda row: is_higher_than_single_threshold(row, max_per_nSIRS_class_dict, 'S+'),
                                                  axis=1)
    merged_df['higher_than_H_nS-'] = merged_df.apply(lambda row: is_higher_than_double_threshold(row, max_H_nS_dict),
                                                     axis=1)
    merged_df['higher_than_nS-_S'] = merged_df.apply(lambda row: is_higher_than_double_threshold(row, max_nS_S_dict),
                                                     axis=1)

    merged_df.to_csv('../../output/02_tables/04_source_data/contaminantFree_bacteria_long.csv', index=False)


