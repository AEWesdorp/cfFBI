import argparse
import numpy as np
import collections
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
from scipy.stats import f
from collections import defaultdict
import os

# ## FIXED PARAMETER
# CID2 = ['H01', 'H02', 'H03', 'S01', 'S02', 'S03', 'NS01', 'NS02', 'S04', 'S05', 'S06', 'H04', 'S07', 'S08-1', 'S08-2',
#        'S09', 'S10', 'S11', 'NS03', 'S12', 'S13', 'NS04', 'S14', 'NS05', 'S15', 'NS06', 'S16', 'S17', 'NS07', 'H05',
#        'H06', 'H07', 'H08', 'NC1MQiso', 'NC1MQlib', 'NTC1', 'NTC2', 'PC1', 'PC2', 'PC3']

# Gram Negative and Gram positive pathogens frequently found in culture of sepsis foals
frequent_pathogen_n = ['Serratia',
                     'Salmonella',
                     'Pseudomonas',
                     'Proteus',
                     'Pasteurella',
                     'Pantoea',
                     'Klebsiella',
                     'Escherichia',
                     'Enterobacter',
                     'Aeromonas',
                     'Actinobacillus',
                     'Acinetobacter',]
frequent_pathogen_p = ['Streptococcus',
                       'Staphylococcus',
                       'Enterococcus',
                       'Bacillus']

## Add labels of culture results + disease type
def disease_class_culture(input_value):
    disease_state = {"NEGATIVE CULTURE - CLINICAL SUSPICION OF SEPSIS": "SepsisCulture(-)",
                     "POSITIVE CULTURE + SEPSIS": "SepsisCulture(+)",
                     "HEALTHY CONTROL": "H",
                     "NEGATIVE CULTURE - NOT SEPSIS SUSPECTED": "NotSepsisHC",
                     "POSITIVE CULTURE - NOT SUSPECTED OF SEPSIS (CONTAMINATION OR BACTEREMIA?)": "NotSepsisCul(+)contamination",
                     "NO CULTURE - NOT SEPSIS SUSPECTED (EHV-1 infection)": "HC_EHV-1",
                     "NO CULTURE - NOT SEPSIS SUSPECTED": "NotSepsisHC", }

    return disease_state[input_value]


## Add labels of disease type
def disease_class2(input_value):
    disease_state = {"NEGATIVE CULTURE - CLINICAL SUSPICION OF SEPSIS": "S",
                     "POSITIVE CULTURE + SEPSIS": "S",
                     "HEALTHY CONTROL": "H",
                     "NEGATIVE CULTURE - NOT SEPSIS SUSPECTED": "NS",
                     "POSITIVE CULTURE - NOT SUSPECTED OF SEPSIS (CONTAMINATION OR BACTEREMIA?)": "NS",
                     "NO CULTURE - NOT SEPSIS SUSPECTED (EHV-1 infection)": "NS",
                     "NO CULTURE - NOT SEPSIS SUSPECTED": "NS", }

    return disease_state[input_value]

def get_FID_from_disease_class(df, special_id_list):
    disease_states = df['nSIRS_class']
    index_numbering = 0
    FID = []
    special_index = 1
    for i, patientID in enumerate(list(df['patientID'])):
        category = list(disease_states)[i]
        if patientID in special_id_list:
            digit2 = (index_numbering+1 ) // 10
            digit1 = (index_numbering+1 ) % 10
            FID.append(f"F{digit2}{digit1}-{special_index}")
            if special_index == len(special_id_list):
                index_numbering += 1
            special_index += 1
        # PCs and NCs has no disease state category
        elif category is np.nan:
            FID.append(list(df['patientID'])[i])
        else:
            digit2 = (index_numbering+1 ) // 10
            digit1 = (index_numbering+1 ) % 10
            FID.append(f"F{digit2}{digit1}")
            index_numbering += 1
    return FID


def get_CID_from_disease_class(df, special_id_list):
    disease_states = df['disease_state']
    index_numbering = collections.Counter()
    patient_CID = []
    special_index = 1
    for i, patientID in enumerate(list(df['patientID'])):
        category = list(disease_states)[i]
        if patientID in special_id_list:
            digit2 = (index_numbering[category]+1 ) // 10
            digit1 = (index_numbering[category]+1 ) % 10
            patient_CID.append(f"{category}{digit2}{digit1}-{special_index}")
            if special_index == len(special_id_list):
                index_numbering[category] += 1
            special_index += 1


        elif category is np.nan:
            patient_CID.append(list(df['patientID'])[i])
        else:
            digit2 = (index_numbering[category]+1 ) // 10
            digit1 = (index_numbering[category]+1 ) % 10
            patient_CID.append(f"{category}{digit2}{digit1}")
            index_numbering[category] += 1
    return patient_CID

def makedir(path1):
    if not os.path.exists(path1):
        os.makedirs(path1)
        return True


def get_stats_table(a_path):
    df_stats_list = []
    for a_file in os.listdir(a_path):
        if a_file.endswith(".txt"):
            with open(a_path + a_file, 'r') as f:
                out = f.readlines()
            # print(out[0])
            pid = a_file.split('_', 1)[0]
            complete_name = a_file.split('_', 1)[-1].split('.txt')[0]
            # print('PID',step)
            count = int(out[0].split(' ')[0].strip())
            mer = complete_name.split('_')[-1]
            step = complete_name.split('_')[1]
            read = complete_name.split('_')[0]
            df_stats_list.append([pid, complete_name, count, step, mer, read])
        else:
            pass
    df_stats_pre = pd.DataFrame(df_stats_list,
                            columns=["patientID", "prefix", "count", "step", "k-mer spike-in count", "read"])
    df_stats_pre = df_stats_pre[df_stats_pre['read'] == 'R1']
    df_stats = df_stats_pre.pivot(index='patientID', columns='prefix', values='count')[
        ['R1_01_raw_fastq',
         'R1_02_50mer',
         'R1_02_150mer',
         'R1_02_100mer',
         'R1_02_spike_in_unmapped',
         'R1_03_uniq_fastq',
         'R1_04_fastp_fastq',
         'R1_05_adapt_remov_fastq',
         'R1_06_host_mapp_fastq',
         ]]

    ## Calculation Reads mapped to host genome is
    # (1) unique reads
    # (2) after fastp quality filter
    # (3) after adapter removal
    # (4) minus unmapped reads
    df_stats['Mapped to host genomes'] = df_stats['R1_05_adapt_remov_fastq'] - df_stats['R1_06_host_mapp_fastq']
    df_stats['Not mapped'] = df_stats['R1_06_host_mapp_fastq'].copy()
    df_stats.drop(columns=[ 'R1_06_host_mapp_fastq'], inplace=True)
    df_stats['R1_02_spike_in_unmapped'] =  df_stats['R1_02_spike_in_unmapped']/ 4
    # kmer_stats = df_stats[df_stats.columns.contains('mer')])]]

    return df_stats


def get_MT(path_MT):
    with open(path_MT, 'r') as f:
        mt_doc = f.readlines()
    sample_names = []
    values_var = []
    for i, line in enumerate(mt_doc):
        line = line.strip()
        if i % 3 == 0:
            sample_names.append(line.split('_')[0])
        elif i % 3 == 1:
            column_names = line.split('\t')
        else:
            values_var.append(line.split('\t'))

    MT_load = pd.DataFrame(columns=column_names, data=values_var, index=sample_names, )
    MT_load = MT_load.astype(float, errors='ignore')
    MT_load.rename(columns={'numreads': 'MT_reads', 'coverage': 'MT_cov'}, inplace=True)
    return MT_load


def parse_combined_kreport_genera(in_file,max_samples=10000):
    df = pd.read_csv(in_file, header=42, sep='\t', skipinitialspace=True)
    # Read in patientID and index table
    with open(in_file, 'r') as f:
        a = f.readlines(max_samples)
    sample_count = int(a[0].split(':')[-1].strip())
    samples_ordered = []
    samples_dict_index = {}
    samples_dict_index2 = {}
    # print(len(a))
    i = 1
    for x in range(2, sample_count + 2):
        # print(x)
        # print(a[x])
        patient_id = a[x].split('/')[-1].split('_')[0].split('#')[-1]
        samples_ordered.append(patient_id)
        samples_dict_index[i] = patient_id
        samples_dict_index2[patient_id] = i
        i += 1

    # Get the amount of all nodes including sub-trees
    count_all = df[df.filter(like='_all').columns].copy()
    count_all.loc[:, ['#perc', 'lvl_type', 'taxid', 'name']] = df.loc[:, ['#perc', 'lvl_type', 'taxid', 'name']].copy()
    ## Alternative: Get the amount of all nodes exactly only each node
    ## count_exact_level = df[df.filter(like='_lvl').columns].copy()

    # Select for genera only
    rows_to_keep = []
    for i, row in count_all.iterrows():
        if row['lvl_type'] == 'G':
            rows_to_keep.append(i)
    all_classified = count_all[count_all['name'] == 'root']
    all_un_classified = count_all[count_all['name'] == 'unclassified']
    all_bacteria = count_all[count_all['name'] == 'Bacteria']

    # Get classified reads, unclassified reads, bacterial read count
    collect_classified_reads = []
    collected_unclassified_reads = []
    collect_bacteria_reads = []
    for x in samples_dict_index.keys():
        collect_classified_reads.append(all_classified[f'S{x}_all'].values[0])
    for x in samples_dict_index.keys():
        collected_unclassified_reads.append(all_un_classified[f'S{x}_all'].values[0])
    for x in samples_dict_index.keys():
        collect_bacteria_reads.append(all_bacteria[f'S{x}_all'].values[0])

    # Summary statistics with correct
    df_summary = pd.DataFrame(columns=samples_ordered,
                              index=['classified_read_count', 'unclassified_read_count', 'bacteria_read_count'],
                              data=[collect_classified_reads, collected_unclassified_reads, collect_bacteria_reads]).T
    df_summary = df_summary.reset_index(names = ['patientID'])

    return count_all, df_summary, samples_ordered


def _get_rows_to_keep(count_all, taxonomy_level = 'G'):
    rows = []
    for i, row in count_all.iterrows():
        if row['lvl_type'] == taxonomy_level:
            rows.append(i)
    return rows

def get_sub_table(count_all, ID, taxonomy_level,):
    if taxonomy_level == 'all':
        df = count_all.copy()
    else:
        rows_to_keep = _get_rows_to_keep(count_all, taxonomy_level=taxonomy_level)
        df = count_all.iloc[rows_to_keep]
    df.columns = ['tot_all'] + ID + ['#perc', 'lvl_type', 'taxid', 'name']

    df.index = df['name']
    df = df.T
    df['FID'] = df.index
    df.reset_index(drop=True)
    # df = df.drop(columns=['tot_all', '#perc', 'lvl_type', 'taxid', 'name']).T
    return df


def get_sub_tree(count_all, ID, taxonomy_level):
    # Step 1: get rows in domain X
    # Step 2: get rows with taxonomy level Y
    subtrees = {}

    rows_to_keep = _get_rows_to_keep(count_all, taxonomy_level='D')
    last_row = count_all.shape[0]
    rows_to_keep.append(last_row)
    row_range_indexes = [(rows_to_keep[i], rows_to_keep[i + 1]) for i in range(len(rows_to_keep) - 1)]
    names = count_all.iloc[rows_to_keep[:-1]]['name']

    for name, range_x in zip(names, row_range_indexes):
        df_per_domain = count_all.iloc[range_x[0]:range_x[1]]
        txt_level_rows = _get_rows_to_keep(df_per_domain, taxonomy_level=taxonomy_level)
        df = count_all.iloc[txt_level_rows]
        df.columns = ['tot_all'] + ID + ['#perc', 'lvl_type', 'taxid', 'name']
        # interchange taxid -> name for columns in decontam table
        df.index = df['taxid']
        df = df.T
        df['FID'] = df.index
        df.reset_index(drop=True, inplace=True)
        subtrees[name] = df

    return subtrees


def prep_decontam(dfs,
                  clinical_stats_wetlab3,
                  domain='Bacteria',
                  metadata_columns= ['FID', 'Yield', 'isolation_batch', 'prep_n', 'strip_n', 'total_filtered_effective_reads', ],
                  subset_max_index=None,
                  ):
    if domain == 'all':
        df = dfs
    else:
       df = dfs[domain]
    temp2 = clinical_stats_wetlab3[metadata_columns].merge(df, on='FID', how='left')
    temp2 = temp2[:subset_max_index]
    for column in metadata_columns:
        temp2[column].to_csv(f'{args.output}/EquAllRS_08_{column}_for_decontam.csv', index=False)
    output = temp2.drop(columns=metadata_columns)
    row_sums = output.sum(axis=1)
    # Get remaining reads
    output['all_filtered_reads_exc_listed_here'] = temp2['total_filtered_effective_reads'] - row_sums
    return output


# parse argument
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--clinical', required=True, help='clinical table')
    parser.add_argument('--kreport', required=True, help='Combined kreport from kraken')
    parser.add_argument('--wetlab-isolation', required=True, help='wetlab table')
    parser.add_argument('--wetlab-dna', required=True, help='wetlab table')
    parser.add_argument('--wetlab-yield', required=True, help='wetlab table')

    parser.add_argument('--stats', required=True, help='Path to folder containing statistics read count table')
    parser.add_argument('--MT', required=True, help='Path to mitochondrial read count file')

    parser.add_argument('--output', required=True, help='Path to folder for exporting merged table and sorted csv')
    args = parser.parse_args()
    return args



if __name__ == '__main__':
    # TODO: Add nSIRS score
    # Healthy = H
    # 3+ = S+
    # 1-2 = nS-
    # 0 = sS-

#    SPECIAL_ID_LIST = ['21048121', '21048122',]
    SPECIAL_ID_LIST = []
    TOTAL_NUMBER_OF_FOAL_SAMPLES_plus1 = 32
    TOTAL_NUMBER_OF_FOAL_SEQUENCING_RUNS_plus1 = 32

    args = parse_arguments()
    makedir(args.output)
    # Read clinical, add disease classes
    clinical = pd.read_csv(args.clinical, skipinitialspace=True)
    clinical.drop(columns = ['Clinical Assessment'], inplace = True)
    # clinical['disease_state_culture'] = clinical['Clinical Assessment'].apply(disease_class_culture)
    # clinical['disease_state'] = clinical['Clinical Assessment'].apply(disease_class2)
    # temp = (clinical.T)
    # temp.insert(0, 'subtable', ['clinical']* clinical.shape[1])
    # clinical = temp.T

    # Read stats (Quality filter, spike-in count, Mapped/Unmapped)
    df_stats = get_stats_table(args.stats)

    # Get Mitochondrial reads and merge to stats table
    MT_load = get_MT(args.MT)
    df_stats = df_stats.join(MT_load['MT_reads'])
    # Exclude foal 21048122 becuase it was a second timepoint after treatment
    df_stats = df_stats[df_stats.index != '21048122']
    # Create CID
    clinical_stats = pd.merge(df_stats, clinical, on='patientID', how='left', )
#    CID = get_CID_from_disease_class(clinical_stats, SPECIAL_ID_LIST)
    FID = get_FID_from_disease_class(clinical_stats, SPECIAL_ID_LIST)
#    clinical_stats.insert(0, 'CID', CID)
    clinical_stats.insert(0, 'FID', FID)

    # Generate Conversion table
#    patientID_to_CID = clinical_stats[['patientID','CID']].set_index('patientID')
#    patientID_to_CID.to_csv(f'{args.output}/patientID_CID.csv')
    patientID_to_FID = clinical_stats[['patientID','FID']].set_index('patientID')
    patientID_to_FID.to_csv(f'{args.output}/patientID_FID.csv')

#    dict_patientID_to_CID = patientID_to_CID.to_dict()['CID']
#    dict_patientID_to_CID = defaultdict(int, dict_patientID_to_CID)
    dict_patientID_to_FID = patientID_to_FID.to_dict()['FID']
    dict_patientID_to_FID = defaultdict(int, dict_patientID_to_FID)

    ## Export clinical table
    clinical2 = clinical_stats.drop(columns = df_stats.columns)
    clinical2.to_csv(f'{args.output}/clinical.csv', index=False)

    # Import wetlab tables
    ## (1) Wetlab isolation
    wetlab_isolation = pd.read_csv(args.wetlab_isolation, sep=',',)

    # (2) Plasma and DNA (saved as UTF-8 (comma separated))
    wetlab_dna = pd.read_csv(args.wetlab_dna, sep=',',
                             skiprows=0)
#    df_exp['CID'] = df_exp.index.to_series().apply(lambda x: dict_patientID_to_CID[x])
    # df_exp.reset_index(drop = True, inplace = True)
    wetlab_dna.rename(mapper = {'DNA yield per ml plasma (if elution in 28 or 35 ul) (ng/ml)':'cfDNA yield (ng/ml)'},axis=1, inplace = True)
    wetlab_dna.drop(columns = ['DNA yield','DNA yield (if elution in 28 or 35 ul)'], inplace = True)
    ## (3) SRSLY PREP DNA yield table saved as UTF-8
    wetlab_PREP = pd.read_csv(args.wetlab_yield, sep=',',
                         skiprows=0)

    clinical_stats_wetlab3 = clinical_stats.merge(wetlab_isolation, on = 'patientID', how='left')
    clinical_stats_wetlab3 = clinical_stats_wetlab3.merge(wetlab_dna, on ='patientID', how='left')
    clinical_stats_wetlab3 = clinical_stats_wetlab3.merge(wetlab_PREP, on = 'patientID', how='left')

    # Note
    # For cfDNA and MT load, we can analyse all foals. For clinical results, exclude S08-2 because they are the same.
    count_all, df_summary, samples_ordered = parse_combined_kreport_genera(args.kreport, )
    # Create
    taxid_dict = count_all[['taxid','name']].set_index('taxid', inplace = False)

    clinical_stats_wetlab3 = clinical_stats_wetlab3.merge(df_summary, on = 'patientID',how='left')
    clinical_stats_wetlab3['total_filtered_effective_reads'] = clinical_stats_wetlab3[['Mapped to host genomes', 'unclassified_read_count', 'classified_read_count']].sum(axis=1)
    clinical_stats_wetlab3.to_csv(f'{args.output}/all_basic_stats.csv', index=False)
    summary_metadata = clinical_stats_wetlab3[['FID','patientID','nSIRS_class', 'isolation_batch', 'prep_n', 'strip_n']]
    summary_metadata.to_csv(f'{args.output}/summary_basic_stats.csv', index=False)
 #   CID_ordered = [dict_patientID_to_CID[x] for x in samples_ordered]
    FID_ordered = [dict_patientID_to_FID[x] for x in samples_ordered]

    # Tables for merging with other metadata
    df_species = get_sub_table(count_all, FID_ordered, 'S')
    df_genera = get_sub_table(count_all, FID_ordered, 'G',)
    df_domains = get_sub_table(count_all, FID_ordered, 'D',)
    df_kindoms = get_sub_table(count_all, FID_ordered, 'K',)
    count_all2 = get_sub_table(count_all, FID_ordered, 'all')

    # Tables for summary statistics
    species = clinical_stats_wetlab3.merge(df_species, on='FID', how='left')
    genera = clinical_stats_wetlab3.merge(df_genera, on = 'FID',how='left')
    domain = clinical_stats_wetlab3.merge(df_domains, on='FID', how='left')
    domain = domain.merge(df_kindoms, on='FID', how='left')

    all = clinical_stats_wetlab3.merge(count_all2, on = 'FID',how='left')

    df_genera.to_csv(f'{args.output}/EquAllRS_08_genera.csv', index=False)
    df_domains.to_csv(f'{args.output}/EquAllRS_08_domains.csv', index=False)
    df_kindoms.to_csv(f'{args.output}/EquAllRS_08_kindoms.csv', index=False)
    count_all2.to_csv(f'{args.output}/EquAllRS_08_all.csv', index=False)

    species.to_csv(f'{args.output}/EquAllRS_08_merged_species.csv', index=False)
    genera.to_csv(f'{args.output}/EquAllRS_08_merged_genera.csv', index=False)
    domain.to_csv(f'{args.output}/EquAllRS_08_merged_domains_kindom.csv', index=False)
    all.to_csv(f'{args.output}/EquAllRS_08_merged_all.csv', index=False)



    # Get df for decontam
    # Taxonomy level G, S
    dfs_per_genera = get_sub_tree(count_all, FID_ordered, taxonomy_level='G', )
    print(dfs_per_genera['Bacteria'].shape)
    decontam_bac_genera = prep_decontam(dfs_per_genera, clinical_stats_wetlab3, domain='Bacteria', subset_max_index=TOTAL_NUMBER_OF_FOAL_SEQUENCING_RUNS_plus1)
    decontam_bac_genera.to_csv(f'{args.output}/EquAllRS_08_bacteria_genera_for_decontam_taxid.csv', index=False)
    dfs_per_species = get_sub_tree(count_all, FID_ordered, taxonomy_level='S', )
    print(dfs_per_species['Bacteria'].shape)


    decontam_bac_species = prep_decontam(dfs_per_species, clinical_stats_wetlab3, domain='Bacteria', subset_max_index=TOTAL_NUMBER_OF_FOAL_SEQUENCING_RUNS_plus1)
    # interchange taxid -> name for columns in decontam table
    decontam_bac_species.to_csv(f'{args.output}/EquAllRS_08_bacteria_species_for_decontam_taxid.csv', index=False)
    #metadata_columns =  ['CID', 'Yield', 'isolation_batch', 'prep_n', 'strip_n', 'total_filtered_effective_reads', ]
    ## All

#    decontam_all_species = prep_decontam(df_species, clinical_stats_wetlab3, domain='all')
    df_species_taxid = df_species.copy()
    # interchange taxid -> name for columns in decontam table
    df_species_taxid.columns = df_species_taxid.loc['taxid']
    df_species_taxid['FID'] = df_species_taxid.index
    df_species_taxid.reset_index(drop=True, inplace= True)
    # interchange taxid -> name for columns in decontam table
    df_species_taxid.drop(columns = ['taxid'], inplace=True)
    print(df_species_taxid.shape)

    decontam_all_species = prep_decontam(df_species_taxid, clinical_stats_wetlab3, domain='all', subset_max_index=TOTAL_NUMBER_OF_FOAL_SEQUENCING_RUNS_plus1)
    # interchange taxid -> name for columns in decontam table
    decontam_all_species.to_csv(f'{args.output}/EquAllRS_08_all_species_for_decontam_taxid.csv', index=False)

