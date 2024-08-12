# Get decontam input tables
import pandas as pd


def get_batches_exclude_small_batch(value_counts, at_least = 4):
    per_batch_count_dict = value_counts.to_dict()
    batches_to_exclude = []
    for key in per_batch_count_dict.keys():
        if per_batch_count_dict[key] < at_least:
            batches_to_exclude.append(key)
        else:
            pass
    return batches_to_exclude


def get_subtable_based_on_values_and_columns(df, key, batch_numbers):
    sub_df = df[~df[key].isin(batch_numbers)]
    return sub_df


def get_samples_for_batched_analysis(df_stats, all_taxid, key, at_least):
    # Get Type A batches
    isolation_batch = df_stats[key]
    batches_to_exclude = get_batches_exclude_small_batch(isolation_batch.value_counts(), at_least)
    # Subset samples that fit the criteria
    df_stats_batched_analysis = df_stats[~df_stats[key].isin(batches_to_exclude)]
    batch_index_remaining = df_stats_batched_analysis[key]
    all_taxid_batched_analysis = all_taxid.loc[df_stats_batched_analysis.index]
    return all_taxid_batched_analysis, df_stats_batched_analysis, batch_index_remaining


if __name__ == '__main__':
    # Path for output:
    outpath = '../../output/02_tables/03_intermediate/'
    # Global parameters
    TOTAL_NUMBER_OF_FOAL_SEQUENCING_RUNS_plus1 = 32
    # Input variables
    df_stats = pd.read_csv("../../output/02_tables/02_data_merged/all_basic_stats.csv")
    df_stats = df_stats[:TOTAL_NUMBER_OF_FOAL_SEQUENCING_RUNS_plus1]
    # Input matrix
    all_taxid = pd.read_csv("../../output/02_tables/02_data_merged/EquAllRS_08_all_species_for_decontam_taxid.csv")
    # Subset samples excluding positive and negative controls because they go through different steps in prep.
    all_taxid = all_taxid[:TOTAL_NUMBER_OF_FOAL_SEQUENCING_RUNS_plus1]


    ## Get batch numbers. Exclude batches that contains less than or equal to 5 samples.
    all_taxid_batched_analysis, df_stats_batched_contamA, batch_index_remaining = get_samples_for_batched_analysis(df_stats, all_taxid, 'isolation_batch', at_least=5)
    all_taxid_batched_analysis.to_csv(f"{outpath}/contamA_iso_all_taxid_BATCHED.csv")
    batch_index_remaining.to_csv(f"{outpath}/contamA_iso_BATCHES.csv", index=False)

    # Detect contaminants that inverse correlate with input volume after DNA isolation. (Taken into account only one
    # term instead of multiplication with Yield after SRSLY because then the influence of second step cannot dominate.
    ALT_contamA = 1 / df_stats_batched_contamA['Input volume Srsly (µl)']
    ## Correction: The correct column name should be the following '1 / Input volume Srsly (µl)' instead of
    ## 'Input volume Srsly (µl)' if not changed accordingly.
    ALT_contamA.columns = ['1 / Input volume Srsly (µl)']
    ALT_contamA.to_csv(f"{outpath}/ALT_contamA_input_volume_BATCHED.csv", index=False)


    ## Get batch numbers. Exclude batches that contains less than or equal to 5 samples.
    all_taxid_batched_contamB, df_stats_batched_contamB, batch_index_contamB = get_samples_for_batched_analysis(df_stats, all_taxid, 'prep_n', at_least=5)
    all_taxid_batched_contamB.to_csv(f"{outpath}/contamB_prep_all_taxid_BATCHED.csv")
    batch_index_contamB.to_csv(f"{outpath}/contamB_prep_BATCHES.csv", index=False)

    ## ConType B: contaminant related to library preparation. They inverse correlate with final library volume pooled.
    ## Volume is inversely correlated with Yield.
    contamB = 1 / (1/df_stats_batched_contamB['Yield'])
    contamB.to_csv(f"{outpath}/contamB_srsly_input_volume_BATCHED.csv", index=False)

    ## Create dataframe for regardless of batches
    all_taxid.to_csv(f"{outpath}/contam_all_taxid_NO_batch.csv")
    contamB_no_batch =  1 /  (1/df_stats['Yield'])
    contamA_alt_no_batch = 1 /   df_stats['Input volume Srsly (µl)']
    contamB_no_batch.to_csv(f"{outpath}/contamB_srsly_input_volume_NO_batch.csv", index=False)
    contamA_alt_no_batch.to_csv(f"{outpath}/ALT_contamA_input_volume_NO_batch.csv", index=False)
