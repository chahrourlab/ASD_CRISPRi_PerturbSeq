#!/usr/bin/env python
import os
import re
import sys
import collections
import argparse
import itertools
import matplotlib
import glob
import math
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import scipy.stats as stats
import scipy.sparse as sp_sparse
import scanpy as sc 
import scanpy.external as sce
from cnmf import cNMF
from collections import defaultdict
from scipy import sparse, io
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import pickle
import time
import gseapy as gp
from gseapy import Msigdb
import requests

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

output_directory = '/output_directory'
run_name = 'Neuron_Prod1_cNMF_run'

cnmf_obj = cNMF(output_dir=output_directory, name=run_name)

density_threshold = 0.1
K_values = [50,100,200,250,500]

region_df = pd.read_pickle('Prod1_neuronal_cell_region_df.pkl')

def load_results(self, K, density_threshold, n_top_genes=300):
        scorefn = self.paths['gene_spectra_score__txt'] % (K, str(density_threshold).replace('.', '_'))
        tpmfn = self.paths['gene_spectra_tpm__txt'] % (K, str(density_threshold).replace('.', '_'))
        usagefn = self.paths['consensus_usages__txt'] % (K, str(density_threshold).replace('.', '_'))
        spectra_scores = pd.read_csv(scorefn, sep='\t', index_col=0).T
        spectra_tpm = pd.read_csv(tpmfn, sep='\t', index_col=0).T

        usage = pd.read_csv(usagefn, sep='\t', index_col=0)
        usage = usage.div(usage.sum(axis=1), axis=0)

        try:
            usage.columns = [int(x) for x in usage.columns]
        except:
            print('Usage matrix columns include non integer values')

        top_genes = []
        for gep in spectra_scores.columns:
            top_genes.append(list(spectra_scores.sort_values(by=gep, ascending=False).index[:n_top_genes]))

        top_genes = pd.DataFrame(top_genes, index=spectra_scores.columns).T
        return(usage, spectra_scores, spectra_tpm, top_genes)

cNMF_output_dict_all_Ks = {}

#Iterate over each K value
for K in K_values:
    usage_norm, gep_scores, gep_tpm, topgenes = load_results(cnmf_obj, K=K, density_threshold=density_threshold)
    # Rename the columns of usage_norm
    usage_norm.columns = ['Usage_%d' % i for i in usage_norm.columns]
    # Store the results in the dictionary with dynamic keys
    cNMF_output_dict_all_Ks[f'usage_norm_k_{K}'] = usage_norm
    cNMF_output_dict_all_Ks[f'gep_scores_k_{K}'] = gep_scores
    cNMF_output_dict_all_Ks[f'gep_tpm_k_{K}'] = gep_tpm
    cNMF_output_dict_all_Ks[f'topgenes_k_{K}'] = topgenes

## PERTURBATION SENSITIVITY
def compare_scores(usage_norm_k, cell_gRNA):
    results = []
    non_targeting_rows = cell_gRNA.index[cell_gRNA['non'] == 1]
    # Loop over all score columns in usage_norm_k
    for score_col in usage_norm_k.columns:
        # Compare non-targeting with itself
        scores_non_targeting_self = usage_norm_k.loc[non_targeting_rows, score_col]
        stat_self, p_value_self = mannwhitneyu(scores_non_targeting_self, scores_non_targeting_self, alternative='two-sided')
        # Store self-comparison results
        results.append({
            'Gene': 'non-targeting (self)',
            'Score Column': score_col,
            'Statistic': stat_self,
            'p-value': p_value_self,
            'Log2 Fold Change': np.nan  # No change when comparing with itself
        })
        # Loop over each gene column in cell_gRNA
        for gene_col in cell_gRNA.columns:
            if gene_col != 'non':
                comparison_rows = cell_gRNA.index[cell_gRNA[gene_col] == 1]
                if len(comparison_rows) == 0:
                    continue
                scores_comparison = usage_norm_k.loc[comparison_rows, score_col]
                # Perform Mann-Whitney U test
                stat, p_value = mannwhitneyu(scores_non_targeting_self, scores_comparison, alternative='two-sided')
                # Calculate log2 fold change
                mean_non_targeting = scores_non_targeting_self.mean()
                mean_comparison = scores_comparison.mean()
                log2fc = np.log2((mean_comparison + 1e-9) / (mean_non_targeting + 1e-9))
                results.append({
                    'Gene': gene_col,
                    'Score Column': score_col,
                    'Statistic': stat,
                    'p-value': p_value,
                    'Log2 Fold Change': log2fc
                })

    # Create a DataFrame from results
    results_df = pd.DataFrame(results)

    # Apply Benjamini-Hochberg correction
    reject, corrected_p_values = multipletests(results_df['p-value'], method='fdr_bh')[:2]

    # Add corrected p-values to the results DataFrame
    results_df['Corrected p-value (BH)'] = corrected_p_values

    return results_df

# Initialize a dictionary to store results for each K
usage_per_perturbation = {}

# Iterate over each K value
for K in K_values:
    # Get the usage_norm DataFrame for the current K
    usage_norm_k = cNMF_output_dict_all_Ks[f'usage_norm_k_{K}']
    
    # Initialize an empty DataFrame to store the results
    result_k = pd.DataFrame(columns=usage_norm_k.columns)

    # Assuming adata_subset is defined and cell_gRNA is extracted
   # cell_gRNA = adata_subset.obs.iloc[:, 4:]
    cell_gRNA = region_df.copy()
    for col in cell_gRNA.columns:
        # Find rows in cell_gRNA where the value is 1
        rows_with_1 = cell_gRNA.index[cell_gRNA[col] == 1]
        
        # Extract those rows from the usage_norm DataFrame
        usage_subset = usage_norm_k.loc[rows_with_1]
        
        # Calculate the average for each column in the usage_subset
        averages = usage_subset.mean()
        
        # Append the averages as a new row in the result DataFrame
        result_k.loc[col] = averages
    
    # Store the result DataFrame in the results dictionary
    usage_per_perturbation[f'result_k_{K}'] = result_k

compare_score_pert_all_ks = {}

for k, df in cNMF_output_dict_all_Ks.items():
    if isinstance(df, pd.DataFrame) and 'usage_norm_' in k:  # Check for the format
        compare_score_pert_all_ks[k] = compare_scores(df, cell_gRNA)

with open('Prod1_neuron_compare_score_pert_all_K.pkl', 'wb') as file:
    pickle.dump(compare_score_pert_all_ks, file)

with open('Prod1_neuron_usage_per_perturbation.pkl', 'wb') as file:
    pickle.dump(usage_per_perturbation, file)

### Doing GO for all Ks
# Assume cNMF_output_dict_all_Ks is already defined and contains your data
top_genes_dict = {}

# Loop through each K value
for K in K_values:
    gep_zscore = cNMF_output_dict_all_Ks[f"gep_scores_k_{K}"]  # Access the DataFrame directly
    top_genes = []
    ngenes = 300
    
    for gep in gep_zscore.columns:
        # Sort and get the top genes
        top_genes.append(list(gep_zscore.sort_values(by=gep, ascending=False).index[:ngenes]))
    
    # Store the result in a dictionary with K as the key
    top_genes_dict[K] = pd.DataFrame(top_genes, index=gep_zscore.columns).T

# Now top_genes_dict contains top genes for each K

def enrich_for_columns_batch(top_genes_df, batch_size=10):
    all_top_terms = []
    total_columns = top_genes_df.shape[1]
    
    for start in range(0, total_columns, batch_size):
        end = min(start + batch_size, total_columns)
        columns = top_genes_df.columns[start:end]

        # Retry logic until the batch is processed successfully
        while True:
            try:
                batch_results = []

                # Process each column in the batch
                for column in columns:
                    # Convert column to a list, dropping NaN values
                    gene_list = top_genes_df[column].dropna().tolist()

                    # Perform enrichment analysis (assumed you have 'gp' imported and ready)
                    enr_res = gp.enrichr(
                        gene_list=gene_list,
                        organism='Human',
                        gene_sets=['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023'],
                        cutoff=0.05
                    )

                    # Get top 10 terms for the current column
                    top_terms = enr_res.results.sort_values(by='Adjusted P-value')
                    top_terms['Source'] = column  # Add the source column for tracking
                    batch_results.append(top_terms)

                # Combine batch results
                all_top_terms.extend(batch_results)
                break  # Exit retry loop if successful

            except Exception as e:
                print(f"Error processing columns {start} to {end}: {e}")
                print("Retrying the batch...")
                time.sleep(10)  # Sleep for 10 seconds before retrying (adjust as needed)

        # Sleep after each batch to avoid hitting the API rate limit
        time.sleep(30)  # Adjust sleep time as needed between batches

    # Combine all top terms into a single DataFrame
    combined_top_terms = pd.concat(all_top_terms, ignore_index=True)

    return combined_top_terms


all_enriched_results = {}

# Loop through the dictionary and apply the enrich_for_columns_batch function to each DataFrame
for key, top_genes_df in top_genes_dict.items():
    print(f"Processing DataFrame with {key} columns.")
    
    # Call the enrich_for_columns_batch function with the key (as it represents the number of columns to process at once)
    enriched_df = enrich_for_columns_batch(top_genes_df, batch_size=10)  # The function still processes in batches of 10
    
    # Store the enriched results corresponding to the key
    all_enriched_results[key] = enriched_df

with open("Prod1_neuron_GO_all_enriched_results.pkl", "wb") as pkl_file:
    pickle.dump(all_enriched_results, pkl_file)
