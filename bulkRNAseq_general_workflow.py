#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: gouthamvasam
"""

# Import necessary libraries
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

# Quality Control of Raw Reads
def run_fastqc(fastq_files, output_dir):
    fastqc_cmd = ['fastqc'] + fastq_files + ['-o', output_dir]
    subprocess.run(fastqc_cmd, check=True)

# Read Alignment
def run_star(genome_dir, read_files, output_dir):
    star_cmd = [
        'STAR', 
        '--genomeDir', genome_dir,
        '--readFilesIn'] + read_files + [
        '--outFileNamePrefix', output_dir,
        '--outSAMtype', 'BAM', 'SortedByCoordinate'
    ]
    subprocess.run(star_cmd, check=True)

# Quantification
def run_featurecounts(input_bam, annotation_file, output_file):
    featurecounts_cmd = [
        'featureCounts',
        '-a', annotation_file,
        '-o', output_file,
        input_bam
    ]
    subprocess.run(featurecounts_cmd, check=True)

# Normalization and Differential Expression Analysis
def run_deseq2(counts, col_data):
    # Enable automatic conversion between Pandas dataframes and R dataframes
    pandas2ri.activate()

    # Import R's "DESeq2" package
    deseq = importr('DESeq2')

    # Convert the Pandas DataFrames into R DataFrames
    counts_r = pandas2ri.py2rpy(counts)
    col_data_r = pandas2ri.py2rpy(col_data)

    # Run DESeq2
    dds = deseq.DESeqDataSetFromMatrix(countData=counts_r, colData=col_data_r, design=ro.Formula('~ condition'))
    dds = deseq.DESeq(dds)
    res = deseq.results(dds)

    # Convert the results to a Pandas DataFrame
    results = pandas2ri.rpy2py(res)

    return results

# MA Plot
def plot_ma(de_results):
    de_results['log2FoldChange'] = np.log2(de_results['foldChange'])
    de_results['baseMean'] = np.log2(de_results['baseMean'] + 1)

    plt.figure(figsize=(10, 6))
    plt.scatter(de_results['baseMean'], de_results['log2FoldChange'], alpha=0.5)
    plt.title('MA Plot')
    plt.xlabel('Average Expression (log scale)')
    plt.ylabel('Log Fold Change')
    plt.axhline(y=0, linestyle='--', color='red')
    plt.show()

# Volcano Plot
def plot_volcano(de_results):
    plt.figure(figsize=(10, 6))
    plt.scatter(de_results['log2FoldChange'], -np.log10(de_results['pvalue']), alpha=0.5)
    plt.title('Volcano Plot')
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10 p-value')
    plt.axhline(y=-np.log10(0.05), linestyle='--', color='red')  # P-value threshold
    plt.axvline(x=1, linestyle='--', color='blue')  # Fold change threshold
    plt.axvline(x=-1, linestyle='--', color='blue')  # Fold change threshold
    plt.show()

# Heatmap
def plot_heatmap(normalized_data, de_results):
    # Select top differentially expressed genes for the heatmap
    top_genes = de_results.nsmallest(50, 'pvalue')['gene']

    # Subset the normalized data for the top genes
    data_subset = normalized_data.loc[top_genes]

    # Create the heatmap
    plt.figure(figsize=(10, 10))
    sns.heatmap(
        data_subset,
        cmap='viridis',
        xticklabels=True,
        yticklabels=True,
        linewidths=0.1,
        linecolor='gray',
        cbar_kws={"shrink": 0.5}
    )

    plt.title('Heatmap of Top Differentially Expressed Genes')
    plt.xlabel('Samples')
    plt.ylabel('Genes')
    plt.tight_layout()
    plt.show()

def main():
    # Example usage:
    fastq_files = ['sample1.fastq', 'sample2.fastq']  # List your FASTQ files here
    output_dir = 'fastqc_results'
    run_fastqc(fastq_files, output_dir)

    genome_dir = 'path/to/genome_indices'
    read_files = ['sample1_R1.fastq', 'sample1_R2.fastq']  # Replace with your read files
    output_dir = 'star_output/'
    run_star(genome_dir, read_files, output_dir)

    input_bam = 'sorted.bam'
    annotation_file = 'genome_annotation.gtf'
    output_file = 'gene_counts.txt'
    run_featurecounts(input_bam, annotation_file, output_file)

    # Load count data and the sample information
    counts = pd.read_csv('gene_counts.csv', index_col=0)
    col_data = pd.read_csv('sample_conditions.csv', index_col=0)

    de_results = run_deseq2(counts, col_data)

    # Save the results
    de_results.to_csv('differential_expression_results.csv')

    plot_ma(de_results)
    plot_volcano(de_results)

    # Load normalized expression data
    normalized_data = pd.read_csv('normalized_expression.csv', index_col=0)

    plot_heatmap(normalized_data, de_results)

if __name__ == "__main__":
    main()