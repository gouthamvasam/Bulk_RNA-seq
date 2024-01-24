#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 12:09:46 2024

@author: gouthamvasam
"""

# Import necessary libraries
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Quality Control of Raw Reads
def run_fastqc(fastq_files, output_dir):
    fastqc_cmd = ['fastqc'] + fastq_files + ['-o', output_dir]
    subprocess.run(fastqc_cmd, check=True)

# Read Alignment
def run_star(genome_dir, read_files, output_prefix):
    star_cmd = [
        'STAR', 
        '--genomeDir', genome_dir,
        '--readFilesIn'] + read_files + [
        '--outFileNamePrefix', output_prefix,
        '--outSAMtype', 'BAM', 'SortedByCoordinate',
        '--runThreadN', '8'  # Adjust based on your available CPU cores
    ]
    subprocess.run(star_cmd, check=True)
    
# Quantification
def run_featurecounts(input_bam, annotation_file, output_file):
    featurecounts_cmd = [
        'featureCounts',
        '-T', '8',  # Adjust based on your available CPU cores
        '-a', annotation_file,
        '-o', output_file,
        input_bam
    ]
    subprocess.run(featurecounts_cmd, check=True)
    
# Normalization and Differential Expression Analysis
def run_deseq2(counts_file, col_data_file, output_file):
    # Write the R script to a file
    r_script = """
    library(DESeq2)
    counts <- as.matrix(read.csv('{counts_file}', row.names=1, header=TRUE))
    col_data <- read.csv('{col_data_file}', row.names=1)
    dds <- DESeqDataSetFromMatrix(countData=counts, colData=DataFrame(col_data), design=~ condition)
    dds <- DESeq(dds)
    res <- results(dds)
    write.csv(as.data.frame(res), file='{output_file}')
    
    # Save the normalized counts
    norm_counts <- counts(dds, normalized=TRUE)
    write.csv(as.data.frame(norm_counts), file='normalized_counts.csv')
    """.format(counts_file=counts_file, col_data_file=col_data_file, output_file=output_file)
    
    with open('deseq2_script.R', 'w') as file:
        file.write(r_script)
    
    # Run the R script
    subprocess.run(['Rscript', 'deseq2_script.R'], check=True)

# MA Plot
def plot_ma(de_results):
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=np.log2(de_results['baseMean']), y=de_results['log2FoldChange'], alpha=0.5)
    plt.title('MA Plot')
    plt.xlabel('Average Expression (log scale)')
    plt.ylabel('Log2 Fold Change')
    plt.axhline(y=0, linestyle='--', color='red')
    plt.show()

# Volcano Plot
def plot_volcano(de_results):
    # Handle potential NA or extreme values in 'padj'
    de_results = de_results.replace([np.inf, -np.inf], np.nan).dropna(subset=["padj"], how="all")
    
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=de_results['log2FoldChange'], y=-np.log10(de_results['padj']), alpha=0.5)
    plt.title('Volcano Plot')
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10 Adjusted p-value')
    plt.axhline(y=-np.log10(0.05), linestyle='--', color='red')  # Adjusted p-value threshold
    plt.axvline(x=1, linestyle='--', color='blue')  # Log2 Fold Change positive threshold
    plt.axvline(x=-1, linestyle='--', color='blue')  # Log2 Fold Change negative threshold
    plt.show()

# Heatmap
def plot_heatmap(normalized_data, de_results):
    # Select top differentially expressed genes for the heatmap
    top_genes = de_results.nsmallest(50, 'padj').index

    # Subset the normalized data for the top genes
    data_subset = normalized_data.loc[top_genes]

    # Create the heatmap
    plt.figure(figsize=(10, 10))
    sns.heatmap(
        data_subset,
        cmap='viridis',
        annot=False,  # Turn this to True if you want to see the values in the heatmap
        xticklabels=True,
        yticklabels=True,
        linewidths=0.1,
        linecolor='gray',
        cbar_kws={"shrink": 0.5}
    )
    plt.title('Heatmap of Top Differentially Expressed Genes')
    plt.xlabel('Samples')
    plt.ylabel('Genes')
    plt.tight_layout()  # Adjust the plot to ensure everything fits without overlapping
    plt.show()
    
def main():
    # Example usage:
    fastq_files = ['sample1_R1.fastq', 'sample1_R2.fastq', 'sample2_R1.fastq', 'sample2_R2.fastq']  # List your FASTQ files here
    output_dir = 'fastqc_results'
    run_fastqc(fastq_files, output_dir)

    genome_dir = 'path/to/genome_indices'
    read_files = ['sample1_R1.fastq', 'sample1_R2.fastq']  # Replace with your read files
    output_prefix = 'star_output/sample1_'
    run_star(genome_dir, read_files, output_prefix)

    input_bam = 'star_output/sample1_Aligned.sortedByCoord.out.bam'  # Adjust this path as needed
    annotation_file = 'path/to/genome_annotation.gtf'
    output_file = 'gene_counts.txt'
    run_featurecounts(input_bam, annotation_file, output_file)

    counts_file = 'gene_counts.txt'
    col_data_file = 'sample_conditions.csv'
    de_output_file = 'differential_expression_results.csv'
    run_deseq2(counts_file, col_data_file, de_output_file)

    # Read the DESeq2 results
    de_results = pd.read_csv(de_output_file, index_col=0)

    plot_ma(de_results)
    plot_volcano(de_results)

    # Load normalized expression data (this would typically be obtained from DESeq2 after running the analysis)
    normalized_data = pd.read_csv('normalized_counts.csv', index_col=0)

    # Assuming the normalized data and DESeq2 results are properly formatted and indexed
    plot_heatmap(normalized_data, de_results)

if __name__ == "__main__":
    main()
