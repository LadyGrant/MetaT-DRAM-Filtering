# MetaT-DRAM-Filtering
How to add expression from metaT to DRAM
DRAM is designed to tell you genomic potential but not expression. We can run DRAM on a filtered annotation file that is matched with the MetaT to get a product output with expression.

# 1.) Filtering DRAM annotations to "on" genes.
First we need to filter the DRAM annotations.tsv file down to the non-zero ("on" genes) gene IDs in the metatranscriptomics expression table (geTMM values). This will generate two new annotations file for each treatment (CT and ST).    
The metatranscriptomics expression table values are averaged abundance values across each sample for each treatment.
```
# Author: Reed Woyda: reed.woyda@colostate.edu

import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import re
import time

# Load the data files
annotations_df = pd.read_csv("annotations.tsv", sep='\t', low_memory=False)
mag_genes_df = pd.read_excel("mag_getmms_ge5_1X.xlsx", sheet_name="mag_genes_average")

# Function to filter annotations based on exact partial match and non-zero average values
def filter_annotations(annotations_chunk, non_zero_genes, gene_map):
    print(f"Processing chunk with {len(annotations_chunk)} rows...")
    start_time = time.time()
    pattern = '|'.join([re.escape(gene) + r'(?=$|_)' for gene in non_zero_genes])
    filtered_chunk = annotations_chunk[annotations_chunk.iloc[:, 0].str.contains(pattern, regex=True)]

    # Add gene_metaT column with corresponding gene values
    filtered_chunk.loc[:, 'gene_metaT'] = filtered_chunk.iloc[:, 0].apply(lambda x: ';'.join(gene_map[gene] for gene in non_zero_genes if re.s$

    # Rearrange columns to make gene_metaT the second column
    cols = list(filtered_chunk.columns)
    cols.insert(1, cols.pop(cols.index('gene_metaT')))
    filtered_chunk = filtered_chunk[cols]

    end_time = time.time()
    print(f"Chunk processed in {end_time - start_time:.2f} seconds.")
    return filtered_chunk

# Get gene_short values with non-zero averages for CT and ST
ct_non_zero_genes = set(mag_genes_df[mag_genes_df['CT_average'] != 0]['gene_short'])
st_non_zero_genes = set(mag_genes_df[mag_genes_df['ST_average'] != 0]['gene_short'])

# Create a map from gene_short to gene
gene_map = mag_genes_df.set_index('gene_short')['gene'].to_dict()

# Parallel processing to filter annotations in chunks
def parallel_filtering(annotations_df, non_zero_genes, gene_map, max_workers=5, chunk_size=1000):
    chunks = [annotations_df[i:i + chunk_size] for i in range(0, len(annotations_df), chunk_size)]
    print(f"Total chunks to process: {len(chunks)}")
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        filtered_annotations = pd.concat(executor.map(filter_annotations, chunks, [non_zero_genes]*len(chunks), [gene_map]*len(chunks)))
    return filtered_annotations

# Filter annotations for CT_average and ST_average using multiple CPUs with chunk processing
print("Starting CT_average filtering...")
filtered_ct_annotations = parallel_filtering(annotations_df, ct_non_zero_genes, gene_map, max_workers=5)
print("CT_average filtering complete.")

print("Starting ST_average filtering...")
filtered_st_annotations = parallel_filtering(annotations_df, st_non_zero_genes, gene_map, max_workers=5)
print("ST_average filtering complete.")

# Save the filtered annotations to output files
filtered_ct_annotations.to_csv("CT_bins_average.annotations.tsv", sep='\t', index=False)
filtered_st_annotations.to_csv("ST_bins_average.annotations.tsv", sep='\t', index=False)

print("Filtering complete. Output files generated.")
```
Slurm scrip to run filter-annotations.py :
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --time=14-00:00:00
#SBATCH --job-name=filter_annotations
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=eryn.grant@colostate.edu
#SBATCH --mem=50gb
#SBATCH --partition=wrighton-hi,wrighton-low

python filter-annotations.py
```
Your output files include:
    CT_bins_average.annotations.tsv    
    ST_bins_average.annotations.tsv
    
