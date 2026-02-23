#!/usr/bin/env python3

"""
Performs comprehensive annotation and filtering of genomic peaks.

This script takes preliminary peak summaries, annotates them with genomic features
(TSS, TES, Genebody, Intergenic), filters them based on various criteria (blacklist,
fold change, peak origin), and generates final summary reports in both CSV and
BED formats.
"""

import argparse
import pandas as pd
import numpy as np
import os

def load_ref_data(ref_file):
    """Loads and processes the reference gene data from a BED file."""
    print("Loading and processing reference gene data...")
    df = pd.read_csv(ref_file, sep='\t', header=None, 
                     names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    df = df[~df['chrom'].str.contains('_|chrUn')]
    df['gene_size'] = df['end'] - df['start']
    df['name'] = df['name'] + '|' + df['gene_size'].astype(str)

    # Create dictionaries for fast lookups
    tss_regions, tes_regions, gene_bodies = {}, {}, {}
    for _, row in df.iterrows():
        # Define TSS, TES, and Gene Body regions
        if row['strand'] == '+':
            tss_pos = row['start']
            tes_pos = row['end']
        else:
            tss_pos = row['end']
            tes_pos = row['start']
        
        # Using a wider range for TSS/TES matching
        for i in range(tss_pos - 1000, tss_pos + 1000, 50):
            key = f"{row['chrom']}:{i // 50 * 50}"
            tss_regions.setdefault(key, []).append(row['name'])
        for i in range(tes_pos - 1000, tes_pos + 1000, 50):
            key = f"{row['chrom']}:{i // 50 * 50}"
            tes_regions.setdefault(key, []).append(row['name'])
        for i in range(row['start'], row['end'], 100):
            key = f"{row['chrom']}:{i // 100 * 100}"
            gene_bodies.setdefault(key, []).append(row['name'])
            
    print(f"Reference data loaded: {len(df)} genes.")
    return {'TSS': tss_regions, 'TES': tes_regions, 'Genebody': gene_bodies}

def annotate_peak(peak_coord, ref_data):
    """Annotates a single peak based on its overlap with genomic features."""
    try:
        chrom, coords = peak_coord.split(':')
        start, end = map(int, coords.split('-'))
    except ValueError:
        return "Intergenic", "NA"

    # Check for overlap in order of priority: TSS, TES, Genebody
    for i in range(start, end, 50):
        key = f"{chrom}:{i // 50 * 50}"
        if key in ref_data['TSS']:
            return "TSS-proximal", ";".join(set(ref_data['TSS'][key]))
        if key in ref_data['TES']:
            return "TES", ";".join(set(ref_data['TES'][key]))
        if key in ref_data['Genebody']:
            return "Genebody", ";".join(set(ref_data['Genebody'][key]))

    return "Intergenic", "NA"

def process_mark(mark_name, input_file, ref_data, blacklist_df, peak_origin_df):
    """Processes and annotates data for a single histone mark."""
    print(f"\nProcessing mark: {mark_name}...")
    try:
        df = pd.read_csv(input_file)
    except FileNotFoundError:
        print(f"  Input file not found: {input_file}")
        return None

    # Annotate peaks
    annotations = df['Coordinates'].apply(lambda x: annotate_peak(x, ref_data))
    df['Annotation'] = [a[0] for a in annotations]
    df['Associated_Genes'] = [a[1] for a in annotations]

    # Calculate Fold Change
    young_cols = [c for c in df.columns if 'Young' in c or 'Sample_1' in c or 'Sample_2' in c]
    old_cols = [c for c in df.columns if 'Old' in c or 'Sample_3' in c or 'Sample_4' in c]
    df['Young_Mean'] = df[young_cols].mean(axis=1)
    df['Old_Mean'] = df[old_cols].mean(axis=1)
    df['log2FC'] = np.log2(df['Old_Mean'] / df['Young_Mean']).replace([np.inf, -np.inf], 0)

    # Add peak origin and blacklist info
    df = df.merge(peak_origin_df, on='Coordinates', how='left')
    df['Peak_Size'] = df['Coordinates'].apply(lambda x: int(x.split('-')[1]) - int(x.split('-')[0]))
    
    # A simple blacklist check (can be improved)
    df['Blacklisted'] = df['Coordinates'].isin(blacklist_df['Coordinates'])

    return df

def main():
    """Main function to run the summary report generation."""
    parser = argparse.ArgumentParser(description="Annotate and filter genomic peaks.")
    parser.add_argument("ref_gene_file", help="UCSC RefSeq gene file in BED format.")
    parser.add_argument("blacklist_file", help="Blacklist regions in BED format.")
    parser.add_argument("--input_dir", default="./preliminary_reports", help="Directory with preliminary reports.")
    parser.add_argument("--origin_dir", default="./peak_origins", help="Directory with peak origin files.")
    parser.add_argument("--output_dir", default="./summary_reports", help="Directory to save final reports.")
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Load reference data
    ref_data = load_ref_data(args.ref_gene_file)
    blacklist_df = pd.read_csv(args.blacklist_file, sep='\t', header=None, names=['chrom', 'start', 'end'])
    blacklist_df['Coordinates'] = blacklist_df['chrom'] + ':' + blacklist_df['start'].astype(str) + '-' + blacklist_df['end'].astype(str)

    marks = ["h3k4me3", "h3k27ac", "h3k27me3", "atac"]
    for mark in marks:
        input_file = os.path.join(args.input_dir, f"summary_{mark}.csv")
        origin_file = os.path.join(args.origin_dir, f"{mark}_peak_origins.bed") # Assuming this format
        
        try:
            peak_origin_df = pd.read_csv(origin_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'Coordinates', 'Origin'])
        except FileNotFoundError:
            print(f"Warning: Peak origin file not found for {mark}. Skipping origin info.")
            peak_origin_df = pd.DataFrame(columns=['Coordinates', 'Origin'])

        annotated_df = process_mark(mark, input_file, ref_data, blacklist_df, peak_origin_df[['Coordinates', 'Origin']])

        if annotated_df is not None:
            # --- Final Filtering and Output ---
            final_df = annotated_df[~annotated_df['Blacklisted']]
            
            # Separate by annotation type
            df_all = final_df
            df_tss = final_df[final_df['Annotation'] == 'TSS-proximal']
            
            # Filter by fold change
            df_fc = final_df[abs(final_df['log2FC']) > np.log2(1.5)]

            # Save files
            df_all.to_csv(os.path.join(args.output_dir, f"summary_{mark}_annotation.All.final.csv"), index=False)
            df_tss.to_csv(os.path.join(args.output_dir, f"summary_{mark}_annotation.TSS.final.csv"), index=False)
            df_fc.to_csv(os.path.join(args.output_dir, f"summary_{mark}_annotationFC1.5All.final.csv"), index=False)
            
            print(f"  Successfully generated final reports for {mark}.")

if __name__ == "__main__":
    main()
