#!/usr/bin/env python3

"""
Generates a preliminary summary report from parsed AUC peak files.

This script reads a list of AUC peak files, aggregates the results for different
histone marks (H3K4me3, H3K27me3, H3K27ac, ATAC-seq), and generates
a summary CSV and BED file for each mark.
"""

import argparse
import pandas as pd

def create_preliminary_report(file_list, input_dir, output_dir):
    """
    Aggregates AUC peak data and creates summary reports.

    Args:
        file_list (str): A file containing the list of AUC peak files to process.
        input_dir (str): The directory where AUC peak files are located.
        output_dir (str): The directory to save summary reports.
    """
    try:
        with open(file_list, 'r') as f:
            files = [line.strip() for line in f.readlines()]
    except FileNotFoundError:
        print(f"Error: File list '{file_list}' not found.")
        return

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Group files by histone mark
    file_groups = {
        "h3k4me3": [f for f in files if "h3k4me3" in f.lower()],
        "h3k27me3": [f for f in files if "h3k27me3" in f.lower()],
        "h3k27ac": [f for f in files if "h3k27ac" in f.lower()],
        "atac": [f for f in files if "atac" in f.lower()],
    }

    for mark, file_group in file_groups.items():
        if not file_group:
            print(f"No files found for mark: {mark}")
            continue

        print(f"Processing mark: {mark}...")
        all_peaks = {}

        for filename in file_group:
            filepath = os.path.join(input_dir, filename)
            try:
                df = pd.read_csv(filepath, sep='\t')
                # Assuming columns: Coordinates, AUC_per_bp
                for _, row in df.iterrows():
                    peak_coord = row['Coordinates']
                    auc_val = row['AUC_per_bp']
                    if peak_coord not in all_peaks:
                        all_peaks[peak_coord] = []
                    all_peaks[peak_coord].append(auc_val)
            except (FileNotFoundError, KeyError, pd.errors.ParserError) as e:
                print(f"  Warning: Could not process file {filepath}. Error: {e}")
                continue

        # Create summary DataFrame
        summary_data = []
        for peak, values in all_peaks.items():
            avg_auc = sum(values) / len(values) if values else 0
            summary_data.append([peak] + values + [avg_auc])
        
        num_samples = max(len(v) for v in all_peaks.values()) if all_peaks else 0
        columns = ['Coordinates'] + [f'Sample_{i+1}' for i in range(num_samples)] + ['Average_AUC']
        summary_df = pd.DataFrame(summary_data, columns=columns)
        summary_df = summary_df.sort_values(by='Average_AUC', ascending=False)

        # Write CSV and BED files
        csv_path = os.path.join(output_dir, f"summary_{mark}.csv")
        bed_path = os.path.join(output_dir, f"summary_{mark}.bed")
        summary_df.to_csv(csv_path, index=False)

        with open(bed_path, 'w') as f_bed:
            for _, row in summary_df.iterrows():
                try:
                    chrom, coords = row['Coordinates'].split(':')
                    start, end = coords.split('-')
                    f_bed.write(f"{chrom}\t{start}\t{end}\t{row['Coordinates']}\n")
                except ValueError:
                    print(f"  Warning: Malformed coordinate string '{row['Coordinates']}'. Skipping BED entry.")

        print(f"  Generated {csv_path} and {bed_path}")

def main():
    """
    Command-line argument parsing.
    """
    parser = argparse.ArgumentParser(description="Generate a preliminary summary report from AUC peak files.")
    parser.add_argument("file_list", help="A text file listing the AUC peak files to process.")
    parser.add_argument("-i", "--input_dir", default="./", help="Directory containing the AUC peak files.")
    parser.add_argument("-o", "--output_dir", default="./preliminary_reports", help="Directory to save the summary reports.")
    args = parser.parse_args()

    create_preliminary_report(args.file_list, args.input_dir, args.output_dir)

if __name__ == "__main__":
    main()
