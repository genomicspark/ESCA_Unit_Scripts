#!/usr/bin/env python3

"""
Calculates the Area Under the Curve (AUC) for peaks from a BigWig file.

This script reads a BED file of genomic regions (peaks), and for each region,
calculates the sum of signal values (AUC) from a given BigWig file.
"""

import pyBigWig
import sys
import argparse
import math

def calculate_auc(bed_file, bigwig_file, output_file):
    """
    Calculates AUC for peaks in a BED file from a BigWig file.

    Args:
        bed_file (str): Path to the input BED file.
        bigwig_file (str): Path to the input BigWig file.
        output_file (str): Path to the output file.
    """
    try:
        bw = pyBigWig.open(bigwig_file)
    except RuntimeError as e:
        print(f"Error opening BigWig file {bigwig_file}: {e}", file=sys.stderr)
        return

    with open(bed_file, 'r') as f_bed, open(output_file, 'w') as f_out:
        f_out.write("GeneSymbol\tRefSeq\tPeakSize\tCoordinates\tAUC\tAUC_per_bp\n")
        
        for line in f_bed:
            data = line.strip().split('\t')
            if len(data) < 4:
                continue

            chrom, start, end = data[0], int(data[1]), int(data[2])
            peak_name = data[3]
            peak_size = end - start

            try:
                values = bw.values(chrom, start, end)
                auc = 0.0
                for val in values:
                    if val is not None and not math.isnan(val):
                        auc += val
                
                auc_per_bp = auc / peak_size if peak_size > 0 else 0.0

                # Attempt to parse gene symbol and refseq from peak name
                try:
                    gene_symbol, refseq = peak_name.split('|')
                except ValueError:
                    gene_symbol, refseq = peak_name, '.'

                f_out.write(f"{gene_symbol}\t{refseq}\t{peak_size}\t{chrom}:{start}-{end}\t{auc:.4f}\t{auc_per_bp:.4f}\n")

            except RuntimeError as e:
                print(f"Error processing region {chrom}:{start}-{end}: {e}", file=sys.stderr)
                continue

    bw.close()
    print(f"AUC calculation complete. Results written to {output_file}")

def main():
    """Command-line argument parsing."""
    parser = argparse.ArgumentParser(description="Calculate AUC for peaks from a BigWig file.")
    parser.add_argument("bed_file", help="Input BED file with peak coordinates.")
    parser.add_argument("bigwig_file", help="Input BigWig file with signal data.")
    parser.add_argument("-o", "--output_file", required=True, help="Output file for AUC results.")
    args = parser.parse_args()

    calculate_auc(args.bed_file, args.bigwig_file, args.output_file)

if __name__ == "__main__":
    main()
