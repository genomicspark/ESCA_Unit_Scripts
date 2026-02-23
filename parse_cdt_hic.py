#!/usr/bin/env python3

"""
Generates a CDT file for heatmaps around gene bodies from a BigWig file.

This script is similar to parse_cdt_bin.py but focuses on the entire gene body,
flanked by a region proportional to the gene's length. It is designed for
visualizing signals across gene bodies (e.g., for Hi-C data).
"""

import operator
import pyBigWig
import sys
import argparse
import math

def generate_cdt_hic(ref_gene_file, bigwig_file, output_file):
    """
    Generates a CDT file for gene body regions.

    Args:
        ref_gene_file (str): Path to the UCSC RefSeq gene file (BED format).
        bigwig_file (str): Path to the input BigWig file.
        output_file (str): Path to the output CDT file.
    """
    regions = {}
    region_metadata = {}

    with open(ref_gene_file, 'r') as f_ref:
        for line in f_ref:
            data = line.strip().split('\t')
            if len(data) < 6:
                continue

            chrom, start, end, strand = data[0], int(data[1]), int(data[2]), data[5]
            gene_length = end - start
            if gene_length <= 0:
                continue

            # Flanking region is 1/5th of the gene length
            heatmap_size = gene_length // 5
            region_start = start - heatmap_size
            region_end = end + heatmap_size
            total_region_size = region_end - region_start

            # Bin size aims for ~2800 bins, as in original script
            bin_size = total_region_size // 2800 if total_region_size > 2800 else 1

            the_key = f"{chrom}:{region_start}-{region_end}"
            regions[the_key] = [chrom, region_start, region_end, str(total_region_size), strand]
            region_metadata[the_key] = {'total_size': total_region_size, 'bin_size': bin_size}

    print(f"Loaded {len(regions)} regions to process.")

    try:
        bw = pyBigWig.open(bigwig_file)
    except RuntimeError as e:
        print(f"Error opening BigWig file {bigwig_file}: {e}", file=sys.stderr)
        return

    auc_values = {}
    binned_signals = {}

    for i, (the_key, region_info) in enumerate(regions.items()):
        if i % 1000 == 0:
            print(f"  ... processed {i} regions")

        chrom, start, end, _, strand = region_info
        bin_size = region_metadata[the_key]['bin_size']

        try:
            raw_values = bw.values(chrom, start, end)
            processed_values = []
            auc = 0.0
            for val in raw_values:
                if val is not None and not math.isnan(val):
                    # Normalization by bin_size from original script
                    norm_val = val / bin_size if bin_size > 0 else 0
                    processed_values.append(f"{norm_val:.2f}")
                    auc += norm_val
                else:
                    processed_values.append("0.00")

            if strand == '-':
                processed_values.reverse()

            # Binning logic
            binned_line = f"{the_key}\t1"
            if bin_size > 1:
                for j in range(0, len(processed_values), bin_size):
                    chunk = processed_values[j:j + bin_size]
                    bin_sum = sum(float(v) for v in chunk)
                    binned_line += f"\t{bin_sum / len(chunk):.4f}"
            else:
                binned_line += "\t" + "\t".join(processed_values)

            binned_signals[the_key] = binned_line
            auc_values[the_key] = auc

        except (RuntimeError, ZeroDivisionError) as e:
            print(f"  Warning: Could not process {the_key}: {e}", file=sys.stderr)
            continue

    # Sort by AUC and write to CDT file
    sorted_regions = sorted(auc_values.items(), key=operator.itemgetter(1), reverse=True)

    with open(output_file, 'w') as f_out:
        # Header from original script, fixed range
        header = '\t'.join(str(i) for i in range(-1399, 1401))
        f_out.write(f"Unique\tWeight\t{header}\n")
        for the_key, _ in sorted_regions:
            if the_key in binned_signals:
                f_out.write(binned_signals[the_key] + '\n')

    print(f"Finished writing {output_file}")
    bw.close()

def main():
    """Command-line argument parsing."""
    parser = argparse.ArgumentParser(description="Generate CDT files for heatmaps across gene bodies.")
    parser.add_argument("ref_gene_file", help="UCSC RefSeq gene file in BED format.")
    parser.add_argument("bigwig_file", help="Input BigWig file.")
    parser.add_argument("-o", "--output_file", required=True, help="Output CDT file path.")
    args = parser.parse_args()

    generate_cdt_hic(args.ref_gene_file, args.bigwig_file, args.output_file)

if __name__ == "__main__":
    main()
