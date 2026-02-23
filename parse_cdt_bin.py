#!/usr/bin/env python3

"""
Generates a CDT file for heatmaps based on signal from a BigWig file.

This script takes a reference gene file (BED format) and a BigWig file.
It extracts signal values around the TSS and TES of each gene, bins them,
and outputs a sorted CDT file ready for visualization.
"""

import operator
import pyBigWig
import sys
import argparse
import math

def generate_cdt(ref_gene_file, bigwig_file, output_prefix, heatmap_size=5000, bin_size=10):
    """
    Generates CDT files for TSS and TES regions.

    Args:
        ref_gene_file (str): Path to the UCSC RefSeq gene file (BED format).
        bigwig_file (str): Path to the input BigWig file.
        output_prefix (str): Prefix for the output CDT files.
        heatmap_size (int): Flanking region size around TSS/TES.
        bin_size (int): Size of bins for averaging signal.
    """
    tss_regions = {}
    tes_regions = {}

    with open(ref_gene_file, 'r') as f_ref:
        for line in f_ref:
            data = line.strip().split('\t')
            if len(data) < 6:
                continue
            
            chrom, strand, gene_name = data[0], data[5], data[3]

            if strand == "+":
                tss_start, tss_end = int(data[1]) - heatmap_size, int(data[1]) + heatmap_size
                tes_start, tes_end = int(data[2]) - heatmap_size, int(data[2]) + heatmap_size
            else:
                tss_start, tss_end = int(data[2]) - heatmap_size, int(data[2]) + heatmap_size
                tes_start, tes_end = int(data[1]) - heatmap_size, int(data[1]) + heatmap_size

            tss_regions[gene_name] = [chrom, tss_start, tss_end, str(tss_end - tss_start), strand]
            tes_regions[gene_name] = [chrom, tes_start, tes_end, str(tes_end - tes_start), strand]

    print(f"Loaded {len(tss_regions)} TSS and {len(tes_regions)} TES regions.")

    try:
        bw = pyBigWig.open(bigwig_file)
    except RuntimeError as e:
        print(f"Error opening BigWig file {bigwig_file}: {e}", file=sys.stderr)
        return

    for region_type, regions in [('tss', tss_regions), ('tes', tes_regions)]:
        output_file = f"{output_prefix}_{region_type}.cdt"
        print(f"Processing {region_type.upper()} regions...")

        auc_values = {}
        binned_signals = {}

        for i, (gene_name, region_info) in enumerate(regions.items()):
            if i % 1000 == 0:
                print(f"  ... processed {i} genes")

            chrom, start, end, _, strand = region_info
            try:
                raw_values = bw.values(chrom, start, end)
                processed_values = []
                auc = 0.0
                for val in raw_values:
                    if val is not None and not math.isnan(val):
                        # Normalization factor from original script
                        norm_val = val / 40.0
                        processed_values.append(f"{norm_val:.2f}")
                        auc += norm_val
                    else:
                        processed_values.append("0.00")

                if strand == '-':
                    processed_values.reverse()

                # Binning logic
                binned_line = f"{gene_name}\t1"
                if bin_size > 1:
                    for j in range(0, len(processed_values), bin_size):
                        chunk = processed_values[j:j + bin_size]
                        bin_sum = sum(float(v) for v in chunk)
                        binned_line += f"\t{bin_sum / len(chunk):.4f}"
                else:
                    binned_line += "\t" + "\t".join(processed_values)
                
                binned_signals[gene_name] = binned_line
                auc_values[gene_name] = auc

            except RuntimeError as e:
                print(f"  Warning: Could not process {gene_name} ({chrom}:{start}-{end}): {e}", file=sys.stderr)
                continue

        # Sort by AUC and write to CDT file
        sorted_genes = sorted(auc_values.items(), key=operator.itemgetter(1), reverse=True)
        
        with open(output_file, 'w') as f_out:
            header = '\t'.join(str(i) for i in range(-heatmap_size + bin_size, heatmap_size, bin_size))
            f_out.write(f"Unique\tWeight\t{header}\n")
            for gene_name, _ in sorted_genes:
                if gene_name in binned_signals:
                    f_out.write(binned_signals[gene_name] + '\n')
        
        print(f"Finished writing {output_file}")

    bw.close()

def main():
    """Command-line argument parsing."""
    parser = argparse.ArgumentParser(description="Generate CDT files for heatmaps from BigWig data.")
    parser.add_argument("ref_gene_file", help="UCSC RefSeq gene file in BED format.")
    parser.add_argument("bigwig_file", help="Input BigWig file.")
    parser.add_argument("-o", "--output_prefix", required=True, help="Prefix for output CDT files (e.g., my_experiment).")
    parser.add_argument("--heatmap_size", type=int, default=5000, help="Size of the flanking region around TSS/TES.")
    parser.add_argument("--bin_size", type=int, default=10, help="Size of bins for averaging signal.")
    args = parser.parse_args()

    generate_cdt(args.ref_gene_file, args.bigwig_file, args.output_prefix, args.heatmap_size, args.bin_size)

if __name__ == "__main__":
    main()
