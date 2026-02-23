#!/usr/bin/env python3

"""
Generates a CDT file for heatmaps centered on peak midpoints from a BigWig file.

This script reads a BED file of peaks (e.g., Differentially Accessible Regions),
finds the midpoint of each peak, and extracts signal from a BigWig file in a
fixed-size window around that midpoint. The result is a CDT file suitable for
creating composite plots or heatmaps.
"""

import operator
import pyBigWig
import sys
import argparse
import math

def generate_cdt_peak(peak_file, bigwig_file, output_file, heatmap_size=5000, bin_size=10):
    """
    Generates a CDT file centered on peak midpoints.

    Args:
        peak_file (str): Path to the input BED file of peaks.
        bigwig_file (str): Path to the input BigWig file.
        output_file (str): Path to the output CDT file.
        heatmap_size (int): Flanking region size around the peak midpoint.
        bin_size (int): Size of bins for averaging signal.
    """
    regions = {}
    with open(peak_file, 'r') as f_peaks:
        for line in f_peaks:
            data = line.strip().split('\t')
            if len(data) < 3:
                continue

            chrom, start, end = data[0], int(data[1]), int(data[2])
            mid_point = start + (end - start) // 2
            region_start = mid_point - heatmap_size
            region_end = mid_point + heatmap_size

            # Use original peak name or construct one if not present
            peak_name = data[3] if len(data) > 3 else f"{chrom}:{start}-{end}"

            regions[peak_name] = [chrom, region_start, region_end, str(region_end - region_start)]

    print(f"Loaded {len(regions)} regions from {peak_file}.")

    try:
        bw = pyBigWig.open(bigwig_file)
    except RuntimeError as e:
        print(f"Error opening BigWig file {bigwig_file}: {e}", file=sys.stderr)
        return

    auc_values = {}
    binned_signals = {}

    for i, (peak_name, region_info) in enumerate(regions.items()):
        if i % 1000 == 0:
            print(f"  ... processed {i} peaks")

        chrom, start, end, _ = region_info

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

            # Binning logic
            binned_line = f"{peak_name}\t1"
            if bin_size > 1:
                for j in range(0, len(processed_values), bin_size):
                    chunk = processed_values[j:j + bin_size]
                    bin_sum = sum(float(v) for v in chunk)
                    # Original script had a bug here, not dividing by bin_size
                    binned_line += f"\t{bin_sum / len(chunk):.4f}"
            else:
                binned_line += "\t" + "\t".join(processed_values)

            binned_signals[peak_name] = binned_line
            auc_values[peak_name] = auc

        except RuntimeError as e:
            print(f"  Warning: Could not process {peak_name} ({chrom}:{start}-{end}): {e}", file=sys.stderr)
            continue

    # Sort by AUC and write to CDT file
    sorted_peaks = sorted(auc_values.items(), key=operator.itemgetter(1), reverse=True)

    with open(output_file, 'w') as f_out:
        header = '\t'.join(str(i) for i in range(-heatmap_size + bin_size, heatmap_size, bin_size))
        f_out.write(f"Unique\tWeight\t{header}\n")
        for peak_name, _ in sorted_peaks:
            if peak_name in binned_signals:
                f_out.write(binned_signals[peak_name] + '\n')

    print(f"Finished writing {output_file}")
    bw.close()

def main():
    """Command-line argument parsing."""
    parser = argparse.ArgumentParser(description="Generate CDT files for heatmaps centered on peak midpoints.")
    parser.add_argument("peak_file", help="Input BED file of peaks.")
    parser.add_argument("bigwig_file", help="Input BigWig file.")
    parser.add_argument("-o", "--output_file", required=True, help="Output CDT file path.")
    parser.add_argument("--heatmap_size", type=int, default=5000, help="Size of the flanking region around the peak midpoint.")
    parser.add_argument("--bin_size", type=int, default=10, help="Size of bins for averaging signal.")
    args = parser.parse_args()

    generate_cdt_peak(args.peak_file, args.bigwig_file, args.output_file, args.heatmap_size, args.bin_size)

if __name__ == "__main__":
    main()
