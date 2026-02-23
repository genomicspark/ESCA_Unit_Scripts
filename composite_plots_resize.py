#!/usr/bin/env python3

"""
Generates composite plots from CDT files.

This script reads one or more CDT files from a directory, calculates the average
signal profile, and plots them on a single graph. It supports moving average
smoothing, normalization, and creating shaded plots.
"""

import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from operator import add

# Expanded color palette for more plots
COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", 
    "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", 
    "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", 
    "#dbdb8d", "#9edae5"
]

def moving_average(interval, window_size):
    """Calculates the moving average of an interval."""
    window = np.ones(int(window_size)) / float(window_size)
    return np.convolve(list(interval), list(window), "valid")

def process_file(infile):
    """Reads a CDT file and calculates the average signal profile."""
    print(f"Processing {infile}...")
    x_coords, y_sum = [], []
    num_lines = 0
    with open(infile, 'r') as f_in:
        for i, line in enumerate(f_in):
            if i == 0: # Header line
                x_coords = [float(x) for x in line.strip().split('\t')[2:]]
                y_sum = [0.0] * len(x_coords)
                continue
            
            try:
                tmplist = [float(x) for x in line.strip().split('\t')[2:]]
                y_sum = list(map(add, y_sum, tmplist))
                num_lines += 1
            except (ValueError, IndexError) as e:
                print(f"  Warning: Skipping malformed line {i+1} in {infile}: {e}", file=sys.stderr)

    if num_lines > 0:
        y_avg = [val / num_lines for val in y_sum]
    else:
        y_avg = y_sum

    return x_coords, y_avg, num_lines

def plot_graph(ax, x, y_avg, window_size, y_limit, shaded, normalize, label, color, num_lines):
    """Plots a single signal profile on the given axes."""
    if not x or not y_avg:
        print(f"  Skipping plot for {label} due to empty data.")
        return

    x_smooth = moving_average(x, window_size)
    y_smooth = moving_average(y_avg, window_size)

    if normalize and max(y_smooth) > 0:
        y_smooth = [val / max(y_smooth) for val in y_smooth]

    plot_label = f"{label} (n={num_lines})"
    ax.plot(x_smooth, y_smooth, color=color, label=plot_label, lw=2.5, zorder=2)
    if shaded:
        ax.fill_between(x_smooth, y_smooth, 0, facecolor=color, alpha=0.3)

    if y_limit is not None:
        ax.set_ylim(0, y_limit)
    else:
        current_max = ax.get_ylim()[1]
        new_max = max(current_max, max(y_smooth) * 1.2)
        ax.set_ylim(0, new_max)

    ax.set_xlim(min(x_smooth), max(x_smooth))

def main():
    """Main function to run the plotting script."""
    parser = argparse.ArgumentParser(description="Generate composite plots from CDT files.")
    parser.add_argument("input_directory", help="Path to directory containing *.cdt files.")
    parser.add_argument("-w", "--window_size", type=int, default=10, help="Window size for moving average. Default=10.")
    parser.add_argument("-y", "--y_limit", type=float, help="Set a fixed y-axis limit.")
    parser.add_argument("--shaded", action="store_true", help="Create shaded plots instead of line plots.")
    parser.add_argument("--normalize", action="store_true", help="Normalize each profile to its max value (0 to 1).")
    parser.add_argument("-o", "--output_file", help="Output file name (e.g., plot.svg). Default is composite_plot.svg in the input directory.")
    args = parser.parse_args()

    if not os.path.isdir(args.input_directory):
        parser.error(f"Input directory not found: {args.input_directory}")

    cdt_files = sorted([os.path.join(args.input_directory, f) for f in os.listdir(args.input_directory) if f.endswith(".cdt")])

    if not cdt_files:
        print("No .cdt files found in the specified directory.")
        return

    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    for i, cdt_file in enumerate(cdt_files):
        x, y_avg, num_lines = process_file(cdt_file)
        label = os.path.basename(cdt_file).replace('.cdt', '')
        color = COLORS[i % len(COLORS)]
        plot_graph(ax, x, y_avg, args.window_size, args.y_limit, args.shaded, args.normalize, label, color, num_lines)

    ax.set_title("Composite Signal Profile", fontsize=16)
    ax.set_xlabel("Genomic Position (bp)", fontsize=12)
    ax.set_ylabel("Average Signal", fontsize=12)
    ax.axhline(0, color='black', lw=0.8)
    ax.legend(loc='upper right')
    fig.tight_layout()

    if args.output_file:
        outfile = args.output_file
    else:
        outfile = os.path.join(args.input_directory, "composite_plot.svg")

    try:
        fig.savefig(outfile)
        print(f"\nPlot saved to {outfile}")
    except (IOError, PermissionError) as e:
        print(f"Error saving plot: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()
