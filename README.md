# Composite Epigenetic Profiling of Aged Hematopoietic Stem Cells
This repository contains the code and resources for the HSC Aging Epigenome Project (2025), which investigates dysregulated chromatin organization in aged hematopoietic stem cells (HSCs) leading to aberrant transcription. This work is conducted by Le Zong, Bongsoo Park, and Isabel Beerman at the Translational Gerontology Branch, National Institute on Aging (NIA), as part of the Epigenetics and Stem Cell Aging research program.

# Bioinformatics Pipeline: Bivalent Domain Detection
Developed by: Bongsoo Park (2025)

This pipeline processes histone modification ChIP-seq data to generate epigenome profiles and detect bivalent domains. 
The following Python scripts are included:

1. generate_consensus_peaks.py: Generates consensus peaks from ChIP-seq data.
2. parse_auc_peak.py: Parses peak areas under the curve (AUC, aka: peak height) for downstream analysis.
3. create_preliminary_report.py: Creates preliminary reports from processed data.
4. create_summary_report.py: Generates summary reports for peak calling results.
5. create_summary_report_final.py: Produces the final summary report.
6. research_bivalent_domain.py: Identifies and analyzes bivalent domains.

Input Requirements:  Reference genomic coordinates in BED format.  
ChIP-seq intensity retrieved using the pyBigWig library (e.g., parse_cdt_bin.py BigWigName, parse_cdt_hic.py BigWigName).

# Visualization Scripts
The following scripts facilitate visualization of epigenomic data:

Composite Plot Generation: Creates composite plots for visualizing epigenetic profiles.
CDT File Generation: Produces CDT files for generating heatmaps and composite plots.

# Usage Instructions
Setup: Ensure all dependencies (e.g., pyBigWig, matplotlib) are installed. tested in Python 3.9.6

Input Files: Provide genomic coordinates in BED format for peak calling and BigWig files for visualization.
Running Scripts: Execute the scripts in the order listed above for the respective pipelines. 
Parse_cdt scripts will generate a CDT file (tab-delimited signal matrix), then you can create a composite plot with it. 

Example commands:
python parse_cdt_bin.py YHSC_H3K4me3.chr19.bigWig
python parse_cdt_hic.py YHSC_H3K4me3.chr19.bigWig
python composite_plots_resize.py -w 10 <cdt_folder_name>

Output: Results include consensus peaks, reports, composite plots, and CDT files for further analysis.

# Single-Cell Analysis Pipeline
Developed by: Bongsoo Park (2022)
This pipeline supports single-cell RNA-seq analysis for breast cancer research using Seurat v3. 
The key script is:NIH-NIA-single-cell-breast-scTransform.R: Implements single-cell analysis with scTransform normalization.

# Contributions and Feedback
We welcome contributions, bug reports, and suggestions to improve this pipeline. Please submit issues or pull requests via GitHub.Contact:
Translational Gerontology Branch, National Institute on Aging

