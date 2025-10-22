# Composite Epigenetic Profiling of Aged Hematopoietic Stem Cells
This repository contains the code and resources for the HSC Aging Epigenome Project (2025), which investigates dysregulated chromatin organization in aged hematopoietic stem cells (HSCs) leading to aberrant transcription. This work is conducted by Le Zong, Bongsoo Park, and Isabel Beerman at the Translational Gerontology Branch, National Institute on Aging (NIA), as part of the Epigenetics and Stem Cell Aging research program.

# Bioinformatics Pipeline: Bivalent Domain Detection
Developed by: Bongsoo Park (2025)

This pipeline processes histone modification ChIP-seq data to generate epigenome profiles and detect bivalent domains. 
The following Python scripts are included:

generate_consensus_peaks.py: Generates consensus peaks from ChIP-seq data.
parse_auc_peak.py: Parses peak areas under the curve (AUC) for analysis.
create_preliminary_report.py: Creates preliminary reports from processed data.
create_summary_report.py: Generates summary reports for peak calling results.
create_summary_report_final.py: Produces the final summary report.
research_bivalent_domain.py: Identifies and analyzes bivalent domains.

Input Requirements:  Reference genomic coordinates in BED format.  
ChIP-seq intensity retrieved using the pyBigWig library (e.g., parse_cdt_bin.py BigWigName, parse_cdt_hic.py BigWigName).

# Visualization Scripts
The following scripts facilitate visualization of epigenomic data:

Composite Plot Generation: Creates composite plots for visualizing epigenetic profiles.
CDT File Generation: Produces CDT files for heatmap and composite plot generation.

# Usage Instructions
Setup: Ensure all dependencies (e.g., pyBigWig, matplotlib) are installed.
Input Files: Provide genomic coordinates in BED format for peak calling and BigWig files for visualization.
Running Scripts: Execute the scripts in the order listed above for the respective pipelines. 

Example commands:
python parse_cdt_bin.py <BigWigName>
python parse_cdt_hic.py <BigWigName>

Output: Results include consensus peaks, reports, composite plots, and CDT files for further analysis.

# Single-Cell Analysis Pipeline
Developed by: Bongsoo Park (2022)
This pipeline supports single-cell RNA-seq analysis for breast cancer research using Seurat v3. 
The key script is:NIH-NIA-single-cell-breast-scTransform.R: Implements single-cell analysis with scTransform normalization.

# Contributions and Feedback
We welcome contributions, bug reports, and suggestions to improve this pipeline. Please submit issues or pull requests via GitHub.Contact:
Translational Gerontology Branch, National Institute on Aging

