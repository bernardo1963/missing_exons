# **Drosophila Y-Linked Exon Sequencing Bias Analysis**

This repository contains supplementary files for the study **"Strong sequencing bias in Nanopore and PacBio prevents assembly of Drosophila melanogaster Y-linked genes."** The study investigates sequencing biases affecting the assembly of Y-linked exons in *Drosophila melanogaster* using **Nanopore, PacBio, and Illumina** technologies.

## **Repository Structure**
The repository is organized into three main directories:

### **Fasta/** – Assembled Contigs and Reference Sequences  
This directory contains FASTA files with assembled "missing exons" and reference sequences used in the study. The assembled contigs were generated from Nanopore long reads using **Minimap/Miniasm** followed by **Racon polishing**. Reference sequences include Y-linked exons and satellite DNA.

### **Scripts/** – Data Processing and Analysis Scripts  
This directory includes **Python, AWK, and Shell scripts** used for read coverage computation, statistical tests, and visualization. Key scripts include `read_coverage_CDS_v4.py` for computing read coverage per base, `missingExon_stat_1jan2025.py` for statistical analysis, and `visualize_repeats_censor.py` for visualizing repetitive sequences.

### **Data/** – Processed Data for Figures and Tables  
This directory contains **tabulated data** used in figures and statistical tests presented in the manuscript. Each file corresponds to data from specific analyses, such as exon coverage profiles (`Fig4_data.txt`), statistical distributions (`Fig5_data.txt`), and read start positions (`Fig6_data.txt`).

## **Installation & Requirements**
The scripts require **Python 3+** and dependencies that can be installed using:

```bash
pip install numpy scipy matplotlib
