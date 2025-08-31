# Scripts for Data Processing, Analysis, and Figure Generation

This folder contains scripts used for data processing, statistical analysis, and figure generation in the study.

## Files

### Data Processing & Analysis
- [`missingExon_stat_1jan2025.py`](missingExon_stat_1jan2025.py) - Statistical analysis of missing exons, including binomial tests, variance analysis, Kolmogorov-Smirnov tests, and Anderson-Darling tests for uniformity. 
- [`missingExon_binomial_stat_v1.py`](missingExon_binomial_stat_v1.py) - Enhanced version that accepts custom expected probabilities (instead of fixed 50:50) and performs binomial tests with Fisher's combined probability analysis.
- [`read_coverage_CDS_v6.py`](read_coverage_CDS_v6.py) - Computes read coverage across coding sequences with visualization capabilities, exon overlay support, and smoothing filters.
- [`find_tandem_repeats_v2.py`](find_tandem_repeats_v2.py) - Detects and summarizes tandem repeats in DNA sequences, reporting their motifs, copy numbers, and positions
- [`AnDarl.c`](AnDarl.c) - C implementation of the Anderson-Darling test for uniformity (original version by Marsaglia).
- [`AnDarl_modified.c`](AnDarl_modified.c) - Modified version of `AnDarl.c` that allows command-line parameter input for easier integration with Python scripts.

### Visualization Scripts
- [`box_violin_plot.py`](box_violin_plot.py) - Generates box and violin plots for read size distributions across different genes using matplotlib and seaborn.
- [`visualize_repeats_censor5.py`](visualize_repeats_censor4.py) - Visualizes repetitive elements detected using `censor` with support for custom gene/repeat classification, exon overlays, and detailed composition statistics.

### Figure Generation
- [`Fig5_graph.syc`](Fig5_graph.syc) - Script for plotting data related to Figure 5.

### File Processing & Utilities
- [`fasta_size.awk`](fasta_size.awk) - AWK script to count the length of sequences in multiple FASTA files and produce a summary table.
- [`fplot_reads.awk`](fplot_reads.awk) - AWK script for plotting read distributions and generating gnuplot scripts to visualize reads as arrows.
- [`process_nonB_gfa.awk`](process_nonB_gfa.awk) - Processes non-B DNA structure analysis output files and generates summary tables or raw data output.

### Assembly Quality Assessment
- [`Ycompleteness_v3.sh`](Ycompleteness_v3.sh) - Tests completeness of genes (usually Y-linked) in assemblies using BLAST searches with customizable parameters and produces tabular output for comparison.

## Complete Analysis Workflow

### Reproducible Research Pipeline
- **[`Supplemental_Code.sh`](Supplemental_Code.sh)** - Complete computational workflow to reproduce all figures, tables, and results from the manuscript. This shell script provides exact command lines, parameters, and data processing steps used in the study "Strong bias in long-read sequencing prevents assembly of Drosophila melanogaster Y-linked genes."
  
## Usage

Each script includes internal documentation or comments. Run scripts using **Python**, **AWK**, **Shell**, or **C** as appropriate.

### Python Scripts
Most Python scripts accept command-line arguments. Use `--help` to see available options:
```bash
python script_name.py --help
```

### AWK Scripts
Run AWK scripts directly on input files:
```bash
./script_name.awk input_file
# or
gawk -f script_name.awk input_file
```

### Shell Scripts
Make executable and run with appropriate parameters:
```bash
chmod +x script_name.sh
./script_name.sh [options] input_files
```

## Anderson-Darling Test Implementation

The Anderson-Darling (AD) test was implemented in Python following these resources:

- [Six Sigma AD Test Explanation](https://www.6sigma.us/six-sigma-in-focus/anderson-darling-normality-test/)
- [SPC for Excel AD Test Guide](https://www.spcforexcel.com/knowledge/basic-statistics/anderson-darling-test-for-normality)

To compute p-values, we utilized the algorithm provided in Marsaglia and Marsaglia (2004). The original implementation (`AnDarl.c`) required interactive input, making it inconvenient for automation. We modified the program (`AnDarl_modified.c`) to accept command-line arguments for seamless integration with Python.

## Compilation & Usage

To compile and use the original version:

```bash
gcc AnDarl.c -lm -o AnDarl
./AnDarl  # Interactive mode, requires manual input of n and z
```

Example interactive session:
```bash
# Enter n and z: 100 0.9
Prob(A_100 <  0.9000)=  0.58603
```

To compile and use the modified version:
```bash
gcc AnDarl_modified.c -lm -o AnDarl_modified
./AnDarl_modified 100 0.9  # Example with n=100, z=0.9
```

Example output:
```bash
Prob(A_100 <  0.9000)=0.586025560369086
```

## Dependencies

### Python Scripts
- **Required packages**: `pandas`, `matplotlib`, `seaborn`, `scipy`, `numpy`, `Bio` (Biopython)
- **Python version**: 3.7+ (some scripts specify 3.9+)

### External Tools
- **BLAST**: Required for assembly completeness analysis ([`Ycompleteness_v3.sh`](Ycompleteness_v3.sh))
- **gnuplot**: Required for read visualization ([`fplot_reads.awk`](fplot_reads.awk))
- **gcc**: Required for compiling Anderson-Darling test programs


