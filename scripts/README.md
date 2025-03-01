# Scripts for Data Processing, Analysis, and Figure Generation

This folder contains scripts used for data processing, statistical analysis, and figure generation in the study.

## Files

### Data Processing & Analysis
- **`missingExon_stat_1jan2025.py`**  - Statistical analysis of missing exons, including variance tests and the Anderson-Darling test.
- **`read_coverage_CDS_v4.py`** - Computes read coverage across coding sequences.
- **`AnDarl.c`** - C implementation of the Anderson-Darling test for uniformity (original version by Marsaglia).
- **`AnDarl_modified.c`** - Modified version of `AnDarl.c` that allows command-line parameter input for easier integration with Python scripts.

### Visualization Scripts
- **`box_violin_plot.py`**  - Generates box and violin plots for data visualization.
- **`visualize_repeats_censor.py`**  - Visualizes repetitive elements detected using `censor`.

### Figure Generation
- **`Fig5_graph.syc`**  - Script for plotting data related to Figure 5.
- **`figures_tables_scripts.sh`**  - Shell script for processing figures and tables.

### Read Coverage Processing
- **`fplot_reads.awk`**  - AWK script for plotting read distributions.

## Usage

Each script includes internal documentation or comments. Run scripts using **Python**, **AWK**, or **Shell** as appropriate.

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

Example output
```bash
Prob(A_100 <  0.9000)=0.586025560369086
```
