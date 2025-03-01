# Scripts for Data Processing, Analysis, and Figure Generation

This folder contains scripts used for data processing, statistical analysis, and figure generation in the study.

## Files

### Data Processing & Analysis
- **`missingExon_stat_1jan2025.py`**  
- **`read_coverage_CDS_v4.py`**
- **`AnDarl.c`** - C implementation of the Anderson-Darling test for uniformity (original version by Marsaglia).
- **`AnDarl_modified.c`** - Modified version of `AnDarl.c` that allows command-line parameter input for easier integration with Python scripts.

### Visualization Scripts
- **`box_violin_plot.py`**  
- **`visualize_repeats_censor.py`**  

### Figure Generation
- **`Fig5_graph.syc`**  
- **`figures_tables_scripts.sh`**  

### Read Coverage Processing
- **`fplot_reads.awk`**  

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

