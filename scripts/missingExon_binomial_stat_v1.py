#!/usr/bin/python3.9 -u
#
# missingExon_binomial_stat_v1.py  based on  missingExon_stat.py v1.0    Fabiana dez2024 / Bernardo
# new in v2: accepts column 4 w/ expected probability (insetad of using always 50:50); changed order of fields in input file
# usage: missingExon_binomial_stat_v1.py    lowcov_start_seq.data 
#
# This script reads an input file with data on count1 and count2 , and the expected probability, performs a binomial test on each gene to see if the observed count1 count2 ratio deviates significantly from the expected probability. Then combines the resulting p-values from all genes using Fisher’s combined probability test
#
# Input file format: gene name (column 1), count1 (column 2), count2  (column 3) and expected probability (in %; column 4)
# CCY            0 13      69.2
# kl5_exon_3-9   0 12      37.5
# kl5_exon13     0 12      84.3
# kl3_region1    0 16      85.0
# ORY_exon_1-2   1 15      86.5
# PprY           0 17      18.7



import sys
import math
from scipy.stats import binomtest, chi2

def main():
    if len(sys.argv) != 2:
        print("Usage: python missingExon_stat.py <input_file> ")
        sys.exit(1)
        
    input_file = sys.argv[1]
    #store the p-values from the binomial tests (one per gene)
    p_values = []
    results = []

    # Read input file
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split() #splits the line by whitespace
            if len(parts) != 4:
                print("Input file format error. Expected 4 columns: Gene count1 count2 prob")
                sys.exit(1)
            
            # Parse columns
            gene = parts[0]
            # Convert the  columns 2 and 3 to integers (counts)
            count1 = int(parts[1])
            count2 = int(parts[2])
            # Convert  columns 4  to float (expected_frequency)
            expected_frequency = float(parts[3])    
            # STEP 1. Perform binomial test
            # null hypothesis: p= expected_frequency
            # Number of trials is the sum of forward and reverse reads (n=count1+count2 trials)
            # k=count1: number of "successes" (forward reads)
            # n=n: total trials
            n = count1 + count2
            observed_frequency = 100*count1/n

            # binomtest returns a result object (test_res), which includes the pvalue.
            test_res = binomtest(k=count1, n=n, p=expected_frequency/100)
            # extracting the p-value
            p_val = test_res.pvalue
            # appends the p-value to p_values
            p_values.append(p_val)
            results.append((gene, count1, count2,observed_frequency,expected_frequency, p_val))

    # STEP 2. Fisher’s combined probability test
    # Combine all p-values
    # k = number of genes/p-values
    # X^2 = -2 * sum(ln(p_i))
    # df = 2k, where k = number of p-values
    # p_value_combined = 1 - chi2.cdf(X^2, df=2k)
    k = len(p_values)
    if k > 0:
        chi_square_stat = -2 * sum(math.log(p) for p in p_values)
        df = 2 * k
        combined_p = 1 - chi2.cdf(chi_square_stat, df)
    else:
        # If k=0 (no data), just set defaults
        chi_square_stat = float('nan')
        df = 0
        combined_p = float('nan')

    # print results to stdout
    
    print(f'{"gene":15s}{"count1":>10s}{"count2":>10s}{"obs_freq":>10s}{"exp_freq":>10s}{"P":>10s}')
    for (gene, count1, count2,observed_frequency,expected_frequency, p_val) in results:
         print(f'{gene:15s}{count1:>10d}{count2:>10d}{observed_frequency:>10.1f}{expected_frequency:>10.1f}{p_val:>10.4f}')
    
    print(f'\n{"Fisher method:     chi square=":s}{chi_square_stat:.2f}{"  df=":s}{df:d}{"    combined P=":s}{combined_p:.4f}')
    # print("\nFisher's combined probability test results:\n")
    # print(f"Chi-square statistic: {chi_square_stat}\n")
    # print(f"Degrees of freedom: {df}\n")
    # print(f"Combined p-value: {combined_p}\n")


if __name__ == "__main__":
    main()

