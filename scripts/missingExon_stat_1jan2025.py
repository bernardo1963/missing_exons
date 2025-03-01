#!/usr/bin/env python3
#
# missingExon_stat_26dez2024.py  based on missingExon_stat.py v1.0  Bernardo/Fabiana dez2024

# new in v  1jan2025: implements Implements Anderson-Darling test with a composite null hypothesis (sum of uniform distributions); Corrected the Kolmogorov-Smirnov function to `composite_cdf += stats.uniform.cdf(curr_start, 0, size)` (previously, it incorrectly used `composite_cdf += stats.uniform.cdf(curr_start, 1, size)`). 
# new in v 30dez2024: Implements Kolmogorov-Smirnov test with a composite null hypothesis (sum of uniform distributions).
# new in v 29dez2024: Implements MAD (Median Absolute Deviation) as a robust alternative to variance.
# new in v 28dez2024: Default offset is now 0 (previously hardcoded as 20 before Dec 27); Implemented `--m8mod` option.
# new in v 27dez2024: Added argparse support. Implemented control for `max_readsize`.
# new in v 26dez2024: Implements start variance analysis in the R reads.
# new in v 25dez2024: Implements Fisher's method for variance instead of summing SS and df. Reason: Some genes require a two-tailed variance test, and individual genes may fall into opposite tails of the F-test. Their effects could cancel each other out, concealing significant variance patterns.

# This script reads an input file with data on forward (F) and reverse (R) read counts for different genes, performs a binomial test on each gene to see if the observed F vs. R ratio deviates significantly from 50:50, and then combines the resulting p-values from all genes using Fisher’s combined probability test
#
# Input file format: target$ type$ gi read$ size strand$ start
#
# 1   6 CCY_exon2_30bp_ONT.txt
# 0   4 kl3_region1_30bp_ONT.txt
# 3   4 kl5_exon13_30bp_ONT.txt
#

import sys
import math
from scipy.stats import binomtest, chi2
from collections import defaultdict # collections.defaultdict: https://stackoverflow.com/questions/29348345/declaring-a-multi-dimensional-dictionary-in-python
# new_dic = defaultdict(dict)
# new_dic[1][2] = 5
from statistics import variance
from statistics import fmean

import scipy.stats as stats
import argparse
import numpy as np
import re
from scipy.stats import kstwo
import subprocess


def get_parameters():
    parser = argparse.ArgumentParser(description='processes a table derived from m8 blast output  to make statistics of non-random redas (eg, F/R usage, start withing read)', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--data', type=str, default="",  help='custom abular  format for this program. Fields:  target$ type$ gi read$ size strand$ start')
    parser.add_argument('--m8mod', type=str, default="", help='provide m8mod file. It is m8 with fileds 9 and 10 swuitched ion +/- alns, and complemented w/  additional fields:  strand$ type$  read$ size')
    parser.add_argument('--max_read_size', type=int, default=0,  help='maximum read size allowed for start comparisons (will not affect  F/R usage). Larger reads will be ignored')
    parser.add_argument('--max_start', type=int, default=0,  help='maximum start allowed for start comparisons (will not affect  F/R usage). Reads with start larger thn limit will be ignored.')    
    parser.add_argument('--n_simulations', type=int, default=1000,  help='number of start point simulations used in function simulation_var_unif_dist')
    parser.add_argument('--offset', type=int, default=0,  help='discout appllied to read size attemptingvto compensate that a F macth cannot start at trhe last position. It is a bad idea')
    '''
    parser.add_argument('--systat_output_file', type=str, default="CDS_coverage_systat.txt", help='final result  file, ready for SYSTAT.')
    parser.add_argument('--identity_cutoff', type=float, default=0,  help='minimum % identity that will be considered for coverage. Hits  with identity below this threshold will be ignored')   
    parser.add_argument('--verbose', type=int, default=1,  help='log verbosity')
    parser.add_argument('--rosetta', default="", type=str, help='optional file containing gene - gene _type table ()to build the  dictionary .')
    '''
    args = parser.parse_args()
    if (args.m8mod != "" and args.data != "") or  (args.m8mod == "" and args.data == "") :
        print("ERROR: either --m8mod OR --data parameters must be set. Setting both or none is not allowed" , flush=True)
        exit(1)
    
    return args


def Fisher_method(p_values_list,replace_zero=0):
    # Implements Fisher’s method to combine P-values
    # k = number of genes/p-values
    # X^2 = -2 * sum(ln(p_i))
    # df = 2k, where k = number of p-values
    # p_value_combined = 1 - chi2.cdf(X^2, df=2k)
    if replace_zero != 0:
        for i in range(len(p_values_list)):
            if p_values_list[i] == 0:
                p_values_list[i] = replace_zero
    
    k = len(p_values_list)
    
    if k > 0:
        chi_square_stat = -2 * sum(math.log(p) for p in p_values_list)
        df = 2 * k
        combined_p = 1 - chi2.cdf(chi_square_stat, df)
    else:
        # If k=0 (no data), just set defaults
        chi_square_stat = float('nan')
        df = 0
        combined_p = float('nan')
    return(chi_square_stat,df,combined_p)
       


def F_R_usage_stats(FR_dict):
    # the F / R  usage  statistics
    results = [] 
    for gene_type2 in gene_type_list:
        p_values = []
        for gene in gene_type_dict:
            if (gene_type_dict[gene] != gene_type2):
                continue   # skips genes that belong to other gene types
            F_reads = FR_dict[gene,"F"]
            R_reads = FR_dict[gene,"R"]
            n_reads = F_reads + R_reads
            # STEP 1. Perform binomial test for each gene    
            # null hypothesis: p=0.5               
            # binomtest returns a result object (test_res), which includes the pvalue.
            test_res = binomtest(k=F_reads, n=n_reads, p=0.5)
            # extracting the p-value
            p_val = test_res.pvalue
            # appends the p-value to p_values
            p_values.append(p_val)
            results.append((gene, F_reads, R_reads, p_val))                
          
        # STEP 2. Fisher’s combined probability test
        (chi_square_stat,df,combined_p) = Fisher_method(p_values)

        # Print  individual and combined results of a given gene_type
        print("gene_type:" , gene_type2)
        # print("Gene\t\tF\tR\tBinomial_pvalue")
        print(f'{"gene":20s}{"F":>6s}{"R":>6s}{"P":>12s}')
        for gene, F, R, p_val in results:
            if (gene_type_dict[gene] == gene_type2):
                print(f'{gene:20s}{F:>6d}{R:>6d}{p_val:>12.4f}')  

                # print(f"{gene}\t{F}\t{R}\t{pval}")           
        print("\nFisher combined probability test results:")
        print(f'Chi-square statistic: {chi_square_stat:.2f}  Degrees of freedom: {df:6d}   Combined p-value: {combined_p:2.3f}\n')
        print("")

def read_start_stats():
    print("read start - variance analytical:")
    results_dict = defaultdict(dict)
    for gene_type2 in gene_type_list:
        curr_gene_list = []
        for gene in gene_type_dict:
            if gene_type_dict[gene] == gene_type2:
                curr_gene_list.append(gene)
        p_values = []
        for gene in curr_gene_list:                         
            results_dict[gene,"F"] = (0, 0, 0, 0, 0) 
            results_dict[gene,"R"] = (0, 0, 0, 0, 0)           
            for strand in ["F","R"]:
            
                start_list = read_start_dict[gene,strand]
                size_list = read_size_dict[gene,strand]
                
                start_list_trimmed = []
                size_list_trimmed = []             
                for i in range(0, len(size_list), 1):
                    if size_list[i] <= 10000:
                        size_list_trimmed.append(size_list[i])
                        start_list_trimmed.append(start_list[i])
                if len(start_list_trimmed) < 2:              
                    continue                       
                var_start = variance(start_list_trimmed)
                n_reads = len(start_list_trimmed)
                df_start = n_reads -1
                max_size = max(size_list_trimmed)
                # print(gene, "size_list_trimmed:", size_list_trimmed)
                # print(gene, "start_list_trimmed:", start_list_trimmed)
                
                # expected_unif_variance = ((max_size_F -1)**2) / 12  # https://en.wikipedia.org/wiki/Continuous_uniform_distribution
                # Below, I removed 20 because this is approximately the minimum alignment length required to be detected in BLAST. 
                # A read that terminates at the first base of the CDS will not be detected. 
                # expected_unif_variance = ((max_size -1 -20)**2) / 12  # https://en.wikipedia.org/wiki/Continuous_uniform_distribution
                expected_unif_variance = ((10000 -1 -20)**2) / 12  # https://en.wikipedia.org/wiki/Continuous_uniform_distribution
                
                # F two-tailed test https://www.geeksforgeeks.org/how-to-perform-an-f-test-in-python/
                f_value = var_start  / expected_unif_variance
                if var_start  >= expected_unif_variance:
                    p_val = ( 1 - stats.f.cdf(f_value , df_start , 100000) ) *2  # I used 100000 to emmulate infinite df  I think this is equivalent to a  chi-square
                else:
                    p_val = (stats.f.cdf(f_value , df_start, 100000) ) *2           
                p_values.append(p_val)
                results_dict[gene,strand] = (n_reads, var_start, expected_unif_variance, f_value, p_val) 
        
        # Print  individual and combined results of a given gene_type
        print("gene_type:" , gene_type2)
        print(f'{"gene":20s}{"n_reads_F":>10s}{"var_start_F":>15s}{"expected_var":>15s}{"F":>10s}{"P":>8s}{"n_reads_R":>15s}{"var_start_R":>15s}{"expected_var":>15s}{"F":>10s}{"P":>8s}')
        
        for gene in curr_gene_list:          
            (n_reads, var_start,expected_unif_variance,f_value, p_val) = results_dict[gene,"F"]
            if n_reads > 0:
                print(f'{gene:20s}{n_reads:>10d}{var_start:>15.1f}{expected_unif_variance:>15.1f}{f_value:>10.2f}{p_val:>8.4f}', end='     ')
            else:
                print(f'{gene:20s}{n_reads:>10d}{"-":>15s}{"-":>15s}{"-":>10s}{"-":>8s}', end='     ')

            (n_reads, var_start,expected_unif_variance,f_value, p_val) = results_dict[gene,"R"]
            print(f'{n_reads:>10d}{var_start:>15.1f}{expected_unif_variance:>15.1f}{f_value:>10.2f}{p_val:>8.4f}')
        # now the combined p_vales:
        (chi_square_stat,df,combined_p) = Fisher_method(p_values)
        print(f'{"Fisher method:     chi square=":s}{chi_square_stat:.2f}{"  df=":s}{df:d}{"    combined P=":s}{combined_p:.4f}')  
        print("")


def Kolmogorov_Smirnov():  # Kolmogorov-Smirnov test with composit null hypothesis (sum of uniform distribuitions)
    print("read start - KS test:")
    results_dict = defaultdict(dict)
    for gene_type2 in gene_type_list:
        curr_gene_list = []
        for gene in gene_type_dict:
            if gene_type_dict[gene] == gene_type2:
                curr_gene_list.append(gene)
        p_values = []
        for gene in curr_gene_list:                         
            results_dict[gene,"F"] = (0, 0, 0, 0, 0) 
            results_dict[gene,"R"] = (0, 0, 0, 0, 0)           
            for strand in ["F","R"]:
                size_list = read_size_dict[gene,strand]
                start_list = read_start_dict[gene,strand]
                if len(start_list) < 2:              
                    continue                      
                start_list_sorted = sorted(start_list) 
                n_reads = len(start_list)                
                KS_diff_list = []
                for i in range( n_reads):
                    curr_start = start_list_sorted[i]
                    composite_cdf = 0
                    for size in size_list:
                        composite_cdf += stats.uniform.cdf(curr_start,0,size)   
                    KS_diff = abs(i - composite_cdf) / n_reads   
                    KS_diff_list.append(KS_diff) 
                KS_stats = max(KS_diff_list)
                # print("KS_stats:", KS_stats , "n_reads:",n_reads, flush=True)
                KS_pvalue = 1 - (kstwo.cdf(KS_stats,n_reads))
                p_values.append(KS_pvalue)
                results_dict[gene,strand] = (n_reads, KS_stats, KS_pvalue) 
        
        # Print  individual and combined results of a given gene_type
        print("gene_type:" , gene_type2)
        print(f'{"gene":20s}{"n_reads_F":>10s}{"KS_stats":>10s}{"P":>8s}{"n_reads_R":>15s}{"KS_stats":>10s}{"P":>8s}')       
        for gene in curr_gene_list:
            # print(gene,"F:", results_dict[gene,"F"], flush=True) 
            # print(gene,"R:", results_dict[gene,"R"], flush=True)             
            (n_reads, KS_stats, p_val) = results_dict[gene,"F"]
            if n_reads > 0:
                print(f'{gene:20s}{n_reads:>10d}{KS_stats:>10.2f}{p_val:>8.4f}', end='     ')
            else:
                print(f'{gene:20s}{n_reads:>10d}{"-":>10s}{"-":>8s}', end='     ')

            (n_reads, KS_stats, p_val) = results_dict[gene,"R"]
            if n_reads > 0:            
                print(f'{n_reads:>10d}{KS_stats:>10.2f}{p_val:>8.4f}')
            else:
                print(f'{n_reads:>10d}{"-":>10s}{"-":>8s}')
                
        # now the combined p_vales:
        # print("p_values:", p_values)
        (chi_square_stat,df,combined_p) = Fisher_method(p_values,0.0001)
        print(f'{"Fisher method:     chi square=":s}{chi_square_stat:.2f}{"  df=":s}{df:d}{"    combined P=":s}{combined_p:.4f}')  
        print("")



def Anderson_Darling_2():  # This is the standard Anderson-Darling test, with p value obtained with Marsaglia and Marsaglia C-code
    # https://www.spcforexcel.com/knowledge/basic-statistics/anderson-darling-test-for-normality/
    # from https://www.6sigma.us/six-sigma-in-focus/anderson-darling-normality-test/
    # A^2 = -n – (1/n) * Sum[(2i – 1) * (ln(Φ(x(i))) + ln(1 – Φ(x(n+1-i)))]
    # I will follow the above formula and the Wikipedia formula (instead of the formula based on Excel )
    # will also replace composite_cdf += stats.uniform.cdf(curr_start,0,size) / n_reads  by composite_cdf += (curr_start/size) / n_reads
    
    print("read start - standard Anderson-Darling test:")
    results_dict = defaultdict(dict)
    for gene_type2 in gene_type_list:
        curr_gene_list = []
        for gene in gene_type_dict:
            if gene_type_dict[gene] == gene_type2:
                curr_gene_list.append(gene)
        p_values = []
        for gene in curr_gene_list:                         
            results_dict[gene,"F"] = (0, 0, 0, 0, 0) 
            results_dict[gene,"R"] = (0, 0, 0, 0, 0)           
            for strand in ["F","R"]:
                size_list = read_size_dict[gene,strand]
                start_list = read_start_dict[gene,strand]
                if len(start_list) < 2:              
                    continue                      
                start_list_sorted = sorted(start_list) 
                n_reads = len(start_list)                
                
                S = 0
                composite_cdf_list = []
                complementary_composite_cdf_list = []
                for i in range( n_reads):
                    curr_start = start_list_sorted[i]
                    composite_cdf = 0
                    for size in size_list:
                        # composite_cdf += (curr_start/size) / n_reads     # replaces stats.uniform.cdf(curr_start,0,size) / n_reads
                        composite_cdf +=  min((curr_start/size),1)  / n_reads
                        # composite_cdf += stats.uniform.cdf(curr_start,0,size) / n_reads
                        if  stats.uniform.cdf(curr_start,0,size) != min((curr_start/size),1):
                            print("cdf discrepancy:" , curr_start, size,   stats.uniform.cdf(curr_start,0,size) ,min((curr_start/size),1) ) 
                    
                    composite_cdf_list.append(composite_cdf) 

                
                # converting the 0-based composite_cdf_list to 1-based :
                composite_cdf_list = ["dummy"] + composite_cdf_list
                for idx in range(1, n_reads+1):                
                    try:
                        S += ((2*idx -1)/n_reads)*(  math.log(composite_cdf_list[idx]) + math.log( 1 -  composite_cdf_list[n_reads +1 -idx] ) )  
                        # S += (2*idx -1)*(  math.log(composite_cdf_list[i])  +  math.log(complementary_composite_cdf_list[i]) )
                    except:
                        print("idx:",idx,"  composite_cdf_list:", composite_cdf_list , flush=True)
                        print(gene, strand, "idx:",idx, "  composite_cdf_list[idx]:",composite_cdf_list[idx], "     composite_cdf_list[n_reads +1  -idx]:",composite_cdf_list[n_reads +1 -idx] )
                        print("start_list_sorted:", start_list_sorted)
                        raise
                # AD_stats = -n_reads -S/n_reads
                AD_stats = -n_reads -S

                cmd = "AnDarl_modified " + str(n_reads) + " " + str(AD_stats) 
                raw_AD = (subprocess.check_output(cmd,shell=True).decode("utf-8")).split("=")
                # print(raw_AD, flush=True)
                AD_pvalue = 1 - float(raw_AD[1])
                
                p_values.append(AD_pvalue)
                results_dict[gene,strand] = (n_reads, AD_stats, AD_pvalue) 
        
        # Print  individual and combined results of a given gene_type
        print("gene_type:" , gene_type2)
        print(f'{"gene":20s}{"n_reads_F":>10s}{"AD_stats":>10s}{"P":>8s}{"n_reads_R":>15s}{"AD_stats":>10s}{"P":>8s}')       
        for gene in curr_gene_list:
            # print(gene,"F:", results_dict[gene,"F"], flush=True) 
            # print(gene,"R:", results_dict[gene,"R"], flush=True)             
            (n_reads, AD_stats, p_val) = results_dict[gene,"F"]
            if n_reads > 0:
                print(f'{gene:20s}{n_reads:>10d}{AD_stats:>10.2f}{p_val:>8.4f}', end='     ')
            else:
                print(f'{gene:20s}{n_reads:>10d}{"-":>10s}{"-":>8s}', end='     ')

            (n_reads, AD_stats, p_val) = results_dict[gene,"R"]
            if n_reads > 0:            
                print(f'{n_reads:>10d}{AD_stats:>10.2f}{p_val:>8.4f}')
            else:
                print(f'{n_reads:>10d}{"-":>10s}{"-":>8s}')
                
        # now the combined p_vales:
        # print("p_values:", p_values)
        (chi_square_stat,df,combined_p) = Fisher_method(p_values,0.0001)
        print(f'{"Fisher method:     chi square=":s}{chi_square_stat:.2f}{"  df=":s}{df:d}{"    combined P=":s}{combined_p:.4f}')  
        print("")



def Anderson_Darling_2_old():  # This is the standard Anderson-Darling test, with p-values obtained using the Marsaglia and Marsaglia C code.
    # https://www.spcforexcel.com/knowledge/basic-statistics/anderson-darling-test-for-normality/
    # from https://www.6sigma.us/six-sigma-in-focus/anderson-darling-normality-test/
    # A^2 = -n – (1/n) * Sum[(2i – 1) * (ln(Φ(x(i))) + ln(1 – Φ(x(n+1-i)))]
    
    print("read start - standard Anderson-Darling test:")
    results_dict = defaultdict(dict)
    for gene_type2 in gene_type_list:
        curr_gene_list = []
        for gene in gene_type_dict:
            if gene_type_dict[gene] == gene_type2:
                curr_gene_list.append(gene)
        p_values = []
        for gene in curr_gene_list:                         
            results_dict[gene,"F"] = (0, 0, 0, 0, 0) 
            results_dict[gene,"R"] = (0, 0, 0, 0, 0)           
            for strand in ["F","R"]:
                size_list = read_size_dict[gene,strand]
                start_list = read_start_dict[gene,strand]
                if len(start_list) < 2:              
                    continue                      
                start_list_sorted = sorted(start_list) 
                n_reads = len(start_list)                
                
                S = 0
                composite_cdf_list = []
                complementary_composite_cdf_list = []
                for i in range( n_reads):
                    curr_start = start_list_sorted[i]
                    composite_cdf = 0
                    for size in size_list:
                        composite_cdf += stats.uniform.cdf(curr_start,0,size) / n_reads   
                    composite_cdf_list.append(composite_cdf)
                    complementary_composite_cdf_list.append(1 -composite_cdf ) 
                # sorting complementary_composite_cdf_list
                complementary_composite_cdf_list.sort()
                for i in range( n_reads):                
                    idx = i+1 
                    try:
                        S += (2*idx -1)*(  math.log(composite_cdf_list[i])  +  math.log(complementary_composite_cdf_list[i]) )
                    except:
                        print(gene, strand, i, "composite_cdf_list[i]:",composite_cdf_list[i], "     complementary_composite_cdf_list[i]:",complementary_composite_cdf_list[i] )
                        print("start_list_sorted:", start_list_sorted)
                        raise
                    # S += (2*idx -1)*(  math.log(composite_cdf_list[i])  +  math.log(complementary_composite_cdf_list[i]) )
                AD_stats = -n_reads -S/n_reads

                cmd = "AnDarl_modified " + str(n_reads) + " " + str(AD_stats) 
                raw_AD = (subprocess.check_output(cmd,shell=True).decode("utf-8")).split("=")
                # print(raw_AD, flush=True)
                AD_pvalue = 1 - float(raw_AD[1])
                
                p_values.append(AD_pvalue)
                results_dict[gene,strand] = (n_reads, AD_stats, AD_pvalue) 
        
        # Print  individual and combined results of a given gene_type
        print("gene_type:" , gene_type2)
        print(f'{"gene":20s}{"n_reads_F":>10s}{"AD_stats":>10s}{"P":>8s}{"n_reads_R":>15s}{"AD_stats":>10s}{"P":>8s}')       
        for gene in curr_gene_list:
            # print(gene,"F:", results_dict[gene,"F"], flush=True) 
            # print(gene,"R:", results_dict[gene,"R"], flush=True)             
            (n_reads, AD_stats, p_val) = results_dict[gene,"F"]
            if n_reads > 0:
                print(f'{gene:20s}{n_reads:>10d}{AD_stats:>10.2f}{p_val:>8.4f}', end='     ')
            else:
                print(f'{gene:20s}{n_reads:>10d}{"-":>10s}{"-":>8s}', end='     ')

            (n_reads, AD_stats, p_val) = results_dict[gene,"R"]
            if n_reads > 0:            
                print(f'{n_reads:>10d}{AD_stats:>10.2f}{p_val:>8.4f}')
            else:
                print(f'{n_reads:>10d}{"-":>10s}{"-":>8s}')
                
        # now the combined p_vales:
        # print("p_values:", p_values)
        (chi_square_stat,df,combined_p) = Fisher_method(p_values,0.0001)
        print(f'{"Fisher method:     chi square=":s}{chi_square_stat:.2f}{"  df=":s}{df:d}{"    combined P=":s}{combined_p:.4f}')  
        print("")




def mad_read_start_sampling_simulations():  # analogous to def read_start_stats_simulation():
    print("read start - mad resampling:")
    results_dict = defaultdict(dict)
    for gene_type2 in gene_type_list:
        curr_gene_list = []
        for gene in gene_type_dict:
            if gene_type_dict[gene] == gene_type2:
                curr_gene_list.append(gene)
        p_values = []
        for gene in curr_gene_list:                         
            results_dict[gene,"F"] = (0, 0, 0, 0, 0) 
            results_dict[gene,"R"] = (0, 0, 0, 0, 0)           
            for strand in ["F","R"]:
            
                start_list = read_start_dict[gene,strand]
                size_list = read_size_dict[gene,strand]
                if len(start_list) < 2:              
                    continue                                                   
                mad_start = stats.median_abs_deviation(start_list)
                n_reads = len(start_list)
                # now the P-value by sampling simulations:
                (expected_unif_mad,two_tailed_prob,rank) = mad_sampling_simulations_engine(mad_start,size_list,args.offset,args.n_simulations)
                
                mad_ratio = mad_start  / expected_unif_mad
                p_values.append(two_tailed_prob)
                # results_dict[gene,strand].append((n_reads, var_start, expected_unif_variance, f_value, p_val))
                results_dict[gene,strand] = (n_reads, mad_start, expected_unif_mad,mad_ratio, two_tailed_prob) 
                # print("gene:",gene,"   cum_SS_start_F:", cum_SS_start_F, "   cum_df_start_F:",cum_df_start_F)
        
        # Print  individual and combined results of a given gene_type
        print("gene_type:" , gene_type2)
        # print("gene\t\tn_reads_F\tvar_start_F\texpected_unif_variance\tf_value\tp_val")
        print(f'{"gene":20s}{"n_reads_F":>10s}{"mad_start_F":>15s}{"expected_mad":>15s}{"mad_ratio":>10s}{"P":>8s}{"n_reads_R":>15s}{"mad_start_R":>15s}{"expected_mad":>15s}{"mad_ratio":>10s}{"P":>8s}')

        
        for gene in curr_gene_list:
            # print(gene,"F:", results_dict[gene,"F"], flush=True) 
            # print(gene,"R:", results_dict[gene,"R"], flush=True) 
            
            (n_reads, mad_start,expected_unif_mad,f_value, p_val) = results_dict[gene,"F"]
            if n_reads > 0:
                print(f'{gene:20s}{n_reads:>10d}{mad_start:>15.1f}{expected_unif_mad:>15.1f}{f_value:>10.2f}{p_val:>8.4f}', end='     ')
            else:
                print(f'{gene:20s}{n_reads:>10d}{"-":>15s}{"-":>15s}{"-":>10s}{"-":>8s}', end='     ')

            (n_reads, mad_start,expected_unif_mad,f_value, p_val) = results_dict[gene,"R"]
            if n_reads > 0:            
                print(f'{n_reads:>10d}{mad_start:>15.1f}{expected_unif_mad:>15.1f}{f_value:>10.2f}{p_val:>8.4f}')
            else:
                print(f'{n_reads:>10d}{"-":>15s}{"-":>15s}{"-":>10s}{"-":>8s}')
                
        # now the combined p_vales:
        # print("p_values:", p_values)
        (chi_square_stat,df,combined_p) = Fisher_method(p_values,0.0001)
        print(f'{"Fisher method:     chi square=":s}{chi_square_stat:.2f}{"  df=":s}{df:d}{"    combined P=":s}{combined_p:.4f}')  
        print("")

def mad_sampling_simulations_engine(mad_obs,size_list_tmp,size_offset,n_sim):    # analogous to def simulation_var_unif_dist
    size_list = []
    mad_sim_list=[]
    mad_obs_bigger = 0
    mad_obs_smaller = 0
    
    for size in size_list_tmp:
        size_list.append(size - size_offset) 
        
    for i in range(0, n_sim, 1):
        start_sim_list = []
        for read_size in size_list:
            start_sim = stats.randint.rvs(1, read_size, size=1)[0]
            start_sim_list.append(start_sim)
        mad_sim = stats.median_abs_deviation(start_sim_list)
        mad_sim_list.append(mad_sim)
        # print(i,"start_sim_list:",start_sim_list)  
    
    
    avg_mad_sim = fmean(mad_sim_list)
    # print("mad_sim_list:", mad_sim_list)
    # print("mad_obs:",mad_obs, "    avg_mad_sim:",avg_mad_sim) 
    for  mad_sim in mad_sim_list:
        if mad_obs >= mad_sim:
            mad_obs_bigger += 1
        else:
            mad_obs_smaller += 1    
    
    # print("mad_obs_smaller:",mad_obs_smaller, "mad_obs_bigger:",mad_obs_bigger)  
    rank = min(mad_obs_bigger,mad_obs_smaller)
    two_tailed_probablility = 2*rank / n_sim
    
    
    return(avg_mad_sim,two_tailed_probablility,rank)




def read_start_stats_simulation():
    print("read start - variance simulation:")
    results_dict = defaultdict(dict)
    for gene_type2 in gene_type_list:
        curr_gene_list = []
        for gene in gene_type_dict:
            if gene_type_dict[gene] == gene_type2:
                curr_gene_list.append(gene)
        p_values = []
        for gene in curr_gene_list:                         
            results_dict[gene,"F"] = (0, 0, 0, 0, 0) 
            results_dict[gene,"R"] = (0, 0, 0, 0, 0)           
            for strand in ["F","R"]:
            
                start_list = read_start_dict[gene,strand]
                size_list = read_size_dict[gene,strand]
                if len(start_list) < 2:              
                    continue                                                   
                var_start = variance(start_list)
                n_reads = len(start_list)
                # now the P-value by simulation:
                (expected_unif_variance,two_tailed_prob,rank) = simulation_var_unif_dist(var_start,size_list,args.offset,args.n_simulations)
                f_value = var_start  / expected_unif_variance
                p_values.append(two_tailed_prob)
                # results_dict[gene,strand].append((n_reads, var_start, expected_unif_variance, f_value, p_val))
                results_dict[gene,strand] = (n_reads, var_start, expected_unif_variance, f_value, two_tailed_prob) 
                # print("gene:",gene,"   cum_SS_start_F:", cum_SS_start_F, "   cum_df_start_F:",cum_df_start_F)
        
        # Print  individual and combined results of a given gene_type
        print("gene_type:" , gene_type2)
        # print("gene\t\tn_reads_F\tvar_start_F\texpected_unif_variance\tf_value\tp_val")
        print(f'{"gene":20s}{"n_reads_F":>10s}{"var_start_F":>15s}{"expected_var":>15s}{"F":>10s}{"P":>8s}{"n_reads_R":>15s}{"var_start_R":>15s}{"expected_var":>15s}{"F":>10s}{"P":>8s}')

        
        for gene in curr_gene_list:
            # print(gene,"F:", results_dict[gene,"F"], flush=True) 
            # print(gene,"R:", results_dict[gene,"R"], flush=True) 
            
            (n_reads, var_start,expected_unif_variance,f_value, p_val) = results_dict[gene,"F"]
            if n_reads > 0:
                print(f'{gene:20s}{n_reads:>10d}{var_start:>15.1f}{expected_unif_variance:>15.1f}{f_value:>10.2f}{p_val:>8.4f}', end='     ')
            else:
                print(f'{gene:20s}{n_reads:>10d}{"-":>15s}{"-":>15s}{"-":>10s}{"-":>8s}', end='     ')

            (n_reads, var_start,expected_unif_variance,f_value, p_val) = results_dict[gene,"R"]
            if n_reads > 0:            
                print(f'{n_reads:>10d}{var_start:>15.1f}{expected_unif_variance:>15.1f}{f_value:>10.2f}{p_val:>8.4f}')
            else:
                print(f'{n_reads:>10d}{"-":>15s}{"-":>15s}{"-":>10s}{"-":>8s}')
                
        # now the combined p_vales:
        # print("p_values:", p_values)
        (chi_square_stat,df,combined_p) = Fisher_method(p_values,0.0001)
        print(f'{"Fisher method:     chi square=":s}{chi_square_stat:.2f}{"  df=":s}{df:d}{"    combined P=":s}{combined_p:.4f}')  
        print("")




def read_start_Anderson_Darling():
    print("read start - Anderson_Darling:")
    results_dict = defaultdict(dict)
    for gene_type2 in gene_type_list:
        curr_gene_list = []
        for gene in gene_type_dict:
            if gene_type_dict[gene] == gene_type2:
                curr_gene_list.append(gene)
        p_values = []
        for gene in curr_gene_list:                         
            results_dict[gene,"F"] = (0, 0, 0, 0, 0) 
            results_dict[gene,"R"] = (0, 0, 0, 0, 0)           
            for strand in ["F","R"]:
            
                start_list = read_start_dict[gene,strand]
                size_list = read_size_dict[gene,strand]
                if len(start_list) < 2:              
                    continue                                                   
                
                
                n_reads = len(start_list)
                # now the P-value by Anderso-Darling + simulation:
                (AD_stats, AD_pvalue, critical_values) = Anderson_Darling_ksamples(start_list,size_list,args.offset,args.n_simulations)
                p_values.append(AD_pvalue)
                # results_dict[gene,strand].append((n_reads, var_start, expected_unif_variance, f_value, p_val))
                results_dict[gene,strand] = (n_reads, AD_stats, AD_pvalue) 
        
        # Print  individual and combined results of a given gene_type
        print("gene_type:" , gene_type2)
        print(f'{"gene":20s}{"n_reads_F":>10s}{"AD_stats":>10s}{"P":>8s}{"n_reads_R":>15s}{"AD_stats":>10s}{"P":>8s}')

        
        for gene in curr_gene_list:
            # print(gene,"F:", results_dict[gene,"F"], flush=True) 
            # print(gene,"R:", results_dict[gene,"R"], flush=True) 
            
            (n_reads, AD_stats, p_val) = results_dict[gene,"F"]
            if n_reads > 0:
                print(f'{gene:20s}{n_reads:>10d}{AD_stats:>10.2f}{p_val:>8.4f}', end='     ')
            else:
                print(f'{gene:20s}{n_reads:>10d}{"-":>10s}{"-":>8s}', end='     ')

            (n_reads, AD_stats, p_val) = results_dict[gene,"R"]
            if n_reads > 0:            
                print(f'{n_reads:>10d}{AD_stats:>10.2f}{p_val:>8.4f}')
            else:
                print(f'{n_reads:>10d}{"-":>10s}{"-":>8s}')
                
        # now the combined p_vales:
        # print("p_values:", p_values)
        (chi_square_stat,df,combined_p) = Fisher_method(p_values,0.0001)
        print(f'{"Fisher method:     chi square=":s}{chi_square_stat:.2f}{"  df=":s}{df:d}{"    combined P=":s}{combined_p:.4f}')  
        print("")






def Anderson_Darling_ksamples(start_obs_list,size_list_tmp,size_offset,n_sim):
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.anderson_ksamp.html#scipy.stats.anderson_ksamp
    size_list = []
    
    for size in size_list_tmp:
        size_list.append(size - size_offset) 
     
    for i in range(0, n_sim, 1):
        start_sim_list = []
        for read_size in size_list:
            start_sim_list_temp = stats.randint.rvs(1, read_size, size=n_sim)
            start_sim_list = [*start_sim_list, *start_sim_list_temp]        # [*l1, *l2]   unpack both iterables in a list literal
        # print(i,"start_sim_list:",start_sim_list)  
    
    print("variance in " , n_sim , "replicas:", variance(start_sim_list))
    start_sim_array = np.asarray(start_sim_list)
    start_obs_array = np.asarray(start_obs_list)
    
    res = stats.anderson_ksamp([start_sim_array,start_obs_array])
    
    return(res.statistic, res.pvalue, res.critical_values)

   



def simulation_var_unif_dist(var_obs,size_list_tmp,size_offset,n_sim):
    size_list = []
    var_sim_list=[]
    var_obs_bigger = 0
    var_obs_smaller = 0
    
    for size in size_list_tmp:
        size_list.append(size - size_offset) 
        
    for i in range(0, n_sim, 1):
        start_sim_list = []
        for read_size in size_list:
            start_sim = stats.randint.rvs(1, read_size, size=1)[0]
            start_sim_list.append(start_sim)
        var_sim = variance(start_sim_list)
        var_sim_list.append(var_sim)
        # print(i,"start_sim_list:",start_sim_list)  
    
    
    avg_var_sim = fmean(var_sim_list)
    # print("var_sim_list:", var_sim_list)
    # print("var_obs:",var_obs, "    avg_var_sim:",avg_var_sim) 
    for  var_sim in var_sim_list:
        if var_obs >= var_sim:
            var_obs_bigger += 1
        else:
            var_obs_smaller += 1    
    
    # print("var_obs_smaller:",var_obs_smaller, "var_obs_bigger:",var_obs_bigger)  
    rank = min(var_obs_bigger,var_obs_smaller)
    two_tailed_probablility = 2*rank / n_sim
    
    
    return(avg_var_sim,two_tailed_probablility,rank)
    

def main():
    if len(sys.argv) < 2:
        print("Usage: python missingExon_stat.py <input_file> ")
        sys.exit(1)
        
    
    # output_file = sys.argv[2]

    #store the p-values from the binomial tests (one per gene)
    # p_values = []
    # results = []
    global read_size_dict, gene_type_dict, gene_type_list, read_start_dict
    read_size_dict = defaultdict(dict) # no need to declare index  in advance
    read_start_dict = defaultdict(dict)  
    read_FR_dict = defaultdict(dict)  
    gene_type_dict = defaultdict(dict)
    gene_type_set = set()
    gene_type_list = []
    # read_size_dict = defaultdict(dict)    # read_size_dict["CCY_exon2"]["F"] = [list of read sizes]
    # I will use the tuple as a key: read_size_dict["CCS","F"]

    
    if args.data == "" and args.m8mod != "":
        input_type = "m8mod"
        input_file = args.m8mod
    if args.data != "" and args.m8mod == "":
        input_type = "data"
        input_file = args.data        
    
    
    # Read input file
    with open(input_file, 'r') as f:
        for line in f:
            # print(line)
            line = line.strip()
            if not line:
                continue
            if re.search("\$", line): # header
                continue           
            
            fields = line.split() #splits the line by whitespace
            # if fields[0] == "target$":  
                # continue
            
            if input_type == "data":
                if len(fields) != 7:   # target$ type$ gi read$ size strand$ start
                    print("Input file format error. --data expected 7 columns: target$ type$ gi read$ size strand$ start")
                    sys.exit(1)
                # Parse columns
                gene = fields[0]
                gene_type = fields[1]
                read_size = int(fields[4])
                read_strand = fields[5]
                read_start = int(fields[6])                
                
            if input_type == "m8mod":
                if len(fields) < 16:   #query$ gi  identity  aln_len mismatches  gap_op  q_start  q_end  s_start  s_end   e_value bit_score  strand$ type$       read$ size
                    print("Input file format error. --m8mod expected 16 or 17 columns: query$   gi  identity  aln_len mismatches  gap_op  q_start  q_end  s_start  s_end   e_value bit_score  strand$ type$ read$ size  cum_count optional]")
                    sys.exit(1)                 
                # Parse columns
                gene = fields[0]
                gene_type = fields[13]
                read_size = int(fields[15])
                read_strand = fields[12]
                read_start = int(fields[8])             
            
            #populate dictionaries            
            if not(gene in gene_type_dict):
               read_size_dict[gene,"F"] = []   # empty list
               read_size_dict[gene,"R"] = []    
               read_start_dict[gene,"F"] = []    
               read_start_dict[gene,"R"] = []  
               read_FR_dict[gene,"F"] = 0 
               read_FR_dict[gene,"R"] = 0 

            read_FR_dict[gene,read_strand] += 1
            if args.max_read_size != 0 and read_size >= args.max_read_size   :
                continue
            if args.max_start != 0 and read_start >= args.max_start   :
                continue            
            
            read_size_dict[gene,read_strand].append(read_size)
            read_start_dict[gene,read_strand].append(read_start)
            
            # if read_start <= 10000:               
                # read_start_dict[gene,read_strand].append(read_start)
            
            gene_type_dict[gene] = gene_type
            gene_type_set.add(gene_type)
            
        gene_type_list = list(gene_type_set) 
        gene_type_list.sort()
        F_R_usage_stats(read_FR_dict)
        # read_start_stats()
        Kolmogorov_Smirnov()
        Anderson_Darling_2()
        # mad_read_start_sampling_simulations()
        # read_start_stats_simulation()
        # read_start_Anderson_Darling()
        
if __name__ == "__main__":
    args = get_parameters()
    main()

