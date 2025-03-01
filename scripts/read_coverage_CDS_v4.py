#!/usr/local/bin/python3.9 -u
# read_coverage_CDS_v3.py   based on digital_kennison_v9.py    Bernardo   moj_Y and nanopore_missing projects   v. 5feb2025
# new in v2: serious bug correction;  implement option --m8_file
# -F "m S"
# new in v4: Generate a depth vs. position plot for each gene using matplotlib; Apply a rolling median filter to smooth out spikes
import os
import subprocess
import logging
from Bio import SeqIO
import argparse
import sys
import re
import statistics
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def get_parameters():
    global regex_info    
    # logger.info("Starting argument parsing")
    parser = argparse.ArgumentParser(description='detects satellite DNA sequences (and any other sequence) in reads, using regex (requires 100% match)',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--fasta_file', type=str, required=True, help='fasta sequence file that will be analysed. Usually a CDS multi-fasta ')
    parser.add_argument('--db',  default="", type=str, help='formatted blast db (usualy reads) ')  
    parser.add_argument('--blast_output_file', type=str, default="blast_output.m8", help='blast output, in m8 format.')
    parser.add_argument('--systat_output_file', type=str, default="CDS_coverage_systat.txt", help='final result  file, ready for SYSTAT.')
    parser.add_argument('--identity_cutoff', type=float, default=0,  help='minimum % identity that will be considered for coverage. Hits  with identity below this threshold will be ignored')   
    parser.add_argument('--ncbi_blast_parameters', type=str, default=' -e 0.001  -F "m D" -b 1000000 -v 1000000 -a 100 ')
    parser.add_argument('--wu_blast_parameters', type=str, default=' E=0.001  -wordmask dust  V=1000000 B=1000000   hspmax=0  cpus=100  ')  
    parser.add_argument('--blast_pgm', type=str, choices=["ncbi","wu"], default="ncbi",help='blast program that will be used. ncbi blast / wu-blast')    
    parser.add_argument('--m8_file', type=str, default="", help='provide m8 file, and do not run blast.')
    parser.add_argument('--verbose', type=int, default=1,  help='log verbosity')
    parser.add_argument('--rosetta', default="", type=str, help='optional file containing gene info that will be attached t the output files. Typically, chrom location.')
    parser.add_argument("--window_size", type=int, default=5, help="Size of the moving median window")
    parser.add_argument('--graph_name_suffix', type=str, default='', help='suffix to be added to automatic name (eg, porechop)')
    parser.add_argument("--figsize_W", type=float, default=12, help="fig size WIDTH (in inches)")
    parser.add_argument("--figsize_H", type=float, default=6, help="fig size HEIGHT (in inches)")
    # parser.add_argument("--Xmin", type=int, default=0, help="Xmin value")
    # parser.add_argument("--Xmax", type=int, default=0, help="Xmax value")
    parser.add_argument("--Ymin", type=float, default=0, help="Xmin value")
    parser.add_argument("--Ymax", type=float, default=-1, help="Ymax value")    
    parser.add_argument('--Xlabel', type=str, default='CDS position (bp)', help='label of the X-axis')
    parser.add_argument('--Ylabel', type=str, default='Depth', help='label of the Y-axis')
 

    args = parser.parse_args()
    if (args.m8_file != "" and args.db != "") or  (args.m8_file == "" and args.db == "") :
        print("ERROR: either --m8_file OR --db parameters must be set. Setting both or none is not allowed" , flush=True)
        exit(1)
    
    if args.m8_file != "" :
        args.blast_output_file = args.m8_file 
    else:
        args.blast_output_file = re.sub("\..*$","",args.systat_output_file) + ".blast_output.m8"
    global sat_pattern, min_copy_log_info
    min_copy_log_info = ''
    return args


def get_rosetta_info(rosetta_file):
    global category_set
    category_set = set()
    rosetta_dict = {}
    logging.info(f"creating the gene rosetta_dict from file {rosetta_file}")
    try:
        with open(rosetta_file, 'r') as f:
            for line in f:
                if re.search("gene\$",line):  # skips header, if present
                    continue
                #rosetta_tuple = re.findall("^([^\s\t]+)[\s\t]+(.+$)",line)[0]   # ('w', 'X eye')   
                rosetta_tuple = re.findall("^([^\s\t]+)[\s\t]+([^\s\t]+)",line)[0]   # ('w', 'X')
                gene_key   = rosetta_tuple[0]
                info_value = rosetta_tuple[1]
                rosetta_dict[gene_key] = info_value
                category_set.add(info_value)
    except Exception as e:
        logging.error(f"creating the gene rosetta_dict from file {args.rosetta}: {e}")
        raise
    return rosetta_dict


def extract_sequence_lengths_from_fasta(fasta_file):
    """Extracts sequence lengths directly from the given FASTA file."""
    size_dict = {}
    logging.info(f"Extracting sequence lengths from {fasta_file}")
    try:
        with open(fasta_file, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                seq_id = record.id
                seq_length = len(record.seq)
                size_dict[seq_id] = seq_length
    except Exception as e:
        logging.error(f"Error extracting sequence lengths from {fasta_file}: {e}")
        raise
    logging.info(f"Extracted lengths for {len(size_dict)} sequences.")
    return size_dict
    
    

def run_ncbi_blast(fasta_file, db, blast_output_file):
    """Runs BLAST using blastall against a given database and saves the output in m8 format."""
    blast_cmd_str = "blastall -p blastn -m 8 " + " -i " + fasta_file +  " -d " + db + " -o " + blast_output_file + " " + args.ncbi_blast_parameters  # default args.ncbi_blast_parameters: '-e 0.001  -F "m D" -b 100000 -v 100000 -a 50 '
    # print(blast_cmd_str, flush=True)
    logging.info(f"Running ncbi-blast search for database {db} with {fasta_file}.")
    logging.info(f"blast command: {blast_cmd_str}")
        
    try:
        subprocess.run(blast_cmd_str,shell=True, check=True)
        logging.info(f"BLAST search completed for {db}, results saved in {blast_output_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"BLAST search failed for database {db}: {e}")
        raise
    except FileNotFoundError:
        logging.critical(f"'blastall' executable not found. Ensure BLAST is installed and accessible.")
        raise


def run_wu_blast(fasta_file, db, blast_output_file):
    """Runs wu-blast against a given database and saves the output in m8 format."""
    WU_bls_output = re.sub("\..*$","",args.systat_output_file) + ".wu-bls"
    
    blast_cmd_str = "blastn " + db + " " + fasta_file   + " -o " + WU_bls_output +  "  mformat=2   -novalidctxok -nonnegok -gapall -restest  Q=7 R=2 kap M=1 N=-3 W=11    S2=14 gapS2=19 X=6 gapX=15 gapW=12    -gi  gapL=1.3741 gapK=0.711 gapH=1.3073 " + args.wu_blast_parameters  # default args.wu_blast_parameters: ' E=0.001 filter=dust -wordmask dust  V=100000 B=100000   hspmax=0  cpus=50 '
    print(blast_cmd_str, flush=True)
    
    awk_cmd_string = ' awk \'BEGIN{OFS="\t"};{print $1,$2,$11,$7,$10,$13,$18,$19,$21,$22,$3,$5}\'  ' + WU_bls_output + ' > ' + blast_output_file
    print(awk_cmd_string,flush=True) 
    
    logging.info(f"Running wu-blast search for database {db} with {fasta_file}.")
    logging.info(f"blast command: {blast_cmd_str}")
        
    try:
        subprocess.run(blast_cmd_str,shell=True, check=True)
        subprocess.run(awk_cmd_string,shell=True, check=True)       
        logging.info(f"BLAST search completed for {db}, results saved in {blast_output_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"BLAST search failed for database {db}: {e}")
        raise
    except FileNotFoundError:
        logging.critical(f"'blastall' executable not found. Ensure BLAST is installed and accessible.")
        raise





def parse_blast_m8_with_identity(blast_file):
    """Parses the BLAST output in m8 format"""
    global gene_cov_table_file, category_coverage_dict, category_count_dict
    previous_query_id = ""
    CDS_found_set = set()
    coverage_dict = {}  # per base
    if args.rosetta != "":
        category_coverage_dict = {}  # will store one list of coverages for each gene category (e.g. chrom) 
        category_count_dict = {}
        for gene_category in category_set:
            category_coverage_dict[gene_category] = [] 
            category_count_dict[gene_category] = 0
    zero_m8_hits = 0
    logging.info(f"Parsing BLAST results from {blast_file}")
    
    output_file_tmp = open(args.systat_output_file, 'w')
    if args.rosetta == "":
        output_file_tmp.write("\t".join(map(str,["gene$", "position", "depth"])) + "\n")
    else:
        output_file_tmp.write("\t".join(map(str,["gene$", "position", "depth","chrom$"])) + "\n")
         
    gene_cov_table_file_name = re.sub("\..*$","",args.systat_output_file) + ".gencov_table.txt"
    gene_cov_table_file = open(gene_cov_table_file_name, 'w')
    if args.rosetta == "":
        gene_cov_table_file.write("\t".join(map(str,["gene$", "mean", "median","mode"])) + "\n")
    else:
        gene_cov_table_file.write("\t".join(map(str,["gene$", "mean", "median","mode","chrom$"])) + "\n")
           
    try:
        with open(blast_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue  # Skip comment lines if any
                parts = line.strip().split('\t')
                if len(parts) < 12:
                    continue  # Skip incomplete lines
                identity = float(parts[2])  # The third column is the percentage identity
                if identity < args.identity_cutoff:
                    continue  # Skip hits with identity below cut-off 
                
                query_id = parts[0]
                
                qstart = int(parts[6])
                qend = int(parts[7])
                # Normalize start and end for consistency
                start = min(qstart, qend)
                # new in v2 : range(2,5,1) produces 2 3 4 
                # end = max(qstart, qend) 
                end = max(qstart, qend) +1  
                
                # print("line:",line,flush=True)
                
                if previous_query_id == "" :   # first gene in m8
                    for base in range(1,size_dict[query_id] +1,1):
                        coverage_dict[base] = 0  
                    previous_query_id = query_id
                    CDS_found_set.add(query_id)                
                    if args.rosetta != "" and previous_query_id in rosetta_dict:
                        gene_category = rosetta_dict[previous_query_id]
                    else:
                        gene_category = "not_found1"
                        
                if query_id != previous_query_id :    # new gene in m8
                    if args.rosetta != "" and previous_query_id in rosetta_dict:
                        gene_category = rosetta_dict[previous_query_id]
                    else:
                        gene_category = "not_found2"                    
                    
                    for base in range(1,size_dict[previous_query_id] +1,1):
                        # print(previous_query_id,base,coverage_dict[base])                       
                        if args.rosetta == "":
                            output_file_tmp.write("\t".join(map(str,[previous_query_id,base,coverage_dict[base]])) + "\n")  
                        else:
                            output_file_tmp.write("\t".join(map(str,[previous_query_id,base,coverage_dict[base],gene_category])) + "\n")                             
                    gene_cov_stats(previous_query_id,coverage_dict)                  
                    coverage_dict.clear()
                    previous_query_id = query_id
                    CDS_found_set.add(query_id)
                    for base in range(1,size_dict[query_id] +1,1):
                        coverage_dict[base] = 0
                    # new in v2a to correct serious bug: 
                    for base in range(start,end,1):    
                        coverage_dict[base] = 1                   
                else:
                    for base in range(start,end,1):    
                        coverage_dict[base] +=1                       
                    
            
            # processes the last gene
            if args.rosetta != "" and query_id in rosetta_dict:
                gene_category = rosetta_dict[query_id]
            else:
                gene_category = "not_found5"            
            for base in range(1,size_dict[query_id] +1,1):    
                # print(previous_query_id,base,coverage_dict[base])    
                if args.rosetta == "":
                    output_file_tmp.write("\t".join(map(str,[previous_query_id,base,coverage_dict[base]])) + "\n")  
                else:
                    output_file_tmp.write("\t".join(map(str,[previous_query_id,base,coverage_dict[base],gene_category])) + "\n")               
            gene_cov_stats(query_id,coverage_dict)
 
    except Exception as e:
        logging.error(f"Error parsing BLAST results from {blast_file}: {e}")
        raise
    
    for seq_id in size_dict:   # handle the fasta entries devoid of any hit in blast.m8
        if seq_id not in CDS_found_set:
            zero_m8_hits +=1
            zerocov_dict = {}
            if args.rosetta != "" and seq_id in rosetta_dict:
                gene_category = rosetta_dict[seq_id]
            else:
                gene_category = "not_found3"             
            
            for base in range(1,size_dict[seq_id] +1,1):
                # print(seq_id,base,0) 
                if args.rosetta == "":               
                    output_file_tmp.write("\t".join(map(str,[seq_id,base,0])) + "\n")
                else:                   
                    output_file_tmp.write("\t".join(map(str,[seq_id,base,0,gene_category])) + "\n")                                            
                zerocov_dict[base] = 0
            gene_cov_stats(seq_id,zerocov_dict)
    output_file_tmp.close()    
    logging.info(f"Found {len(CDS_found_set)} query sequences from {blast_file}")
    if zero_m8_hits > 0:
        logging.info(f"  {zero_m8_hits} CDS sequences from {args.fasta_file} did not have hits in the blast m8 file")
    # gene_cov_table_file.close()
    return


def gene_cov_stats(gene,cov_dict):
    global category_coverage_dict , category_count_dict
    cov_list = list(cov_dict.values())
    if args.rosetta == "":    
        # cov_list = list(cov_dict.values())
        print(gene,round(statistics.mean(cov_list),1),round(statistics.median(cov_list)),statistics.mode(cov_list),sep="\t")
        gene_cov_table_file.write("\t".join(map(str,[gene,round(statistics.mean(cov_list),1),round(statistics.median(cov_list)),statistics.mode(cov_list)])) + "\n")
    else:
        if args.rosetta != "" and gene in rosetta_dict:
            gene_category = rosetta_dict[gene]
            category_coverage_dict[gene_category] = category_coverage_dict[gene_category] + cov_list 
            category_count_dict[gene_category] += 1
        else:
            gene_category = "not found4"
        # cov_list = list(cov_dict.values())
        print(gene,round(statistics.mean(cov_list),1),round(statistics.median(cov_list)),statistics.mode(cov_list),gene_category,sep="\t")
        gene_cov_table_file.write("\t".join(map(str,[gene,round(statistics.mean(cov_list),1),round(statistics.median(cov_list)),statistics.mode(cov_list),gene_category])) + "\n")       
    return


def plot_gene_coverage_with_smoothing(systat_output_file, window_size):
    """Reads the systat output file and generates smoothed coverage plots for each gene."""
    logging.info(f"Generating smoothed coverage plots from {systat_output_file} with window size {window_size}")

    try:
        # Read the data into a pandas DataFrame
        df = pd.read_csv(systat_output_file, sep='\t')

        # Check if the required columns exist
        if not {'gene$', 'position', 'depth'}.issubset(df.columns):
            logging.error("Error: The expected columns (gene$, position, depth) are not found in the file.")
            return

        # Apply a moving median filter separately for each gene
        df['smoothed_depth'] = 0
        for gene in df['gene$'].unique():
            gene_df = df[df['gene$'] == gene].copy()
            gene_df['smoothed_depth'] = gene_df['depth'].rolling(window=window_size, center=True, min_periods=1).median()
            df.loc[df['gene$'] == gene, 'smoothed_depth'] = gene_df['smoothed_depth']

        # Ensure proper formatting of the output file
        df.to_csv(systat_output_file, sep='\t', index=False, float_format='%.1f')
        logging.info(f"Updated {systat_output_file} with smoothed depth values.")

        # Set Seaborn style for better visuals
        sns.set_style("white")
        sns.set_style("ticks")
        # Process each gene separately
        for gene in df['gene$'].unique():
            gene_df = df[df['gene$'] == gene]

            # Create the plot
            # plt.figure(figsize=(12, 6))
            plt.figure(figsize=(args.figsize_W, args.figsize_H))

#           plt.plot(gene_df['position'], gene_df['depth'], marker='o', linestyle='-', color='lightgray', alpha=0.6, label="Original")
            plt.plot(gene_df['position'], gene_df['smoothed_depth'], marker='', linestyle='-', linewidth=2, color='#005f99', label="Smoothed (Moving Median)")
            plt.xlabel(args.Xlabel, fontsize=12)
            plt.ylabel(args.Ylabel, fontsize=12)
#           plt.title(f"Coverage Depth for {gene} (Smoothed)", fontsize=14, fontweight='bold')
#           plt.legend(fontsize=10, frameon=True, loc='upper right')
            plt.grid(False)

            # Set x and y axis limits to start at 0
            plt.xlim(left=0)
            if args.Ymax == -1: 
                plt.ylim(bottom=args.Ymin)
            else:
                plt.ylim(args.Ymin, args.Ymax)
            
            # Adjust left and bottom spines to align with the border
            plt.gca().spines['left'].set_position(('outward', 0))
            plt.gca().spines['bottom'].set_position(('outward', 0))

            # Save the plot
            # output_filename = f"{gene}_coverage_smoothed_plot.png"           
            if args.graph_name_suffix == '':
                output_filename = f"{gene}_coverage_smoothed{args.window_size}.png"
            else:
                output_filename = f"{gene}_coverage_smoothed{args.window_size}_{args.graph_name_suffix}.png"
                
            plt.savefig(output_filename, dpi=300, bbox_inches='tight')
            plt.close()
            logging.info(f"Saved smoothed coverage plot: {output_filename}")


    except Exception as e:
        logging.error(f"Error generating smoothed coverage plots: {e}")
        raise


def main():
    global size_dict, rosetta_dict
    # Initialize dictionary to hold results per database
    logging.info(f"Starting occupancy analysis for {args.fasta_file}")
    # loading gnew info that will be attached to output files (typically chromosome location)
    if args.rosetta != "":
        rosetta_dict = get_rosetta_info(args.rosetta)
    # Extract sequence lengths from the FASTA file
    try:
        size_dict = extract_sequence_lengths_from_fasta(args.fasta_file)
    except Exception as e:
        logging.critical(f"Failed to extract sequence lengths: {e}")
        return

    # Run BLAST and parse all hits
    if args.m8_file == "" :
        try:
            logging.debug(f"Running BLAST against {args.db}")
            if args.blast_pgm == "ncbi":
                run_ncbi_blast(args.fasta_file, args.db, args.blast_output_file)
            if args.blast_pgm == "wu":
                run_wu_blast(args.fasta_file, args.db, args.blast_output_file)                  
        except Exception as e:
            logging.critical(f"Failed to run BLAST for database {args.db}: {e}")        
            return
        
    # Parse BLAST results with identity percentages and produce the final files   
    try:
        logging.debug(f"parsing blast m8 file ")
        parse_blast_m8_with_identity(args.blast_output_file)
    except Exception as e:
        logging.critical(f"Failed to  parse blast m8 file: {e}")
        return
    if args.rosetta != "":
        gene_cov_table_file.write("\n\nclass\tgenes\tbp\tmean\tmedian\tmode\n")
        print("\n\nclass","genes","bp","mean","median","mode",sep="\t")
        for gene_category in category_set:
            data = category_coverage_dict[gene_category]
            # print("gene_category:" , gene_category, flush=True)
            print(gene_category,category_count_dict[gene_category],len(data), round(statistics.mean(data),1),round(statistics.median(data)),statistics.mode(data),sep="\t")
            gene_cov_table_file.write("\t".join(map(str,[gene_category,category_count_dict[gene_category],len(data),round(statistics.mean(data),1),round(statistics.median(data)),statistics.mode(data)])) + "\n")   
    gene_cov_table_file.close()
    logging.info("coverage analysis completed successfully.")

if __name__ == "__main__":
    args = get_parameters()   
    # Set up logging
    log_file = re.sub("\..*$","",args.systat_output_file) + "_coverage_log.txt"
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(levelname)s: %(message)s',                        
                        handlers=[
                            logging.FileHandler(log_file),
                            logging.StreamHandler(sys.stdout)
                        ]
                        )
    main()
    plot_gene_coverage_with_smoothing(args.systat_output_file, args.window_size)
