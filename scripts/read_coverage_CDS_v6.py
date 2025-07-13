#!/usr/bin/python3.9 -u
# #!/usr/local/bin/python3.9 -u
# read_coverage_CDS_v6.py   based on digital_kennison_v9.py    Bernardo   moj_Y and nanopore_missing projects   v. 15may2025
# new in v2: serious bug correction;  implement option --m8_file
# -F "m S"
# new in v4: Generate a depth vs. position plot for each gene using matplotlib; Apply a rolling median filter to smooth out spikes
# new in v5:  Added --cds_file and --cds_Yaxis parameters for exon visualization
# new in v6    proportional Y-axis position where exon lines will be displayed (default: 5.0, meaning Y=5 in a graph with Ymax=100). <mantains exons at vthe same relatve position in the graphs 
# 15may2025 version: new option --Xlabel gene_name   new parameter --select_fasta
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
    global effective_fasta_file
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
    parser.add_argument('--rosetta', default="", type=str, help='optional file containing gene info that will be attached to the output text files. Typically, chrom location.')
    parser.add_argument("--window_size", type=int, default=5, help="Size of the moving median window")
    parser.add_argument('--graph_name_suffix', type=str, default='', help='suffix to be added to automatic name (eg, porechop)')
    parser.add_argument("--figsize_W", type=float, default=12, help="fig size WIDTH (in inches)")
    parser.add_argument("--figsize_H", type=float, default=6, help="fig size HEIGHT (in inches)")
    # parser.add_argument("--Xmin", type=int, default=0, help="Xmin value")
    # parser.add_argument("--Xmax", type=int, default=0, help="Xmax value")
    parser.add_argument("--Ymin", type=float, default=0, help="Xmin value")
    parser.add_argument("--Ymax", type=float, default=-1, help="Ymax value")
    parser.add_argument('--Xlabel', type=str, default='CDS position (bp)', help='label of the X-axis. deafult is "CDS position (bp)" . gene_name will put something like "position in AMPdean-RC HiFi" ')    
    parser.add_argument('--Ylabel', type=str, default='Coverage', help='label of the Y-axis')
    parser.add_argument('--select_fasta', type=str, default='', help='only process a subset of fasta ("Myo81F-RB Mitf-RA" ) ')
    # NEW PARAMETERS FOR EXON VISUALIZATION
    parser.add_argument('--cds_file', type=str, default='', help='optional file containing exon positions (start end exon_id) to display on the plot')
    parser.add_argument('--cds_Yaxis', type=float, default=5.0, help='proportional Y-axis position where exon lines will be displayed (default: 5.0, meaning Y=5 in a graph with Ymax=100)')
    parser.add_argument("--exon_fontsize", type=int, default=8, help="font size of the exon labels")
    parser.add_argument("--exon_fontweight", type=str, default='normal', help="font weight of the exon labels   normal/bold/heavy/light/ultrabold/ultralight")
    parser.add_argument("--graph_format", type=str, default='png' ,  help="output graph  format: png/pdf/svg/eps")
    parser.add_argument("--linewidth", type=float, default=1.5, help="linewidth of the coverage data")
    parser.add_argument("--output", default="", type=str, help="graph output file for visualization. File extension dtermines formnat PNG/PDF/SVG. Overrides --graph_name_suffix and --graph_format ")



    args = parser.parse_args()
    if (args.m8_file != "" and args.db != "") or  (args.m8_file == "" and args.db == "") :
        print("ERROR: either --m8_file OR --db parameters must be set. Setting both or none is not allowed" , flush=True)
        exit(1)

    if args.m8_file != "" :
        args.blast_output_file = args.m8_file
    else:
        args.blast_output_file = re.sub("\..*$","",args.systat_output_file) + ".blast_output.m8"

    if args.select_fasta != "" :
        if os.path.exists("temp.fasta"):
            os.remove("temp.fasta")
        selected_fasta_list = args.select_fasta.split()
        print(args.select_fasta, selected_fasta_list)
        for record in SeqIO.parse(args.fasta_file, "fasta"):
            if record.id in selected_fasta_list:
                print(record.id, " found")
                with open("temp.fasta", "a") as handle:
                    SeqIO.write(record, handle, "fasta")
        effective_fasta_file = "temp.fasta"
    else:
        effective_fasta_file = args.fasta_file

    return args


def parse_cds_file(cds_file):
    """Parse the CDS/exon file and return a list of exon information."""
    exons = []
    if not cds_file or not os.path.exists(cds_file):
        return exons
    
    logging.info(f"Parsing CDS/exon file: {cds_file}")
    try:
        with open(cds_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):  # Skip empty lines and comments
                    continue
                
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        start = int(parts[0])
                        end = int(parts[1])
                        exon_id = parts[2]
                        exons.append({'start': start, 'end': end, 'id': exon_id})
                    except ValueError as e:
                        logging.warning(f"Invalid data in CDS file line {line_num}: {line} - {e}")
                else:
                    logging.warning(f"Insufficient columns in CDS file line {line_num}: {line}")
    
    except Exception as e:
        logging.error(f"Error parsing CDS file {cds_file}: {e}")
        return []
    
    logging.info(f"Successfully parsed {len(exons)} exons from {cds_file}")
    return exons


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
        category_coverage_dict = {}  # will store one list of coverages for each gene category (eg, chromn)
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
        logging.info(f"  {zero_m8_hits} CDS sequences from {effective_fasta_file} did not have hits in the blast m8 file")
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

    # Parse CDS/exon file if provided
    exons = []
    if args.cds_file:
        exons = parse_cds_file(args.cds_file)

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
            plt.figure(figsize=(args.figsize_W, args.figsize_H))

            # Plot the main coverage data
            plt.plot(gene_df['position'], gene_df['smoothed_depth'], marker='', linestyle='-', linewidth=args.linewidth, color='#005f99', label="Smoothed (Moving Median)")
            
            '''
            # Add exon visualization if CDS file is provided
            if exons:
                logging.info(f"Adding {len(exons)} exons to plot for gene {gene}")
                
                # Calculate vertical offset for each exon to prevent overlap
                exon_height = 2.0  # Height of each exon line
                base_y = args.cds_Yaxis
                
                for i, exon in enumerate(exons):
                    # Calculate Y position with slight vertical offset for each exon
                    # Odd-numbered exons (1,3,5...) go to higher position
                    y_pos = base_y + ((i + 1) % 2) * exon_height  # Cycle through 2 levels
                    
                    # Draw the exon as a horizontal line
                    plt.hlines(y=y_pos, xmin=exon['start']+4, xmax=exon['end']-4, 
                              colors='black', linewidth=2, alpha=0.8)
                    
                    # Add exon label - all labels at the higher position
                    label_y_pos = base_y + exon_height + 1  # All labels at higher position
                    mid_point = (exon['start'] + exon['end']) / 2
                    plt.text(mid_point, label_y_pos + 0.5, f"{exon['id']}", 
                            ha='center', va='bottom', fontsize=8, fontweight='bold')
            '''
            # Set labels and formatting
            if args.Xlabel != "gene_name":
                plt.xlabel(args.Xlabel, fontsize=12)
            else:
                xlabel_2 = gene + " " + args.graph_name_suffix + "  (bp)"
                plt.xlabel(xlabel_2, fontsize=12)
            plt.ylabel(args.Ylabel, fontsize=12)
            plt.grid(False)

            # Set x and y axis limits
            plt.xlim(left=0)
            if args.Ymax == -1:
                plt.ylim(bottom=args.Ymin)
            else:
                plt.ylim(args.Ymin, args.Ymax)

            
              
            
            # Add exon visualization if CDS file is provided
            if exons:
                logging.info(f"Adding {len(exons)} exons to plot for gene {gene}")
                curr_ymin, curr_ymax = plt.ylim()
                # print("curr_ymax=", curr_ymax)
                Yexon_scale = curr_ymax/100  # new in 4jun2025. To put the exons at the same relative position, irrespective of coverage (ie, Ymax)
                # Calculate vertical offset for each exon to prevent overlap
                exon_height = 2.0  # Height of each exon line
                
                base_y = args.cds_Yaxis   
                
                for i, exon in enumerate(exons):
                    # Calculate Y position with slight vertical offset for each exon
                    # Odd-numbered exons (1,3,5...) go to higher position
                    # y_pos =  base_y + ((i + 1) % 2) * exon_height     # Cycle through 2 levels
                    y_pos = ( base_y + ((i + 1) % 2) * exon_height ) * Yexon_scale  # Cycle through 2 levels
                    
                    # Draw the exon as a horizontal line
                    plt.hlines(y=y_pos, xmin=exon['start']+4, xmax=exon['end']-4, 
                              colors='black', linewidth=2, alpha=0.8)
                    
                    # Add exon label - all labels at the higher position
                    # label_y_pos = base_y + exon_height + 1  # All labels at higher position
                    # mid_point = (exon['start'] + exon['end']) / 2
                    # plt.text(mid_point, label_y_pos + 0.5, f"{exon['id']}", 
                            # ha='center', va='bottom', fontsize=8, fontweight='bold')

                    # Add exon label - all labels at the higher position
                    label_y_pos = ( base_y + exon_height + 1 ) * Yexon_scale  # All labels at higher position
                    
                    mid_point = (exon['start'] + exon['end']) / 2
                    plt.text(mid_point, label_y_pos, f"{exon['id']}", 
                            ha='center', va='bottom', fontsize=args.exon_fontsize, fontweight=args.exon_fontweight)
            
            # Adjust left and bottom spines to align with the border
            plt.gca().spines['left'].set_position(('outward', 0))
            plt.gca().spines['bottom'].set_position(('outward', 0))

            # Save the plot
            if args.output == '':      # new in 10jun2025
                if args.graph_name_suffix == '':
                    output_filename = f"{gene}_coverage_smoothed{args.window_size}.{args.graph_format}"
                else:
                    output_filename = f"{gene}_coverage_smoothed{args.window_size}_{args.graph_name_suffix}.{args.graph_format}" # new in 7jun2025: chose png, svg, etc
            else:
                output_filename = args.output
            plt.savefig(output_filename, dpi=600, bbox_inches='tight')   # new in 7jun2025: 600dpi instead of 300dpi
            plt.close()
            logging.info(f"Saved smoothed coverage plot with exons: {output_filename}")

    except Exception as e:
        logging.error(f"Error generating smoothed coverage plots: {e}")
        raise


def main():
    global size_dict, rosetta_dict
    # Initialize dictionary to hold results per database
    logging.info(f"Starting occupancy analysis for {effective_fasta_file}")
    # loading gnew info that will be attached to output files (typically chromosome location)
    if args.rosetta != "":
        rosetta_dict = get_rosetta_info(args.rosetta)
    # Extract sequence lengths from the FASTA file
    try:
        # size_dict = extract_sequence_lengths_from_fasta(args.fasta_file)
        size_dict = extract_sequence_lengths_from_fasta(effective_fasta_file)
    except Exception as e:
        logging.critical(f"Failed to extract sequence lengths: {e}")
        return

    # Run BLAST and parse all hits
    if args.m8_file == "" :
        try:
            logging.debug(f"Running BLAST against {args.db}")
            if args.blast_pgm == "ncbi":
                run_ncbi_blast(effective_fasta_file, args.db, args.blast_output_file)
            if args.blast_pgm == "wu":
                run_wu_blast(effective_fasta_file, args.db, args.blast_output_file)
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
    if args.select_fasta != "" and os.path.exists("temp.fasta"):
        os.remove("temp.fasta")
