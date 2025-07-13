#!/usr/bin/env python3.9

#############################################################################################################################################
# visualize_repeats_censor4.py  :  v.  8jun2025
# streamlined repeat classification. Sat DNA is now identificed by a regex; TEs are now the uncclassified stuff, after  aseries of elif
#                                                                                                                                        
# Fabiana                                                                                                         
#                                                                                                                                           
# The script processes a Censor output file (.map), classifies sequence regions (e.g., simple repeats, TEs, and Ycds), generates a          
# graphical visualization, and outputs sequence composition statistics.                                                                    
#                                                                                                                                           
# new in v2.5: Adds "Details" column with the most common satellite DNA in the summary output.
# new in v2.6:  new  parameter -len_graph: used to make all graphs proportional (eh, 100kb), w/o  disturbing the "Details" column             
# new in v2.7: new parameters --gene_regex (instead Ycds fixed))  and --gene_color (instead of fixed red)  17may2025
# New in v3.1: Modified to read exon/CDS information from external --cds_file   with exon ID labels.
# new in v4:   added parameters    --Xmin  --Xmax  --exon_fontsize  --exon_fontweight   
#############################################################################################################################################


import argparse
import re
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import numpy as np

def calculate_at_content(sequence):
    """Calculate AT% for a simple repeat sequence."""
    match = re.match(r'^\((.*?)\)n$', sequence)
    if not match:
        return None  # Not a simple repeat
    repeat_bases = match.group(1)
    at_count = repeat_bases.count('A') + repeat_bases.count('T')
    total_count = len(repeat_bases)
    return (at_count / total_count) * 100 if total_count > 0 else None


def parse_cds_file(cds_file):
    """Parse the CDS/exon file and return a list of exon information."""
    exons = []
    if not cds_file or not os.path.exists(cds_file):
        return exons
    
    print(f"Parsing CDS/exon file: {cds_file}")
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
                        if start >= args.Xmin and end <= args.Xmax:  # new in 6jun2025
                            exons.append({'start': start, 'end': end, 'id': exon_id})
                    except ValueError as e:
                        print(f"Warning: Invalid data in CDS file line {line_num}: {line} - {e}")

                elif len(parts) == 2:   # new in 10jun2025  allow unlabeled exons
                    try:
                        start = int(parts[0])
                        end = int(parts[1])
                        exon_id = " "
                        if start >= args.Xmin and end <= args.Xmax:  # new in 6jun2025
                            exons.append({'start': start, 'end': end, 'id': exon_id})
                    except ValueError as e:
                        print(f"Warning: Invalid data in CDS file line {line_num}: {line} - {e}")    
                        
                else:
                    print(f"Warning: Insufficient columns in CDS file line {line_num}: {line}")
    
    except Exception as e:
        print(f"Error parsing CDS file {cds_file}: {e}")
        return []
    
    print(f"Successfully parsed {len(exons)} exons from {cds_file}")
    return exons
    
    


def process_censor_output(censor_file, total_length):
    """Parses the Censor output file and classifies sequence regions."""
    classifications = []
    satellite_counts = {}  # Tracks satellite counts for summaries
    debug_log = []

    with open(censor_file, 'r') as file:
        for line in file:
            if line.strip():
                try:
                    columns = line.split()
                    start = int(columns[1])
                    end = int(columns[2])
                    repeat = columns[3]  # 4th column
                    bases = end - start + 1
                    
                    if re.match(fr"({args.exclude_regex})", repeat):   # skips undesired entries in censor, eg, anaother gene 
                        continue
                    if  start < args.Xmin or end > args.Xmax:  # new in 6jun2025
                        continue
                                                            
                    # if re.match(r'^Ycds_', repeat):  # Handle Ycds sequences      
                    if re.match(fr"({args.gene_regex})", repeat):  # Handle gene names  w/ dynamic regex     
                        classifications.append({
                            'type': 'gene',  # Rename for clarity
                            'start': start,
                            'end': end
                        })
                        debug_log.append(f"gene: {repeat} Start: {start} End: {end} Bases: {bases}")

                    elif bases < args.min_repeat_len :   # skips small sat or TE  
                        continue
                        
                    elif re.match(fr"({args.repeat_regex})", repeat):  # Handle a specific repeat name  w/ dynamic regex    e.g.,  1.688_GK_42F7/2-X
                        classifications.append({
                            'type': 'repeat_regex',   
                            'start': start,
                            'end': end
                        })
                        debug_log.append(f"repeat_regex: {repeat} Start: {start} End: {end} Bases: {bases}")                    
                    
                    elif re.match("\([ATGC]+\)n", repeat):    #new in v3. simple repeat
                        at_content = calculate_at_content(repeat)
                        classifications.append({
                            'type': 'simple_repeat',
                            'start': start,
                            'end': end,
                            'at_content': at_content
                        })
                        if repeat in satellite_counts:
                            satellite_counts[repeat] += bases
                        else:
                            satellite_counts[repeat] = bases
                        debug_log.append(f"Simple Repeat: {repeat} Start: {start} End: {end} Bases: {bases}")                    
                    
                    else:  # Transposable elements
                        # if re.match(r'^[A-Za-z0-9_-]+$', repeat):
                        classifications.append({
                            'type': 'TE',
                            'start': start,
                            'end': end
                        })
                        debug_log.append(f"TE: {repeat} Start: {start} End: {end} Bases: {bases}")
                        
                    # else:  # Unclassified  will never reach here, will be cathed before , as a TE
                        # debug_log.append(f"Unclassified: {repeat} Start: {start} End: {end} Bases: {bases}")

                except (ValueError, IndexError) as e:
                    debug_log.append(f"Skipping malformed line: {line.strip()} Error: {e}")

    # Output debug log
    with open("debug_log.txt", "w") as debug_file:
        debug_file.write("\n".join(debug_log))

    # Prepare satellite summary
    satellite_summary = [
        (repeat, bases, (bases / total_length) * 100)
        for repeat, bases in satellite_counts.items()
    ]
    satellite_summary = sorted(satellite_summary, key=lambda x: x[2], reverse=True)
    print("\ninput file:", args.censor)    
    # Print Top 10 satellite
    print("Top 10 Satellite DNA Summary:")
    print(f"{'Repeat':<15} {'Bases':<10} {'Percentage':<10}")
    for repeat, bases, percentage in satellite_summary[:10]:
        print(f"{repeat:<15} {bases:<10} {percentage:.1f}%")

    return classifications, satellite_summary


def summarize_categories(classifications, total_length, satellite_summary):
    """Summarize the sequence composition by category (single-copy, TE, satellite)."""
    total_sat_bases = 0
    total_te_bases = 0
    total_ycds_bases = 0
    total_repeat_regex_bases = 0    
    covered_bases = 0  # Tracks bases covered by satellites, TEs, and gene

    for entry in classifications:
        start = entry['start']
        end = entry['end']
        bases = end - start + 1

        if entry['type'] == 'simple_repeat':
            total_sat_bases += bases
        elif entry['type'] == 'gene':
            total_ycds_bases += bases
        elif entry['type'] == 'repeat_regex':
            total_repeat_regex_bases += bases
        elif entry['type'] == 'TE':
            total_te_bases += bases
        covered_bases += bases

    # Calculate single-copy bases as the remainder
    single_copy_bases = total_length - covered_bases

    # Get the most common satellite DNA
    if satellite_summary:
        most_common_satellite, _, most_common_percentage = satellite_summary[0]
        details = f"{most_common_satellite} {most_common_percentage:.1f}%"
    else:
        details = "None"

    # Print the summary
    print("\nSequence Composition Summary:")
    print(f"{'Category':<15} {'Bases':<10} {'Percentage':<15} {'Details':<20}")
    print(f"{'single-copy':<15} {single_copy_bases + total_ycds_bases:<10} {100 * (single_copy_bases + total_ycds_bases) / total_length:.1f}%")
    print(f"{'TE':<15} {total_te_bases:<10} {100 * total_te_bases / total_length:.1f}%")
    print(f"{'sat':<15} {total_sat_bases:<10} {100 * total_sat_bases / total_length:.1f}%           {details}")
    if args.repeat_regex != 'not_used':
        print(f"{'repeat_regex':<15} {total_repeat_regex_bases:<10} {100 * total_repeat_regex_bases / total_length:.1f}%")



def create_visualization(classifications, total_length, output_file, exons, cds_yaxis):
    """Creates the graphical output using Matplotlib."""
    fig, ax = plt.subplots(figsize=(10, 2))
    # ax.set_xlim(0, total_length)
    ax.set_xlim(args.Xmin, total_length)  # new in 6jun2025
    ax.set_ylim(0, 1)

    # Plot classified regions
    for entry in classifications:
        start = entry['start']
        end = entry['end']
        width = end - start

        if entry['type'] == 'simple_repeat':
            at_content = entry['at_content']
            color = (
                0 + (0.7 - 0) * (1 - at_content / 100),  # Red component
                0 + (0.9 - 0) * (1 - at_content / 100),  # Green component
                0.5 + (1 - 0.5) * (1 - at_content / 100)  # Blue component
            )
            ax.add_patch(
                patches.Rectangle(
                    (start, 0.375), width, 0.25, facecolor=color, zorder=2
                )
            )
        elif entry['type'] == 'TE':
            ax.add_patch(
                patches.Rectangle(
                    (start, 0.375), width, 0.25,
                    facecolor="lightcoral", edgecolor=None, linewidth=0, zorder=2
                )
            )

        
        elif entry['type'] == 'repeat_regex':
            ax.add_patch(
                patches.Rectangle(
                    (start, 0.375), width, 0.25,
                    facecolor=args.repeat_color, edgecolor=None, linewidth=0, zorder=2
                )
            )

        elif entry['type'] == 'gene' and len(exons)==0 :  # no external exon location file. Using censor coordinates (less accurate)
            ax.add_patch(
                patches.Rectangle(
                    (start, 0.65), width, 0.1,
                    facecolor=args.gene_color, edgecolor=None, linewidth=0, zorder=3
                )
            )


    # Plot exons from external CDS file, if available
    if exons:
        print(f"Adding {len(exons)} exons to visualization")
        for exon in exons:
            start = exon['start']
            end = exon['end']
            width = end - start
            exon_id = exon['id']
            
            # Draw exon rectangle
            ax.add_patch(
                patches.Rectangle(
                    (start, cds_yaxis), width, 0.1,
                    facecolor=args.gene_color, edgecolor=None, linewidth=0, zorder=3
                )
            )
            
            # Add exon label above the rectangle
            mid_point = (start + end) / 2
            ax.text(mid_point, cds_yaxis + 0.15, f"{exon_id}", 
                   ha='center', va='bottom', fontsize=args.exon_fontsize, fontweight=args.exon_fontweight)


    # Draw baseline  (flat line for unmatched  sequences (presumably single-copy)  in the middle of rectangles)
    # ax.plot([0, total_length], [0.5, 0.5], color="black", lw=0.5, zorder=1)
    ax.plot([args.Xmin, total_length], [0.5, 0.5], color="black", lw=0.5, zorder=1)  # new in 6jun2025
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    # Set x-axis formatting       3jun2025:  copied from previous version
    tick_interval_minor = args.tick_minor  # Minor ticks every 5 kb
    tick_interval_major = args.tick_major  # Major labels every 10 kb
    # minor_ticks = np.arange(0, total_length + tick_interval_minor, tick_interval_minor)
    # major_ticks = np.arange(0, total_length + tick_interval_major, tick_interval_major)
    # new in v4d: 
    # minor_ticks = np.arange(0, total_length , tick_interval_minor)
    # major_ticks = np.arange(0, total_length , tick_interval_major) 
    minor_ticks = np.arange(args.Xmin, total_length , tick_interval_minor)    # new in 6jun2025
    major_ticks = np.arange(args.Xmin, total_length , tick_interval_major)   # new in 6jun2025
    major_labels = [f"{int(tick/1000)} kb" for tick in major_ticks]
    ax.set_xticks(minor_ticks, minor=True)  # Set minor ticks
    ax.set_xticks(major_ticks)  # Set major ticks
    # ax.set_xticklabels(major_labels, fontsize=8)  # Reduce font size for readability
    ax.set_xticklabels(major_labels, fontsize=args.scale_fontsize)
    # Customize tick appearance
    ax.tick_params(axis="x", which="minor", length=4, width=0.5)  # Smaller minor ticks
    ax.tick_params(axis="x", which="major", length=6, width=1.0)  # Larger major ticks


    ax.set_yticks([])
    # ax.set_xlabel("Position")
    plt.tight_layout()
    plt.savefig(output_file, dpi=600)  # new in 7jun2025  was 300dpi
    plt.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process Censor output and classify repeats.")
    parser.add_argument("--censor", required=True, help="Path to the Censor output file.")
    parser.add_argument("--len", type=int, required=True, help="Total length of the FASTA sequence.")
    parser.add_argument("--output", required=True, help="graph output file for visualization. File extension dtermines formnat PNG/PDF/SVG")
    parser.add_argument("--len_graph", type=int, default=0,  help="Total length of the FASTA sequence for the graph.")
    parser.add_argument("--tick_minor", type=int, default=5000,  help="tick_interval_minor (unlabelled marks) for the graph.")
    parser.add_argument("--tick_major", type=int, default=10000,  help="tick_interval_major (labelled marks) for the graph.")   
    parser.add_argument("--gene_regex", type=str, default='Ycds',  help="gene names regex to detect gene hits.")   
    parser.add_argument("--gene_color", type=str, default='red',  help="color that marks gene sequences")   
    parser.add_argument("--repeat_regex", type=str, default='not_used',  help="repeat name  regex to detect specific  repeats. eg, 1.688_GK.")   
    parser.add_argument("--repeat_color", type=str, default='orange',  help="color that marks specific  repeats. eg, 1.688_GK.")   
    parser.add_argument("--exclude_regex", type=str, default='not_used',  help="  regex to ignore specific  repeats")   
    parser.add_argument("--min_repeat_len", type=int, default=0,  help="repeats smaller than this value will be ignored") 
    # NEW PARAMETERS FOR EXTERNAL CDS FILE
    parser.add_argument('--cds_file', type=str, default='', help='optional file containing exon positions (start end exon_id) to display on the plot')
    parser.add_argument('--cds_Yaxis', type=float, default=0.65, help='Y-axis position where exon rectangles will be displayed (default: 0.65)') 
    parser.add_argument('--Xmin', type=int, default=0,  help="start of the graph (bp)") 
    parser.add_argument('--Xmax', type=int, default=0,  help="end of the graph (bp)") 
    parser.add_argument("--exon_fontsize", type=float, default=10, help="font size of the exon labels")
    parser.add_argument("--exon_fontweight", type=str, default='normal', help="font weight of the exon labels   normal/bold/heavy/light/ultrabold/ultralight")
    parser.add_argument("--scale_fontsize", type=float, default=10, help="font size of the X-axis scale")


    
    args = parser.parse_args()
    if args.len_graph == 0:
        args.len_graph = args.len
    if args.Xmax == 0:
        args.Xmax = args.len_graph        

    # Parse external CDS file if provided
    exons = []
    if args.cds_file:
        exons = parse_cds_file(args.cds_file)


    classifications, satellite_summary = process_censor_output(args.censor, args.len)
    summarize_categories(classifications, args.len, satellite_summary)  # FU added "exons", but I will not use them here
    # create_visualization(classifications, args.len_graph, args.output, exons, args.cds_Yaxis)  # FU used args.len here, but should be args.len_graph
    create_visualization(classifications, args.Xmax, args.output, exons, args.cds_Yaxis)  
