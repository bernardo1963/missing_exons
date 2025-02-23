#!/usr/bin/env python3.9

#############################################################################################################################################
#                                                                                                                                           #
# Fabiana v2.5        v2.6 cosmetic changes + new parameter -len_graph                                                                                                         #
#                                                                                                                                           #
# The script processes a Censor output file (.map), classifies sequence regions (e.g., simple repeats, TEs, and Ycds), generates a          #
# graphical visualization, and outputs sequence composition statistics.                                                                    #
#                                                                                                                                           #
# NEW(2.5): Adds "Details" column with the most common satellite DNA in the summary output.       
# new in 2.6  new  parameter -len_graph: used to make all graphs proportional (eh, 100kb), w/o  disturbing the "Details" column                                                                           #
#############################################################################################################################################

import argparse
import re
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np  # Import numpy for tick spacing

def calculate_at_content(sequence):
    """Calculate AT% for a simple repeat sequence."""
    match = re.match(r'^\((.*?)\)n$', sequence)
    if not match:
        return None  # Not a simple repeat
    repeat_bases = match.group(1)
    at_count = repeat_bases.count('A') + repeat_bases.count('T')
    total_count = len(repeat_bases)
    return (at_count / total_count) * 100 if total_count > 0 else None


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

                    if re.match(r'^Ycds_', repeat):  # Handle Ycds sequences
                        classifications.append({
                            'type': 'Ycds',  # Rename for clarity
                            'start': start,
                            'end': end
                        })
                        debug_log.append(f"Ycds: {repeat} Start: {start} End: {end} Bases: {bases}")
                    else:
                        at_content = calculate_at_content(repeat)
                        if at_content is not None:  # Satellite DNA
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
                            if re.match(r'^[A-Za-z0-9_-]+$', repeat):
                                classifications.append({
                                    'type': 'TE',
                                    'start': start,
                                    'end': end
                                })
                                debug_log.append(f"TE: {repeat} Start: {start} End: {end} Bases: {bases}")
                            else:  # Unclassified
                                debug_log.append(f"Unclassified: {repeat} Start: {start} End: {end} Bases: {bases}")

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

    covered_bases = 0  # Tracks bases covered by satellites, TEs, and Ycds

    for entry in classifications:
        start = entry['start']
        end = entry['end']
        bases = end - start + 1

        if entry['type'] == 'simple_repeat':
            total_sat_bases += bases
        elif entry['type'] == 'TE':
            total_te_bases += bases
        elif entry['type'] == 'Ycds':
            total_ycds_bases += bases

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


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def create_visualization(classifications, total_length, output_file):
    """Creates the graphical output using Matplotlib with improved x-axis formatting."""
    fig, ax = plt.subplots(figsize=(10, 2))
    ax.set_xlim(0, total_length)
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
        elif entry['type'] == 'Ycds':
            ax.add_patch(
                patches.Rectangle(
                    (start, 0.65), width, 0.1,
                    facecolor="red", edgecolor=None, linewidth=0, zorder=3
                )
            )

    # Draw flat line for unclassified sequences in the middle of rectangles
    ax.plot([0, total_length], [0.5, 0.5], color="black", lw=0.5, zorder=1)

    # Set x-axis formatting
    tick_interval_minor = args.tick_minor  # Minor ticks every 5 kb
    tick_interval_major = args.tick_major  # Major labels every 10 kb

    minor_ticks = np.arange(0, total_length + tick_interval_minor, tick_interval_minor)
    major_ticks = np.arange(0, total_length + tick_interval_major, tick_interval_major)

    major_labels = [f"{int(tick/1000)} kb" for tick in major_ticks]

    ax.set_xticks(minor_ticks, minor=True)  # Set minor ticks
    ax.set_xticks(major_ticks)  # Set major ticks
    ax.set_xticklabels(major_labels, fontsize=8)  # Reduce font size for readability

    # Customize tick appearance
    ax.tick_params(axis="x", which="minor", length=4, width=0.5)  # Smaller minor ticks
    ax.tick_params(axis="x", which="major", length=6, width=1.0)  # Larger major ticks

    # Remove unnecessary spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    ax.set_yticks([])
    # ax.set_xlabel("Position")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process Censor output and classify repeats.")
    parser.add_argument("--censor", required=True, help="Path to the Censor output file.")
    parser.add_argument("--len", type=int, required=True, help="Total length of the FASTA sequence.")
    parser.add_argument("--output", required=True, help="Output PNG/PDF file for visualization.")
    parser.add_argument("--len_graph", type=int, default=0,  help="Total length of the FASTA sequence for the graph.")
    parser.add_argument("--tick_minor", type=int, default=5000,  help="tick_interval_minor (unlabelled marks) for the graph.")
    parser.add_argument("--tick_major", type=int, default=10000,  help="tick_interval_major (labelled marks) for the graph.")    
    args = parser.parse_args()
    if args.len_graph == 0:
        args.len_graph = args.len

    classifications, satellite_summary = process_censor_output(args.censor, args.len)
    summarize_categories(classifications, args.len, satellite_summary)
    create_visualization(classifications, args.len_graph, args.output)


