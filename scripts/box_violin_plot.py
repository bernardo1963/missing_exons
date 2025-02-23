#!/usr/bin/env python
#   Fabiana Uno   v. 22feb2025  
# This script generates box plots or violin plots for read size distributions across different genes using matplotlib
#

import argparse
import pandas as pd

import matplotlib as mpl
# https://stackoverflow.com/questions/50506076/is-there-an-efficient-way-to-store-2d-plots-as-a-vector-graphic-in-python
mpl.use('svg')
new_rc_params = {"svg.fonttype": 'none'} #to store text as text, not as path
mpl.rcParams.update(new_rc_params)
import matplotlib.pyplot as plt
import seaborn as sns

def plot_read_size(input_file, plot_type, output_file, gene_order, show_outliers, Ymax):
    # Load data
    df = pd.read_csv(input_file, sep=r"\s+", engine="python")
    
    # Order genes if provided (otherwise use alphabetical order)
    if gene_order:
        gene_order = gene_order.split(',')
        df['gene$'] = pd.Categorical(df['gene$'], categories=gene_order, ordered=True)
    else:
        df['gene$'] = pd.Categorical(df['gene$'], ordered=True)
    
    # Filter for specific types
    # df = df[(df['type$'] == "lowcov") | (df['type$'] == "Ychrom")]
    
    # Plot
    plt.figure(figsize=(14, 7))  # Slightly larger for better readability
    
    # seaborn style (aesthetics)
    # sns.set_style("whitegrid")
    sns.set_style("ticks")
    set2_colors = sns.color_palette("Set2") 
    blue = set2_colors[2] 
    green = set2_colors[0]

    # Determine if outliers should be shown
    showfliers = True if show_outliers.lower() == "y" else False
    
    if plot_type == "box":
        # sns.boxplot(x='gene$', y='size', data=df, order=gene_order if gene_order else None, showfliers=showfliers, linewidth=1.2, width=0.6, color=blue)
        sns.boxplot(x='gene$', y='size', data=df, order=gene_order if gene_order else None, showfliers=showfliers, linewidth=1.2, width=0.6,  fill=False)
    elif plot_type == "violin":
        sns.violinplot(x='gene$', y='size', data=df, order=gene_order if gene_order else None, inner="quartile", linewidth=1.2, color=green)
    else:
        print("Invalid plot type. Use 'box' or 'violin'.")
        return
    
    # setting labels and titles
    plt.xticks(rotation=45, ha='right', fontsize=args.Xfontsize)
    plt.yticks(fontsize=args.Yfontsize)
    plt.xlabel("Gene", fontsize=(args.Yfontsize +2), fontweight='bold')
    plt.ylabel("Read Size", fontsize=(args.Yfontsize +2), fontweight='bold')
    
    # Set Y-axis limit
    if Ymax:
        plt.ylim(0, Ymax)
    
    plt.tight_layout()
    
    # Save plot
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    # print(f"Plot saved as {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate box or violin plot for read sizes.")
    parser.add_argument("--input_file", help="Path to input tab-separated file")
    parser.add_argument("--plot_type", choices=["box", "violin"], help="Type of plot: 'box' or 'violin'")
    parser.add_argument("--output_file", help="Path to save the output image")
    parser.add_argument("--gene_order", help="Comma-separated list of gene names in desired order", default=None)
    parser.add_argument("--outliers", choices=["y", "n"], default="y", help="Show outliers in box plots (y/n)")
    parser.add_argument("--Ymax", type=float, help="Maximum value for Y-axis", default=None)
    parser.add_argument("--Xfontsize", type=float, default=12, help="font size X axis")
    parser.add_argument("--Yfontsize", type=float, default=12, help="font size Y axis")    
    args = parser.parse_args()
    plot_read_size(args.input_file, args.plot_type, args.output_file, args.gene_order, args.outliers, args.Ymax)
