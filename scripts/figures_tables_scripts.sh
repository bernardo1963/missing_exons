# figures_tables_scripts.sh     v. 22feb2025
# command lines used to produce the figures and tables for the manuscript "Strong sequencing bias in Nanopore and PacBio prevents assembly of Drosophila melanogaster Y-linked gene"

***************************************************   Figure 1 *****************************************************
# Fig 1a: PacBio CLR (Kim et al 2014)
read_coverage_CDS_v4.py    --blast_pgm wu --db iso1_PacBio_CLR_WU  --fasta_file mel_kl3_CDS.fasta  --systat_output_file kl3_PacBioCLR_covWU_i75_26jan2025.txt --identity_cutoff 75 --window_size 12 --graph_name_suffix PacBioCLR --figsize_W 6 --figsize_H 3

# Fig 1b: ONT Q20  (Kim et al 2024)
read_coverage_CDS_v4.py    --blast_pgm wu --db 282_929_930_1k_porechop_WU  --fasta_file mel_kl3_CDS.fasta  --systat_output_file kl3_ONT_covWU_i90_26jan2025.txt --identity_cutoff 90 --window_size 12 --graph_name_suffix ONT --figsize_W 6 --figsize_H 3

# Fig 1c:  PacBio HiFi (Shukla et al 2024)
read_coverage_CDS_v4.py    --blast_pgm wu --db SRR29479668_WU  --fasta_file mel_kl3_CDS.fasta  --systat_output_file kl3_HiFi_covWU_i90_26jan2025.txt --identity_cutoff 90 --window_size 12 --graph_name_suffix HiFi --figsize_W 6 --figsize_H 3

# Fig 1d:  Illumina 2024 (ABC)
read_coverage_CDS_v4.py    --blast_pgm wu --db iso1_m_2024_WU  --fasta_file mel_kl3_CDS.fasta  --systat_output_file kl3_Illumina_covWU_i90_26jan2025.txt --identity_cutoff 90 --window_size 12 --graph_name_suffix Illumina --figsize_W 6 --figsize_H 3



***************************************************   Figure 2 *****************************************************
read_coverage_CDS_v4.py    --blast_pgm wu --db 282_929_930_1k_porechop_WU  --fasta_file Ycds2.fasta  --systat_output_file Ycds2_ONT_covWU_i90_26jan2025.txt --identity_cutoff 90 --window_size 12 --graph_name_suffix ONT  --figsize_W 6 --figsize_H 3


***************************************************   Figure 3 *****************************************************
# shown for one region:
censor CCYexon2_13reads_mmplain1_stage1.racon1.fasta -lib Ycds -lib dro -lib inv   -mode sens   -bprm 'cpus=100'  -bprm '-filter=none'  -show_simple  -nomasked   -nofound
visualize_repeats_censor.py --censor ./CCYexon2_13reads_mmplain1_stage1.racon1.fasta.map --len 98042 --output CCYexon2_3jan2024.racon.7feb2025.png 

***************************************************   Figure 4 *****************************************************
box_violin_plot.py  --input_file Fig4_data.txt --plot_type box --gene_order CCY_exon_2,kl-5_exon_13,kl-5_exons_3-9,kl-3_exons_12-13,ORY_exons_1-2,Ppr-Y_exon_3,Pp1-Y1,Pp1-Y2,Ppr-Y_exon_4,kl-5_exons_10-12 --outliers n --Ymax 100000  --Xfontsize 16  --Yfontsize 16 --output_file Fig4_22feb2025.svg 

***************************************************   Figure 5 *****************************************************
# produced with SYSTAT command file Fig5_graph.syc , reading data file Fig5_data.txt

***************************************************   Figure 6 *****************************************************
fplot_reads.awk   Fig6_data.txt  # produces file  Fig6_data.svg 


***************************************************   Table 1 *****************************************************
missingExon_stat_1jan2025.py  --m8mod    Fig5_data.txt > Table1_raw.txt