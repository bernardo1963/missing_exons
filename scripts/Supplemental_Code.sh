# Supplemental_Code.sh			version 30Aug2025
# Essential computing codes to generate Figures, Tables, and Results from the ms.      
# Strong bias in long-read sequencing prevents assembly of Drosophila melanogaster Y-linked genes
# A. Bernardo Carvalho, Bernard Y. Kim, and Fabiana Uno

# Related programs, scripts and data files are available at GitHub ( https://github.com/bernardo1963/missing_exons )

step=Figure 1 24aug2025
{

ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/intron_locations/kl3/kl3_RC_exon_boundaries.txt
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/intron_locations/kl3/kl3_RC_cds.fasta

# PacBio CLR
read_coverage_CDS_v6.py    --blast_pgm wu --db iso1_PacBio_CLR_WU  --fasta_file kl3_RC_cds.fasta --identity_cutoff 75 --window_size 12   --figsize_W 3.5 --figsize_H 1.75 --Xlabel "kl-3 CDS position (bp)"   --output  Fig_1A_25aug2025.svg  --exon_fontsize 8  --axes_fontsize 8 --label_fontsize 10 --exon_linewd 1.2 --skip_even_exons_labels yes  --Ylabel "PacBio CLR coverage"   

# ONT (> 1kb)
read_coverage_CDS_v6.py --blast_pgm wu --db 282_929_930_1k_porechop_WU --fasta_file kl3_RC_cds.fasta --identity_cutoff 90 --window_size 12  --figsize_W 3.5 --figsize_H 1.75 --Xlabel "kl-3 CDS position (bp)"   --output  Fig_1B_25aug2025.svg  --exon_fontsize 8  --axes_fontsize 8 --label_fontsize 10 --exon_linewd 1.2 --skip_even_exons_labels yes  --Ylabel "ONT Q20+ coverage"   

# HiFi
read_coverage_CDS_v6.py --blast_pgm wu --db SRR29479668_WU --fasta_file kl3_RC_cds.fasta --identity_cutoff 90 --window_size 12    --figsize_W 3.5 --figsize_H 1.75 --Xlabel "kl-3 CDS position (bp)"   --output  Fig_1C_25aug2025.svg  --exon_fontsize 8  --axes_fontsize 8 --label_fontsize 10 --exon_linewd 1.2 --skip_even_exons_labels yes  --Ylabel "PacBio HiFi coverage"   

# Illumina 
read_coverage_CDS_v6.py --blast_pgm wu --db iso1_m_2024_WU  --fasta_file kl3_RC_cds.fasta --identity_cutoff 90 --window_size 12   --figsize_W 3.5 --figsize_H 1.75 --Xlabel "kl-3 CDS position (bp)"  --cds_file kl3_RC_exon_boundaries.txt  --output  Fig_1D_25aug2025.svg  --exon_fontsize 8  --axes_fontsize 8 --label_fontsize 10 --exon_linewd 1.2 --skip_even_exons_labels yes  --Ylabel "Illumina coverage"   


rm    CDS_coverage_systat*


}


step=Figure 2 25aug2025
{
# kl-3
read_coverage_CDS_v6.py --blast_pgm wu --db 282_929_930_1k_porechop_WU --fasta_file kl3_RC_cds.fasta --identity_cutoff 90 --window_size 12  --Ymax 140 --cds_Yaxis 85    --figsize_W 3.5 --figsize_H 1.75 --Xlabel "kl-3 CDS position (bp)"  --cds_file kl3_RC_exon_boundaries.txt  --output  Fig2_kl3_25aug2025.svg  --exon_fontsize 8  --axes_fontsize 8 --label_fontsize 10 --exon_linewd 1.2 --skip_even_exons_labels yes  

# kl-5
read_coverage_CDS_v6.py --blast_pgm wu --db 282_929_930_1k_porechop_WU --fasta_file kl5_cds.fasta --identity_cutoff 90 --window_size 12  --Ymax 180 --cds_Yaxis 85    --figsize_W 3.5 --figsize_H 1.75 --Xlabel "kl-5 CDS position (bp)"  --cds_file kl5_exon_boundaries.txt  --output  Fig2_kl5_25aug2025.svg  --exon_fontsize 8  --axes_fontsize 8 --label_fontsize 10 --exon_linewd 1.2 --skip_even_exons_labels yes  


# ORY
read_coverage_CDS_v6.py --blast_pgm wu --db 282_929_930_1k_porechop_WU --fasta_file ORY_cds.fasta --identity_cutoff 90 --window_size 12  --Ymax 160 --cds_Yaxis 85    --figsize_W 3.5 --figsize_H 1.75 --Xlabel "ORY CDS position (bp)"  --cds_file ORY_exon_boundaries.txt  --output  Fig2_ORY_25aug2025.svg  --exon_fontsize 8  --axes_fontsize 8 --label_fontsize 10 --exon_linewd 1.2   

# PprY
read_coverage_CDS_v6.py --blast_pgm wu --db 282_929_930_1k_porechop_WU --fasta_file PprY_cds.fasta --identity_cutoff 90 --window_size 12  --Ymax 700 --cds_Yaxis 85    --figsize_W 3.5 --figsize_H 1.75 --Xlabel "Ppr-Y CDS position (bp)"  --cds_file PprY_exon_boundaries.txt  --output  Fig2_PprY_25aug2025.svg  --exon_fontsize 8  --axes_fontsize 8 --label_fontsize 10 --exon_linewd 1.2   

# CCY
read_coverage_CDS_v6.py --blast_pgm wu --db 282_929_930_1k_porechop_WU --fasta_file CCY_cds.fasta --identity_cutoff 90 --window_size 12  --Ymax 90 --cds_Yaxis 85    --figsize_W 3.5 --figsize_H 1.75 --Xlabel "CCY CDS position (bp)"  --cds_file CCY_exon_boundaries.txt  --output  Fig2_CCY_25aug2025.svg  --exon_fontsize 8  --axes_fontsize 8 --label_fontsize 10 --exon_linewd 1.2   

# Pp1Y1
read_coverage_CDS_v6.py --blast_pgm wu --db 282_929_930_1k_porechop_WU --fasta_file Pp1-Y1_cds.fasta --identity_cutoff 90 --window_size 12  --Ymax 200    --figsize_W 3.5 --figsize_H 1.75 --Xlabel "Pp1-Y1 CDS position (bp)"    --output  Fig2_Pp1Y1_25aug2025.svg  --exon_fontsize 8  --axes_fontsize 8 --label_fontsize 10 --exon_linewd 1.2   

# Pp1Y2
read_coverage_CDS_v6.py --blast_pgm wu --db 282_929_930_1k_porechop_WU --fasta_file Pp1-Y2_cds.fasta --identity_cutoff 90 --window_size 12  --Ymax 150    --figsize_W 3.5 --figsize_H 1.75 --Xlabel "Pp1-Y2 CDS position (bp)"    --output  Fig2_Pp1Y2_25aug2025.svg  --exon_fontsize 8  --axes_fontsize 8 --label_fontsize 10 --exon_linewd 1.2   


rm    CDS_coverage_systat*
}

step=Figure 3 26aug2025
{
# Release6 is a composite including RT-PCR sequence, so it is unreliable to get exon_boundaores etc
# We obtained the approximate exon boundaries by doing a blastN search of the CDS againd a database of Nanopore reads


# I will use the previously prepared map files (produced by the censor program) .  
# getting the censor map links:

ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/figures/fig3/kl5_exon13_7feb2025.racon.fasta.map
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/figures/fig3/white_ctg72_flyehq2_19216_19316.fasta.map
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/figures/fig3/PprY_exon4_ctg166_flyehq2_130_245.fasta.map
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/figures/fig3/PprY_exon3_10jan2025.racon.fasta.map
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/figures/fig3/Pp1Y2_ctg246_14jun24_128_228.fasta.map
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/figures/fig3/Pp1Y1_ctg239_14jun24_350_450.fasta.map
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/figures/fig3/ORY_exon1_2_10jan2025.racon.fasta.map
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/figures/fig3/kl5_exon3_9_9jan2025.racon.fasta.map
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/figures/fig3/kl5_exon13_8jan2025.racon.fasta.map   #   replaced by a newer one 7feb2025
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/figures/fig3/kl5_exon10_12_ctg33_hq2.fasta.map
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/figures/fig3/kl3_region1_10jan2025.racon.fasta.map
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/figures/fig3/CCYexon2_13reads_mmplain1_stage1.racon1.fasta.map

for gene in CCY kl-3 kl-5 ORY Pp1-Y1 Pp1-Y2 Ppr-Y
do
  gene_simple="${gene//-/}"
  grep -A 1 ${gene} /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/reviewer2/other_genes/dmel-all-CDS-r6.63.edited2.nonredundant.fasta | sed '/>/{s/$/_cds/}' > ${gene_simple}_cds.fasta
done

# white 
grep -A 1 ">w-R"  /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/reviewer2/other_genes/dmel-all-CDS-r6.63.edited2.nonredundant.fasta | sed '/>/{s/$/_cds/}' > white_cds.fasta

# linking the exon_boundaries files produced on 10jun2025
ln -s ~/projects/mel_Y/assembly_evaluation/XY_coverage/intron_locations/*/*_exon_boundaries.txt  ./

# exon boundaries files:
CCY_exon_boundaries.txt
kl2_exon_boundaries.txt
kl3_RC_exon_boundaries.txt
kl5_exon_boundaries.txt
ORY_exon_boundaries.txt
PprY_exon_boundaries.txt


# CCY_exon2
visualize_repeats_censor5.py --censor  CCYexon2_13reads_mmplain1_stage1.racon1.fasta.map --min_repeat_len  0 --gene_color black --cds_file CCY_exon2_locations_gene_region.txt --len_graph 100001 --len 98042 --figsize_W 4.9 --figsize_H_nodepth 1.15  --exon_fontsize 8 --scale_fontsize 8 --no_background --output Fig3_CCY_exon2_27ago2025.svg 

# kl5_exon14
visualize_repeats_censor5.py --censor kl5_exon13_7feb2025.racon.fasta.map --min_repeat_len  0  --gene_color black --cds_file kl5_exon14_locations_gene_region.txt --len_graph 100001 --len 32405 --figsize_W 4.9 --figsize_H_nodepth 1.15  --exon_fontsize 8 --scale_fontsize 8 --no_scale  --no_background  --output Fig3_kl5_exon14_27ago2025.svg 


# kl5_exon3_10
visualize_repeats_censor5.py --censor kl5_exon3_9_9jan2025.racon.fasta.map --min_repeat_len  0  --gene_color black --cds_file kl5_exon_3_10_locations_gene_region.txt --len_graph 100001 --len 12473 --figsize_W 4.9 --figsize_H_nodepth 1.15  --exon_fontsize 8 --scale_fontsize 8 --no_scale  --no_background  --output Fig3_kl5_exon3_10_27ago2025.svg 

# kl3_exon15_16
visualize_repeats_censor5.py --censor kl3_region1_10jan2025.racon.fasta.map --min_repeat_len  0  --gene_color black --cds_file kl3_exon_15_16_locations_gene_region.simple2.txt --len_graph 100001 --len 28613 --figsize_W 4.9 --figsize_H_nodepth 1.15  --exon_fontsize 8 --scale_fontsize 8 --no_scale  --no_background  --output Fig3_kl3_exon15_16_27ago2025.svg 

# ORY_exon1_2
visualize_repeats_censor5.py --censor ORY_exon1_2_10jan2025.racon.fasta.map --min_repeat_len  0  --gene_color black --cds_file ORY_exon1_2_locations_gene_region.simple2.txt --len_graph 100001 --len 53956 --figsize_W 4.9 --figsize_H_nodepth 1.15  --exon_fontsize 8 --scale_fontsize 8 --no_scale  --no_background  --output Fig3_ORY_exon1_2_27ago2025.svg 


# PprY_exon3
visualize_repeats_censor5.py --censor PprY_exon3_10jan2025.racon.fasta.map --min_repeat_len  0  --gene_color black --cds_file PprY_exon3_locations_gene_region.txt   --len 38049 --len_graph 100001   --figsize_W 4.9 --figsize_H_nodepth 1.15  --exon_fontsize 8 --scale_fontsize 8 --no_scale  --no_background  --output Fig3_PprY_exon3_27ago2025.svg 


# Pp1Y1
visualize_repeats_censor5.py --censor Pp1Y1_ctg239_14jun24_350_450.fasta.map --min_repeat_len  0  --gene_color black --cds_file Pp1Y1_exon_locations_gene_region.txt  --len 100001   --figsize_W 4.9 --figsize_H_nodepth 1.15  --exon_fontsize 8 --scale_fontsize 8 --no_scale  --no_background  --output Fig3_Pp1Y1_27ago2025.svg 

# Pp1Y2
visualize_repeats_censor5.py --censor Pp1Y2_ctg246_14jun24_128_228.fasta.map --min_repeat_len  0  --gene_regex Ycds --gene_color black --cds_file Pp1Y2_exon_locations_gene_region.txt --len 100001   --figsize_W 4.9 --figsize_H_nodepth 1.15  --exon_fontsize 8 --scale_fontsize 8 --no_scale  --no_background  --output Fig3_Pp1Y2_27ago2025.svg 

# PprY_exon4
visualize_repeats_censor5.py --censor PprY_exon4_ctg166_flyehq2_130_245.fasta.map --min_repeat_len  0  --gene_color black --cds_file PprY_exon4_locations_gene_region.simple2.txt   --len 115001   --figsize_W 4.9 --figsize_H_nodepth 1.15  --exon_fontsize 8 --scale_fontsize 8 --no_scale  --no_background  --output Fig3_PprY_exon4_27ago2025.svg 


# kl5_exon11_13
visualize_repeats_censor5.py --censor kl5_exon10_12_ctg33_hq2.fasta.map --min_repeat_len  0  --gene_color black --cds_file kl5_exon_11_13_locations_gene_region.txt   --len 132946   --figsize_W 4.9 --figsize_H_nodepth 1.15  --exon_fontsize 8 --scale_fontsize 8 --no_scale  --no_background  --output Fig3_kl5_exon11_13_27ago2025.svg 

# white
visualize_repeats_censor5.py --censor white_ctg72_flyehq2_19216_19316.fasta.map --min_repeat_len  0  --gene_color black --cds_file white_exon_locations_gene_region.txt   --len 100001   --figsize_W 4.9 --figsize_H_nodepth 1.15  --exon_fontsize 8 --scale_fontsize 8 --no_scale  --no_background  --output Fig3_white_27ago2025.svg 


}


step=Figure 4 18jun2025
{
conda activate /home6/bernardo/my_pandas_env
box_violin_plot.py  --input_file Fig4_11jun2025.data --plot_type box --gene_order CCY_exon_2,kl-5_exon_14,kl-5_exons_3-10,kl-3_exons_15-16,ORY_exons_1-2,Ppr-Y_exon_3,Pp1-Y1,Pp1-Y2,Ppr-Y_exon_4,kl-5_exons_11-13,white --outliers n --Ymax 100000  --Xfontsize 16  --Yfontsize 16 --output_file Fig4_18jun2025.svg   

}


step=Figure 5  20jun2025
{
# Figure 5 was produced in SYSTAT, running the script   Fig5_graph_20jun2025.syc in SYSTAT, using  data file  Fig5_data.txt

}

step=Figure 6  22feb2025
{
fplot_reads.awk   Fig6_data.txt  # produced file  Fig6_data.svg . This will print the arrows, which represent the reads
# the central part showing the satellites, exon, and TEs:
visualize_repeats_censor5.py --censor PprY_exon3_10jan2025.racon.fasta.map --min_repeat_len  0  --gene_color red --cds_file PprY_exon3_locations_gene_region.txt   --len 38049 --len_graph 40001   --figsize_W 7 --figsize_H_nodepth 1.6  --exon_fontsize 10 --scale_fontsize 10 --no_background  --output Fig6_PprY_exon3_7inches.svg 

}


step=Figure 7  14jun2025
{

cd ~/projects/mel_Y/assembly_evaluation/XY_coverage/reviewer1/non-B/Y_regions

ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/assembly/CCY_exon2/CCYexon2_13reads_mmplain1_stage1.racon1.fasta ./CCY_exon2_region
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/assembly/kl5_exon13/kl5_exon13_7feb2025.racon.fasta     ./kl5_exon14_region
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/assembly/kl5_exon10_12/ctg33_hq2.fasta                  ./kl5_exon11_13_region
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/assembly/PprY_exon4/ctg166_flyehq2_130_245.fasta        ./PprY_exon4_region   
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/assembly/Pp1Y2/ctg246_14jun24_128_228.fasta             ./Pp1Y2_region
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/assembly/Pp1Y1/ctg239_14jun24_350_450.fasta             ./Pp1Y1_region  
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/assembly/PprY_exon3/PprY_exon3_10jan2025.racon.fasta    ./PprY_exon3_region
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/assembly/ORY_exon1_2/ORY_exon1_2_10jan2025.racon.fasta  ./ORY_exon1_2_region 
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/assembly/kl3_region1/kl3_region1_10jan2025.racon.fasta  ./kl3_exon15_16_region
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/assembly/kl5_exon3_9/kl5_exon3_9_9jan2025.racon.fasta   ./kl5_exon3_10_region
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/assembly/white/white_ctg72_flyehq2_19216_19316.fasta    ./white_region

# editing deflines
awk '/^>/{print ">" FILENAME}; !/^>/{print $0}' *_region  > Y_regions.fasta
rm *_region

# gfa is the command line version of the nBMST program
gfa    -seq  Y_regions.fasta  -out  Y_regions.out   

process_nonB_gfa.awk  output=raw  Y_regions.out_*.tsv > nonB_Yregions.txt


}



step= Table 1 
{ 
missingExon_stat_1jan2025.py  --m8mod    Fig5_data.txt > Table1_raw.txt
}

step= Table 2 24jun2025 
{ 
cat  lowcov_start_seq_24jun2025.data
# CCY   0 13       69.2
# kl5_exon_3-10   0 12      37.5
# kl5_exon14   0 12        75.0
# kl3_region1   0 16       85.0
# ORY_exon_1-2   1 15      86.5
# PprY   0 18      18.7

missingExon_binomial_stat_v1.py    lowcov_start_seq_24jun2025.data

# CCY                     0        13       0.0      69.2    0.0000
# kl5_exon_3-10           0        12       0.0      37.5    0.0051
# kl5_exon14              0        12       0.0      75.0    0.0000
# kl3_region1             0        16       0.0      85.0    0.0000
# ORY_exon_1-2            1        15       6.2      86.5    0.0000
# PprY                    0        18       0.0      18.7    0.0350

# count1 i the number of reads starting with satellite DNA; count2 is the number of reads NOT starting with satellite DNA
}


step= Table 3   
{ 

find_tandem_repeats_v2.py    allreads.ORY_exon_1-2.fasta  --min_tandem_copies_block 4   --min_tandem_copies_summary 50 --print_mode summary

# repeat the above script for each fasta file

}




step=Supplemental Table S1  4feb2025 
{
for db in mel_fail_400   mel_fail_260    mel_fail_Pro
do 
  blastall -p blastn -i Fig4_targets_11jun2025.fasta -d /draft10/bernardo/BKim_nanopore_data/db/$db -a 100 -e 0.00001 -F "m D" -m 8 >  Fig4targets_$db.m8
  echo $db " done"
done

# kl-5_exons_3-10 gi|169348       97.37   38      1       0       507     544     378     415     4e-08   67.9
# kl-3_exons_15-16        gi|237448       84.35   313     27      16      2474    2781    169     464     8e-23    117
# ORY_exons_1-2   gi|313488       99.15   236     1       1       277     511     24968   24733   5e-122   444

# combining the three serches, removing redundat reads (when the blast aln is split in two)
awk '{array[$1,$2]++} ; (array[$1,$2]==1){hit_array[$1]++} ; END{ for(hit in hit_array){print hit, hit_array[hit]}}'  Fig4targets_*.m8

# arranging , addiing the genes w/ 0 hits,  making a table
CCY_exon_2	0
kl-5_exon_14 2
kl-5_exons_3-10 1
kl-3_exons_15-16 1
ORY_exons_1-2 2
Ppr-Y_exon_3 2
Pp1-Y1 12
Pp1-Y2 4
Ppr-Y_exon_4 24
kl-5_exons_11-13 9
white 39

}

step=Supplemental Table S2 5feb2025
{ 
missingExon_stat_1jan2025.py  --m8mod     FigS6_data_simulated.txt > TableS2_raw.txt

# part 1 of Table 1: F/R reads usage:  [not useful]
# part 2 of Table 1: Anderson-Dalring test of read starts (simulated)

}

step=Supplemental Table S3 
{ 
# same code of Table S1
}

step=Supplemental Table S4  
{
# We then ran the shell script Ycompleteness_v3.sh, which builds a blast database for the target assembly, runs megablast of the CDS against the formatted database, and outputs for each CDS how much of its sequence is present in the assembly. The command line we used was:

Ycompleteness_v3.sh   -p " -W 20 -F F " -P 98 -a 28 -t 100  -s yes -o no -r dmel-all-CDS-r6.63.edited2.nonredundant.fasta   HiFi_hifiasmL0_p.fasta > completeness_HiFi_hifiasmL0p_a28_FF.txt 

# We then extract only the cases missing more than 50 bp, to avoid too many false positives. In the intial runs we found several spurious hits which upon close examination were found to be either mitochondrial genes (which assemble poorly possibly due to too high coverage: 2,000 × for HiFi; 7,000 × for Nanopore), or genes with repetitive exons that  confound blast (several mucins, Sgs1, etc.). We blacklisted them as follows: 

awk '(NF==5)&&(($3-$4)>50){print $0 "\t\t" ($3-$4) }'   completeness_HiFi_hifiasmL0p_a28_FF.txt | grep  -v  -f completeness_FBallgenes_a28_FF_blacklist.txt

# The above scripts were run for each of the four assemblies, and the outputs were combined in Supplemental Table S4.
}



step=Supplemental Figure S3
{
read_coverage_CDS_v6.py --blast_pgm wu --db 282_929_930_1k_porechop_WU --fasta_file 10Xgenes.fasta --identity_cutoff 98 --window_size 12 --graph_name_suffix ONT_1kb --figsize_W 6 --figsize_H 3 --Xlabel gene_name  --Ymax 300

rm CDS_coverage*
}




step=Supplemental Figure S21
{
# below is the example fopr one gene (Syt12)
# download the Syt12 gene region (plus 2 kb), and Syt12 mRNA  from FlyBase.


# adding gene mRNA to censor library [in not already there] 
sed  '/^[ATGC]/{y/ATGCN/atgcn/}'  Syt12_mRNA.fasta  >>  /home6/tools/miniconda3/envs/censor/share/censor-4.2.31/biolib/genes_with_sat_mRNArep.ref
 
#getting exon locations
bl2seq -p blastn -i Syt12_mRNA.fasta  -j Syt12_gene.fasta -F F -e 0.001 -W 28 -D 2 | sort -k 7,7n | sed '/#/d' | awk '{print $0 "\t" NR}'
Syt12-RB_mRNA   Syt12_gene      100.00  190     0       0       1       190     1       190     1e-106   377    1
Syt12-RB_mRNA   Syt12_gene      100.00  240     0       0       188     427     3994    4233    2e-136   476    2
Syt12-RB_mRNA   Syt12_gene      100.00  1290    0       0       428     1717    4337    5626    0.0     2557    3
Syt12-RB_mRNA   Syt12_gene      100.00  148     0       0       1717    1864    5684    5831    1e-81    293    4
Syt12-RB_mRNA   Syt12_gene      100.00  358     0       0       1860    2217    5891    6248    0.0      710    5
Syt12-RB_mRNA   Syt12_gene      100.00  2609    0       0       2212    4820    6324    8932    0.0     5172    6

#getting exon boundaries file for read_coverage_CDS_v6.py  
bl2seq -p blastn -i Syt12_mRNA.fasta  -j Syt12_gene.fasta -F F -e 0.001 -W 28 -D 2 | sort -k 7,7n | sed '/#/d' | awk '{print $7 "\t" $8 "\t" NR}' > Syt12_exon_boundaries.txt
#getting exon locations in the gene region, for   visualize_repeats_censor4_abc.py     
bl2seq -p blastn -i Syt12_mRNA.fasta  -j Syt12_gene.fasta -F F -e 0.001 -W 28 -D 2 | sort -k 7,7n | sed '/#/d' | awk '{print $9 "\t" $10 "\t" NR}' > Syt12_exon_locations_gene_region.txt	

conda activate censor
censor  Syt12_gene.fasta -lib GK1688genes -lib GK1688 -lib dro  -mode norm   -bprm 'cpus=100'  -bprm '-filter=none'  -show_simple -nomasked   -nofound
rm *.idx  *.log *.aln

awk '($4 ~/Syt12|1\.688_GK/) || (($3-$2+1)>=100){print $0 "\t" $3-$2+1 " bp"}'  Syt12_gene.fasta.map   
Syt12_gene         1 190         Syt12-RB       1 190     d   1.0000 99.0000 1791       2.13    3.94    190 bp
Syt12_gene      1500 1709    1.688_GK_27O05-B2/4-X-11F-12A       1 212     c   0.8957 1.9091 1429       2.35    58.89   210 bp
Syt12_gene      1710 2069    1.688_GK_27O05-B2/5-X-11F-12A       1 360     c   1.0000 99.0000 3192      4.03    100.00  360 bp
Syt12_gene      2070 2429    1.688_GK_27O05-B2/4-X-11F-12A       1 360     c   1.0000 99.0000 3195      4.03    100.00  360 bp
Syt12_gene      2430 2789    1.688_GK_27O05-B2/3-X-11F-12A       1 360     c   1.0000 99.0000 3207      4.03    100.00  360 bp
Syt12_gene      2790 3149    1.688_GK_27O05-B2/2-X-11F-12A       1 360     c   1.0000 99.0000 3198      4.03    100.00  360 bp
Syt12_gene      3152 3318    1.688_GK_02B12/2-X-13E-13E     192 359     c   0.8242 1.8462 759   1.87    46.41   167 bp
Syt12_gene      3994 4236        Syt12-RB     188 430     d   0.9959 1.0000 2294        2.72    5.04    243 bp
Syt12_gene      4337 5831        Syt12-RB     428 1864    d   0.9993 99.0000 13455      16.74   29.81   1495 bp
Syt12_gene      5891 8932        Syt12-RB    1860 4820    d   0.9997 99.0000 27502      34.06   61.43   3042 bp
# largest 1.6888 block: 1817

# below generates the figure with repetitive elements and exon locations
visualize_repeats_censor4.py --censor  Syt12_gene.fasta.map   --repeat_regex "1.688_GK" --gene_regex Syt12  --gene_color black    --len 8932 --output  Syt12_gene.png   --min_repeat_len 100  --cds_file Syt12_exon_locations_gene_region.txt 



# read  coverage
# HiFi 
read_coverage_CDS_v6.py    --blast_pgm wu --db SRR29479668_WU  --fasta_file Syt12_mRNA.fasta  --identity_cutoff 90 --window_size 12 --graph_name_suffix HiFi  --figsize_W 6 --figsize_H 3  --Xlabel gene_name --cds_file Syt12_exon_boundaries.txt
# ONT (> 1kb)
read_coverage_CDS_v6.py    --blast_pgm wu --db 282_929_930_1k_porechop_WU  --fasta_file Syt12_mRNA.fasta   --identity_cutoff 90 --window_size 12 --graph_name_suffix ONT_1kb  --figsize_W 6 --figsize_H 3   --Xlabel gene_name --cds_file Syt12_exon_boundaries.txt
rm *.m8  CDS_coverage_systat*  debug_log.txt 


}

step=Supplemental Figure S22
{
# same code of Figure 7
}
