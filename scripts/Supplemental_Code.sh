# Supplemental_Code.sh			version 24jun2025
# Essential computing codes to generate Figures, Tables, and Results from the ms.      
# Strong bias in long-read sequencing prevents assembly of Drosophila melanogaster Y-linked genes
# A. Bernardo Carvalho, Bernard Y. Kim, and Fabiana Uno

# Related programs, scripts and data files are available at GitHub ( https://github.com/bernardo1963/missing_exons )

step=Figure 1 19jun2025
{

ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/intron_locations/kl3/kl3_RC_exon_boundaries.txt
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/intron_locations/kl3/kl3_RC_cds.fasta


# PacBio CLR
read_coverage_CDS_v6.py    --blast_pgm wu --db iso1_PacBio_CLR_WU  --fasta_file kl3_RC_cds.fasta --identity_cutoff 75 --window_size 12 --graph_name_suffix PacBioCLR_19jun2025 --figsize_W 6 --figsize_H 3 --Xlabel "CDS position (bp)"   
# ONT (> 1kb)
read_coverage_CDS_v6.py --blast_pgm wu --db 282_929_930_1k_porechop_WU --fasta_file kl3_RC_cds.fasta --identity_cutoff 90 --window_size 12 --graph_name_suffix ONT_19jun2025 --figsize_W 6 --figsize_H 3 --Xlabel "CDS position (bp)"
# HiFi
read_coverage_CDS_v6.py --blast_pgm wu --db SRR29479668_WU --fasta_file kl3_RC_cds.fasta --identity_cutoff 90 --window_size 12 --graph_name_suffix HiFi_19jun2025 --figsize_W 6 --figsize_H 3 --Xlabel "CDS position (bp)"
# Illumina  # will need to edit the lables of some exons (too close)
read_coverage_CDS_v6.py --blast_pgm wu --db iso1_m_2024_WU --fasta_file kl3_RC_cds.fasta --identity_cutoff 90 --window_size 12 --graph_name_suffix Illumina_19jun2025 --figsize_W 6 --figsize_H 3 --Xlabel "CDS position (bp)"  --cds_file kl3_RC_exon_boundaries.txt  --output  kl-3-RC_cds_coverage_smoothed12_Illumina_19jun2025_edited.svg 

rm *.idx  *.log *.aln  debug_log.txt CDS_coverage_systat*

}


step=Figure 2 19jun2025
{
#######################   Figure 2 kl3 ##################

cd /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/reviewer3/point0_Fig2/
ln -s ~/projects/mel_Y/assembly_evaluation/XY_coverage/intron_locations/kl3/kl3_RC_cds.fasta
ln -s ~/projects/mel_Y/assembly_evaluation/XY_coverage/intron_locations/kl3/kl3_RC_exon_boundaries.txt

# attempt2: increase Ymax, --cds_Yaxis to put exons at the top of the figure , svg output (to add blue/black arrows)
read_coverage_CDS_v6.py --blast_pgm wu --db 282_929_930_1k_porechop_WU --fasta_file kl3_RC_cds.fasta --identity_cutoff 90 --window_size 12 --graph_name_suffix ONT --figsize_W 6 --figsize_H 3 --Xlabel "CDS position (bp)"  --cds_file kl3_RC_exon_boundaries.txt  --Ymax 130 --cds_Yaxis 90  --output kl3-RC_cds_coverage_smoothed12_ONT_19jun2025_edited.svg
rm *.idx  *.log *.aln  debug_log.txt CDS_coverage_systat*


#######################   Figure 2 kl5 ##################

ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/intron_locations/kl5/kl5_exon_boundaries.txt
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/reviewer1/point2_Fig3/kl5_cds.fasta

# attempt2: increase Ymax, --cds_Yaxis to put exons at the top of the figure ,will need Illkustrrator edition (save in svg) 
read_coverage_CDS_v6.py --blast_pgm wu --db 282_929_930_1k_porechop_WU --fasta_file kl5_cds.fasta --identity_cutoff 90 --window_size 12 --graph_name_suffix ONT --figsize_W 6 --figsize_H 3 --Xlabel "CDS position (bp)"  --cds_file kl5_exon_boundaries.txt  --Ymax 170 --cds_Yaxis 90  --output kl5_cds_coverage_smoothed12_ONT_19jun2025_edited.svg
rm *.idx  *.log *.aln  debug_log.txt CDS_coverage_systat*


#######################   Figure 2 ORY ##################

ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/intron_locations/ORY/ORY_exon_boundaries.txt
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/reviewer1/point2_Fig3/ORY_cds.fasta


# I will make two variations: exons in the top (as above) and exons in the bottom. 
# top
read_coverage_CDS_v6.py --blast_pgm wu --db 282_929_930_1k_porechop_WU --fasta_file ORY_cds.fasta --identity_cutoff 90 --window_size 12 --graph_name_suffix ONT --figsize_W 6 --figsize_H 3 --Xlabel "CDS position (bp)"  --cds_file ORY_exon_boundaries.txt  --Ymax 150 --cds_Yaxis 90  --output ORY_cds_coverage_smoothed12_ONT_19jun2025_top_edited.svg
# bottom
read_coverage_CDS_v6.py --blast_pgm wu --db 282_929_930_1k_porechop_WU --fasta_file ORY_cds.fasta --identity_cutoff 90 --window_size 12 --graph_name_suffix ONT --figsize_W 6 --figsize_H 3 --Xlabel "CDS position (bp)"  --cds_file ORY_exon_boundaries.txt     --output ORY_cds_coverage_smoothed12_ONT_19jun2025_bottom.png
rm *.idx  *.log *.aln  debug_log.txt CDS_coverage_systat*
# ANSWER: perhaps top is better , to leave all figures with the smae visual lgg

#######################   Figure 2 PprY ##################

ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/intron_locations/PprY/PprY_exon_boundaries.txt
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/reviewer1/point2_Fig3/PprY_cds.fasta
# top
read_coverage_CDS_v6.py --blast_pgm wu --db 282_929_930_1k_porechop_WU --fasta_file PprY_cds.fasta --identity_cutoff 90 --window_size 12 --graph_name_suffix ONT --figsize_W 6 --figsize_H 3 --Xlabel "CDS position (bp)"  --cds_file PprY_exon_boundaries.txt  --Ymax 700 --cds_Yaxis 90  --output PprY_cds_coverage_smoothed12_ONT_19jun2025_edited.svg

#######################   Figure 2 CCY ##################

ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/intron_locations/CCY/CCY_exon_boundaries.txt
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/reviewer1/point2_Fig3/CCY_cds.fasta
# top
read_coverage_CDS_v6.py --blast_pgm wu --db 282_929_930_1k_porechop_WU --fasta_file CCY_cds.fasta --identity_cutoff 90 --window_size 12 --graph_name_suffix ONT --figsize_W 6 --figsize_H 3 --Xlabel "CDS position (bp)"  --cds_file CCY_exon_boundaries.txt  --Ymax 80 --cds_Yaxis 90  --output CCY_cds_coverage_smoothed12_ONT_19jun2025_edited.svg


#######################   Figure 2 Pp1Y1 ##################
ln -s ~/projects/mel_Y/assembly_evaluation/XY_coverage/Pp1Y1/Pp1-Y1_cds.fasta

read_coverage_CDS_v6.py --blast_pgm wu --db 282_929_930_1k_porechop_WU --fasta_file Pp1-Y1_cds.fasta --identity_cutoff 90 --window_size 12 --graph_name_suffix ONT_19jun2025 --figsize_W 6 --figsize_H 3 --Xlabel "CDS position (bp)"    --Ymax 200     


#######################   Figure 2 Pp1Y2 ##################
ln -s /home6/bernardo/projects/mel_Y/assembly_evaluation/XY_coverage/assembly/Pp1Y2/Pp1-Y2_cds.fasta
read_coverage_CDS_v6.py --blast_pgm wu --db 282_929_930_1k_porechop_WU --fasta_file Pp1-Y2_cds.fasta --identity_cutoff 90 --window_size 12 --graph_name_suffix ONT_19jun2025 --figsize_W 6 --figsize_H 3 --Xlabel "CDS position (bp)"    --Ymax 150     


rm *.idx  *.log *.aln  debug_log.txt CDS_coverage_systat*

}

step=Figure 3 10jun2025
{
# Release6 is a composite including RT-PCR sequence, so it is unreliable to get exon_boundaores etc
# We obtained the approximate exon boundaries by doing a blastN search of the CDS againd a database of Nanopore reads


cd ~/projects/mel_Y/assembly_evaluation/XY_coverage/reviewer1/point2_Fig3
# getting the fasta links:







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




##################### CCY  10jun2025  #####################


bl2seq -p blastn -i CCY_cds.fasta  -j CCYexon2_13reads_mmplain1_stage1.racon1.fasta -F F -e 0.001 -W 16 -D 2 | sort -k 7,7n | sed '/#/d' | awk '{print $0 "\t" NR}' 
# CCY-RB_cds      CCYexon2_assembly_racon 100.00  158     0       0       1072    1229    21404   21561   1e-86    313    1
cat  CCY_exon_boundaries.txt
# 1 1071 1
# 1072 1229 2
# 1227 3880 3
# 3880 4002 4

echo "21404   21561  2" > CCY_exon2_locations_gene_region.txt


# run visualize_repeats_censor4.py 

fasta_size.awk CCYexon2_13reads_mmplain1_stage1.racon1.fasta # 98042

visualize_repeats_censor4.py --censor ./CCYexon2_13reads_mmplain1_stage1.racon1.fasta.map --min_repeat_len  0  --gene_regex Ycds --gene_color black --cds_file CCY_exon2_locations_gene_region.txt --len_graph 100001 --len 98042  --output CCY_exon2_gene_region_10jun2025.png  --scale_fontsize 10 --exon_fontsize 10





##################### kl5 exon 3-10    10jun2025  #####################     

fasta_size.awk kl5_exon3_9_9jan2025.racon.fasta #  12473

bl2seq -p blastn -i kl5_cds.fasta -j kl5_exon3_9_9jan2025.racon.fasta -F F -e 0.001 -W 16 -D 2 | sort -k 7,7n | sed /#/d | awk '{print $0 "\t" NR+2}' > kl5cds_exon3_10region.bls
cat kl5cds_exon3_10region.bls
kl-5-RC_cds     kl5_exon3_9_9jan2025_racon      100.00  175     0       0       129     303     5314    5488    4e-97    347    3
kl-5-RC_cds     kl5_exon3_9_9jan2025_racon      100.00  261     0       0       304     564     5544    5804    2e-148   517    4
kl-5-RC_cds     kl5_exon3_9_9jan2025_racon      100.00  341     0       0       563     903     5860    6200    0.0      676    5
kl-5-RC_cds     kl5_exon3_9_9jan2025_racon      100.00  267     0       0       903     1169    6257    6523    6e-152   529    6
kl-5-RC_cds     kl5_exon3_9_9jan2025_racon      99.82   570     0       1       1170    1738    6580    7149    0.0     1114    7
kl-5-RC_cds     kl5_exon3_9_9jan2025_racon      100.00  113     0       0       1738    1850    7206    7318    4e-60    224    8
kl-5-RC_cds     kl5_exon3_9_9jan2025_racon      100.00  292     0       0       1851    2142    7370    7661    7e-167   579    9
kl-5-RC_cds     kl5_exon3_9_9jan2025_racon      100.00  156     0       0       2142    2297    7719    7874    1e-85    309    10


awk '{print $9 "\t" $10 "\t" $13}' kl5cds_exon3_10region.bls > kl5_exon_3_10_locations_gene_region.txt

# I will use SVG format output because we will need to edit (exon numbers are superposed))
 
 visualize_repeats_censor4.py --censor kl5_exon3_9_9jan2025.racon.fasta.map --min_repeat_len  0  --gene_regex Ycds --gene_color black --cds_file kl5_exon_3_10_locations_gene_region.txt --len_graph 100001 --len 12473   --output kl5_exon3_10_gene_region_10jun2025.svg   --scale_fontsize 10 --exon_fontsize 10

##################### kl5  exon_14   10jun2025 / 12jun2025  #####################     

fasta_size.awk kl5_exon13_7feb2025.racon.fasta #  32405

bl2seq -p blastn -i kl5_cds.fasta -j kl5_exon13_7feb2025.racon.fasta -F F -e 0.001 -W 16 -D 2 | sort -k 7,7n | sed /#/d | awk '{print $0 "\t" NR}' 
# kl-5-RC_cds     kl5_exon13_7feb2025_racon       99.94   1674    1       0       6387    8060    18904   17231   0.0     3311    1

echo -e "18904   17231 14"  > kl5_exon14_locations_gene_region.txt


visualize_repeats_censor4.py --censor kl5_exon13_7feb2025.racon.fasta.map --min_repeat_len  0  --gene_regex Ycds --gene_color black --cds_file kl5_exon14_locations_gene_region.txt --len_graph 100001 --len 32405   --output kl5_exon14_gene_region_12jun2025.png --scale_fontsize 10 --exon_fontsize 10

############ kl5 exon  exon   11_13    10jun2025  #####################   
fasta_size.awk ctg33_hq2.fasta #  132946


bl2seq -p blastn -i kl5_cds.fasta -j ctg33_hq2.fasta  -F F -e 0.001 -W 16 -D 2 | sort -k 7,7n | sed /#/d | awk '{print $0 "\t" NR+10}' > kl5cds_exon11_13_region.bls
cat kl5cds_exon11_13_region.bls
kl-5-RC_cds     gi|33   100.00  831     0       0       2297    3127    103284  104114  0.0     1647    11
kl-5-RC_cds     gi|33   100.00  454     0       0       3127    3580    104170  104623  0.0      900    12
kl-5-RC_cds     gi|33   100.00  2810    0       0       3578    6387    104667  107476  0.0     5570    13

awk '{print $9 "\t" $10 "\t" $13}' kl5cds_exon11_13_region.bls > kl5_exon_11_13_locations_gene_region.txt
# svg output (needs editing)

visualize_repeats_censor4.py --censor kl5_exon10_12_ctg33_hq2.fasta.map --min_repeat_len  0  --gene_regex Ycds --gene_color black --cds_file kl5_exon_11_13_locations_gene_region.txt   --len 132946   --output kl5_exon11_13_gene_region_10jun2025_edited.svg --scale_fontsize 10 --exon_fontsize 10

 
 ############   kl3 exon 15-16     10jun2025  #####################   
cat  kl3_RC_exon_boundaries.txt

10056 12931 15
12929 13110 16
13110 13782 17

fasta_size.awk  kl3_region1_10jan2025.racon.fasta  #   28613


bl2seq -p blastn -i kl3_cds.fasta -j kl3_region1_10jan2025.racon.fasta   -F F -e 0.001 -W 16 -D 2 | sort -k 7,7n | sed /#/d | awk '{print $0 "\t" NR+14}' > kl3_exon15_16_region.bls
cat kl3_exon15_16_region.bls

kl-3-RC_cds     kl3_region1_10jan2025_racon     100.00  2876    0       0       10056   12931   3570    695     0.0     5701    15
kl-3-RC_cds     kl3_region1_10jan2025_racon     100.00  66      0       0       10056   10121   10085   10020   1e-31    131    p15
kl-3-RC_cds     kl3_region1_10jan2025_racon     98.48   66      0       1       10056   10121   23642   23578   7e-27    115    p15
kl-3-RC_cds     kl3_region1_10jan2025_racon     98.48   66      1       0       10056   10121   7797    7732    3e-29    123    p15
kl-3-RC_cds     kl3_region1_10jan2025_racon     98.90   182     0       1       12929   13110   642     463     3e-94    339    16

awk '{print $9 "\t" $10 "\t" $13}' kl3_exon15_16_region.bls > kl3_exon_15_16_locations_gene_region.txt

sed 's/p.*/\./'  kl3_exon_15_16_locations_gene_region.txt >  kl3_exon_15_16_locations_gene_region.simple2.txt

visualize_repeats_censor4.py --censor kl3_region1_10jan2025.racon.fasta.map --min_repeat_len  0  --gene_regex Ycds --gene_color black --cds_file kl3_exon_15_16_locations_gene_region.simple2.txt   --len 28613  --len_graph 100001   --output kl3_exon15_16_gene_region_10jun2025_edited.svg --scale_fontsize 10 --exon_fontsize 10



##################### ORY exon 1 2 10jun2025  #####################  
cat ORY_exon_boundaries.txt
1 165 1
162 512 2
513 1136 3
1137 1453 4
1452 1632 5
1633 2276 6
2277 2415 7
2414 2745 8

fasta_size.awk ORY_exon1_2_10jan2025.racon.fasta  #  53956


bl2seq -p blastn -i ORY_cds.fasta -j ORY_exon1_2_10jan2025.racon.fasta  -F F -e 0.001 -W 16 -D 2 | sort -k 7,7n  -k 9,9n | sed /#/d | awk '{print $0 "\t" NR}' > ORY_exon1_2_region.bls
cat ORY_exon1_2_region.bls

ORY-RD_cds      ORY_exon1_2_10jan2025_racon     99.39   165     1       0       1       165     3541    3377    9e-89    319    1
ORY-RD_cds      ORY_exon1_2_10jan2025_racon     97.37   38      0       1       1       38      37885   37849   1e-10   60.0    p1
ORY-RD_cds      ORY_exon1_2_10jan2025_racon     100.00  165     0       0       1       165     42733   42569   4e-91    327    1
ORY-RD_cds      ORY_exon1_2_10jan2025_racon     95.86   169     3       3       1       165     46974   46806   3e-70    258    p1
ORY-RD_cds      ORY_exon1_2_10jan2025_racon     100.00  68      0       0       98      165     37806   37739   3e-33    135    p1
ORY-RD_cds      ORY_exon1_2_10jan2025_racon     99.43   351     2       0       162     512     3313    2963    0.0      680    2
ORY-RD_cds      ORY_exon1_2_10jan2025_racon     99.72   351     1       0       162     512     42505   42155   0.0      688    2
ORY-RD_cds      ORY_exon1_2_10jan2025_racon     99.43   352     1       1       162     512     46742   46391   0.0      674    p2
ORY-RD_cds      ORY_exon1_2_10jan2025_racon     100.00  62      0       0       179     240     37671   37610   1e-29    123    p2
ORY-RD_cds      ORY_exon1_2_10jan2025_racon     100.00  34      0       0       233     266     4854    4821    5e-13   67.9    p2
ORY-RD_cds      ORY_exon1_2_10jan2025_racon     100.00  22      0       0       348     369     4752    4731    8e-06   44.1    p2
ORY-RD_cds      ORY_exon1_2_10jan2025_racon     98.15   54      1       0       375     428     37531   37478   2e-22   99.6    p2
# I used vim to glue above


awk '{print $9 "\t" $10 "\t" $13}' ORY_exon1_2_region.bls > ORY_exon1_2_locations_gene_region.txt
# I will simplify it, got  too clutered 

sed 's/p.//g'  ORY_exon1_2_locations_gene_region.txt  > ORY_exon1_2_locations_gene_region.simple2.txt

visualize_repeats_censor4.py --censor ORY_exon1_2_10jan2025.racon.fasta.map --min_repeat_len  0  --gene_regex Ycds --gene_color black --cds_file ORY_exon1_2_locations_gene_region.simple2.txt   --len 53956 --len_graph 100001  --output ORY_exon1_2_gene_region_10jun2025_edited.svg --scale_fontsize 10 --exon_fontsize 10



##################### PprY exon 3  10jun2025  #####################
 cat PprY_exon_boundaries.txt
1 208 1
209 492 2
493 671 3
672 1249 4
1246 1583 5
1582 1710 6

fasta_size.awk PprY_exon3_10jan2025.racon.fasta  #  38049

bl2seq -p blastn -i PprY_cds.fasta -j PprY_exon3_10jan2025.racon.fasta  -F F -e 0.001 -W 16 -D 2 | sort -k 7,7n  -k 9,9n | sed /#/d | awk '{print $0 "\t" NR}' > PprY_exon3_region.bls
cat PprY_exon3_region.bls
# Ppr-Y-RB_cds    PprY_exon3_10jan2025_racon      100.00  179     0       0       493     671     7236    7414    7e-100   355    1


echo "7236    7414 3" > PprY_exon3_locations_gene_region.txt


visualize_repeats_censor4.py --censor PprY_exon3_10jan2025.racon.fasta.map --min_repeat_len  0  --gene_regex Ycds --gene_color black --cds_file PprY_exon3_locations_gene_region.txt   --len 38049 --len_graph 100001  --output PprY_exon3_gene_region_10jun2025.png  --scale_fontsize 10 --exon_fontsize 10

##################### PprY exon 4  10jun2025  #####################
 
fasta_size.awk  ctg166_flyehq2_130_245.fasta  #  115001

bl2seq -p blastn -i PprY_cds.fasta -j  ctg166_flyehq2_130_245.fasta  -F F -e 0.001 -W 16 -D 2 | sort -k 7,7n  -k 9,9n | sed /#/d | awk '{print $0 "\t" NR}' > PprY_exon4_region.bls
cat PprY_exon4_region.bls
# Ppr-Y-RB_cds    gi|166:130000-245000    99.83   579     0       1       672     1249    31333   30755   0.0     1132    p4
# Ppr-Y-RB_cds    gi|166:130000-245000    100.00  578     0       0       672     1249    86977   87554   0.0     1146    4

echo -e "31333   30755  \n86977   87554 4" > PprY_exon4_locations_gene_region.simple2.txt
visualize_repeats_censor4.py --censor PprY_exon4_ctg166_flyehq2_130_245.fasta.map --min_repeat_len  0  --gene_regex Ycds --gene_color black --cds_file PprY_exon4_locations_gene_region.simple2.txt   --len 115001  --output PprY_exon4_gene_region_10jun2025.png  --scale_fontsize 10 --exon_fontsize 10

##################### Pp1Y1  10jun2025  #####################
fasta_size.awk  ctg239_14jun24_350_450.fasta  #  100001

bl2seq -p blastn -i Pp1Y1_cds.fasta -j  ctg239_14jun24_350_450.fasta  -F F -e 0.001 -W 16 -D 2 | sort -k 7,7n  -k 9,9n | sed /#/d | awk '{print $0 "\t" NR}' > Pp1Y1_region.bls
cat Pp1Y1_region.bls
# Pp1-Y1-RC_cds   gi|239:350000-450000    100.00  954     0       0       1       954     49417   48464   0.0     1891    1

echo   "49417   48464 1" > Pp1Y1_exon_locations_gene_region.txt

visualize_repeats_censor4.py --censor Pp1Y1_ctg239_14jun24_350_450.fasta.map --min_repeat_len  0  --gene_regex Ycds --gene_color black --cds_file Pp1Y1_exon_locations_gene_region.txt   --len 100001  --output Pp1Y1_gene_region_10jun2025.png  --scale_fontsize 10 --exon_fontsize 10

##################### Pp1Y2  10jun2025  #####################
fasta_size.awk  ctg246_14jun24_128_228.fasta  #  100001

bl2seq -p blastn -i Pp1Y2_cds.fasta -j ctg246_14jun24_128_228.fasta  -F F -e 0.001 -W 16 -D 2 | sort -k 7,7n  -k 9,9n | sed /#/d | awk '{print $0 "\t" NR}' > Pp1Y2_region.bls
cat Pp1Y2_region.bls
# Pp1-Y2-RC_cds   gi|246:128000-228000    100.00  942     0       0       1       942     50870   51811   0.0     1867    1
echo   "50870   51811 1" > Pp1Y2_exon_locations_gene_region.txt

visualize_repeats_censor4.py --censor Pp1Y2_ctg246_14jun24_128_228.fasta.map --min_repeat_len  0  --gene_regex Ycds --gene_color black --cds_file Pp1Y2_exon_locations_gene_region.txt   --len 100001  --output Pp1Y2_gene_region_10jun2025.png  --scale_fontsize 10 --exon_fontsize 10


##################### white  10jun2025  #####################
fasta_size.awk white_ctg72_flyehq2_19216_19316.fasta  #  100001

bl2seq -p blastn -i white_cds.fasta -j white_ctg72_flyehq2_19216_19316.fasta  -F F -e 0.001 -W 16 -D 2 | sort -k 7,7n  -k 9,9n | sed /#/d | awk '{print $0 "\t" NR}' > white_region.bls
cat white_region.bls
w-RA_cds        gi|72:19216000-19316000 100.00  73      0       0       1       73      50007   50079   4e-36    145    1
w-RA_cds        gi|72:19216000-19316000 100.00  274     0       0       73      346     53187   53460   5e-156   543    2
w-RA_cds        gi|72:19216000-19316000 100.00  656     0       0       346     1001    53534   54189   0.0     1300    3
w-RA_cds        gi|72:19216000-19316000 100.00  316     0       0       1002    1317    54251   54566   0.0      626    4
w-RA_cds        gi|72:19216000-19316000 100.00  135     0       0       1317    1451    54769   54903   4e-73    268    5
w-RA_cds        gi|72:19216000-19316000 100.00  615     0       0       1450    2064    54972   55586   0.0     1219    6

#getting exon boundaries   
awk '($13 !~/p/){print $7 "\t" $8 "\t" $13}' white_region.bls > white_exon_boundaries.txt
#getting exon locations in the gene region, for   visualize_repeats_censor3.py     
awk '{print $9 "\t" $10 "\t" $13}' white_region.bls > white_exon_locations_gene_region.txt	

# saving in svg to edit

visualize_repeats_censor4.py --censor white_ctg72_flyehq2_19216_19316.fasta.map --min_repeat_len  0  --gene_regex Ycds --gene_color black --cds_file white_exon_locations_gene_region.txt   --len 100001  --output white2_gene_region_10jun2025.svg  --scale_fontsize 10 --exon_fontsize 10


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
fplot_reads.awk   Fig6_data.txt  # produced file  Fig6_data.svg . This will print the arrows, whcih represent the reads
# the central part showing the satellites, exon, and TEs came from Fig. 3
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

simple_satellite_detector.py --input_file CCY_exon2_13reads.fasta --satellite AAAC --min_copies 4 --output CCY_exon2_13reads_AAAC_mincopy1.txt
# repeat the above script for each combination of fasta sequence and main satellite

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
