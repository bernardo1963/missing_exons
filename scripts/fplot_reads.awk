#! /usr/bin/awk -f
# fplot_reads.awk    Bernardo   16feb20255 for missing exon projewct   
# usage: fplot_reads.awk 18reads_data_strand_edited.txt 

# reads association file produced by web dgenies (https://dgenies.toulouse.inra.fr)  and produce gnuplot script  to plot the reads as arrows

($1 == "Query"){next}  # header, if present

(flag_first_read == 0){	
	ctg_size = $7
	flag_first_read ++
	}

# F reads:
($3 == "+"){
	read_count ++
	start_array[read_count] = $8 
	end_array[read_count] = $9 
	}
	
# R reads:
($3 == "-"){
	read_count ++
	start_array[read_count] = $9 
	end_array[read_count] = $8 
	}


END{
	gnuplot_script = gensub(/\..+/,"",1,FILENAME) ".gp"
	system("rm -fr " gnuplot_script)
	print "set terminal svg dynamic size 1036,300" > gnuplot_script
	print "set size ratio 0.29" > gnuplot_script
	# print "set border 4"  > gnuplot_script
	print "unset border" > gnuplot_script
	print "set output \x27" gensub(/\..+/,"",1,FILENAME) ".svg\x27" > gnuplot_script
	print "unset grid" > gnuplot_script
	print "unset ytics" > gnuplot_script
	print "unset xtics" > gnuplot_script	
	# print "set xtics out" > gnuplot_script
	print "set style arrow 1 head  filled size screen 0.02,4,90   lw 1 lc rgb \x22" "black\x22" > gnuplot_script
	# print "set style arrow 2 backhead lw 1 lc rgb \x22black\x22" > gnuplot_script
	print "set xrange [-1000:"  ctg_size +1000  "]" > gnuplot_script
	print "set yrange [0:" (read_count +1) "]" > gnuplot_script
	print "set arrow 1 from 1,0 to " ctg_size ",0 nohead lw 4 lc rgb \x22" "black\x22" > gnuplot_script
	for (i=1; i<= read_count; i++){
		print "set arrow " (i+1) " arrowstyle 1 from " start_array[i]"," i " to " end_array[i]","i > gnuplot_script
		}
	print "plot NaN t \x27\x27" > gnuplot_script
	print "unset output" > gnuplot_script
	system("gnuplot " gnuplot_script)
	}


# desired gnuplot commands

# set terminal svg
# set output 'arrows.svg'
# unset grid
#set title "arrows"
#F reads:
# set style arrow 1 head     lw 1 lc rgb "black"
#R reads 
# set style arrow 2 backhead lw 1 lc rgb "black"

# set xrange [0:38000]
# set yrange [0:5]
# set arrow 1 from 0,1    to 38000,1 nohead lw 4 lc rgb "black"
# set arrow 2 arrowstyle 1 from 1000,2 to 15000,2 
# set arrow 3 arrowstyle 2 from 2000,3 to 20000,3 
# set arrow 4 arrowstyle 1 from 5000,4 to 10000,4
# plot NaN t ''

# unset output




# Query	Target	Strand	Q-len	Q-start	Q-stop	T-len	T-start	T-stop
# SRR22822929.277518	PprY_exon3_10jan2025_racon	+	12122	33	12082	38049	7290	19341
# SRR26246282.503891	PprY_exon3_10jan2025_racon	+	34185	37	32643	38049	7230	38048
# SRR26246282.1589819	PprY_exon3_10jan2025_racon	+	35354	3427	35342	38049	4864	38048   quimeric?
# SRR26246282.13567	PprY_exon3_10jan2025_racon	-	24282	15	24269	38049	2433	27232
# SRR26246282.998787	PprY_exon3_10jan2025_racon	-	25044	24	24912	38049	2	25068
# SRR26246282.1364601	PprY_exon3_10jan2025_racon	-	25054	35	25021	38049	5750	30777
# SRR26246282.1416290	PprY_exon3_10jan2025_racon	-	43828	198	43583	38049	6810	38048
# SRR26246282.1452085	PprY_exon3_10jan2025_racon	-	10796	39	10786	38049	5786	16539
# SRR26246282.1456379	PprY_exon3_10jan2025_racon	-	30248	22	30105	38049	3474	34158
# SRR26246282.1589807	PprY_exon3_10jan2025_racon	-	15931	5	15814	38049	3421	21217
# SRR26246282.1615639	PprY_exon3_10jan2025_racon	-	32120	9	30562	38049	2	24152
# SRR26246282.1673132	PprY_exon3_10jan2025_racon	-	29964	61	29961	38049	6057	36321
# SRR22822929.19709	PprY_exon3_10jan2025_racon	-	16230	37	16176	38049	393	22224
# SRR22822930.148478	PprY_exon3_10jan2025_racon	-	22085	34	21602	38049	6810	28452
# SRR26246282.512264	PprY_exon3_10jan2025_racon	-	13511	34	13271	38049	7315	20613
# SRR26246282.1077839	PprY_exon3_10jan2025_racon	-	3730	44	3715	38049	7315	11028
# SRR22822930.1058108	PprY_exon3_10jan2025_racon	-	4214	31	4200	38049	7363	11523
# SRR26246282.308816	PprY_exon3_10jan2025_racon	+-	1250	73	1173	38049	7256	7841    artifact
