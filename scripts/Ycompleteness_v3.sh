#!/bin/bash
# Ycompleteness_v3.sh    27mar2024  t test completeness of  genes (usually Y-linked) in   assemblies   v. 18may2025
# new in v3: accests any CDS file, so can be used with other scpeciwes; customizable and more stringent blast parameters ; ref cCDS can be anywghere (uses locate); no longer uses blasthits_seq_occupancy2.awk [it assumes that all hits doto a cds came together in the m8, which is risky: megablast -D 3 does not do this. ; produces tabular output in the end
 
# differnce from v1: full paths in the files and folders, so the script can run anywhere. Also, all files written in Y_completeness folder; b) formatdb uses -o F, so contigs can have the same bame 

set -e

function PrintUsage() {
   echo "Usage: Ycompleteness_v3.sh  -s yes -r dmel_Y_CDS_2.fasta -P 98 -a 100   assembly_1.fasta  assembly_2.fasta  ... assembly_N.fasta "
   exit 1
}

if [ "$#" == "0" ]; then
	PrintUsage
	exit 1
fi

while getopts "hr:s:d:t:P:a:D:p:T:o:" OPTION
do
   case $OPTION in
        h) PrintUsage
         ;;
        r) ref_fasta=$OPTARG
         ;;
        p) blast_parameters=$OPTARG
         ;;      
        s) save_files=$OPTARG    # do not erase temporary  folder 
         ;;          
        t) threads=$OPTARG
         ;;    
        d) fasta_directory="$OPTARG/"
         ;;        
        P) min_perc_id=$OPTARG
         ;;   
        a) min_aln_len=$OPTARG
         ;;	
        D) run_dir=$OPTARG
         ;;   
		T) tab_space=$OPTARG
         ;;   
		o) final_tab=$OPTARG   # prints the side by side tab (usefiul for comparing more than one assembly)
         ;;   		
   esac
done


# setting default values if not specified by the user:
[[ $threads = "" ]] && threads=50
[[ $fasta_directory = "" ]] && fasta_directory="$PWD/"
[[ $blast_parameters = "" ]] && blast_parameters=" -W 20 "
[[ $min_perc_id = "" ]] && min_perc_id=97
[[ $min_aln_len = "" ]] && min_aln_len=50
[[ $save_files = "" ]]  && save_files="no"
[[ $tab_space = "" ]] && tab_space=20
[[ $final_tab = "" ]] && final_tab="yes"


init_dir="$PWD/"


# empty or create $run_dir
[[ $run_dir  = "" ]] && run_dir="/tmp/temp_Ycompleteness"
# [ -d $run_dir  ]  &&  rm $run_dir/*    # gives error when $run_dir exists but is empty
[ -d $run_dir  ]  &&  rm -f $run_dir/*.m8  $run_dir/CDS.size $run_dir/results.txt
[ ! -d $run_dir ] &&  mkdir  $run_dir


# [[ "${ref_fasta:0:1}" == '/' ]] && full_path_ref_fasta=${ref_fasta}
# [[ "$1" =~ "/" ]]    && full_path_ref_fasta=${ref_fasta}
# ! [[ "$1" =~ "/" ]]  && full_path_ref_fasta=`locate $ref_fasta | head -1 `

if [[ "$ref_fasta" =~ "/" ]]
then
  full_path_ref_fasta=${ref_fasta}
else
  full_path_ref_fasta=`locate $ref_fasta | head -1` 
fi

[[ "$full_path_ref_fasta" = "" ]] && echo "file " $ref_fasta " not found. Please check name and give full path" && exit 1



echo "using as  ref_fasta: "  $full_path_ref_fasta


# creating CDS size file
fasta_size2.awk   $full_path_ref_fasta | sed '/ATCGN/d' > $run_dir/CDS.size
 
while [[ "$1" =~ ^- ]]
   do
      shift ; shift   # skips for example  "-e"  "fasta"  , which were already taken care by getopts
   done



# while [ "$1" != "" ]
while (( "$#" ))  # $# in bash is the number of paramaters 
do	
	assembly=${1/".fasta"/""}  # removes ".fasta"   but keeps the full path, if present (/home6/xxx)
	echo $assembly "     @" | tee -a $run_dir/results.txt
	gawk 'BEGIN{printf("%-10s %8s %8s %8s %8s\n" ,"gene","duplic","size","occup","% occup")}'
	assembly_short_name=${assembly//*\//""}  # removes full path, for the formatdb  etc  [will remove all characters until (inclusive) the last "/"
	# echo "assembly_short_name: " $assembly_short_name	
	formatdb -i $assembly.fasta   -p F -o F -n $run_dir/$assembly_short_name -l $run_dir/formatdb.log  2>$run_dir/formatdb.error
	# megablast -d $rundir/$assembly_short_name -i $rundir/dmel_Y_CDS_2.fasta -p 98 -m 8 -a 20 -F "m D" > $rundir/$assembly_short_name.bls
	megablast -d $run_dir/$assembly_short_name -i $full_path_ref_fasta -F "m D"  $blast_parameters  -a $threads  -m 8 -b 100000 -v 100000 2>$run_dir/$assembly_short_name.error  >  $run_dir/$assembly_short_name.unfiltered.m8
	gawk -v min_perc_id=$min_perc_id   -v min_aln_len=$min_aln_len '(($3>=min_perc_id)&&($4>=min_aln_len)){print $0}'   $run_dir/$assembly_short_name.unfiltered.m8  > $run_dir/$assembly_short_name.m8
	# cd $run_dir
	gawk -v run_dir=$run_dir 'BEGIN{while ((getline line < (run_dir "/CDS.size")) > 0){split(line,a);size[a[1]]=a[2]; sizetot +=a[2]}}; !/^#/{g_array[$1]++; for (j=$7;j<= $8;j++){gb_array[$1][j]++}};END{bptot=0;duptot=0;for(g in size){bp=0;dup=0; for(b in gb_array[g]){bp ++; bptot ++; if(gb_array[g][b] >1){dup ++; duptot++}};printf("%-10s %8s %8s %8s %8.1f\n" ,g,dup,size[g],bp,100*(bp/size[g]))}; printf("%-10s %8s %8s %8s %8.1f\n" ,"TOTAL",duptot,sizetot,bptot,100*(bptot/sizetot))}'  $run_dir/$assembly_short_name.m8 | tee -a $run_dir/results.txt	
	# cd $init_dir	
	echo ""
	rm -f  $run_dir/*.nhr   $run_dir/*.nsd  $run_dir/*.nsi $run_dir/*.nsq   $run_dir/*.nin   # 2> /dev/null
	# [ "$save_files" = "no" ] && rm -f $run_dir/*.m8  $run_dir/CDS.size &&  rmdir  $run_dir  # 2> /dev/null
	shift   # Shift to next assembly
done

[[ $final_tab = "yes" ]] && gawk -v tab_space=$tab_space 'BEGIN{ts="%" tab_space "s"; PROCINFO["sorted_in"]="@val_num_asc"};/@/{asm=$1;asm_array[asm]=NR}; (NF==5){b[$1,asm]=$5;gene_array[$1]=NR}; END{printf("          ");for(a in asm_array){printf(ts,a)};for(g in gene_array){printf("\n%-10s",g);for(a in asm_array){printf(ts,b[g,a]) }}print ""}' $run_dir/results.txt

[ "$save_files" = "no" ] && rm -f $run_dir/*.m8  $run_dir/CDS.size $run_dir/*.txt $run_dir/*.log $run_dir/*.error   &&  rmdir  $run_dir  # 2> /dev/null



# kl-2    contig_514      100.00  667     0       0       1       667     60725   60059   0.0     1322
# kl-3    contig_514      100.00  1026    0       0       1       1026    283085  284110  0.0     2034
# kl-3    contig_514      99.90   1026    1       0       1       1026    294319  295344  0.0     2026
# kl-3    contig_177      100.00  673     0       0       13110   13782   13830   14502   0.0     1334
# WDY     contig_73       100.00  1080    0       0       612     1691    360569  359490  0.0     2141
# WDY     contig_73       99.85   669     1       0       2542    3210    57661   56993   0.0     1318
# WDY     contig_73       99.85   669     1       0       2542    3210    61373   60705   0.0     1318
