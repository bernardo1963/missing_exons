#! /bin/gawk -f
# process_nonB_gfa.awk   Bernardp for missing exons p[roject  13jun2025
# usage: process_nonB_gfa.awk  output=raw [default: table]  Y_regions.out_*.tsv
    
       
BEGIN{} 

(NR==1){
	if (output==""){output="table"}
	if (output=="raw"){
		printf("%-15s\t%-20s\t%7s\n","region$","motif$","bp")
		}
	}
 
/Sequence_name/{next} 

 
(flag_start==0){curr_region=$1; curr_motif=$3}

(($1!=curr_region) || ($3 != curr_motif)){
	# if (output=="raw"){printf("%-15s\t%-20s\t%7s\n",gensub(/_region/,"",1,curr_region),curr_motif,0+length(base_array)) }
	result_array[curr_region,curr_motif] = 0+length(base_array)
	region_array[$1]++ 
	motif_array[$3]++
	curr_region=$1; curr_motif=$3
	delete base_array
	}
	
(($1==curr_region) && ($3 == curr_motif)){
	for(idx=$4; idx<=$5; idx++){base_array[idx]++}
	flag_start ++
	}

END{
	result_array[curr_region,curr_motif] = 0+length(base_array)
	if (output=="raw"){
		for (m in motif_array){
			for(r in region_array){			
				# printf("%20s",0+result_array[r,m])				
				printf("%-15s\t%-20s\t%7s\n",gensub(/_region/,"",1,r),m,0+result_array[r,m]) }				
				} 
			}
		# printf("%-15s\t%-20s\t%7s\n",gensub(/_region/,"",1,curr_region),curr_motif,0+length(base_array)) }
	# {print gensub(/_region/,"",1,curr_region) ,curr_motif, 0+length(base_array)}
	else{
		printf("%15s",".")
		for(m in motif_array){printf("%20s",m)}
		print ""

		for(r in region_array){
			printf("%15s",gensub(/_region/,"",1,r))
			for (m in motif_array){
				printf("%20s",0+result_array[r,m])
				} 
			print ""
			} 
		}
	}