#! /bin/gawk -f
#counts the lenght of multiple fasta files and produces a table in the output     v. 25jul2018
#usage: fasta_size.awk    mel_Fplain_h35h60.fasta  
# fasta_size.awk    defline_type=1 mel_Fplain_h35h60.fasta
 
BEGIN {processed_fasta=0
	total_size_fasta = 0
      }
/^>/ {
	if (NR != 1){print name_fasta, size_fasta}
    ++processed_fasta
    size_fasta = 0
    if ( (defline_type=="simple")||(defline_type=="") ) {name_fasta=substr($1,2)}  			# >78031    (gets 78031)
    if   (defline_type=="gi")                       {split($1,a,"|") ;name_fasta=a[2]} 	# gi|55845883|gb|AAFS01000001.1|    (gets 5845883)
    if   (defline_type=="acc")                       {split($1,a,"|") ;name_fasta=a[4]}		# >gi|7289045|gb|AE003219.1|AE003219    (gets AE003219)     # >gi|18118649|gb|Tr18118648|MP18576395 46306 1 10000 R 79347
    if   (defline_type=="gnl")                       {split($1,a,"|") ;name_fasta=a[3]}		# >gnl|ti|156548644 41850699    (gets 156548644)
    }


!/^>/{
	size_fasta = size_fasta + length($0)
	total_size_fasta = total_size_fasta + length($0)
	}

END {   print name_fasta, size_fasta
	#print "fasta sequences found: ",processed_fasta
	#print "total seq length: ", total_size_fasta
	}

