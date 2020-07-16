
######TO EDIT#######

##put prefix to be used for output, usually just genome name
prefix="AAR.Sdn_fl_r3_m2"

##put path to genome for analysis
genome="AAR.Sdn_fl_r3_m2.fa"

##size to be selected at ends of contigs, to be searched for telomeric and subtelomeric sequences
end_size="1500"

#######DO NOT EDIT#######
###CURRENTLY THIS SCRIPT REQUIRES THAT THE TELOMERE REPEAT TG{1-3} or AC{1-3} IS REPEATED AT LEAST 6 TIMES IN A ROW
#to change this use the below command and change 'NUMBER' to desired value
#cat SCRIPT.sh | sed 's/{6,}/{NUMBER,}/g' > SCRIPT_NEW.sh

##First take any genome and make sure that each contig is placed on a single line, i.e. some may come wrapped
cat $genome | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' | grep "\S" > ${prefix}.UNWRAPED.fa

##get all contig names
cat ${prefix}.UNWRAPED.fa | grep '>' | while read contig
	do 
		##get contig size to make sure is large enough to trim 
		contig_name=$(echo $contig | awk '{print $1}' | sed 's/>//g')
		contig_size=$(grep -A 1 "^$contig$" ${prefix}.UNWRAPED.fa | tail -n 1 | awk '{print length}')
		
		echo $contig_name
		echo $contig_size >> ${prefix}_contig_length.txt
			##use contig names to get the actual contig and then take only the first #kb as selected above
			#check if the telomeric sequence is present
			grep -A 1 "^$contig$" ${prefix}.UNWRAPED.fa | tail -n 1 | cut -c1-$end_size |\
				grep -oE '(T{1,2}G{1,5}){6,}|(A{1,2}C{1,5}){6,}|(G{1,5}T{1,2}){6,}|(C{1,5}A{1,2}){6,}' | awk '{print length"\t"$0}' | sort -nr | cut -f2 | head -n1 | awk -v contig="$contig_name" 'END{if(NR < 1) {print contig"\tNA\tNA"} else {print contig"\t"$0"\t"length}}' >> ${prefix}_Telomere_L_presence.txt
			#create variable containing the fragment labelled as the telomeric fragment
			telo=$(grep -A 1 "^$contig$" ${prefix}.UNWRAPED.fa | tail -n 1 | cut -c1-$end_size |\
				grep -oE '(T{1,2}G{1,5}){6,}|(A{1,2}C{1,5}){6,}|(G{1,5}T{1,2}){6,}|(C{1,5}A{1,2}){6,}' | awk '{print length"\t"$0}' | sort -nr | cut -f2 | head -n1)
			#create variable to see if a telomere end was found
			found=$(grep -A 1 "^$contig$" ${prefix}.UNWRAPED.fa | tail -n 1 | cut -c1-$end_size |\
				grep -oE '(T{1,2}G{1,5}){6,}|(A{1,2}C{1,5}){6,}|(G{1,5}T{1,2}){6,}|(C{1,5}A{1,2}){6,}' | awk '{print length"\t"$0}' | sort -nr | cut -f2 | head -n1 | awk -v contig="$contig_name" 'END{if(NR < 1) {print "NA"} else {print "FOUND"}}')
			##get the distance of the first telomeric repeat to the start
			if [ "$found" = "FOUND" ]
			then
			grep -A 1 "^$contig$" ${prefix}.UNWRAPED.fa | tail -n 1 | cut -c1-$end_size |\
				sed "s/$telo/;/g" |\
				cut -f1 -d ';' | awk '{print length}' >> ${prefix}_Telomere_L_distance.txt
			else
				echo "NA" >> ${prefix}_Telomere_L_distance.txt
			fi
			##same but using the last #kb, so first get the whole length
			contig_end=$(($contig_size-$end_size))
			grep -A 1 "^$contig$" ${prefix}.UNWRAPED.fa | tail -n 1 | cut -c$contig_end-$contig_size |\
				grep -oE '(T{1,2}G{1,5}){6,}|(A{1,2}C{1,5}){6,}|(G{1,5}T{1,2}){6,}|(C{1,5}A{1,2}){6,}' | awk '{print length"\t"$0}' | sort -nr | cut -f2 | head -n1 | awk -v contig="$contig_name" 'END{if(NR < 1) {print contig"\tNA\tNA"} else {print contig"\t"$0"\t"length}}' >> ${prefix}_Telomere_R_presence.txt
			telo2=$(grep -A 1 "^$contig$" ${prefix}.UNWRAPED.fa | tail -n 1 | cut -c$contig_end-$contig_size |\
				grep -oE '(T{1,2}G{1,5}){6,}|(A{1,2}C{1,5}){6,}|(G{1,5}T{1,2}){6,}|(C{1,5}A{1,2}){6,}' | awk '{print length"\t"$0}' | sort -nr | cut -f2 | head -n1)
			found2=$(grep -A 1 "^$contig$" ${prefix}.UNWRAPED.fa | tail -n 1 | cut -c$contig_end-$contig_size |\
				grep -oE '(T{1,2}G{1,5}){6,}|(A{1,2}C{1,5}){6,}|(G{1,5}T{1,2}){6,}|(C{1,5}A{1,2}){6,}' | awk '{print length"\t"$0}' | sort -nr | cut -f2 | head -n1 | awk -v contig="$contig_name" 'END{if(NR < 1) {print "NA"} else {print "FOUND"}}')
			if [ "$found2" = "FOUND" ]
			then
			##get the distance of the first telomeric repeat to the end
			grep -A 1 "^$contig$" ${prefix}.UNWRAPED.fa | tail -n 1 | cut -c$contig_end-$contig_size |\
				sed "s/$telo2/;/g" |\
				awk -F";" '{print $NF}' | awk '{print length}' >> ${prefix}_Telomere_R_distance.txt
			else
				echo "NA" >> ${prefix}_Telomere_R_distance.txt
			fi
	done

paste ${prefix}_Telomere_L_presence.txt ${prefix}_Telomere_L_distance.txt ${prefix}_Telomere_R_presence.txt ${prefix}_Telomere_R_distance.txt ${prefix}_contig_length.txt |\
	awk 'BEGIN{print "Left_boarder\ttelomere\tsize\tdistance_to_L\tRight_boarder\ttelomere\tsize\tdistance_to_R\tcontig_size"} {print $0}'> ${prefix}_summary.csv

rm ${prefix}_contig_length.txt
rm ${prefix}_Telomere_L_presence.txt
rm ${prefix}_Telomere_L_distance.txt
rm ${prefix}_Telomere_R_presence.txt
rm ${prefix}_Telomere_R_distance.txt
#rm ${prefix}.UNWRAPED.fa





####GREP EXAPLANTION
# () groups the string
# {1} says find this letter one time
# {1,3} says one to three times
# {6,} says find all in () 6 times in a row or more
# | says or
##i.e. look for TG1-3 or AC1-3 repeated six times in a row minimum
# this finds every line and prints the whole line if regular expression is found
# grep -E '(T{1,2}G{1,3}){6,}|(A{1,2}C{1,3}){6,}' genome.fa
# using the -o option prints each match on an individual line
# grep -oE '(T{1,2}G{1,3}){6,}|(A{1,2}C{1,3}){6,}' genome.fa