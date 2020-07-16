import sys
from Bio.Seq import Seq
from Bio import SeqIO
import os.path


## !!! MAKE SURE TO DELETE THE OUTPUT FILE "telom_length.csv" FROM YOUR DIRECTORY BEFORE RUNNING THE SCRIPT !!!

## command to launch the script: $python 00_analyze_telom_length.py file.fasta

## This script looks if telomere sequences were present at the junction between merged contigs during the scaffolding step
## and caculates the size of these sequences and of a potential offset sequences between the scaffolding mark (50 Ns)
## and the telomere sequence

## define the genme fasta file
if len(sys.argv) != 2:
    sys.exit("ERROR: you must specify the name of genome fasta file or its path as the first argument")

## function#1: return the strain name, each chromosome name and the size of putative telomere internally merged on each side
## of the scaffolding sequence (50 Ns)
def chr_start_end(genome_fasta):

    ## get strain name from fasta files either from a path or in the src directory
    ## from a path
    if "/" in genome_fasta:
    	file_name = genome_fasta.split("/")
    	assembly = file_name[-1].split(".")
    	strain = assembly[0]
    ## from the src directory
    else:
    	file_name = genome_fasta.split(".")
    	strain = file_name[0]
    
    print("\n", "=================================",
        "\n", strain,
        "\n", "=================================",
        "\n")

    ## parse the multifasta file and reverse complement each sequence
    for seq_record in SeqIO.parse(genome_fasta, "fasta"):
        revcomp = seq_record.reverse_complement()

        ## stdout
        chrom = seq_record.id
        print("---------------------")
        print(chrom)

        ## Look for 50 Ns:
        n_seq = "N"*50
        nb_n = 0
        starts_N = [] # list of start coordinates of the firsts of the 50 N stretches (strand +)
        ends_N = [] # list of end coordinates of the lasts of the 50 N stretches  (strand +)
        revends_N = [] # list of end coordinates of the lasts of the 50 N stretches  (strand -)

        ## look for telomere sequence on the RIGHT side of the scaffolding mark (50 Ns)
        if n_seq in seq_record:
            ## define coordinates of the 50 Ns
            for i in range(0, len(seq_record)-49):
                if seq_record.seq[i] == "N" and seq_record.seq[i+49] == "N" and seq_record.seq[i-1] != "N":
                    nb_n += 1
                    starts_N.append(i) ## append the coordinate of the first N for each stretch of 50 Ns in the fasta seq
                    ends_N.append(i+50) ## append the coordinate of the last N for each stretch of 50 Ns in the fasta seq
            ## stdout
            #print("the scaffold contains {} stretch(es) of 50 Ns starting at coordinates {} and ending at coordinates {}:"\
            #    .format(nb_n, starts_N, ends_N))

            ## calculate internal telomere and offset size
            for coord in ends_N:
                count = 0
                #revcount = 0
                offset = 0
                #revoffset = 0

                ## set up a limit of 1000 nt at the end of N stretch to browse until we find a telomere sequence
                if len(seq_record.seq) - coord < 1500:
                    limit = len(seq_record.seq) - coord - 20
                else:
                    limit = 1500
                
                ## estimate the size of the offset sequence before the telomere sequence
                for j in range(coord, coord+limit):
                    mot = str(seq_record.seq[j:j+20] + "\n")
                    if mot.count("C") + mot.count("A") < 19:
                        offset += 1
                    else:
                        break
                ## estimate the size of the telomere sequence
                for k in range(coord+offset, len(seq_record.seq) -20):
                    mot = str(seq_record.seq[k:k+20] + "\n")
                    if mot.count("C") >= 8 and mot.count("A") >= 3 and mot.count("C") + mot.count("A") >= 17 and "AAAA" not in mot:
                        count += 1
                    else:
                        break

                ## definition of telomere and offset lengths    
                if count == 0:
                    rmerge_tel = 0
                    offset = 0
                elif count != 0 and offset != 0:
                    rmerge_tel = 20+count-3
                    offset = offset +1
                elif count != 0 and offset == 0:
                    rmerge_tel = 20+count-3


        ## look for telomere sequence on the LEFT side of the scaffolding mark (50 Ns) - same code as above but on revcomp sequence
        if n_seq in revcomp:
            ## define coordinates of the 50 Ns
            for i in range(0, len(seq_record)):
                if revcomp.seq[i] == "N" and revcomp.seq[i-1] != "N":
                    revends_N.append(i+50) ## append the coordinate of the last N for each stretch of 50 Ns in the fasta seq

            ## calculate internal telomere and offset size
            for coord in revends_N:
                revcount = 0
                revoffset = 0

                ## set up a limit of 1000 nt at the end of N stretch to browse until we find a telomere sequence
                if len(seq_record.seq) - coord < 1500:
                    revlimit = len(seq_record.seq) - coord - 20
                else:
                    revlimit = 1500
                
                ## estimate the size of the offset sequence before the telomere sequence
                for j in range(coord, coord+revlimit):
                    mot = str(revcomp.seq[j:j+20] + "\n")
                    if mot.count("C") + mot.count("A") < 19:
                        revoffset += 1
                    else:
                        break
                ## estimate the size of the telomere sequence
                for k in range(coord+revoffset, len(seq_record.seq) -20):
                    mot = str(revcomp.seq[k:k+20] + "\n")
                    if mot.count("C") >= 8 and mot.count("A") >= 3 and mot.count("C") + mot.count("A") >= 17 and "AAAA" not in mot:
                        revcount += 1
                    else:
                        break

                ## definition of telomere and offset lengths    
                if revcount == 0:
                    lmerge_tel = 0
                    revoffset = 0
                elif revcount != 0 and revoffset != 0:
                    lmerge_tel = 20+revcount-3
                    revoffset = revoffset +1
                elif revcount != 0 and revoffset == 0:
                    lmerge_tel = 20+revcount-3

                # stdout
                if lmerge_tel != 0 or rmerge_tel != 0:
                    print("\t\tright merged telom length = ", rmerge_tel)
                    print("\t\tright offset = ", offset)
                    print("\t\tleft merged telom length = ", lmerge_tel)
                    print("\t\tleft offset = ", revoffset)    

        
                ## save the results in a csv file
                file_exists = os.path.isfile("telom_merged_scaffold.csv")
                with open("telom_merged_scaffold.csv", "a") as filout:
                    if file_exists:
                        filout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(strain, chrom, nb_n, starts_N, ends_N, lmerge_tel, offset, rmerge_tel, revoffset))
                    else:
                        filout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\n"\
                            .format("Strain", "Chromosome", "nb_50Ns", "start_50Ns", "end_50Ns", "size_tel_merged_left", "offset_left", "size_tel_merged_right", "offset_right", \
                            strain, chrom, nb_n, starts_N, ends_N, lmerge_tel, offset, rmerge_tel, revoffset)) # file doesn't exist yet, write a header


        


if __name__ == "__main__":
    chr_start_end(sys.argv[1])
