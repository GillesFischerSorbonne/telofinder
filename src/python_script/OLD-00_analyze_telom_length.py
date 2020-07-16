import sys
from Bio.Seq import Seq
from Bio import SeqIO
import os.path

## to launch the script: $python 00_analyze_telom_length.py file.fasta

## define the genme fasta file
if len(sys.argv) != 2:
    print("ERROR: you must specify the genome fasta file (genome.fasta) as argument")
    sys.exit()

## function#1: return the strain name, each chromosome name and their telomere length on each side
def chr_start_end(genome_fasta):
    file_name = genome_fasta.split(".")
    strain = file_name[-2]
    print("\n", "=================================",
        "\n", strain,
        "\n", "=================================",
        "\n")

    for seq_record in SeqIO.parse(genome_fasta, "fasta"):
        count = 0
        revcount = 0
        revcomp = seq_record.reverse_complement()
        for i in range(0, len(seq_record.seq)-9):
            mot = str(seq_record.seq[i:i+10] + "\n")
            if mot.count("C") + mot.count("A") >= 9:
                count += 1
            else:
                break
        for i in range(0, len(seq_record.seq)-9):
            revmot = str(revcomp.seq[i:i+10] + "\n")
            if revmot.count("C") + revmot.count("A") >= 9:
                revcount += 1
            else:
                break

        if count == 0:
            left_tel = 0
        else:
            left_tel = 10+count-1
            
        if  revcount == 0:
            right_tel = 0
        else:
            right_tel = 10+revcount-1

        chrom = seq_record.id
        print(chrom)
        print("left telom length = ", left_tel)
        print("right telom length = ", right_tel)
        print("\n",
            "-------------------------------",
            "\n")
        
        ## save the results in a csv file
        file_exists = os.path.isfile("telom_length.csv")
        with open("telom_length.csv", "a") as filout:
            if file_exists:
                filout.write("{0}\t{1}\tL\t{2}\n{0}\t{1}\tR\t{3}\n".format(strain, chrom, left_tel, right_tel))
            else:
                filout.write("{0}\t{1}\t{2}\t{3}\n{4}\t{5}\tL\t{6}\n{4}\t{5}\tR\t{7}\n".format("Strain", "Chromosome", "Contig_side", "Telom_length", strain, chrom, left_tel, right_tel)) # file doesn't exist yet, write a header

# main program
chr_start_end(sys.argv[1])






### PRELIMINARY TRIALS
# from io import StringIO
## function#1: returning strain name , each chromosome name and their first 1000 bp on each side
# def chr_start_end(genome_fasta):
#     file_name = genome_fasta.split(".", maxsplit = 1)
#     strain = file_name[0]
#     seq = ""
#     with open(genome_fasta, "r") as filin:
#         for line in filin:
#             if line.startswith(">"):
#                 chr_info = line
#             else:
#                 seq += line.strip()
#                 start_left = (">{} left_end\n{}\n".format(strain, seq[:1000]))
#                 sequ_right = Seq(seq[-1000:])
#                 start_right = (">{} right_end\n{}\n".format(strain, sequ_right.reverse_complement()))
#     return strain, start_left, start_right

# ## function#2: sliding window counting the number of C+A bases
# def CA_window(sequence):
#     for record in SeqIO.parse(StringIO(sequence), "fasta"):
#         count = 0
#         for i in range(0, len(record.seq)-9):
#             #print(record.id + "\n")
#             #print(str(record.seq[i:i+10] + "\n"))
#             mot = str(record.seq[i:i+10] + "\n")
#             if mot.count("C") + mot.count("A") >= 9:
#                 count += 1
#             else:
#                 break
#         break
#     print("telom length = ", 10+count-1)            


# ## main program
# strain, left, right = chr_start_end(sys.argv[1])
# print("-------------------------------",
#     "\n",
#     strain,
#     "\n",
#     "-------------------------------",
#     "\n",
#     left, 
#     "\n",
#      "-------------------------------",
#     "\n",
#     right, 
#     "\n",
#     "==========================================================================",
#      "\n")

# CA_window(left)

###################################
###################################
## OPTION WITH LIST OF FASTA FILES

#import glob
#import re

## Make a list of all fasta files in the current directory
# fasta_files = []
# for file in glob.glob("*.fasta"):
#     fasta_files.append(file)

## définition of functions:
## function_1: returning a dictionary with strain name , each chromosome name and their extreme 1000 bp on each strand
# def chr_start_end(genome_fasta):
#     file_name = genome_fasta.split(".", maxsplit = 1)
#     strain = file_name[0]
#     seq = ""
#     dico_seq = {}
#     count_contig = 0
#     with open(genome_fasta, "r") as filin:
#         for line in filin:
#             if line.startswith(">"):
#                 count_contig += 1
#                 dico_seq["chr_info: {}".format(count_contig)] = line
#             else:
#                 seq += line.strip()
#                 seq_start = Seq(seq[:1000])
#                 dico_seq["start_contig: {}".format(count_contig)] = str(seq_start)
#                 seq_end = Seq(seq[-1000:])
#                 revcomp_end = seq_end.reverse_complement()
#                 dico_seq["end_contig: {}".format(count_contig)] = str(revcomp_end)
#     return strain, dico_seq

# ## main program
# for i in range(len(fasta_files)):
#     strain, chr_ends = chr_start_end(fasta_files[i])
#     print(chr_ends.keys())
#     # print("-------------------------------",
#     # "\n",
#     # strain,
#     # "\n",
#     # "-------------------------------",
#     # "\n",
#     # chr_ends, 
#     # "\n",
#     # "==========================================================================",
#     #  "\n")
#     print("_+_+_+_+_+_+_+_+_+_+_")
#     chr_ends[re.compile('e.*')] = 'word that starts with e'
#     # def find_matching_regexen(word, dicts=chr_ends):
#     #     return [description for regex, description in dicts.items() if regex.match(word)]
#     # print(find_matching_regexen('end'))
