import sys
from Bio.Seq import Seq
from Bio import SeqIO
import os.path

## !!! MAKE SURE TO DELETE THE OUTPUT FILE "telom_length.csv" FROM YOUR DIRECTORY BEFORE RUNNING THE SCRIPT !!!
# TODO: Add a function testing if the file exists and exits program if it does.

## command to launch the script: $python 00_analyze_telom_length.py file.fasta
# TODO: To be added in the documentation of the command line interface (CLI)
# (using argparse for example)

## define the genme fasta file
if len(sys.argv) != 2:
    sys.exit(
        "ERROR: you must specify the name of genome fasta file or its path as the first argument"
    )
# TODO: To be added in doc/setup of argparse

## function#1: return the strain name, each chromosome name and their telomere length on each side
def chr_start_end(genome_fasta):

    ## get strain name from fasta files either from a path or in the src directory
    ## from a path
    # TODO: This can be a function get_strain_name, which will only extract the
    # strain name from the filename and return it. You can also use python
    # module pathlib to improve this part

    if "/" in genome_fasta:
        file_name = genome_fasta.split("/")
        assembly = file_name[-1].split(".")
        strain = assembly[0]
    ## from the src directory
    else:
        file_name = genome_fasta.split(".")
        strain = file_name[0]

    print(
        "\n",
        "=================================",
        "\n",
        strain,
        "\n",
        "=================================",
        "\n",
    )

    ## parse the multifasta file and reverse complement each sequence
    for seq_record in SeqIO.parse(genome_fasta, "fasta"):
        revcomp = seq_record.reverse_complement()

        count = 0
        revcount = 0
        offset = 0
        revoffset = 0

        ## set up a limit of 1000 nt at contig ends to browse until we find a telomere sequence
        if len(seq_record.seq) < 1500:
            limit = len(seq_record.seq) - 9
        else:
            limit = 1500

        ### LEFT TELOMERE
        ## estimate the size of the offset sequence before the telomere sequence
        for i in range(0, limit):
            mot = str(seq_record.seq[i : i + 20] + "\n")
            if mot.count("C") + mot.count("A") < 19:
                offset += 1
            else:
                break
        ## estimate the size of the telomere sequence
        for j in range(offset, len(seq_record.seq) - 19):
            mot = str(seq_record.seq[j : j + 20] + "\n")
            if (
                mot.count("C") >= 8
                and mot.count("A") >= 3
                and mot.count("C") + mot.count("A") >= 17
                and "AAAA" not in mot
            ):
                count += 1
            else:
                break

        ### RIGHT TELOMERE
        ## estimate the size of the offset sequence before the telomere sequence
        for i in range(0, limit):
            revmot = str(revcomp.seq[i : i + 20] + "\n")
            if revmot.count("C") + revmot.count("A") < 19:
                revoffset += 1
            else:
                break
        ## estimate the size of the telomere sequence
        for j in range(revoffset, len(seq_record.seq) - 19):
            revmot = str(revcomp.seq[j : j + 20] + "\n")
            if (
                revmot.count("C") >= 8
                and revmot.count("A") >= 3
                and revmot.count("C") + revmot.count("A") >= 17
                and "AAAA" not in revmot
            ):
                revcount += 1
            else:
                break

        ## definition of telomere and offset lengths
        if count == 0:
            left_tel = 0
            offset = 0
        elif count != 0 and offset != 0:
            left_tel = 20 + count - 3
            offset = offset + 1
        elif count != 0 and offset == 0:
            left_tel = 20 + count - 3

        if revcount == 0:
            right_tel = 0
            revoffset = 0
        elif revcount != 0 and revoffset != 0:
            right_tel = 20 + revcount - 3
            revoffset = revoffset + 1
        elif revcount != 0 and revoffset == 0:
            right_tel = 20 + revcount - 3

        ## stdout
        chrom = seq_record.id
        print(chrom)
        print("---------------------")
        print("left telom length = ", left_tel)
        print("left offset = ", offset)
        print("right telom length = ", right_tel)
        print("right offset = ", revoffset)
        print("\n", "-------------------------------", "\n")

        ## save the results in a csv file
        file_exists = os.path.isfile("telom_length.csv")
        with open("telom_length.csv", "a") as filout:
            if file_exists:
                filout.write(
                    "{0}\t{1}\tL\t{2}\t{4}\n{0}\t{1}\tR\t{3}\t{5}\n".format(
                        strain, chrom, left_tel, right_tel, offset, revoffset
                    )
                )
            else:
                filout.write(
                    "{0}\t{1}\t{2}\t{3}\t{4}\n{5}\t{6}\tL\t{7}\t{9}\n{5}\t{6}\tR\t{8}\t{10}\n".format(
                        "Strain",
                        "Chromosome",
                        "Contig_side",
                        "Telom_length",
                        "Offset",
                        strain,
                        chrom,
                        left_tel,
                        right_tel,
                        offset,
                        revoffset,
                    )
                )  # file doesn't exist yet, write a header


if __name__ == "__main__":
    chr_start_end(sys.argv[1])
