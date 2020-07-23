import sys
from Bio import SeqIO
import os.path
import argparse
from pathlib import Path
import pandas as pd


def output_exists(force):
    """Function to test if the output file already exists, force overwriting the output file or
    exit program.
    """
    file_exists = os.path.isfile("telom_length.csv")
    if file_exists:
        if force:
            print("Replacing the existing file.")
        else:
            sys.exit(
                print(
                    "\n",
                    "Warning!!! A file called 'telom_length.csv' already exists.",
                    "\n",
                    "delete this file before running the script or use the option --force",
                    "\n",
                )
            )


# TODO: Add a parameter for output file name
def parse_arguments():
    """Function to parse and reuse the arguments of the command line"""
    parser = argparse.ArgumentParser(
        description="This program determines telomere length at all sequence ends from a single or multiple (multi)fasta file(s)\
        and outputs a csv file called 'telom_length.csv'"
    )
    parser.add_argument(
        "fasta_path",
        help="Path to a single fasta file or to a directory containing multiple fasta files.",
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Automatically replace the output csv file if present.",
    )
    return parser.parse_args()


def get_strain_name(filename):
    """Function to get the strain name from the path to the fasta file"""
    filepath = Path(filename)
    return filepath.stem


# FIXME: verify that offset can be equal to i+1
def get_offset(sequence):
    """estimate the size of the offset sequence before the telomere sequence"""
    offset = 0
    limit = min(1500, len(sequence) - 19)
    for i in range(0, limit):
        mot = str(sequence[i : i + 20])
        # TODO: Should this be equal to set(seq) = {"A", "T"} ? ie all As or Cs ?
        if mot.count("C") + mot.count("A") < 19:
            print("OKKK")
            offset = i + 1
            if offset == limit:
                offset = 0

        else:
            break
    return offset


# TODO: Maybe integrate offset directly in serach for telomere
def get_telom_size(sequence):
    """Function to calculate telomere length at all contig ends"""
    tel_size = 0
    offset = 0
    limit = min(1500, len(sequence) - 20)
    for i in range(0, len(sequence) - 19):
        if i == limit and tel_size == 0:
            offset = 0
            tel_size = 0
            return offset, tel_size

        mot = str(sequence[i : i + 20])
        if (
            mot.count("C") >= 8
            and mot.count("A") >= 3
            and mot.count("C") + mot.count("A") >= 17
            and "AAAA" not in mot
        ):
            tel_size += 1
        else:
            if tel_size != 0:
                tel_size = 20 + tel_size - 3
                return offset, tel_size

            offset += 1


def get_nucleotide_proportions(sequence):
    """ Get nucleotide proportion along the sliding window
    """
    from collections import Counter

    nucleotide_proportions = {}

    for i in range(0, len(sequence.seq) - 19):
        mot = str(sequence[i : i + 20])
        if (
            mot.count("C") >= 8
            and mot.count("A") >= 3
            and mot.count("C") + mot.count("A") >= 17
            and "AAAA" not in mot
        ):
            nucleotide_proportions[i] = 1
        else:
            nucleotide_proportions[i] = 0

    # for i in range(0, len(sequence.seq) - 1):

    #     mot = str(sequence.seq[i : i + 2])
    #     if mot in ["AC", "CA", "CC"]:
    #         nucleotide_proportions[i] = 1
    #     else:
    #         nucleotide_proportions[i] = 0

    nucleotide_df = pd.DataFrame(nucleotide_proportions, index=[sequence.name])
    return nucleotide_df


# FIXME: missing docstring
def generate_output(
    strain, chrom, left_tel, right_tel, left_offset, right_offset, output_file
):
    print(
        "\n",
        "=================================",
        "\n",
        strain,
        "\n",
        "=================================",
        "\n",
    )
    print(chrom)
    print("---------------------")
    print("left telom length = ", left_tel)
    print("left offset = ", left_offset)
    print("right telom length = ", right_tel)
    print("right offset = ", right_offset)
    print("\n", "-------------------------------", "\n")

    # TODO: Could use pandas DataFrames here
    ## save the results in a csv file
    file_exists = os.path.isfile(output_file)
    with open("telom_length.csv", "a") as filout:
        if file_exists:
            filout.write(
                "{0}\t{1}\tL\t{2}\t{4}\n{0}\t{1}\tR\t{3}\t{5}\n".format(
                    strain, chrom, left_tel, right_tel, left_offset, right_offset
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
                    left_offset,
                    right_offset,
                )
            )  # file doesn't exist yet, write a header


def run_single_or_iterative(fasta_path):
    """Check if fasta_path is a single file or a directory and run the telomere
    detection sequence accordingly."""

    if Path(fasta_path).is_file():
        print(f"'{fasta_path}' is a file. Running in single mode.")
        strain = get_strain_name(fasta_path)

        dinucleotide_df_list = []
        for seq_record in SeqIO.parse(fasta_path, "fasta"):
            # TODO: This should be incorporated in a function (1)
            revcomp = seq_record.reverse_complement()
            dinucleotide_df_list.append(get_nucleotide_proportions(seq_record))
            # left_offset = get_offset(seq_record.seq)
            # left_offset, left_tel = get_telom_size(seq_record.seq)
            # right_offset = get_offset(revcomp)
            # right_offset, right_tel = get_telom_size(revcomp)
            # generate_output(
            #     strain,
            #     seq_record.id,
            #     left_tel,
            #     right_tel,
            #     left_offset,
            #     right_offset,
            #     "telom_length.csv",
            # )
        return pd.concat(dinucleotide_df_list)

    if Path(fasta_path).is_dir():
        print(f"'{fasta_path}' is a directory. Running in iterative mode.")

        for ext in ["*.fasta", "*.fas", "*.fa"]:
            for fasta in Path(fasta_path).glob(ext):

                print(fasta)
                strain = get_strain_name(fasta)

                for seq_record in SeqIO.parse(fasta_path, "fasta"):
                    # TODO: This should be incorporated in a function, see TODO (1)
                    revcomp = seq_record.reverse_complement()
                    left_offset = get_offset(seq_record.seq)
                    left_tel = get_telom_size(seq_record.seq, left_offset)
                    right_offset = get_offset(revcomp)
                    right_tel = get_telom_size(revcomp, right_offset)
                    generate_output(
                        strain,
                        seq_record.id,
                        left_tel,
                        right_tel,
                        left_offset,
                        right_offset,
                        "telom_length.csv",
                    )


# Main program
if __name__ == "__main__":
    args = parse_arguments()
    output_exists(args.force)
    run_single_or_iterative(args.fasta_path)
