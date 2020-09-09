import sys
from Bio import SeqIO
import os.path
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from collections import Counter


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


# TODO: idea for later have a sliding window checking match for polynucleotide
# of various sizes and get a score from that: ie AC: 1 ACC: 1.5, ACCA: 2 etc ...

# def get_polynuc_occurence(sequence, polynucleotide_list):
#     """ Get polynucleotide proportion along the sliding window
#     All polynucleotide should be of the same size.
#     """

#     polynucleotide_occurence = {}

#     for i in range(0, len(sequence) - 1):

#         mot = str(sequence[i : i + 2])
#         if mot in polynucleotide_list:
#             polynucleotide_occurence[i] = 1
#         else:
#             polynucleotide_occurence[i] = 0

#     return polynucleotide_occurence


def sliding_window(sequence, start, end, size):
    if size > len(sequence):
        sys.exit("The window size must be smaller than the sequence")
    for i in range(start, end - (size - 1)):
        window = str(sequence[i : i + size])
        yield window


def get_pattern_occurences(window):
    """Presence/absence of the telomere pattern in a sliding window"""
    if (
        window.count("C") >= 8
        and window.count("A") >= 3
        and window.count("C") + window.count("A") >= 17
        and "AAAA" not in window
    ):
        return 1
    else:
        return 0


def count_polynuc_occurence(window, polynucleotide_list):
    """ Get polynucleotide proportion along the sliding window
    All polynucleotide should be of the same size.
    """
    if window in polynucleotide_list:
        return 1
    else:
        return 0


def get_polynuc(window, dinuc_list):
    sum_dinuc = 0
    # print(window)
    for sub_window in sliding_window(window, 0, len(window), 2):
        sum_dinuc += count_polynuc_occurence(sub_window, dinuc_list)
    freq_dinuc = sum_dinuc / (len(window) - 1)
    return freq_dinuc


def get_skewness(window):
    """ Get AT, GC skewness from a sequence
    """

    base_compos = Counter(window)
    a = base_compos["A"]
    t = base_compos["T"]
    g = base_compos["G"]
    c = base_compos["C"]

    # a = window.count("A")
    # t = window.count("T")
    # g = window.count("G")
    # c = window.count("C")

    if (a + t) == 0:
        at_skew = None
    else:
        at_skew = (a - t) / (a + t)

    if (g + c) == 0:
        gc_skew = None
    else:
        gc_skew = (g - c) / (g + c)

    if at_skew is None or gc_skew is None:
        skewness = None
    else:
        skewness = ((a - t) - (g - c)) / len(window)

    # "at_skew": at_skew,
    # "gc_skew": gc_skew,

    return skewness


def get_entropy(window):
    """ Calculate frequency (probability) of nt in window.
    """
    entropy = 0

    for base in ["A", "T", "G", "C"]:
        if window.count(base) == 0:
            proba_base = 0
        else:
            freq_base = window.count(base) / len(window)
            proba_base = -(freq_base * np.log(freq_base))

        entropy += proba_base

    return entropy


# FIXME return chi2 or skew_norm?
def get_norm_freq_base(window):
    """ Calculate the difference between observed and expected base composition of the window"""

    base_compos = Counter(window)

    fexp_A_T = 0.617
    fexp_G_C = 0.383

    sum_freq = 0

    for base in ["A", "T"]:
        sum_freq += (fexp_A_T - (base_compos[base] / len(window))) ** 2
    for base in ["G", "C"]:
        sum_freq += (fexp_G_C - (base_compos[base] / len(window))) ** 2

    skew_norm = (
        (fexp_A_T - (base_compos["A"] / len(window))) ** 2
        - (fexp_A_T - (base_compos["T"] / len(window))) ** 2
        - (fexp_G_C - (base_compos["G"] / len(window))) ** 2
        + (fexp_G_C - (base_compos["C"] / len(window))) ** 2
    )

    chi2 = sum_freq / 3

    return skew_norm


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
                    strain,
                    chrom,
                    left_tel,
                    right_tel,
                    left_offset,
                    right_offset,
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


def run_on_single_fasta(fasta_path):
    """Run the telomere detection algorithm on a single fasta file.
    """

    strain = get_strain_name(fasta_path)

    polynucleotide_dict = {}

    seq_dict = {}

    for seq_record in SeqIO.parse(fasta_path, "fasta"):
        # TODO: make parameters for window's start, end and size
        # TODO: add a 'start' value if program reads the sequence not from its beginning
        limit_seq = min(20000, len(seq_record.seq))
        for i, window in enumerate(
            sliding_window(seq_record.seq, 0, limit_seq, 20)
        ):

            seq_dict[(strain, seq_record.name, i)] = {
                "pattern": get_pattern_occurences(window),
                "skew": get_skewness(window),
                "entropy": get_entropy(window),
                "polynuc": get_polynuc(window, ["AC", "CA", "CC"]),
                "chi2": get_norm_freq_base(window),
            }

        revcomp = seq_record.reverse_complement()
        left_offset, left_tel = get_telom_size(seq_record.seq)
        right_offset, right_tel = get_telom_size(revcomp)
        generate_output(
            strain,
            seq_record.id,
            left_tel,
            right_tel,
            left_offset,
            right_offset,
            "telom_length.csv",
        )

    df = pd.DataFrame(seq_dict).transpose()
    df.loc[:, "combined"] = (
        df.loc[:, "chi2"]
        + df.loc[:, "skew"]
        + df.loc[:, "polynuc"]
        - df.loc[:, "entropy"]
    )
    df.loc[:, "skew-ent"] = df.loc[:, "skew"] - df.loc[:, "entropy"]

    return df


def run_on_fasta_dir(fasta_dir_path):
    """Run iteratively the telemore detection algorithm on all fasta files in a
    directory
    """
    telom_dfs = []

    for ext in ["*.fasta", "*.fas", "*.fa"]:
        for fasta in fasta_dir_path.glob(ext):

            telom_dfs.append(run_on_single_fasta(fasta))

    total_telom_df = pd.concat(telom_dfs)

    return total_telom_df


def run_telofinder(fasta_path):
    """ Run telofinder on a single fasta file or on a fasta directory
    """

    fasta_path = Path(fasta_path)

    if fasta_path.is_dir():
        print(
            f"Running in iterative mode on all '*.fasta', '*.fas', '*.fa' files in '{fasta_path}'"
        )
        return run_on_fasta_dir(fasta_path)
    elif fasta_path.is_file():
        print(f"Running in single fasta mode on '{fasta_path}'")
        return run_on_single_fasta(fasta_path)
    else:
        raise IOError(f"'{fasta_path}' is not a directory or a file.")


# Main program
if __name__ == "__main__":
    args = parse_arguments()
    output_exists(args.force)
    run_telofinder(args.fasta_path)
