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


# def get_pattern_occurences(sequence):
#     "Presence/absence coding of the telomere pattern in a sliding window"

#     pattern_occurence = {}

#     for i in range(0, len(sequence) - 19):
#         mot = str(sequence[i : i + 20])
#         if (
#             mot.count("C") >= 8
#             and mot.count("A") >= 3
#             and mot.count("C") + mot.count("A") >= 17
#             and "AAAA" not in mot
#         ):
#             pattern_occurence[i] = 1
#         else:
#             pattern_occurence[i] = 0

#     return pattern_occurence


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


def get_polynuc_occurence(window, polynucleotide_list):
    """ Get polynucleotide proportion along the sliding window
    All polynucleotide should be of the same size.
    """
    if window in polynucleotide_list:
        return 1
    else:
        return 0


def get_skewness(window):
    """ Get AT, GC skewness from a sequence
    """
    a = window.count("A")
    t = window.count("T")
    g = window.count("G")
    c = window.count("C")

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

    skewness_stats[i] = skewness
    # "at_skew": at_skew,
    # "gc_skew": gc_skew,


return at_skew, gc_skew, skewness


# def get_skewness(sequence):
#     """ Get AT, GC skewness from a sequence
#     """

#     skewness_stats = {}

#     for i in range(0, len(sequence) - 19):
#         mot = str(sequence[i : i + 20])

#         a = mot.count("A")
#         t = mot.count("T")
#         g = mot.count("G")
#         c = mot.count("C")

#         if (a + t) == 0:
#             at_skew = None
#         else:
#             at_skew = (a - t) / (a + t)

#         if (g + c) == 0:
#             gc_skew = None
#         else:
#             gc_skew = (g - c) / (g + c)

#         if at_skew is None or gc_skew is None:
#             skewness = None
#         else:
#             skewness = ((a - t) - (g - c)) / len(mot)

#         skewness_stats[i] = skewness
#         # "at_skew": at_skew,
#         # "gc_skew": gc_skew,

#     return skewness_stats


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


def run_single_or_iterative(fasta_path):
    """Check if fasta_path is a single file or a directory and run the telomere
    detection sequence accordingly."""

    if Path(fasta_path).is_file():
        print(f"'{fasta_path}' is a file. Running in single mode.")
        strain = get_strain_name(fasta_path)

        polynucleotide_dict = {}
        pattern_dict = {}
        for seq_record in SeqIO.parse(fasta_path, "fasta"):
            # TODO: This should be incorporated in a function (1)
            polynucleotide_dict[seq_record.name] = get_polynuc_occurence(
                seq_record.seq, ["AC", "CA", "CC"]
            )
            pattern_dict[seq_record.name] = get_pattern_occurences(
                seq_record.seq
            )
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

            skewness_stats = get_skewness(seq_record.seq)

        return (
            pd.DataFrame(polynucleotide_dict).transpose(),
            pd.DataFrame(pattern_dict).transpose(),
            pd.DataFrame(skewness_stats).tranpose(),
        )

    if Path(fasta_path).is_dir():
        print(f"'{fasta_path}' is a directory. Running in iterative mode.")

        polynucleotide_dict = {}
        pattern_dict = {}
        skewness_stats = {}

        for ext in ["*.fasta", "*.fas", "*.fa"]:
            for fasta in Path(fasta_path).glob(ext):

                print(fasta)
                strain = get_strain_name(fasta)

                for seq_record in SeqIO.parse(fasta, "fasta"):
                    # TODO: This should be incorporated in a function, see TODO (1)
                    polynucleotide_dict[
                        fasta.stem, seq_record.name
                    ] = get_polynuc_occurence(
                        seq_record.seq, ["AC", "CA", "CC"]
                    )
                    pattern_dict[
                        fasta.stem, seq_record.name
                    ] = get_pattern_occurences(seq_record.seq)
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
                    skewness_stats[fasta.stem, seq_record.name] = get_skewness(
                        seq_record.seq
                    )

        return (
            pd.DataFrame(polynucleotide_dict).transpose(),
            pd.DataFrame(pattern_dict).transpose(),
            pd.DataFrame(skewness_stats).transpose(),
        )


# Main program
if __name__ == "__main__":
    args = parse_arguments()
    output_exists(args.force)
    run_single_or_iterative(args.fasta_path)
