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
        description="This program determines the location and the size of telomeric repeats\
        from genome assemblies. It runs both on single and multiple (multi)fasta file(s)\
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
    parser.add_argument(
        "-e",
        "--entropy_threshold",
        default=1.1,
        type=float,
        help="Entropy threshold for telomere prediction.",
    )
    parser.add_argument(
        "-n",
        "--polynuc_threshold",
        default=0.7,
        type=float,
        help="Poly-nucleotide threshold for telomere prediction.",
    )

    return parser.parse_args()


def get_strain_name(filename):
    """Function to get the strain name from the path to the fasta file"""
    filepath = Path(filename)
    return filepath.stem


def sliding_window(sequence, start, end, size):
    """Apply a sliding window of length = size to a sequence from start to end"""
    if size > len(sequence):
        sys.exit("The window size must be smaller than the sequence")
    for i in range(start, end - (size - 1)):
        window = str(sequence[i : i + size])
        yield i, window


def base_compos(sequence, base):
    """Return the number of base in the sequence"""
    count = Counter(sequence)[base]
    return count


def count_polynuc_occurence(window, polynucleotide_list):
    """Get polynucleotide proportion along the sliding window. 
    All polynucleotide should be of the same size.
    """
    if window in polynucleotide_list:
        return 1
    else:
        return 0


def get_polynuc(window, dinuc_list):
    sum_dinuc = 0
    for _, sub_window in sliding_window(window, 0, len(window), 2):
        sum_dinuc += count_polynuc_occurence(sub_window, dinuc_list)
    freq_dinuc = sum_dinuc / (len(window) - 1)
    return freq_dinuc


def get_skewness(window):
    """Get AT, GC skewness from a sequence"""

    if base_compos(window, "A") + base_compos(window, "T") == 0:
        at_skew = None
    else:
        at_skew = (base_compos(window, "A") - base_compos(window, "T")) / (
            base_compos(window, "A") + base_compos(window, "T")
        )

    if base_compos(window, "G") + base_compos(window, "C") == 0:
        gc_skew = None
    else:
        gc_skew = (base_compos(window, "G") - base_compos(window, "C")) / (
            base_compos(window, "G") + base_compos(window, "C")
        )

    if at_skew is None or gc_skew is None:
        skewness = None
    else:
        skewness = (
            (base_compos(window, "A") - base_compos(window, "T"))
            - (base_compos(window, "G") - base_compos(window, "C"))
        ) / len(window)

    return skewness


def get_cg_skew(window):
    """Get CG skewn from a sequence"""

    if base_compos(window, "G") + base_compos(window, "C") == 0:
        cg_skew = None
    else:
        cg_skew = (base_compos(window, "C") - base_compos(window, "G")) / (
            base_compos(window, "G") + base_compos(window, "C")
        )

    return cg_skew


def get_entropy(window):
    """Calculate frequency (probability) of nt in window"""
    entropy = 0

    for base in ["A", "T", "G", "C"]:
        if window.count(base) == 0:
            proba_base = 0
        else:
            freq_base = window.count(base) / len(window)
            proba_base = -(freq_base * np.log(freq_base))

        entropy += proba_base

    return entropy


def get_norm_freq_base(window):
    """Calculate the difference between observed and expected base composition of the window"""

    fexp_A_T = 0.617
    fexp_G_C = 0.383

    skew_norm = (
        (fexp_A_T - (base_compos(window, "A") / len(window))) ** 2
        - (fexp_A_T - (base_compos(window, "T") / len(window))) ** 2
        - (fexp_G_C - (base_compos(window, "G") / len(window))) ** 2
        + (fexp_G_C - (base_compos(window, "C") / len(window))) ** 2
    )

    return skew_norm


def get_add_freq_diff(window):
    """Calculate the difference between observed and expected base composition of the window"""

    fexp_A_T = 0.617
    fexp_G_C = 0.383

    max_diff = (
        (fexp_A_T - (base_compos(window, "A") / len(window))) ** 2
        + (fexp_A_T - (base_compos(window, "T") / len(window))) ** 2
        + (fexp_G_C - (base_compos(window, "G") / len(window))) ** 2
        - (fexp_G_C - (base_compos(window, "C") / len(window))) ** 2
    )

    return max_diff


def get_freq_norm_C(window):
    """Calculate the difference between observed and expected frequence of C in the window"""

    fexp_G_C = 0.383

    freq_norm_C = fexp_G_C - (base_compos(window, "C") / len(window))

    return freq_norm_C


def get_freq_norm_T(window):
    """Calculate the difference between observed and expected frequence of T in the window"""

    fexp_A_T = 0.617

    freq_norm_T = fexp_A_T - (base_compos(window, "T") / len(window))

    return freq_norm_T


def get_chi2(window):
    """Calculate the difference between observed and expected base composition of the window"""

    fexp_A_T = 0.617
    fexp_G_C = 0.383

    sum_freq = 0

    for base in ["A", "T"]:
        sum_freq += (fexp_A_T - (window.count(base) / len(window))) ** 2
    for base in ["G", "C"]:
        sum_freq += (fexp_G_C - (window.count(base) / len(window))) ** 2

    chi2 = sum_freq / 3

    return chi2


def compute_metrics(window, dinuc_list=["AC", "CA", "CC"]):

    metrics = {
        # "skew": get_skewness(window),
        # "cg_skew": get_cg_skew(window),
        "entropy": get_entropy(window),
        "polynuc": get_polynuc(window, dinuc_list),
        # "skew_norm": get_norm_freq_base(window),
        # "chi2": get_chi2(window),
        # "freq_norm_T": get_freq_norm_T(window),
        # "freq_norm_C": get_freq_norm_C(window),
        # "max_diff": get_add_freq_diff(window),
    }

    return metrics


def get_consecutive_groups(df_chrom):
    """From the raw dataframe get start and end of each telom==1 groups.
    Applied to detect start and end of telomere in nucleotide positions.
    """
    df = df_chrom.reset_index()
    chrom_groups = {}
    for strand in ["W", "C"]:
        nums = list(
            df.query("(level_3==@strand) and (predict_telom==1)").level_2
        )
        nums = sorted(set(nums))
        gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s + 1 < e]
        edges = iter(nums[:1] + sum(gaps, []) + nums[-1:])
        chrom_groups[strand] = list(zip(edges, edges))

        # [{"start":x, "end":y} for x, y in b["W"]]

    return chrom_groups


def classify_telomere(interval_chrom, chrom_len):
    """From a list of tuples obtained from get_consecutive_groups, identify if
    interval corresponds to terminal or interal telomere
    """
    classif_dict_list = []

    interval_W = interval_chrom["W"][:]
    if interval_W == []:
        classif_dict_list.append(
            {"start": None, "end": None, "side": "Left", "type": "term"}
        )
        classif_dict_list.append(
            {"start": None, "end": None, "side": "Left", "type": "intern"}
        )
    elif min(interval_W)[0] == 0:
        classif_dict_list.append(
            {
                "start": 0 + 1,
                "end": min(interval_W)[1] + 1 + 19,
                "side": "Left",
                "type": "term",
            }
        )
        interval_W.remove(min(interval_W))
        for interval in interval_W:
            classif_dict_list.append(
                {
                    "start": interval[0] + 1,
                    "end": interval[1] + 1 + 19,
                    "side": "Left",
                    "type": "intern",
                }
            )
    else:
        for interval in interval_W:
            classif_dict_list.append(
                {
                    "start": interval[0] + 1,
                    "end": interval[1] + 1 + 19,
                    "side": "Left",
                    "type": "intern",
                }
            )

    interval_C = interval_chrom["C"][:]
    if interval_C == []:
        classif_dict_list.append(
            {"start": None, "end": None, "side": "Right", "type": "term"}
        )
        classif_dict_list.append(
            {"start": None, "end": None, "side": "Right", "type": "intern"}
        )
    elif max(interval_C)[1] == (chrom_len - 1):
        classif_dict_list.append(
            {
                "start": max(interval_C)[0] + 1,
                "end": max(interval_C)[1] + 1 + 19,
                "side": "Right",
                "type": "term",
            }
        )
        interval_C.remove(max(interval_C))
        for interval in interval_C:
            classif_dict_list.append(
                {
                    "start": interval[0] + 1,
                    "end": interval[1] + 1 + 19,
                    "side": "Right",
                    "type": "intern",
                }
            )

    else:
        for interval in interval_C:
            classif_dict_list.append(
                {
                    "start": interval[0] + 1,
                    "end": interval[1] + 1 + 19,
                    "side": "Right",
                    "type": "intern",
                }
            )

    return classif_dict_list


def plot_telom(telom_df):
    """Visualization of the watson and crick strand telomere location
    """
    df = telom_df.reset_index()
    for strand in ["W", "C"]:
        ax = (
            df.query("level_3==@strand")
            .loc[:, ["polynuc", "entropy", "predict_telom"]]
            .plot()
            .legend(loc="center left", bbox_to_anchor=(1, 0.5))
        )
        ax.set_title(strand)


def export_results(
    raw_df,
    telom_df,
    raw_outfile="raw_df.csv",
    telom_outfile="telom_df.csv",
    outdir="telofinder_results",
):
    """ Produce output table files 
    """
    outdir = Path(outdir)
    outdir.mkdir()
    raw_df.to_csv(outdir / raw_outfile)
    telom_df.to_csv(outdir / telom_outfile)


def run_on_single_fasta(fasta_path, polynuc_thres, entropy_thres):
    """Run the telomere detection algorithm on a single fasta file"""
    strain = get_strain_name(fasta_path)
    print("\n", "-------------------------------", "\n")
    print(f"file {strain} executed")

    df_list = []
    telo_df_list = []

    for seq_record in SeqIO.parse(fasta_path, "fasta"):
        seqW = str(seq_record.seq)
        revcomp = seq_record.reverse_complement()
        seqC = str(revcomp.seq)

        # TODO: add a 'start' value if program reads the sequence not from its beginning
        limit_seq = min(20000, len(seqW))

        seq_dict_W = {}
        seq_dict_C = {}

        for i, window in sliding_window(seqW, 0, limit_seq, 20):
            seq_dict_W[(strain, seq_record.name, i, "W")] = compute_metrics(
                window
            )

        df_W = pd.DataFrame(seq_dict_W).transpose()

        for i, window in sliding_window(seqC, 0, limit_seq, 20):
            seq_dict_C[
                (strain, seq_record.name, (len(seqC) - i - 1), "C")
            ] = compute_metrics(window)

        df_C = pd.DataFrame(seq_dict_C).transpose()

        ## NOT USED ANYMORE######
        # df_W["entropy_med"] = df_W.rolling(100, min_periods=1).entropy.median()
        # df_W["polynuc_med"] = df_W.rolling(100, min_periods=1).polynuc.median()
        # ##################

        # ## NOT USED ANYMORE######
        # df_C["entropy_med"] = df_C.rolling(100, min_periods=1).entropy.median()
        # df_C["polynuc_med"] = df_C.rolling(100, min_periods=1).polynuc.median()
        ##################

        df_chro = pd.concat([df_W, df_C])

        df_chro.loc[
            (df_chro["entropy"] < entropy_thres)
            & (df_chro["polynuc"] > polynuc_thres),
            "predict_telom",
        ] = 1.0

        df_chro["predict_telom"].fillna(0, inplace=True)

        telo_groups = get_consecutive_groups(df_chro)
        telo_list = classify_telomere(telo_groups, len(seq_record.seq))
        telo_df = pd.DataFrame(telo_list)
        telo_df["strain"] = strain
        telo_df["chrom"] = seq_record.name
        telo_df = telo_df[["strain", "chrom", "side", "type", "start", "end"]]

        df_list.append(df_chro)
        telo_df_list.append(telo_df)

        print(f"chromosome {seq_record.name}")

    df = pd.concat(df_list)
    telo_df = pd.concat(telo_df_list)
    telo_df["len"] = telo_df["end"] - telo_df["start"] + 1

    telo_df = telo_df.astype({"start": "Int64", "end": "Int64", "len": "Int64"})

    return df, telo_df


def run_on_fasta_dir(fasta_dir_path, polynuc_thres, entropy_thres):
    """Run iteratively the telemore detection algorithm on all fasta files in a directory"""
    raw_dfs = []
    telom_dfs = []

    for ext in ["*.fasta", "*.fas", "*.fa"]:
        for fasta in fasta_dir_path.glob(ext):

            raw_df, telom_df = run_on_single_fasta(
                fasta, polynuc_thres, entropy_thres
            )
            raw_dfs.append(raw_df)
            telom_dfs.append(telom_df)

    total_raw_df = pd.concat(raw_dfs)
    total_telom_df = pd.concat(telom_dfs)

    return total_raw_df, total_telom_df


def run_telofinder(fasta_path, polynuc_thres, entropy_thres):
    """Run telofinder on a single fasta file or on a fasta directory"""
    fasta_path = Path(fasta_path)

    if fasta_path.is_dir():
        print(
            f"Running in iterative mode on all '*.fasta', '*.fas', '*.fa' files in '{fasta_path}'"
        )
        raw_df, telom_df = run_on_fasta_dir(
            fasta_path, polynuc_thres, entropy_thres
        )
        export_results(raw_df, telom_df)
        return raw_df, telom_df

    elif fasta_path.is_file():
        print(f"Running in single fasta mode on '{fasta_path}'")

        raw_df, telom_df = run_on_single_fasta(
            fasta_path, polynuc_thres, entropy_thres
        )
        export_results(raw_df, telom_df)
        return raw_df, telom_df
    else:
        raise IOError(f"'{fasta_path}' is not a directory or a file.")


# Main program
if __name__ == "__main__":
    args = parse_arguments()
    output_exists(args.force)
    run_telofinder(
        args.fasta_path, args.polynuc_threshold, args.entropy_threshold
    )
