import sys
from Bio import SeqIO
import os.path
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from collections import Counter
import pybedtools
from multiprocessing import Pool
from functools import partial


def output_dir_exists(force):
    """Function to test if the 'telofinder_results' output diretory already exists.

    :param force: overwrites the output directory otherwise exits program --force or -f
    :return: nothing
    """
    dir_exists = os.path.isdir("telofinder_results")
    if dir_exists:
        if force:
            print("Replacing the existing 'telofinder_results' directory.")
        else:
            sys.exit(
                print(
                    "\n",
                    "Warning!!! A directory called 'telofinder_results' already exists.",
                    "\n",
                    "delete this directory before running the script or use the option --force",
                    "\n",
                )
            )


def parse_arguments():
    """Function to parse and reuse the arguments of the command line

    :param fasta_path: path to a single fasta file or to a directory containing multiple fasta files
    :param force: force optional, overwrites the output directory otherwise exits program, optional
    :param entropy_threshold: optional, default = 0.8 
    :param polynuc_threshold: optional, default = 0.8
    :param nb_scanned_nt: number of scanned nucleotides at each chromosome end, optional, default = 20 000
    :param threads: Number of threads to use. Multithreaded calculations currently occurs at the level of sequences within a fasta file."
    :param raw: Outputs raw_df.csv containing the values of all sliding windows
    :return: parser arguments
    """
    parser = argparse.ArgumentParser(
        description="This program determines the location and the size of telomeric repeats\
        from genome assemblies. It runs both on single and multiple (multi)fasta file(s).\
        It outputs 3 csv and 2 bed files containing the telomere calls and their coordinates\
        either as raw output or after merging consecutive calls"
    )
    parser.add_argument(
        "fasta_path",
        help="Path to a single (multi)fasta file or to a directory containing multiple fasta files.",
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Automatically replace the output 'telofinder_results' directory if present.",
    )
    parser.add_argument(
        "-e",
        "--entropy_threshold",
        default=0.8,
        type=float,
        help="Entropy threshold for telomere prediction. default=0.8",
    )
    parser.add_argument(
        "-n",
        "--polynuc_threshold",
        default=0.8,
        type=float,
        help="Poly-nucleotide threshold for telomere prediction. default=0.8",
    )
    parser.add_argument(
        "-s",
        "--nb_scanned_nt",
        default=20000,
        type=int,
        help="Number of nucleotides scanned in sliding window starting from each sequence\
    extrimity. If set to -1, the whole sequence will be scanned. default=20000",
    )
    parser.add_argument(
        "-t",
        "--threads",
        default=1,
        type=int,
        help="Number of threads to use. Multithreaded calculations currently occurs\
    at the level of sequences within a fasta file.",
    )
    parser.add_argument(
        "-r",
        "--raw",
        action="store_true",
        help="Outputs the raw dataframe (raw_df.csv) containing the values of all sliding windows.",
    )

    return parser.parse_args()


def get_strain_name(filename):
    """Function to get the strain name from the name of the fasta file

    :param filename: path of fasta file
    :return: sequence name
    """
    filepath = Path(filename)
    return filepath.stem


def sliding_window(sequence, start, end, size):
    """Apply a sliding window of length = size to a sequence from start to end

    :param sequence: fasta sequence
    :param start: starting coordinate of the sequence
    :param end: ending coordinate of the sequence
    :param size: size of the sliding window
    :return: the coordinate and the sequence of the window
    """
    if size > len(sequence):
        sys.exit("The window size must be smaller than the sequence")
    for i in range(start, end - (size - 1)):
        window = str(sequence[i : i + size])
        yield i, window


def base_compos(sequence, base):
    """Counts the number of a given base in a sequence

    :param sequence: fasta sequence
    :param base: base to count in the sequence
    :return: the number of that base in the sequence
    """
    count = Counter(sequence)[base]
    return count


def count_polynuc_occurence(sub_window, polynucleotide_list):
    """Define presence of polynucleotide in the window.

    :param sub_window: the sequence of a sub_window
    :param polynuleotide_list: a list of polynucleotides. Note that all polynucleotides must be of the same size
    :return: a boolean for the presence of the sub_window in the polynucleotide list
    """
    if sub_window in polynucleotide_list:
        return 1
    else:
        return 0


def get_polynuc(window, polynucleotide_list):
    """ get the propbortion of polynuceotides in the window

    :param window: sliding window
    :param polynucleotide_list: a list of polynucleotides. Note that all polynucleotides must be of the same size
    :return: total polynucleotide proportion in the sliding window
    """
    sum_dinuc = 0
    for _, sub_window in sliding_window(window, 0, len(window), 2):
        sum_dinuc += count_polynuc_occurence(sub_window, polynucleotide_list)
    freq_dinuc = sum_dinuc / (len(window) - 1)
    return freq_dinuc


def get_entropy(window):
    """Calculate the entropy of the window DNA sequence

    :param window: sliding window
    :return: entropy value of the sequence window
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


def compute_metrics(window, polynucleotide_list=["AC", "CA", "CC"]):
    """Compute entropy and polynucleotide proportion in the sequence window

    :param window: sliding window
    :param polynucleotide_list: a list of polynucleotides, default value is ["AC", "CA", "CC"]
    :return: a dictionary of entropy and polynucleotide proportion of the sequence window
    """

    metrics = {
        "entropy": get_entropy(window),
        "polynuc": get_polynuc(window, polynucleotide_list),
        # "skew": get_skewness(window),
        # "cg_skew": get_cg_skew(window),
        # "skew_norm": get_norm_freq_base(window),
        # "chi2": get_chi2(window),
        # "freq_norm_T": get_freq_norm_T(window),
        # "freq_norm_C": get_freq_norm_C(window),
        # "max_diff": get_add_freq_diff(window),
    }

    return metrics


def get_consecutive_groups(df_chrom):
    """From the raw dataframe get start and end of each telomere window.
    Applied to detect start and end of telomere in nucleotide positions.
    """
    df = df_chrom.reset_index()
    chrom_groups = {}
    for strand in ["W", "C"]:
        nums = list(df.query("(level_3==@strand) and (predict_telom==1)").level_2)
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
                "start": max(interval_C)[0] + 1 - 19,
                "end": max(interval_C)[1] + 1,
                "side": "Right",
                "type": "term",
            }
        )
        interval_C.remove(max(interval_C))
        for interval in interval_C:
            classif_dict_list.append(
                {
                    "start": interval[0] + 1 - 19,
                    "end": interval[1] + 1,
                    "side": "Right",
                    "type": "intern",
                }
            )

    else:
        for interval in interval_C:
            classif_dict_list.append(
                {
                    "start": interval[0] + 1 - 19,
                    "end": interval[1] + 1,
                    "side": "Right",
                    "type": "intern",
                }
            )

    return classif_dict_list


def plot_telom(telom_df):
    """Plotting the telomere detection on both left and right chromosome ends
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
    raw_df, telom_df, merged_telom_df, raw, outdir="telofinder_results",
):
    """ Produce output table files
    """
    outdir = Path(outdir)
    try:
        outdir.mkdir()
    except FileExistsError:
        pass

    telom_df.to_csv(outdir / "telom_df.csv", index=False)
    merged_telom_df.to_csv(outdir / "merged_telom_df.csv", index=False)

    bed_df = telom_df[["chrom", "start", "end", "type"]].copy()
    bed_df.dropna(inplace=True)
    bed_df.to_csv(outdir / "telom.bed", sep="\t", header=None, index=False)

    merged_bed_df = merged_telom_df[["chrom", "start", "end", "type"]].copy()
    merged_bed_df.dropna(inplace=True)
    merged_bed_df.to_csv(outdir / "telom_merged.bed", sep="\t", header=None, index=False)

    if raw:
        raw_df.to_csv(outdir / "raw_df.csv", index=True)


def run_on_single_seq(seq_record, strain, polynuc_thres, entropy_thres, nb_scanned_nt):
    seqW = str(seq_record.seq)
    revcomp = seq_record.reverse_complement()
    seqC = str(revcomp.seq)

    if nb_scanned_nt == -1:
        limit_seq = len(seqW)
    else:
        limit_seq = min(nb_scanned_nt, len(seqW))

    seq_dict_W = {}
    seq_dict_C = {}

    for i, window in sliding_window(seqW, 0, limit_seq, 20):
        seq_dict_W[(strain, seq_record.name, i, "W")] = compute_metrics(window)

    df_W = pd.DataFrame(seq_dict_W).transpose()

    for i, window in sliding_window(seqC, 0, limit_seq, 20):
        seq_dict_C[(strain, seq_record.name, (len(seqC) - i - 1), "C")] = compute_metrics(
            window
        )

    df_C = pd.DataFrame(seq_dict_C).transpose()

    df_chro = pd.concat([df_W, df_C])

    df_chro.loc[
        (df_chro["entropy"] < entropy_thres) & (df_chro["polynuc"] > polynuc_thres),
        "predict_telom",
    ] = 1.0

    df_chro["predict_telom"].fillna(0, inplace=True)

    telo_groups = get_consecutive_groups(df_chro)
    telo_list = classify_telomere(telo_groups, len(seq_record.seq))
    telo_df = pd.DataFrame(telo_list)
    telo_df["chrom"] = seq_record.name
    telo_df["chrom_size"] = len(seq_record.seq)

    if telo_df["start"].isnull().sum() == 4:
        telo_df_merged = telo_df.copy()
    else:
        bed_df = telo_df[["chrom", "start", "end", "type"]].copy()
        bed_df.dropna(inplace=True)
        bed_df = bed_df.astype({"start": int, "end": int})
        bed_file = pybedtools.BedTool().from_dataframe(bed_df)
        bed_sort = bed_file.sort()
        bed_merge = bed_sort.merge(d=20)
        bed_df_merged = bed_merge.to_dataframe()
        telo_df_merged = pd.merge(
            bed_df_merged,
            telo_df.dropna()[["chrom", "side", "type", "start", "chrom_size"]],
            on=["chrom", "start"],
            how="left",
        )
        telo_df_merged.loc[telo_df_merged.end > len(seq_record.seq) - 20, "type"] = "term"
        telo_df_merged.loc[telo_df_merged.start < 20, "type"] = "term"

    telo_df_merged["strain"] = strain
    telo_df_merged = telo_df_merged[
        ["strain", "chrom", "side", "type", "start", "end", "chrom_size"]
    ]

    telo_df["strain"] = strain
    telo_df = telo_df[["strain", "chrom", "side", "type", "start", "end"]]

    print(f"chromosome {seq_record.name} done")

    return (df_chro, telo_df, telo_df_merged)


def run_on_single_fasta(fasta_path, polynuc_thres, entropy_thres, nb_scanned_nt, threads):
    """Run the telomere detection algorithm on a single fasta file

    :param fasta_path: path to fasta file
    :return: a tuple of df, telo_df and telo_df_merged
    """
    strain = get_strain_name(fasta_path)
    print("\n", "-------------------------------", "\n")
    print(f"file {strain} executed")

    partial_ross = partial(
        run_on_single_seq,
        strain=strain,
        polynuc_thres=polynuc_thres,
        entropy_thres=entropy_thres,
        nb_scanned_nt=nb_scanned_nt,
    )

    with Pool(threads) as p:

        results = p.map(partial_ross, SeqIO.parse(fasta_path, "fasta"))

    raw_df = pd.concat([r[0] for r in results])

    telo_df = pd.concat([r[1] for r in results])
    telo_df["len"] = telo_df["end"] - telo_df["start"] + 1
    telo_df = telo_df.astype({"start": "Int64", "end": "Int64", "len": "Int64"})

    telo_df_merged = pd.concat([r[2] for r in results])
    telo_df_merged["len"] = telo_df_merged["end"] - telo_df_merged["start"] + 1
    telo_df_merged = telo_df_merged.astype(
        {"start": "Int64", "end": "Int64", "len": "Int64", "chrom_size": "Int64"}
    )
    telo_df_merged = telo_df_merged[
        ["strain", "chrom", "side", "type", "start", "end", "len", "chrom_size"]
    ]

    return raw_df, telo_df, telo_df_merged


def run_on_fasta_dir(
    fasta_dir_path, polynuc_thres, entropy_thres, nb_scanned_nt, threads
):
    """Run iteratively the telemore detection algorithm on all fasta files in a directory

    :param fasta_dir: path to fasta directory
    :return: a tuple of df, telo_df and telo_df_merged
    """
    raw_dfs = []
    telom_dfs = []
    merged_telom_dfs = []

    for ext in ["*.fasta", "*.fas", "*.fa", "*.fsa"]:
        for fasta in fasta_dir_path.glob(ext):

            raw_df, telom_df, merged_telom_df = run_on_single_fasta(
                fasta, polynuc_thres, entropy_thres, nb_scanned_nt, threads
            )
            raw_dfs.append(raw_df)
            telom_dfs.append(telom_df)
            merged_telom_dfs.append(merged_telom_df)

    total_raw_df = pd.concat(raw_dfs)
    total_telom_df = pd.concat(telom_dfs)
    total_merged_telom_df = pd.concat(merged_telom_dfs)

    return total_raw_df, total_telom_df, total_merged_telom_df


def run_telofinder(fasta_path, polynuc_thres, entropy_thres, nb_scanned_nt, threads, raw):
    """Run telofinder on a single fasta file or on a fasta directory"""
    fasta_path = Path(fasta_path)

    if fasta_path.is_dir():
        print(
            f"Running in iterative mode on all '*.fasta', '*.fas', '*.fa', '*.fsa' files in '{fasta_path}'"
        )
        raw_df, telom_df, merged_telom_df = run_on_fasta_dir(
            fasta_path, polynuc_thres, entropy_thres, nb_scanned_nt, threads
        )
        export_results(raw_df, telom_df, merged_telom_df, raw)
        return raw_df, telom_df, merged_telom_df

    elif fasta_path.is_file():
        print(f"Running in single fasta mode on '{fasta_path}'")

        raw_df, telom_df, merged_telom_df = run_on_single_fasta(
            fasta_path, polynuc_thres, entropy_thres, nb_scanned_nt, threads
        )
        export_results(raw_df, telom_df, merged_telom_df, raw)
        return raw_df, telom_df, merged_telom_df
    else:
        raise IOError(f"'{fasta_path}' is not a directory or a file.")


# Main program
if __name__ == "__main__":
    args = parse_arguments()
    output_dir_exists(args.force)
    run_telofinder(
        args.fasta_path,
        args.polynuc_threshold,
        args.entropy_threshold,
        args.nb_scanned_nt,
        args.threads,
        args.raw,
    )
    # export_raw_df(args.raw)
