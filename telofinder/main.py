import argparse
import os
from pathlib import Path

from telofinder.telofinder import (run_on_single_seq, run_on_fasta_dir, 
    run_on_single_fasta, export_results)


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

def main():
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

# Main program
if __name__ == "__main__":
    main()
