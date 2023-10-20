from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pysam
import random
from itertools import combinations
from difflib import SequenceMatcher
import subprocess

from Bio.Seq import Seq
import sequana
from collections import Counter


def get_read_anchor_position(bam, reference_position, out_fas, reverse_complement=False, max_size=1000):
    """From pysam.AlignedSegment.get_aligned_pairs,
    get the position of the read which are corresponding
    to a position in the reference sequence

    TODO: Have to deal with read not aligning exactly on 'reference_position'
    """
    stats_dict = {"filtered_too_long": 0, "filtered_not_match_ref_pos": 0, "ok": 0}

    with open(out_fas, "w") as fas:
        with pysam.AlignmentFile(bam, "rb") as bam:

            for rd in bam:
                for read_pos, ref_pos in rd.get_aligned_pairs():
                    if ref_pos == reference_position:
                        if not read_pos:
                            stats_dict["filtered_not_match_ref_pos"] += 1
                        elif read_pos > max_size:
                            stats_dict["filtered_too_long"] += 1
                        else:
                            if reverse_complement:
                                seq = str(Seq(rd.query_sequence[:read_pos]).reverse_complement())
                            else:
                                seq = rd.query_sequence[:read_pos]

                            fas.write(f">{rd.query_name}\n{seq}\n")
                            stats_dict["ok"] += 1
                            break
    print(stats_dict)
    return stats_dict


def chunk_fasta(fasta, start=0, window=10, step=1, n_sample=20, do_sample=False):

    if do_sample:
        outdir = Path(f"chunks_sample_{n_sample}") / fasta.stem
    else:
        outdir = Path(f"chunks") / fasta.stem

    outdir.mkdir(parents=True, exist_ok=True)

    fas = sequana.fasta.FastA(fasta)
    # Check how many sequences will be "alignable"
    ali_seqs = [x for x in fas if len(x.sequence) >= start + window + 1]
    while len(ali_seqs) > n_sample:
        # Randomly sample n_sample sequences
        if do_sample:
            sample_seqs = random.sample(ali_seqs, n_sample)
        else:
            sample_seqs = ali_seqs

        out_fas = outdir / f"{start}.fas"

        with open(out_fas, "w") as f:

            for seq in sample_seqs:
                f.write(f">{seq.name}\n{seq.sequence[start:start + window]}\n")

        start += step
        ali_seqs = [x for x in fas if len(x.sequence) > start + window + 1]

    return outdir


def align_fastas(fasta_dir, threads):
    """Align all fastas in fasta_dir using muscle."""

    outdir = Path(fasta_dir).parent / (Path(fasta_dir).stem + "_aligned")
    outdir.mkdir()

    for fas in Path(fasta_dir).glob("*.fas"):
        output = subprocess.run(
            f"mafft --thread {threads} {fas} > {outdir / fas.name}",
            shell=True,
            capture_output=True,
        )
    return outdir


def get_diversity(fasta, col_name):
    """Get nucleotiqe diversity (pi) on the whole fasta"""
    div_dict = {}
    div_dict["start"] = int(fasta.stem)

    fas = sequana.fasta.FastA(fasta)
    seqs = [x.sequence for x in fas]

    seq_counter = Counter(seqs)

    divs = []
    for seq1, seq2 in combinations(set(seqs), 2):
        freq_seq1 = seq_counter[seq1] / len(seqs)
        freq_seq2 = seq_counter[seq2] / len(seqs)
        divs.append(SequenceMatcher(a=seq1, b=seq2).ratio() * freq_seq1 * freq_seq2)

    pi = sum(divs) * (len(seqs) / (len(seqs) + 1))

    div_dict[f"div_{col_name}"] = pi
    div_dict[f"n_{col_name}"] = len(seqs)

    return div_dict


def get_diversity_old(fasta):
    """This depends on current naming of the fastas !!
    ie as test_10.fas with 10 being the start"""

    div_dict = {}
    div_dict["start"] = int(fasta.stem)

    fasta = sequana.fasta.FastA(fasta)
    seqs = [x.sequence for x in fasta]
    div_dict["div"] = len(Counter(seqs)) / len(fasta)
    div_dict["n"] = len(seqs)

    return div_dict


def get_all_diversities(fastas, col_name):
    df_list = [get_diversity(fas, col_name) for fas in fastas]
    df = pd.DataFrame(df_list).set_index("start").sort_index()

    return df


def main(fasta, threads):
    fasta = Path(fasta)
    outdir = chunk_fasta(fasta)
    outdir_ali = align_fastas(outdir, threads)
    df_no_sample = get_all_diversities(outdir.glob("*.fas"), "no_sample")
    df_no_sample_ali = get_all_diversities(outdir_ali.glob("*.fas"), "no_sample_ali")

    outdir_sample = chunk_fasta(fasta, do_sample=True)
    outdir_sample_ali = align_fastas(outdir_sample, threads)
    df_sample = get_all_diversities(outdir_sample.glob("*.fas"), "sample")
    df_sample_ali = get_all_diversities(outdir_sample_ali.glob("*.fas"), "sample_ali")

    df = pd.concat([df_no_sample, df_no_sample_ali, df_sample, df_sample_ali], axis=1)

    sns.lineplot(data=df.loc[:, df.columns.str.startswith("div")], alpha=0.5)
    return df
