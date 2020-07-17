import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="A program for telomere characterization"
    )
    parser.add_argument("fasta", help="Path to the fasta file")
    parser.add_argument("-o", "--output", help="Path to the output file")

    return parser.parse_args()


def get_genome_name(filename):
    filepath = Path(filename)
    test = 0.1235

    print(f"This is the parent: {filepath.parent}")
    print(f"{test:.2f}")
    print(f"{filepath.name}")
    return filepath.stem


if __name__ == "__main__":
    args = parse_args()
    print(get_genome_name(args.fasta))
