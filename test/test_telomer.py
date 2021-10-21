from . import test_dir

import telofinder.telofinder as tf

filename = f"{test_dir}/data/AFH_chrI.fasta"


def test_run_on_single_fasta():
    df = tf.run_on_single_fasta(filename, 0.8, 0.8, 8000, 1)
