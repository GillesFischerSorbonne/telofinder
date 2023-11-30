# telofinder

[![Tests](https://github.com/GillesFischerSorbonne/telofinder/actions/workflows/main.yml/badge.svg)](https://github.com/GillesFischerSorbonne/telofinder/actions/workflows/main.yml)


A python package to determine the location and the size of telomeric repeats (both terminal and internal) from genome assemblies.

Telomere detection is based on calculation in a 20 bp sliding window of the following two metrics:

- DNA sequence entropy < entropy_threshold (default=0.8)  
- proportion of polynucleotides (default_list = ["CC", "CA", "AC"]) >  polynuc_threshold (default=0.8) 

The complete documentation is available at [Read The Docs](https://telofinder.readthedocs.io/en/latest/) 

## Usage

    telofinder fasta_path_to_file(s)


## Help

    telofinder.py --help


## Output

Telofinder outputs a directory called `telofinder_results` including 3 csv and 2 bed files containing the telomere calls and their coordinates, either as raw output or after merging consecutive calls

## Reference

Please cite O’Donnell, S., Yue, JX., Saada, O.A. et al. Telomere-to-telomere assemblies of 142 strains characterize the genome structural landscape in Saccharomyces cerevisiae. Nat Genet 55, 1390–1399 (2023). [DOI](https://doi.org/10.1038/s41588-023-01459-y)

