# dubii_project

## DATA
### Organism: Saccharomyces cerevisiae
### Genome dataset is composed of:
  - 99 de novo genome assemblies from novel Oxford Nanopore sequencing data
  - 18 de novo genome assemblies from published Oxford Nanopore sequencing data
  - 19 de novo genome assemblies from published PacBio sequencing data
  - 1 refernce genome S288C
### SV list as a csv file
### Annotation files (gff3)

## Analyses
### Telomere sequence length
### tRNA gene families
### Retrotransposons (Ty)

## Programming

### Refactoring

#### In progress / Done

- Setting up a good developping environment (SublimeText + plugins)

- Improved CLI:
  - Replacing sys.argv by argparse.
  - Argparse help.
  - See Different types of arguments/options.
  
- Refactoring spaghetti code into specialized functions + docstrings.

- Accepting single fasta file or recursively applied to fasta files in a
  directory.
  
- Usage of pathlib.

#### To do

- shared conda environment 

- algorithm study/improvements
  - GC/AT skewness
  - Entropy/complexity by window
  - window size as a parameter
  - number of consecutive failures in pattern matching

- tests

if time:
- logging in python

## References

### Other Telomere efforts

TelSeq: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4027178/

Telomere to Telomere consortium (human): https://sites.google.com/ucsc.edu/t2tworkinggroup/

Telomere motif identification: https://www.future-science.com/doi/10.2144/btn-2018-0057

TelomereHunter (cancer related telomere variants): https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2851-0

### Text Editor

https://code.visualstudio.com/

https://www.jetbrains.com/pycharm/

### Formating

https://github.com/jgirardet/sublack

