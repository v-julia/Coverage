# Coverage
This folder contains python scripts that calculate coverage of reference nucleotide sequence by some set of sequences.

In paper "The effect of sample bias and experimental artifacts on statistical phylogenetic analysis of picornaviruses" (in review) these scripts were used to show the coverage of most common picornavirus species genomes by sequences availabe in GenBank Nucleotide database.

## Requirements

Python libraries:

* numpy==1.16.4
* pandas==0.24.2
* matplotlib==3.1.0
* Bio==0.1.0

Other programms:

* standalone blast
* mafft


## genbank_coverage.py


Performs standalone blast of sequences from *input_file* against reference sequence, saves blast output table in 'blast_out' folder.
Checks whether sequence from *input_file* overlaps with reference sequence, calculates coverage for each position of reference, plots coverage, 
saves coverage values in the directory of *input_file*. If *input_file* doesn't exist, loads sequences from GenBank Nucleotide database 
using *query*, saves them in fasta-format in current working directory.


### Usage

```
genbank_coverage.py [-h] [-i INPUT_FILE] [-q QUERY] -ref REFERENCE
                       -path_blast PATH_BLAST

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input_file INPUT_FILE
                        File with sequences in fasta or genbank format for
                        calculating coverage
  -q QUERY, --query QUERY
                        Query for GenBank Nucleotide database
  -ref REFERENCE, --reference REFERENCE
                        Path to file with reference sequence in fasta-format
                        or GenBank accession number of reference sequence
  -path_blast PATH_BLAST, --path_blast PATH_BLAST
                        Path to blast program
```


##  merge_coverages.py


Calculates coverage using several blast outputs (blasts against several reference sequences) which 
are saved in *input_dir*. For each query sequence finds the blast table with maximal hit length,
 merges rows with maximal hit lengths from different blast tables, calculates final coverage, 
 draws a coverage plot, saves coverage values to output dir.

### Usage
```
merge_coverages.py [-h] -i INPUT_DIR -o OUTPUT_DIR -al ALIGNMENT -t
                      TITLE

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_DIR, --input_dir INPUT_DIR
                        Input directory with blast output tables
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory
  -al ALIGNMENT, --alignment ALIGNMENT
                        Path to file with alignment of reference sequences
  -t TITLE, --title TITLE
                        Title of output file figure

```
### Requirements:

* argparse
* Biopython
* matplotlib
* subprocess
* pandas
* numpy
* pathlib
* time
* genbank_coverage
* standalone blast
* mafft

                        
##  genbank_coverage_complex.py

For each sequences in *input_file* (or sequences loaded from GenBank nucleotide using query) 
performs standalone blast against reference sequences from *reference*, saves blast output table in 'blast_out' folder.
Then merges results using merge_coverage.py. Draws a coverage plot and saves coverage values to the directory of input file.

### Usage
```
genbank_coverage_complex.py [-h] [-i INPUT_FILE] [-q QUERY] -ref
                               REFERENCE -t TITLE -path_blast PATH_BLAST
                               -path_mafft PATH_MAFFT

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input_file INPUT_FILE
                        File with sequences in fasta or genbank format for
                        calculating coverage
  -q QUERY, --query QUERY
                        Query for GenBank Nucleotide database
  -ref REFERENCE, --reference REFERENCE
                        List with reference sequences'ids (with version!)
                        separated by comma. Example: AB084913.1,AB0834732.2
  -t TITLE, --title TITLE
                        Title for coverage plot
  -path_blast PATH_BLAST, --path_blast PATH_BLAST
                        Path to blast program
  -path_mafft PATH_MAFFT, --path_mafft PATH_MAFFT
                        Path to mafft program
```
