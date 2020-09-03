# Bioinformatics unit coding challenge - Satwant Kaur

## Extract & Translate genes with coordinates from their chromosome sequences.

- Input files (found in *input* folder):
  - sequences.fasta - A FASTA file with two chromosome sequences.
  - intervals.gff - A GFF file with a set of mRNA coordinates for the sequences in    sequences.fasta.
  - standard_code.txt - A tab delimited table mapping codons to single letter amino acids.

Please use this folder to substitute any of the relevant files for testing.

- Expected output (found in *output* folder):
  genes.fasta - The longest translation of each mRNA sequence to standard output in FASTA file format.

## Pre-requisite

- Install [Homebrew](https://brew.sh/)

- Install bedtool using following [link](https://bedtools.readthedocs.io/en/latest/content/installation.html) and choosing the command that matches your operating system.
  - Gotcha :  If using MAC OSX, you can ignore "brew tap homebrew/science" command.

- The current working directory contains **run_task.py**, which is a command line tool, which takes following inputs
  - --input ./input/sequences.fasta (default)
  - --intervals_file ./input/intervals.gff (default)
  - --aa_codons ./input/standard_code.txt (default)
  - --output_location ./output/ (default)
and generates *genes.fasta* in output folder.

- Internally, the script calls workflow.py present in **solution** folder. In case the user wants to use differently named input files, the user can enter the file names using  arguments. You can find more information by executing *python run_task.py -h* from current folder.

## Brief Code Explanation
- The script searches for header in bed file (intervals.gff) and if present re-writes the intervals.gff in a temporary file without the header, which is later deleted. Taking sequences.fasta and intervals.gff (without header) as input the script calls BedTools getfasta module and generates the mRNA sequences, which are saved to the outfile.fasta.
- The mRNA sequences are stored in a list (mRNA_seqs) in forward frame and reversed frame (rev_comp function). Six reading frames are generated  and for each mRNA sequence, the script searches for coding genes using regular expression pattern. Translation function is used to translate the coding genes and the longest of them is appended to the list called protein.
- An output file called genes.fasta is opened in write mode and the longest translated sequence for each mRNA sequence is written to standard output in FASTA file format.
