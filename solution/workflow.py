#! /usr/bin/env python
import os
import sys
import subprocess
import re

def codons(txtfile):
    """
    Extracts the amino acid and their corrosponding codons in a dictionary
    with key as single letter amino acid and value as codon
    :param txtfile: .txt file path as a string
    :type str
    """
    codon_to_aa = {}
    with open(txtfile) as aa_file:
        for line in aa_file:
            if not line.startswith('na'):
                row = line.rstrip("\n").split('\t')
                codon_to_aa[row[0]] = row[1]
        return codon_to_aa

def rev_comp(seq):
    """
    Generates and return reverse complement of the sequnece
    :param seq: a string sequence
    :type str
    """
    # Dictionary with nucleotide complements for generating reverse complement of the given sequence
    nuc = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    rev_seq = (seq[::-1])
    rc = ''.join([nuc.get(rev_seq[i], 'X') for i in range(len(seq))])
    return rc

def translation(sequence, codon2aa):
    """
    Function that performs translation of mRNA sequence using single letter amino acids from codon2aa
    :param sequence: a string sequence
    :type str
    :param codon2aa: amino acid and their corrosponding codons
    :type dictionary
    """
    end = len(sequence) - (len(sequence) % 3) - 1
    translate = ''.join([codon2aa.get(sequence[i:i+3], 'X') for i in range(0, end, 3)])
    return translate

def process_gff(gff_file, input_file, output_genes_file):
    """
    check the intervals.gff if it has a header and creates a temp gff file without header
    :param gff_file: .gff file path as a string
    :type str
    :param input_file: .fasta file path as a string
    :type str
    :param output_genes_file: .fasta file where output from bedtools is saved
    :type str
    """
    seq_header = []

    with open(input_file) as f:
        for line in f.readlines():
            if line.startswith(">"):
                seq_header.append(line.strip(">").strip("\n"))

    with open(gff_file, 'r') as fin:
        first = (fin.readline()).split("\t")[0]
        data = fin.read().splitlines(True)
        if first in seq_header: # GFF file does not contain a header line - using same file for processing
            gene_fasta = "  ".join(["bedtools getfasta  -fi", input_file, "-bed", gff_file, "-name  -fo", output_genes_file])
        else: # GFF file contains a header line - creating a temperory no_header_gff_file file and using it for processing
            no_header_gff_file = gff_file.replace('.gff','_noheader.gff')
            with open(no_header_gff_file, 'w') as fout:
                fout.writelines(data[:])
            gene_fasta = "  ".join(["bedtools getfasta  -fi", input_file, "-bed", no_header_gff_file, "-name  -fo", output_genes_file])

    subprocess.call(gene_fasta, shell=True)

    #GFF file processing finished...Now deleting the temperory no_header_gff_file file if it exists
    if os.path.exists(no_header_gff_file):
        os.remove(no_header_gff_file)

def write_out(output_folder_path, data):
    """
    writing protein sequences with their headers into output file,
    which conatins the longest protein from each mRNA sequence
    :param output_folder_path: folder location to write the output protein to
    :type str
    :param data: protein sequence and its headers
    :type list
    """
    with open(output_folder_path + "genes.fasta", "w") as genes_fa:
        for item in data:
            if item.endswith("*"):
                item = item[:-1]
                genes_fa.write("%s\n" % item)
            else:
                genes_fa.write("%s\n" % item)

def run(args):
    """
    Function to process the mRNA, generate the six translation frames, translate and then save the
    longest proteins in the output file with the gene name.
    :param args
    :type <class 'argparse.Namespace'>
    """
    protein = []

    # this is the file containing gene sequences that will be created by bedtools
    genes_file = args.output_location + "outfile.fasta"

    process_gff(args.intervals_file, args.input, genes_file)

    codon_aa = codons(args.aa_codons)

    with open(genes_file) as f_in:
        for line in f_in.readlines():
            aa_seqs = []
            mRNA_seqs = []
            if line[0] == '>':
                protein.append(line.rstrip('\n').split(":")[0])
            else:
                for i in range(3): # get all the six open reading frames (ORFs)
                    mRNA_seqs.append(line.rstrip("\n")[i:])
                    mRNA_seqs.append(rev_comp(line.rstrip("\n")[:len(line)-(i+1)]))
                for j in range(6): # search for dna sequences from start codon to stop codon in  all the ORFs of each mRNA sequence
                    seq_len = len(mRNA_seqs[j])
                    for k in range(0, seq_len, 3):
                        codon = mRNA_seqs[j][k:k+3]
                        if (codon == "ATG"): # pattern starting with start codon and ending with any of the stop codons
                            gene = re.findall(r'ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)', mRNA_seqs[j][k:])
                            # translation into protein
                            aa_seqs += [translation(g, codon_aa) for g in gene]
                            break
                protein.append(max(aa_seqs, key=len))

    write_out(args.output_location, protein)

    # Performing cleanup
    if os.path.exists(args.input + '.fai'):
        os.remove(args.input + '.fai')
    if os.path.exists(genes_file):
        os.remove(genes_file)
