import sys
import os
from traceback import format_exc
from argparse import ArgumentParser
from solution.workflow import run


def main(argv=None):

    program_name = os.path.basename(sys.argv[0])
    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    try:
        parser = ArgumentParser(description='Generate protein from the given fasta sequences')

        parser.add_argument('--input', required=False, help='A FASTA file with two chromosome sequences', default='./input/sequences.fasta')

        parser.add_argument('--intervals_file', required=False, help='A GFF file with a set of mRNA coordinates', default='./input/intervals.gff')

        parser.add_argument('--aa_codons',required=False, help='A tab delimited table mapping codons to single letter amino acids', default='./input/standard_code.txt')

        parser.add_argument('--output_location',required=False, help='Location where output files are saved', default='./output/')

        args = parser.parse_args()

        print("Translating the genes sequences")

        run(args)

        print("Task accomplished! Output files in {} location".format(args.output_location))

    except Exception as e:
        print(program_name + ": " + repr(e) + '\n' + format_exc() + '\n')
        raise(e)


if __name__ == "__main__":
    sys.exit(main())
