# Annotate_nucleotide_sequences

usage: ./blast_annotate_fasta.py [-h] -s input_sequences -r input_references
                                 -o output_annotations_file

Script for annotating multiple sequences in fasta format using a list of
reference genomes in GenBank format (downloaded from NCBI or in NCBI-compliant
format)

optional arguments:
  -h, --help            show this help message and exit
  -s input_sequences, --sequences input_sequences
                        Input (multi-) fasta file containing nucleotide
                        sequences to be annotated
  -r input_references, --references input_references
                        Input file containing locations to the reference
                        genomes to be used for annotation (one per line)
  -o output_annotations_file, --output output_annotations_file
                        Output file containing the annotated nucleotide
                        sequences

Written by Chrispin Chaguza, Yale School of Public Health, Yale University,
2021
