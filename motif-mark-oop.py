import argparse
import cairo

parser = argparse.ArgumentParser(description='Visualize sequence motifs on DNA sequences.')
parser.add_argument('-f', '--fasta', help='fasta file (required)')
parser.add_argument('-m', '--motifs', help='motif file')

args = parser.parse_args()

print (args.motifs)