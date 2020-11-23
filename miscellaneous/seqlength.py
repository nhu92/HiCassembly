# Fetch seq length form FASTA
# Nan Hu
# Oct. 13, 2020

# import modules
from Bio import SeqIO
import argparse

# argparse claims
Argparse = argparse.ArgumentParser()
Argparse.add_argument("-f", "--IN", required = True, help = "Input FASTA file.")

argv = vars(Argparse.parse_args())

# Read in FASTA and output ID\tLength
                        
seqiter = SeqIO.parse(open(argv["IN"]), 'fasta')
for seq in seqiter:
	print(seq.id+"\t"+str(len(seq.seq)))
