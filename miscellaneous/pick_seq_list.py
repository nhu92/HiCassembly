# Pick sequences in FASTA files based on a list file
# Nan Hu
# Feb. 11, 2020

# import modules
from Bio import SeqIO
import argparse

# argparse claims
Argparse = argparse.ArgumentParser()
Argparse.add_argument("-f", "--IN", required = True, help = "Input FASTA file.")
Argparse.add_argument("-l", "--INLIST", required = True, help = "Input list of sequence name.")
Argparse.add_argument("-o", "--OUTPUT", required = True, help = "Output FASTA file.")
argv = vars(Argparse.parse_args())

# Read FASTA files and write in file
wanted = [line.strip() for line in open(argv["INLIST"])]
seqiter = SeqIO.parse(open(argv["IN"]), 'fasta')
seqID = []
seqSeq = []
for seq in seqiter:
	seqID.append(seq.id)
	seqSeq.append(seq.seq)
with open(argv["OUTPUT"], "wt") as f:
	for item in wanted:
		for pos in range(0,len(seqID)):
			if item == seqID[pos]:
				f.write(">"+seqID[pos]+"\n")
				f.write(str(seqSeq[pos])+"\n")
