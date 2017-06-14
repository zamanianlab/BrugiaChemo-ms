#!/bin/python

import collections
from collections import Counter
from sys import argv 
script, filename1, filename2, filename3 = argv 
import re
from Bio import SeqIO


p = open(filename1) #list of GeneIDs in list/txt format
record_dict = SeqIO.index(filename2, "fasta") #proteome fasta file
output_handle = open(filename3, "w") #output file of exracted sequences

def seqextract(file):
	for line in file:
		ids = []
		line = line.strip()
		ids.append(line)
		print ids
		for entry in ids:
			current = record_dict[entry]
			SeqIO.write(current, output_handle, "fasta")
		else: 
			pass

run = seqextract(p)
output_handle.close()