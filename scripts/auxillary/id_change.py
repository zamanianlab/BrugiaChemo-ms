#!/bin/python

import collections
from collections import Counter
from sys import argv 
script, filename1, filename2 = argv 
import re, os.path
from Bio import SeqIO

# Add species labels to each ID in the FASTA file

def headerchange(file):
	with open(filename1) as original, open(filename2, 'w') as corrected:
		base =  os.path.basename(filename1)
		genus, species = base.split("_", 1)
		label = genus[0] + species[0:4] 
		records = SeqIO.parse(original, 'fasta')
		for record in records:
			#print record.id
			record.id = label + "-" + record.id
			record.description = ""
			#print record.id
			SeqIO.write(record, corrected, 'fasta')

run = headerchange(filename1)

