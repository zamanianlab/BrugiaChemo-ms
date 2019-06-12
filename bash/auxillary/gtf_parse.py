#Author: Mostafa Zamanian 
#!/bin/python
import collections
from collections import Counter
from sys import argv 
script, filename = argv 
import re
import numpy as np
import subprocess
import csv

g = filename #caenorhabditis_elegans.PRJNA13758.WBPS9.canonical_geneset.gtf

#python gtf_parse.py caenorhabditis_elegans.PRJNA13758.WBPS9.canonical_geneset.gtf > GTF_pseudogenes.txt

#find pseudogenes genes: 
#IV	WormBase	gene	64400	66128	.	-	.	gene_id "WBGene00235259"; gene_source "WormBase"; gene_biotype "pseudogene";
#V	WormBase	transcript	2007360	2007709	.	+	.	gene_id "WBGene00023306"; transcript_id "F59A7.10"; gene_source "WormBase"; gene_biotype "pseudogene"; transcript_source "WormBase"; transcript_biotype "pseudogene";


#pull out pseudogenes
def pseudoparse(file):
	with open(file, "r") as f:
		for line in f:
			line = line.strip()
			match = re.search('([I|V|X]+)\s+WormBase\s+transcript\s+(\d+)\s+(\d+).*?([-|+]).*?gene_id\s+\"(WB.*?)\".*?transcript_id\s+\"(.*?)\".*?gene_biotype\s+\"(pseudogene)\"', line)
			if match:
				chrom = match.group(1)
				start = match.group(2)
				stop = match.group(3)
				strand = match.group(4)
				WBid = match.group(5)
				Tid = match.group(6)
				biotype = match.group(7)
				out_line = "{chrom} {start} {stop} {strand} {WBid} {Tid} {biotype}".format(**locals())
				print out_line
			else: 
				pass

pseudoparse(g)




