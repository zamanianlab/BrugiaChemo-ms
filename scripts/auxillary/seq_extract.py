#!/Users/njwheeler/software/miniconda3/bin/python

from sys import argv
from Bio import SeqIO

script, filename1, filename2, filename3 = argv

p = open(filename1)  # list of GeneIDs in list/txt format
record_dict = SeqIO.index(filename2, "fasta")  # proteome fasta file
output_handle = open(filename3, "w")  # output file of extracted sequences


def seqextract(file):
    for line in file:
        ids = []
        line = line.strip()
        ids.append(line)
        # print ids
        for entry in ids:
            # print entry
            try:
                current = record_dict[entry]
                SeqIO.write(current, output_handle, "fasta")
            except KeyError:
                print("Could not find {} in {}.".format(entry, filename2))
            else:
                pass


run = seqextract(p)
output_handle.close()
