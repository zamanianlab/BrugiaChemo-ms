#!/Users/njwheeler/software/miniconda3/bin/python

from sys import argv
import os.path
import pathlib
from Bio import SeqIO

script, filename1 = argv


# Add species labels to each ID in the FASTA file


def headerchange(file):
    with open(filename1) as original:
        path = pathlib.Path(filename1)
        base = path.name
        stem = path.stem
        genus, species = base.split("_", 1)
        label = genus[0] + species[0:4]
        records = SeqIO.parse(original, 'fasta')
        for record in records:
            # print record.id
            record.id = label + "-" + record.id
            record.description = ""
            # print record.id
            with open(pathlib.Path.joinpath(path.parent, stem + "_label.fa"), 'a') as corrected:
                SeqIO.write(record, corrected, 'fasta')
                # print(pathlib.Path.joinpath(path.parent, stem + "_label.fa"))


run = headerchange(filename1)
