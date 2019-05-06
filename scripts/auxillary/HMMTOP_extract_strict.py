#!/Users/njwheeler/software/miniconda3/bin/python

import collections
from collections import Counter
from sys import argv
import re
from Bio import SeqIO

script, filename1, filename2, filename3 = argv

p = open(filename1)  # HMMTOP output file
record_dict = SeqIO.index(filename2, "fasta")  # FASTA file
outfile = open(filename3, "w+")  # Filename for extracted sequences


def seqextract(file):
    for line in file:
        ids = []
        line = line.strip()
        # print line
        match = re.search('>HP:\s+\d+\s+(.*?)\s.*?[IN|OUT]\s+(\d+)', line)
        if match and int(match.group(2)) == 7:
            protein_id = match.group(1)
            print(protein_id + " " + match.group(2))
            current = record_dict[protein_id]
            SeqIO.write(current, outfile, "fasta")
        else:
            pass


def TMseqextract(file):
    for line in file:
        ids = []
        line = line.strip()
        # print line
        match = re.search(
            '>HP:\s+\d+\s+(.*?)\s.*?[IN|OUT]\s+(\d+)\s+(.+)', line)
        if match and int(match.group(2)) > 6 and int(match.group(2)) < 8:
            protein_id = match.group(1)
            TM_count = match.group(2)
            TM_coords = match.group(3)
            current = record_dict[protein_id]
            current_id = str(current.id)
            current_seq = str(current.seq)
            match_coords = re.findall('(\d+)\s+(\d+)', TM_coords)
            TM_seq = str('')
            for i in range(1, len(match_coords)):
                # print match_coords
                coord1 = int(match_coords[i][0]) - 6
                coord2 = int(match_coords[i][1]) + 6
                TM_segment = current_seq[coord1:coord2] + '-'
                # print TM_segment
                TM_seq += TM_segment
            else:
                pass
            output_seq = ">" + current_id + "\n" + TM_seq + "\n"
            # print output_seq
            outfile.write(output_seq)

        else:
            pass

# cuts from beginning of first TM to end of last TM (+10 aa on each end)


def TMDextract(file):
    for line in file:
        ids = []
        line = line.strip()
        # print line
        match = re.search(
            '>HP:\s+\d+\s+(.*?)\s.*?[IN|OUT]\s+(\d+)\s+(\d+)\s+.*\s+(\d+)', line)
        if match and int(match.group(2)) == 7:
            protein_id = match.group(1)
            TM_count = match.group(2)
            TM_coords_start = match.group(3)
            TM_coords_end = match.group(4)
            current = record_dict[protein_id]
            current_id = str(current.id)
            current_seq = str(current.seq)
            coord1 = int(TM_coords_start) - 1
            coord2 = int(TM_coords_end) + 1
            TMD = current_seq[coord1:coord2]
            output_seq = ">" + current_id + "\n" + TMD + "\n"
            # print output_seq
            outfile.write(output_seq)

        else:
            pass


run = seqextract(p)
#run = TMseqextract(p)
# run = TMDextract(p)
outfile.close()
