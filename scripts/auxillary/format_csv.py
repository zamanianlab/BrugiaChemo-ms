#Author: Nic Wheeler
#!/bin/python


import collections
from collections import Counter
from sys import argv 
script, filename1, filename2 = argv 
import re
import numpy as np
import subprocess
import csv


def format(input, output):
	with open(input, "r") as f, open(output, "w") as o:
		for line_no, line in enumerate(f):
			line = line.strip()
			clade_list = []
			clade_list = line.split("\t")
			for gene in clade_list:
				o.write(gene + "\t" + str(line_no + 1) + "\n")


run = format(filename1, filename2)




