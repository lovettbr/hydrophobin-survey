#!/usr/bin/env python

import glob

fastafiles = []
for file in glob.glob("*.faa"):
    fastafiles.append(file)

files=fastafiles

for input in files:

	from Bio import SeqIO
	import pandas as pd 
	import numpy as np 

	for seq_record in SeqIO.parse(input, "fasta"):
		ids = seq_record.id
		C = seq_record.seq.count("C")
		Length = len(seq_record.seq)
		Perc = C/Length
		seq = seq_record.seq
		print(input, ids, C, Length, Perc, seq)
