#!/usr/bin/env python

import sys
from Bio import SeqIO

with open('unique_' +  sys.argv[1], 'w+') as outFile:
    record_ids = list()
    for record in SeqIO.parse(sys.argv[1], 'fasta'):
        if record.id not in record_ids:
            record_ids.append( record.id)
            SeqIO.write(record, outFile, 'fasta')
