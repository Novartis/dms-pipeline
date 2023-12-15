#!/usr/bin/env python
# coding: utf-8

# In[43]:



'''
Creates a table with entries containing unique sequences and their counts per sample.
Optionally checks if the unique sequences are expected or not

'''

import os
import sys
import csv
import numpy as np
from Bio import SeqIO
from collections import defaultdict


log_file = open(snakemake.log.log1,"w")
sys.stdout = log_file

# check for correct file of expected mutants
expected_mutants = False
if os.path.isfile(snakemake.params.expected_mutants):
    if snakemake.params.expected_mutants.endswith('.txt'):
        # check for column 'amplicon'
        with open(snakemake.params.expected_mutants) as csvfile:
            csvrows = csv.reader(csvfile, delimiter='\t')
            fieldnames=np.array(next(csvrows))
            csvfile.close()
            if 'amplicon' in fieldnames:
                expected_mutants = True
            else:
                print('The file containing the expected mutant sequences does not have a column called "amplicon".')
                print("Please specify a tab separated .txt file with at least a column called 'amplicon'.")
                print("The analysis will be performed without counting the expected mutant sequences.")
    else:
        print('The file containing the expected mutant sequences has an incorrect file format.')
        print("Please specify a tab separated .txt file with at least a column called 'amplicon'.")
        print("The analysis will be performed without counting the expected mutant sequences.")
else:
    print("File with expected mutants was not provided or does not exist.")
    print("The analysis will be performed without counting the expected mutant sequences.")



# Count up reads
fasta_sequences = SeqIO.parse(open(snakemake.input.sequences),'fasta')

unique_seq = defaultdict(int)
for fasta in fasta_sequences:
    sequence = str(fasta.seq)
    unique_seq[sequence] += 1

# Read in expected mutants
if expected_mutants:
    #TODO: implement input for different file formats
    with open(snakemake.params.expected_mutants) as csvfile:
        csvrows = csv.reader(csvfile, delimiter='\t')
        fieldnames=np.array(next(csvrows))
        seq_idx = int(np.where(fieldnames == 'amplicon')[0])
        mutants = []
        for row in csvrows:
            mutants.append(row[seq_idx])
        csvfile.close()

    with open(snakemake.output.counts, 'w') as csv_file: 
        writer = csv.writer(csv_file)
        writer.writerow(['seq', 'count', 'expected'])
        for key, value in unique_seq.items():
            if key in mutants:
                writer.writerow([key, value, True])
            else:
                writer.writerow([key, value, False])
        csv_file.close()

else:
    with open(snakemake.output.counts, 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['seq', 'count'])
        for key, value in unique_seq.items():
            writer.writerow([key, value])
        csv_file.close()

log_file.close()

