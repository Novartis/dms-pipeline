#!/usr/bin/env python
# coding: utf-8

# In[ ]:


'''
Creates a n x m count matrix with n=number of unique sequences, m=number of samples.

In case expected mutants are part of the input a second count matrix is generated that contains counts of expected sequences only.

'''

import os
import sys
import csv
import numpy as np
from collections import defaultdict
from collections import OrderedDict


log_file = open(snakemake.log.log1,"w")
sys.stdout = log_file

samples = snakemake.params.samples
OUTPUT_DIR = snakemake.params.output_dir
expected_mutants = False



def initialize_dict(m):
    return lambda: [0] * m


# read counts for each sample into dictionary
counts = defaultdict(initialize_dict(len(samples)))
for i, s in enumerate(samples):
    with open(snakemake.input.counts[i]) as csv_file:
        reader = csv.reader(csv_file)
        header = np.array(next(reader))
        # print(s)
        if 'expected' in header:
            expected_mutants = True
        for row in reader:
            counts[row[int(np.where(header == 'seq')[0])]][i] = int(row[int(np.where(header == 'count')[0])])
    

# order sequences in dictionary
counts = OrderedDict(sorted(counts.items()))
    
    
# print dictionary to csv file
with open(snakemake.output.count_matrix, 'w') as csv_file:
        writer = csv.writer(csv_file)
        samples.insert(0, 'seq')
        writer.writerow(samples)
        for key, value in counts.items():
            value.insert(0, key)
            writer.writerow(value)
       
if expected_mutants:
    with open(snakemake.params.expected_mutants) as csvfile:
        csvrows = csv.reader(csvfile, delimiter='\t')
        fieldnames=np.array(next(csvrows))
        seq_idx = int(np.where(fieldnames == 'amplicon')[0])
        mutants = []
        for row in csvrows:
            mutants.append(row[seq_idx])
    
    with open('%s/count_matrix_expected_mutants.csv' % OUTPUT_DIR, 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(samples)
        for key, value in counts.items():
            if key in mutants:
                writer.writerow(value)


log_file.close()

